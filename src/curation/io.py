"""Ingestão resiliente, escrita atômica e manifesto de execução.

Duas garantias sustentam este módulo.

**Ingestão não posicional.** O formato ``.smi`` não tem convenção única: os arquivos
de referência do ChEMBL usam ``NOME<TAB>SMILES``, o inverso do padrão Daylight
(``SMILES<espaço>id``). Assumir a estrutura na coluna 0 lê ``Tetrafluoroboranuide``
como molécula. A coluna de estrutura é, portanto, determinada por **parseabilidade**
sobre uma amostra, nunca por posição.

**Atomicidade.** Um lote interrompido não pode deixar para trás um arquivo truncado
indistinguível de um completo — num projeto cuja tese é reprodutibilidade, isso é pior
que estourar a memória. A escrita ocorre em ``.tmp``, a promoção é por ``os.replace``,
e ``manifest.json`` é o **marcador de completude**: sua ausência significa execução
incompleta, independentemente do que exista em disco.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Generator, Iterable, Iterator, Optional, Sequence, TextIO, Union

import rdkit
from rdkit import Chem

import chembl_structure_pipeline as csp

from curation.dedup import CONFLICT_COLUMNS, Collision
from curation.engine import PIPELINE_VERSION
from curation.models import CurationRecord

Source = Union[str, Path, TextIO]

#: Amostra usada para inferir separador, cabeçalho e coluna de estrutura.
_SAMPLE_ROWS = 50

#: Fração mínima de valores parseáveis para uma coluna ser aceita como estrutura.
_SMILES_COLUMN_THRESHOLD = 0.6

_DELIMITERS = (",", "\t", ";", "|")

_SMILES_HEADERS = frozenset(
    {"smiles", "canonical_smiles", "smi", "structure", "mol", "molecule"}
)
_ID_HEADERS = frozenset(
    {"id", "input_id", "compound_id", "cmpd_id", "name", "identifier", "chembl_id"}
)

_ENCODINGS = ("utf-8-sig", "cp1252", "latin-1")


# --- Leitura ---------------------------------------------------------------------


def compute_policy_hash(decisions_path: Union[str, Path]) -> str:
    """SHA-256 de ``docs/decisions.md``.

    É o identificador da política sob a qual um lote foi produzido. Dois lotes com
    hashes diferentes não são comparáveis, mesmo que o código seja idêntico.
    """
    return hashlib.sha256(Path(decisions_path).read_bytes()).hexdigest()


def _detect_encoding(path: Path) -> str:
    """Escolhe a codificação lendo um trecho inicial.

    UTF-8 primeiro, em modo estrito, para que a detecção seja real. ``latin-1``
    encerra a lista porque decodifica qualquer byte — é o fallback que nunca falha,
    ao custo de possivelmente produzir mojibake em campos de texto livre. Nomes
    corrompidos são recuperáveis; uma exceção no meio de um lote de milhões, não.
    """
    probe = path.open("rb").read(65536)
    for encoding in _ENCODINGS:
        try:
            probe.decode(encoding)
            return encoding
        except UnicodeDecodeError:
            continue
    return "latin-1"


def _iter_raw_lines(source: Source) -> Iterator[str]:
    """Normaliza qualquer fonte suportada em linhas de texto.

    Um ``str`` é ambíguo entre caminho e texto colado. A resolução é conservadora:
    só é tratado como caminho se contiver uma única linha e o arquivo existir.
    """
    if isinstance(source, Path):
        encoding = _detect_encoding(source)
        with source.open("r", encoding=encoding, errors="replace", newline="") as handle:
            yield from handle
        return

    if isinstance(source, str):
        candidate = source.strip()
        if "\n" not in candidate and candidate:
            try:
                path = Path(candidate)
                if path.is_file():
                    yield from _iter_raw_lines(path)
                    return
            except OSError:
                pass
        yield from source.splitlines()
        return

    yield from source


def _clean(lines: Iterable[str]) -> Iterator[str]:
    """Remove espaços laterais, linhas em branco e comentários ``#``.

    Comentários iniciados por ``#`` são convenção nos ``.smi`` do ecossistema — o
    próprio leitor da referência os descarta.
    """
    for line in lines:
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        yield text


def _split(line: str, delimiter: Optional[str]) -> list[str]:
    if delimiter is None:
        return line.split()
    return [cell.strip() for cell in line.split(delimiter)]


def _detect_delimiter(sample: Sequence[str]) -> Optional[str]:
    """Escolhe o separador que produz o maior número consistente de colunas.

    ``None`` significa divisão por espaços em branco, o caso comum de ``.smi``.
    """
    best: Optional[str] = None
    best_columns = 1
    for delimiter in _DELIMITERS:
        counts = {line.count(delimiter) for line in sample}
        if len(counts) == 1 and counts != {0}:
            columns = counts.pop() + 1
            if columns > best_columns:
                best, best_columns = delimiter, columns

    if best is not None:
        return best

    whitespace_counts = {len(line.split()) for line in sample}
    if len(whitespace_counts) == 1 and whitespace_counts != {1}:
        return None
    return None


def _parses_as_smiles(text: str) -> bool:
    """Verifica se um valor é interpretável como SMILES.

    ``sanitize=False`` deliberadamente: a decisão aqui é sintática, e valências
    anômalas são problema do motor, não do leitor (D-09).
    """
    if not text:
        return False
    try:
        return Chem.MolFromSmiles(text, sanitize=False) is not None
    except Exception:
        return False


def _looks_like_header(row: Sequence[str]) -> bool:
    lowered = [cell.strip().lower() for cell in row]
    if any(cell in _SMILES_HEADERS or cell in _ID_HEADERS for cell in lowered):
        return True
    return not any(_parses_as_smiles(cell) for cell in row)


def _detect_columns(
    header: Optional[Sequence[str]], rows: Sequence[Sequence[str]]
) -> tuple[int, Optional[int]]:
    """Localiza a coluna de estrutura e, se houver, a de identificador.

    A coluna de estrutura é escolhida pela fração de valores que o RDKit consegue
    interpretar. Isso resolve simultaneamente ``SMILES<sep>nome`` e ``nome<sep>SMILES``
    sem que o chamador precise declarar o layout.
    """
    if not rows:
        return 0, None

    width = max(len(row) for row in rows)

    if header is not None:
        lowered = [cell.strip().lower() for cell in header]
        named_smiles = next(
            (i for i, cell in enumerate(lowered) if cell in _SMILES_HEADERS), None
        )
        if named_smiles is not None:
            named_id = next(
                (
                    i
                    for i, cell in enumerate(lowered)
                    if cell in _ID_HEADERS and i != named_smiles
                ),
                None,
            )
            return named_smiles, named_id

    scores: list[tuple[float, int]] = []
    for index in range(width):
        values = [row[index] for row in rows if index < len(row)]
        if not values:
            continue
        rate = sum(_parses_as_smiles(value) for value in values) / len(values)
        scores.append((rate, index))

    if not scores:
        return 0, None

    best_rate, smiles_index = max(scores, key=lambda item: (item[0], -item[1]))
    if best_rate < _SMILES_COLUMN_THRESHOLD:
        smiles_index = 0

    id_index: Optional[int] = None
    if header is not None:
        lowered = [cell.strip().lower() for cell in header]
        id_index = next(
            (
                i
                for i, cell in enumerate(lowered)
                if cell in _ID_HEADERS and i != smiles_index
            ),
            None,
        )
    if id_index is None and width > 1:
        id_index = next((i for i in range(width) if i != smiles_index), None)

    return smiles_index, id_index


def read_input(source: Source) -> Generator[tuple[str, str], None, None]:
    """Lê estruturas de arquivo, handle ou texto colado, em streaming.

    Aceita ``.csv``, ``.tsv``, ``.smi`` e texto bruto. Separador, cabeçalho e papel
    de cada coluna são inferidos de uma amostra inicial; o restante do arquivo é
    consumido de forma preguiçosa, sem carregar o conteúdo em memória.

    Yields:
        Pares ``(input_id, raw_smiles)``. O identificador vem da coluna detectada
        quando existe e é preenchido; caso contrário é gerado deterministicamente
        como ``CMPD_0000001`` a partir da posição do registro, de modo que duas
        leituras da mesma fonte produzem os mesmos identificadores.
    """
    lines = _clean(_iter_raw_lines(source))

    sample: list[str] = []
    for line in lines:
        sample.append(line)
        if len(sample) >= _SAMPLE_ROWS:
            break

    if not sample:
        return

    delimiter = _detect_delimiter(sample)
    sample_rows = [_split(line, delimiter) for line in sample]

    header: Optional[Sequence[str]] = None
    data_rows = sample_rows
    if len(sample_rows) > 1 and _looks_like_header(sample_rows[0]):
        header = sample_rows[0]
        data_rows = sample_rows[1:]

    smiles_index, id_index = _detect_columns(header, data_rows)

    counter = 0
    for row in _chain_rows(data_rows, lines, delimiter):
        if smiles_index >= len(row):
            continue
        smiles = row[smiles_index].strip()
        if not smiles:
            continue

        counter += 1
        identifier = ""
        if id_index is not None and id_index < len(row):
            identifier = row[id_index].strip()
        yield (identifier or f"CMPD_{counter:07d}", smiles)


def _chain_rows(
    buffered: Sequence[Sequence[str]],
    remaining: Iterator[str],
    delimiter: Optional[str],
) -> Iterator[Sequence[str]]:
    yield from buffered
    for line in remaining:
        yield _split(line, delimiter)


# --- Escrita -----------------------------------------------------------------------

_CURATED_COLUMNS = (
    "input_id",
    "raw_smiles",
    "curated_smiles",
    "inchikey",
    "inchikey_block1",
    "removed_fragments",
    "delta_net_charge",
    "delta_formula",
    "delta_stereocenters",
    "n_stereocenters_total",
    "n_undefined_stereocenters",
    "n_components_parent",
    "excluded_flag",
)

_REJECTED_COLUMNS = (
    "input_id",
    "raw_smiles",
    "rejection_code",
    "rejection_stage",
    "rejection_detail",
)

_AUDIT_COLUMNS = (
    "input_id",
    "raw_smiles",
    "status",
    "curated_smiles",
    "inchikey",
    "inchikey_block1",
    "rejection_code",
    "rejection_stage",
    "rejection_detail",
    "removed_fragments",
    "delta_net_charge",
    "delta_formula",
    "delta_stereocenters",
    "n_stereocenters_total",
    "n_undefined_stereocenters",
    "n_components_parent",
    "excluded_flag",
    "transformations",
    "rdkit_version",
    "csp_version",
    "pipeline_version",
    "policy_hash",
)

MANIFEST_NAME = "manifest.json"

#: Fluxos de saída, todos sob a mesma garantia de atomicidade.
_STREAMS: tuple[tuple[str, Sequence[str]], ...] = (
    ("curated", _CURATED_COLUMNS),
    ("rejected", _REJECTED_COLUMNS),
    ("audit", _AUDIT_COLUMNS),
    ("conflicts", CONFLICT_COLUMNS),
)


class BatchWriter:
    """Escreve as saídas de um lote com garantia de atomicidade.

    Uso::

        with BatchWriter(out_dir, policy_hash) as writer:
            for record in records:
                writer.write(record)

    Ao sair sem exceção, os ``.tmp`` são promovidos por ``os.replace`` — atômico no
    mesmo sistema de arquivos — e ``manifest.json`` é gravado por último. Ao sair com
    exceção, os ``.tmp`` são removidos e nenhum arquivo final é criado.

    **Limite honesto:** um crash do sistema operacional impede ``__exit__`` de rodar,
    então nenhum código pode limpar nada naquele instante. Por isso a limpeza de
    resíduos acontece em ``__enter__``, que remove ``.tmp`` de execuções anteriores
    interrompidas, e por isso ``manifest.json`` — escrito por último — é o único
    marcador confiável de completude. Saída sem manifesto é saída inválida.
    """

    def __init__(
        self,
        out_dir: Union[str, Path],
        policy_hash: str,
        pipeline_version: str = PIPELINE_VERSION,
    ) -> None:
        self.out_dir = Path(out_dir)
        self.policy_hash = policy_hash
        self.pipeline_version = pipeline_version

        self.n_total = 0
        self.n_passed = 0
        self.n_rejected = 0
        self.rejection_counts: dict[str, int] = {}
        self.collision_counts: dict[str, int] = {}
        self.n_conflicts = 0
        self.n_unique = 0
        self.policy_hash_mismatches = 0

        self._started_at: Optional[datetime] = None
        self._handles: dict[str, TextIO] = {}
        self._writers: dict[str, csv.DictWriter] = {}
        self._closed = False

    # --- Protocolo de contexto ---------------------------------------------------

    def __enter__(self) -> "BatchWriter":
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self._remove_stale_temporaries()
        self._started_at = datetime.now(timezone.utc)

        for name, columns in _STREAMS:
            handle = self._temp_path(name).open("w", encoding="utf-8", newline="")
            writer = csv.DictWriter(handle, fieldnames=list(columns))
            writer.writeheader()
            self._handles[name] = handle
            self._writers[name] = writer

        return self

    def __exit__(self, exc_type, exc_value, traceback) -> bool:
        if exc_type is not None:
            self._abort()
            return False
        try:
            self._commit()
        except Exception:
            self._abort()
            raise
        return False

    # --- Escrita ------------------------------------------------------------------

    def write(self, record: CurationRecord) -> None:
        """Registra um resultado nas saídas correspondentes.

        Todo registro vai para ``audit``; a segregação em ``curated``/``rejected`` é
        uma projeção. Um registro cujo ``policy_hash`` divirja do lote é contabilizado
        como divergência e reportado no manifesto: misturar saídas produzidas sob
        políticas diferentes invalidaria silenciosamente a comparação entre lotes.
        """
        if self._closed:
            raise RuntimeError("BatchWriter ja foi encerrado")

        payload = record.model_dump(mode="json")
        payload["transformations"] = json.dumps(
            payload.get("transformations", []), ensure_ascii=False
        )

        self.n_total += 1
        if record.policy_hash != self.policy_hash:
            self.policy_hash_mismatches += 1

        self._writers["audit"].writerow(self._project(payload, _AUDIT_COLUMNS))

        if record.status == "PASSED":
            self.n_passed += 1
            self._writers["curated"].writerow(self._project(payload, _CURATED_COLUMNS))
        else:
            self.n_rejected += 1
            code = record.rejection_code.value if record.rejection_code else "UNKNOWN"
            self.rejection_counts[code] = self.rejection_counts.get(code, 0) + 1
            self._writers["rejected"].writerow(
                self._project(payload, _REJECTED_COLUMNS)
            )

    def write_conflict(self, collision: Collision) -> None:
        """Registra uma colisão de identidade em ``conflicts.csv``.

        O relatório é sempre criado, mesmo vazio: ausência de conflitos é um
        resultado, e precisa ser distinguível de "o relatório não foi gerado".
        """
        if self._closed:
            raise RuntimeError("BatchWriter ja foi encerrado")

        self.n_conflicts += 1
        key = collision.collision_type.value
        self.collision_counts[key] = self.collision_counts.get(key, 0) + 1
        self._writers["conflicts"].writerow(collision.as_row())

    def write_all(self, records: Iterable[CurationRecord]) -> None:
        for record in records:
            self.write(record)

    # --- Internos ------------------------------------------------------------------

    @staticmethod
    def _project(payload: dict, columns: Sequence[str]) -> dict:
        return {column: payload.get(column, "") for column in columns}

    def _temp_path(self, name: str) -> Path:
        return self.out_dir / f"{name}.csv.tmp"

    def _final_path(self, name: str) -> Path:
        return self.out_dir / f"{name}.csv"

    def _remove_stale_temporaries(self) -> None:
        """Descarta resíduos de uma execução anterior interrompida.

        Este é o ponto em que a limpeza pós-crash realmente acontece: nada roda no
        instante em que o processo morre.
        """
        for leftover in self.out_dir.glob("*.csv.tmp"):
            leftover.unlink(missing_ok=True)
        (self.out_dir / f"{MANIFEST_NAME}.tmp").unlink(missing_ok=True)

    def _close_handles(self, *, durable: bool) -> None:
        for handle in self._handles.values():
            try:
                handle.flush()
                if durable:
                    os.fsync(handle.fileno())
            except Exception:
                pass
            finally:
                handle.close()
        self._handles.clear()
        self._writers.clear()

    def _abort(self) -> None:
        """Encerra descartando tudo: sem arquivos finais, sem manifesto."""
        self._close_handles(durable=False)
        for name, _ in _STREAMS:
            self._temp_path(name).unlink(missing_ok=True)
        self._closed = True

    def _commit(self) -> None:
        """Promove os temporários e grava o manifesto por último."""
        self._close_handles(durable=True)
        for name, _ in _STREAMS:
            os.replace(self._temp_path(name), self._final_path(name))

        self._write_manifest()
        self._sync_directory()
        self._closed = True

    def _write_manifest(self) -> None:
        finished_at = datetime.now(timezone.utc)
        manifest = {
            "pipeline_version": self.pipeline_version,
            "policy_hash": self.policy_hash,
            "started_at": self._started_at.isoformat() if self._started_at else None,
            "finished_at": finished_at.isoformat(),
            "versions": {
                "python": sys.version.split()[0],
                "rdkit": rdkit.__version__,
                "chembl_structure_pipeline": getattr(csp, "__version__", "unknown"),
            },
            "counts": {
                "total": self.n_total,
                "passed": self.n_passed,
                "rejected": self.n_rejected,
            },
            "rejection_counts": dict(sorted(self.rejection_counts.items())),
            "deduplication": {
                "unique": self.n_unique,
                "conflicts": self.n_conflicts,
                "by_type": dict(sorted(self.collision_counts.items())),
            },
            "policy_hash_mismatches": self.policy_hash_mismatches,
            "outputs": {
                name: self._final_path(name).name for name, _ in _STREAMS
            },
        }

        temp = self.out_dir / f"{MANIFEST_NAME}.tmp"
        with temp.open("w", encoding="utf-8") as handle:
            json.dump(manifest, handle, indent=2, ensure_ascii=False)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp, self.out_dir / MANIFEST_NAME)

    def _sync_directory(self) -> None:
        """Persiste as entradas de diretório criadas pelos renames.

        Sem isso, os renames podem não sobreviver a uma queda de energia em alguns
        sistemas de arquivos. Não suportado em todas as plataformas; a falha é
        tolerada porque a alternativa seria abortar um lote já íntegro.
        """
        try:
            fd = os.open(self.out_dir, os.O_RDONLY)
        except OSError:
            return
        try:
            os.fsync(fd)
        except OSError:
            pass
        finally:
            os.close(fd)


def read_manifest(out_dir: Union[str, Path]) -> Optional[dict]:
    """Lê o manifesto de um lote, ou ``None`` se a execução não foi concluída.

    É a forma correta de verificar se um diretório de saída pode ser consumido.
    """
    path = Path(out_dir) / MANIFEST_NAME
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
