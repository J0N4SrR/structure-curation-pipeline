"""Orquestração de execução: registro único e lote.

A ordem dos estágios químicos — parse, instrumentação, ``standardize_mol``,
``get_parent_mol``, portão de valência, canonicalização — é implementada em
:class:`~curation.engine.EngineWrapper`, que é quem detém o ciclo de vida das
moléculas. Este módulo não reimplementa nenhum estágio: duplicar a sequência em dois
lugares criaria exatamente o tipo de divergência silenciosa que a suíte de regressão
existe para impedir.

O que este módulo acrescenta é orquestração: normalização da entrada, barreira global
de exceção, composição de ingestão e escrita, e o sumário do lote.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Iterator, Optional, Sequence, Union

import time

from curation.dedup import DedupIndex
from curation.engine import PIPELINE_VERSION, EngineWrapper
from curation.filters import EligibilityCriteria
from curation.io import BatchWriter, Source, read_input
from curation.models import CurationRecord, RejectionCode, Stage


@dataclass(frozen=True)
class BatchSummary:
    """Resultado agregado de um lote."""

    total: int
    passed: int
    rejected: int
    rejection_counts: dict[str, int] = field(default_factory=dict)
    unique: int = 0
    conflicts: int = 0
    collision_counts: dict[str, int] = field(default_factory=dict)
    out_dir: Optional[Path] = None

    @property
    def pass_rate(self) -> float:
        return self.passed / self.total if self.total else 0.0

    def top_rejections(self, limit: int = 5) -> list[tuple[str, int]]:
        """Motivos de rejeição mais frequentes, do maior para o menor."""
        return sorted(
            self.rejection_counts.items(), key=lambda item: (-item[1], item[0])
        )[:limit]

    def format_report(self) -> str:
        """Sumário legível para o terminal."""
        width = 22
        lines = [
            "-" * 52,
            "SUMARIO DA EXECUCAO",
            "-" * 52,
            f"{'entradas processadas':<{width}} {self.total:>8,}",
            f"{'aprovados':<{width}} {self.passed:>8,}   {self.pass_rate:6.1%}",
            f"{'rejeitados':<{width}} {self.rejected:>8,}   "
            f"{1 - self.pass_rate:6.1%}",
        ]

        top = self.top_rejections()
        if top:
            lines.append("")
            lines.append("top motivos de rejeicao")
            for code, count in top:
                share = count / self.total if self.total else 0.0
                lines.append(f"    {code:<20s} {count:>8,}   {share:6.1%}")

        if self.unique or self.conflicts:
            lines.append("")
            lines.append("deduplicacao")
            lines.append(f"    {'identidades unicas':<20s} {self.unique:>8,}")
            lines.append(f"    {'colisoes':<20s} {self.conflicts:>8,}")
            for kind, count in sorted(
                self.collision_counts.items(), key=lambda item: (-item[1], item[0])
            ):
                lines.append(f"        {kind:<18s} {count:>8,}")

        if self.out_dir is not None:
            lines.append("")
            lines.append(f"saida: {self.out_dir}")
        lines.append("-" * 52)
        return "\n".join(lines)


class CurationPipeline:
    """Executa a curadoria sobre registros individuais ou lotes inteiros.

    Args:
        policy_hash: SHA-256 de ``docs/decisions.md``, carimbado em cada registro.
        engine: motor a utilizar. O padrão constrói um
            :class:`~curation.engine.EngineWrapper` com as sondas vigentes; injetar
            outro serve a testes e a estudos de ablação (Fase 5).
    """

    def __init__(
        self,
        policy_hash: str,
        engine: Optional[EngineWrapper] = None,
        pipeline_version: str = PIPELINE_VERSION,
        criteria: Optional[EligibilityCriteria] = None,
        deduplicate: bool = True,
    ) -> None:
        self.policy_hash = policy_hash
        self.pipeline_version = pipeline_version
        self.criteria = criteria or EligibilityCriteria()
        self.deduplicate = deduplicate
        self._engine = engine or EngineWrapper(
            policy_hash=policy_hash,
            pipeline_version=pipeline_version,
            criteria=self.criteria,
        )

    # --- Registro único --------------------------------------------------------

    def process_single(self, raw_smiles: str, compound_id: str) -> CurationRecord:
        """Cura um composto, devolvendo sempre um registro.

        Ordem executada pelo motor::

            parse (sanitize=False + UpdatePropertyCache(strict=False))
              -> instrumentacao (sondas, sem rejeitar)
              -> standardize_mol
              -> get_parent_mol
              -> portao de valencia (SanitizeMol estrito)
              -> canonicalizacao (SMILES isomerico + InChIKey)

        O portão de valência é posterior ao motor por decisão de política (D-09):
        amônio quaternário neutro e diazônio neutro são reparados pelo normalizador
        da referência e seriam descartados por uma sanitização antecipada.

        Esta é a **barreira global**: nenhuma exceção escapa. Uma molécula atípica
        que quebrasse o motor de forma imprevista viraria um registro
        ``ERR_INTERNAL``, jamais interromperia o lote.
        """
        try:
            return self._engine.curate(compound_id, raw_smiles)
        except Exception as error:  # noqa: BLE001 - barreira deliberada
            return self._internal_failure(compound_id, raw_smiles, error)

    def process_many(
        self, records: Iterable[tuple[str, str]]
    ) -> Iterator[CurationRecord]:
        """Cura uma sequência de pares ``(input_id, raw_smiles)``, preguiçosamente."""
        for input_id, raw_smiles in records:
            yield self.process_single(raw_smiles, input_id)

    # --- Lote ------------------------------------------------------------------

    def run(self, source: Source, out_dir: Union[str, Path]) -> BatchSummary:
        """Executa um lote completo: ingestão, curadoria e escrita atômica.

        A ingestão é um gerador e os registros são escritos à medida que saem, de
        modo que a memória não cresce com o tamanho da entrada. Se algo falhar no
        meio, o :class:`~curation.io.BatchWriter` descarta os temporários e não grava
        manifesto — a saída fica reconhecidamente incompleta em vez de parecer boa.
        """
        out_path = Path(out_dir)
        index = DedupIndex() if self.deduplicate else None

        with BatchWriter(
            out_path,
            self.policy_hash,
            pipeline_version=self.pipeline_version,
            parameters=self.default_parameters(),
        ) as writer:
            for record in self.process_many(read_input(source)):
                writer.write(record)
                if index is not None:
                    collision = index.add(record)
                    if collision is not None:
                        writer.write_conflict(collision)

            if index is not None:
                writer.n_unique = index.unique_count

            summary = BatchSummary(
                total=writer.n_total,
                passed=writer.n_passed,
                rejected=writer.n_rejected,
                rejection_counts=dict(writer.rejection_counts),
                unique=writer.n_unique,
                conflicts=writer.n_conflicts,
                collision_counts=dict(writer.collision_counts),
                out_dir=out_path,
            )
        return summary

    def run_report(
        self,
        source: Source,
        parameters: Optional[dict] = None,
        input_bytes: Optional[bytes] = None,
        input_name: str = "input",
        policy_path: str = "docs/decisions.md",
        progress: Optional[Callable[[int, int, str], None]] = None,
    ):
        """Executa um lote em memória e devolve um :class:`RunReport` completo.

        Diferente de :meth:`run`, não escreve em disco: serve a clientes interativos
        que precisam do relatório, da proveniência e dos artefatos de exportação sem
        materializar um diretório de saída.
        """
        from curation.provenance import RunProvenance
        from curation.reporting import build_run_report, finalize

        records_in = list(read_input(source))
        raw = input_bytes if input_bytes is not None else b""
        provenance = RunProvenance.start(
            policy_hash=self.policy_hash,
            parameters=parameters or self.default_parameters(),
            input_bytes=raw,
            input_name=input_name,
            policy_path=policy_path,
        )

        self._engine.reset_timings()
        index = DedupIndex() if self.deduplicate else None
        started = time.perf_counter()
        records = []
        total = len(records_in)
        for position, (input_id, raw_smiles) in enumerate(records_in, start=1):
            record = self.process_single(raw_smiles, input_id)
            records.append(record)
            if index is not None:
                index.add(record)
            if progress is not None:
                progress(position, total, input_id)
        elapsed = time.perf_counter() - started

        report = build_run_report(
            records,
            provenance,
            stage_seconds=dict(self._engine.stage_seconds),
            index=index,
        )
        return finalize(report, elapsed)

    def default_parameters(self) -> dict:
        """Parâmetros efetivos desta instância, para registro no manifesto."""
        return {
            "max_mw": self.criteria.max_molecular_weight,
            "max_ha": self.criteria.max_heavy_atoms,
            "require_carbon": self.criteria.require_carbon,
            "deduplicate": self.deduplicate,
            "pipeline_version": self.pipeline_version,
        }

    @staticmethod
    def summarize(records: Sequence[CurationRecord]) -> BatchSummary:
        """Agrega registros já produzidos, sem escrever nada."""
        counts: dict[str, int] = {}
        passed = 0
        for record in records:
            if record.passed:
                passed += 1
            else:
                code = record.rejection_code.value if record.rejection_code else "UNKNOWN"
                counts[code] = counts.get(code, 0) + 1
        return BatchSummary(
            total=len(records),
            passed=passed,
            rejected=len(records) - passed,
            rejection_counts=counts,
        )

    # --- Internos ---------------------------------------------------------------

    def _internal_failure(
        self, compound_id: str, raw_smiles: str, error: Exception
    ) -> CurationRecord:
        return CurationRecord(
            input_id=compound_id,
            raw_smiles=raw_smiles,
            status="REJECTED",
            rejection_code=RejectionCode.ERR_INTERNAL,
            rejection_stage=Stage.PARSE,
            rejection_detail=f"{type(error).__name__}: {error}",
            rdkit_version=self._engine._rdkit_version,
            csp_version=self._engine._csp_version,
            pipeline_version=self.pipeline_version,
            policy_hash=self.policy_hash,
        )
