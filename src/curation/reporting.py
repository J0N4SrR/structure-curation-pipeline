"""Derivação de relatórios a partir dos registros reais de uma execução.

Toda contagem, status e duração aqui vem de :class:`~curation.models.CurationRecord`
e dos cronômetros do motor. Nada é estimado. A interface consome estas estruturas;
ela não as calcula.
"""

from __future__ import annotations

import csv
import io
import json
import zipfile
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Iterable, Optional, Sequence

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from curation.dedup import CollisionType, DedupIndex
from curation.models import CurationRecord, Stage
from curation.provenance import RunProvenance, sha256_of


class StageStatus(str, Enum):
    PENDING = "PENDING"
    RUNNING = "RUNNING"
    SUCCESS = "SUCCESS"
    WARNING = "WARNING"
    FAILED = "FAILED"
    SKIPPED = "SKIPPED"


#: Ordem real de execução. São os estágios que o motor de fato percorre - não uma
#: simplificação didática - para que o funil exibido corresponda ao código.
PIPELINE_STAGES: tuple[tuple[str, str], ...] = (
    (Stage.PARSE.value, "Leitura e parsing não-restritivo do SMILES"),
    (Stage.STANDARDIZE.value, "Normalização e neutralização (ChEMBL)"),
    (Stage.GET_PARENT.value, "Remoção de sais e isolamento da estrutura-mãe"),
    (Stage.VALENCE_GATE.value, "Sanitização estrita, posterior ao motor"),
    (Stage.ELIGIBILITY.value, "Cortes de escopo químico (MW, átomos pesados)"),
    (Stage.CANONICALIZE.value, "SMILES canônico isomérico e InChIKey"),
)

DEDUP_STAGE = "DEDUPLICATION"


@dataclass(frozen=True)
class StageReport:
    """Estado observado de um estágio ao longo do lote."""

    name: str
    description: str
    status: StageStatus
    n_input: int
    n_output: int
    duration_seconds: float = 0.0
    rejection_counts: dict[str, int] = field(default_factory=dict)

    @property
    def n_excluded(self) -> int:
        return self.n_input - self.n_output

    @property
    def has_exclusions(self) -> bool:
        return self.n_excluded > 0


@dataclass(frozen=True)
class ExclusionGroup:
    """Uma categoria de exclusão, com os registros responsáveis.

    Existe para cumprir a regra de que nenhuma estrutura desapareça silenciosamente:
    toda diferença entre entrada e saída de um estágio é atribuível a registros
    nomeados.
    """

    reason: str
    stage: str
    count: int
    record_ids: list[str]


@dataclass
class RunReport:
    """Visão completa de uma execução, derivada dos registros."""

    provenance: RunProvenance
    records: list[CurationRecord]
    stages: list[StageReport]
    dedup: Optional[StageReport] = None
    index: Optional[DedupIndex] = None

    # --- KPIs --------------------------------------------------------------------

    @property
    def total(self) -> int:
        return len(self.records)

    @property
    def approved(self) -> list[CurationRecord]:
        return [record for record in self.records if record.passed]

    @property
    def rejected(self) -> list[CurationRecord]:
        return [record for record in self.records if not record.passed]

    @property
    def n_duplicates(self) -> int:
        """Duplicatas reais - colisões de bloco1 não são duplicatas (D-07)."""
        if self.index is None:
            return 0
        counts = self.index.collision_counts()
        return counts.get(CollisionType.EXACT_DUPLICATE.value, 0) + counts.get(
            CollisionType.ANNOTATION_CONFLICT.value, 0
        )

    @property
    def n_warnings(self) -> int:
        """Registros aprovados que sofreram alguma transformação observada.

        Não é erro: é curadoria que alterou a estrutura e merece inspeção.
        """
        return sum(1 for record in self.approved if record.transformations)

    @property
    def n_unique(self) -> int:
        return self.index.unique_count if self.index else len(self.approved)

    def kpis(self) -> dict[str, int]:
        return {
            "Total Processed": self.total,
            "Approved": len(self.approved),
            "Rejected": len(self.rejected),
            "Duplicates": self.n_duplicates,
            "Warnings": self.n_warnings,
        }

    # --- Exclusões ------------------------------------------------------------------

    def exclusion_groups(self) -> list[ExclusionGroup]:
        """Motivos de exclusão com os registros responsáveis, do maior ao menor."""
        grouped: dict[tuple[str, str], list[str]] = {}
        for record in self.rejected:
            reason = record.rejection_code.value if record.rejection_code else "UNKNOWN"
            stage = record.rejection_stage.value if record.rejection_stage else "UNKNOWN"
            grouped.setdefault((reason, stage), []).append(record.input_id)

        if self.index is not None:
            for collision in self.index.collisions:
                if collision.collision_type is CollisionType.BLOCK1_COLLISION:
                    continue
                key = (collision.collision_type.value, DEDUP_STAGE)
                grouped.setdefault(key, []).append(collision.incoming_id)

        groups = [
            ExclusionGroup(reason=reason, stage=stage, count=len(ids), record_ids=ids)
            for (reason, stage), ids in grouped.items()
        ]
        return sorted(groups, key=lambda group: (-group.count, group.reason))

    def largest_reduction(self) -> Optional[StageReport]:
        """Etapa que mais removeu estruturas, ou ``None`` se nenhuma removeu.

        Existe aqui, e não na interface, porque é uma leitura dos dados do lote -
        a UI apresenta a frase, não a deriva.
        """
        candidates = [
            stage
            for stage in list(self.stages) + ([self.dedup] if self.dedup else [])
            if stage.n_excluded > 0
        ]
        return max(candidates, key=lambda stage: stage.n_excluded, default=None)

    def records_by_id(self, identifiers: Iterable[str]) -> list[CurationRecord]:
        wanted = set(identifiers)
        return [record for record in self.records if record.input_id in wanted]


# --- Construção ---------------------------------------------------------------------


def _stage_status(n_input: int, n_output: int) -> StageStatus:
    if n_input == 0:
        return StageStatus.SKIPPED
    if n_output == 0:
        return StageStatus.FAILED
    if n_output < n_input:
        return StageStatus.WARNING
    return StageStatus.SUCCESS


def build_stage_reports(
    records: Sequence[CurationRecord], stage_seconds: Optional[dict[str, float]] = None
) -> list[StageReport]:
    """Monta o funil a partir do estágio em que cada registro foi rejeitado."""
    timings = stage_seconds or {}
    rejections: dict[str, dict[str, int]] = {}
    for record in records:
        if record.passed or record.rejection_stage is None:
            continue
        stage = record.rejection_stage.value
        code = record.rejection_code.value if record.rejection_code else "UNKNOWN"
        rejections.setdefault(stage, {})
        rejections[stage][code] = rejections[stage].get(code, 0) + 1

    reports: list[StageReport] = []
    carried = len(records)
    for name, description in PIPELINE_STAGES:
        counts = rejections.get(name, {})
        n_output = carried - sum(counts.values())
        reports.append(
            StageReport(
                name=name,
                description=description,
                status=_stage_status(carried, n_output),
                n_input=carried,
                n_output=n_output,
                duration_seconds=timings.get(name, 0.0),
                rejection_counts=dict(counts),
            )
        )
        carried = n_output
    return reports


def build_dedup_report(
    approved: int, index: Optional[DedupIndex]
) -> Optional[StageReport]:
    if index is None:
        return None
    counts = index.collision_counts()
    removed = counts.get(CollisionType.EXACT_DUPLICATE.value, 0) + counts.get(
        CollisionType.ANNOTATION_CONFLICT.value, 0
    )
    return StageReport(
        name=DEDUP_STAGE,
        description="Identidade por InChIKey completo; bloco1 apenas agrupa",
        status=_stage_status(approved, approved - removed),
        n_input=approved,
        n_output=approved - removed,
        rejection_counts=dict(counts),
    )


def build_run_report(
    records: Sequence[CurationRecord],
    provenance: RunProvenance,
    stage_seconds: Optional[dict[str, float]] = None,
    index: Optional[DedupIndex] = None,
) -> RunReport:
    approved = sum(1 for record in records if record.passed)
    return RunReport(
        provenance=provenance,
        records=list(records),
        stages=build_stage_reports(records, stage_seconds),
        dedup=build_dedup_report(approved, index),
        index=index,
    )


# --- Rastreabilidade por estrutura -------------------------------------------------


@dataclass(frozen=True)
class LineageStep:
    """Um passo observado na trajetória de uma estrutura."""

    stage: str
    rule: str
    before_smiles: str
    after_smiles: str
    detail: str


def structure_lineage(record: CurationRecord) -> list[LineageStep]:
    """Trajetória de uma estrutura, estritamente a partir do que foi observado.

    Os passos vêm dos eventos de transformação registrados pelas sondas. Estágios
    que não produziram efeito observável **não** aparecem como passos: o motor é uma
    caixa-preta de oito operações internas e inventar estados intermediários daria
    uma falsa impressão de granularidade.
    """
    return [
        LineageStep(
            stage=event.stage.value,
            rule=event.rule,
            before_smiles=event.before_smiles,
            after_smiles=event.after_smiles,
            detail=event.detail,
        )
        for event in record.transformations
    ]


# --- Exportação -----------------------------------------------------------------------

#: Colunas nativas do registro.
BASE_COLUMNS: tuple[str, ...] = (
    "input_id",
    "raw_smiles",
    "curated_smiles",
    "inchikey",
    "inchikey_block1",
    "status",
    "rejection_code",
    "rejection_stage",
    "rejection_detail",
    "removed_fragments",
    "salt_removed",
    "delta_net_charge",
    "delta_formula",
    "delta_stereocenters",
    "n_stereocenters_total",
    "n_undefined_stereocenters",
    "n_components_parent",
    "parent_mw",
    "parent_heavy_atoms",
    "excluded_flag",
    "transformations",
    "rdkit_version",
    "csp_version",
    "pipeline_version",
    "policy_hash",
)

#: Colunas calculadas sob demanda, porque custam InChI ou fórmula molecular.
DERIVED_COLUMNS: tuple[str, ...] = ("inchi", "molecular_formula")

EXPORTABLE_COLUMNS: tuple[str, ...] = BASE_COLUMNS + DERIVED_COLUMNS

REJECTED_COLUMNS: tuple[str, ...] = (
    "input_id",
    "raw_smiles",
    "status",
    "rejection_stage",
    "rejection_code",
    "rejection_detail",
)


def _derived(record: CurationRecord, column: str) -> str:
    if not record.curated_smiles:
        return ""
    mol = Chem.MolFromSmiles(record.curated_smiles, sanitize=False)
    if mol is None:
        return ""
    try:
        mol.UpdatePropertyCache(strict=False)
        if column == "inchi":
            return Chem.MolToInchi(mol) or ""
        if column == "molecular_formula":
            return rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        return ""
    return ""


def record_row(record: CurationRecord, columns: Sequence[str]) -> dict[str, Any]:
    """Uma linha de saída, com ``None`` preservado nos campos ausentes.

    ``None`` não é convertido em string vazia de propósito. ``csv.DictWriter`` já
    escreve célula vazia para ``None``, e preservá-lo mantém as colunas numéricas
    homogêneas - coagir para ``""`` produz uma coluna mista de ``float`` e ``str``
    que o Arrow recusa converter ao renderizar a tabela.
    """
    payload = record.model_dump(mode="json")
    payload["transformations"] = json.dumps(
        payload.get("transformations", []), ensure_ascii=False
    )
    payload["salt_removed"] = bool(record.removed_fragments)

    row: dict[str, Any] = {}
    for column in columns:
        if column in DERIVED_COLUMNS:
            row[column] = _derived(record, column)
        else:
            row[column] = payload.get(column)
    return row


def to_csv(
    records: Sequence[CurationRecord], columns: Sequence[str] = BASE_COLUMNS
) -> str:
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(columns))
    writer.writeheader()
    for record in records:
        writer.writerow(record_row(record, columns))
    return buffer.getvalue()


def full_csv(report: RunReport) -> str:
    """Todos os registros, aprovados e rejeitados, com status e motivo."""
    return to_csv(report.records, BASE_COLUMNS)


def rejected_csv(report: RunReport) -> str:
    return to_csv(report.rejected, REJECTED_COLUMNS)


def run_manifest(report: RunReport) -> dict[str, Any]:
    """Manifesto auditável da execução."""
    return {
        "run_id": report.provenance.run_id,
        "status": "SUCCESS" if report.total else "EMPTY",
        "provenance": report.provenance.to_dict(),
        "metrics": report.kpis(),
        "stages": [
            {
                "name": stage.name,
                "status": stage.status.value,
                "input": stage.n_input,
                "output": stage.n_output,
                "excluded": stage.n_excluded,
                "duration_seconds": round(stage.duration_seconds, 4),
                "rejections": stage.rejection_counts,
            }
            for stage in report.stages
            + ([report.dedup] if report.dedup is not None else [])
        ],
        "exclusions": [
            {
                "reason": group.reason,
                "stage": group.stage,
                "count": group.count,
            }
            for group in report.exclusion_groups()
        ],
    }


def parameters_yaml(report: RunReport) -> str:
    """Parâmetros da execução em YAML, sem dependência externa."""
    lines = ["# Parametros da execucao", f"run_id: {report.provenance.run_id}"]
    for key, value in sorted(report.provenance.parameters.items()):
        lines.append(f"{key}: {value}")
    lines += [
        f"policy_hash: {report.provenance.policy_hash}",
        f"git_commit: {report.provenance.git.commit}",
        f"git_dirty: {report.provenance.git.dirty}",
        "versions:",
    ]
    for key, value in sorted(report.provenance.versions.items()):
        lines.append(f"  {key}: {value}")
    return "\n".join(lines) + "\n"


def package_readme(report: RunReport) -> str:
    provenance = report.provenance
    blockers = provenance.reproduction_blockers()
    verdict = (
        "Esta execucao e reproduzivel a partir dos metadados abaixo."
        if not blockers
        else "**Esta execucao NAO e exatamente reproduzivel.** Motivos:\n"
        + "\n".join(f"- {item}" for item in blockers)
    )
    kpis = "\n".join(f"| {name} | {value} |" for name, value in report.kpis().items())
    parameters = "\n".join(
        f"| `{key}` | {value} |" for key, value in sorted(provenance.parameters.items())
    )

    return f"""# Execucao {provenance.run_id}

Curadoria estrutural com RDKit e `chembl_structure_pipeline`.

- Inicio: {provenance.started_at}
- Termino: {provenance.finished_at}
- Duracao: {provenance.duration_seconds:.2f} s
- Entrada: `{provenance.input_name}` (sha256 `{provenance.input_hash[:16]}...`)

## Parametros

| parametro | valor |
| --- | --- |
{parameters}

## Resultados

| metrica | valor |
| --- | ---: |
{kpis}

## Arquivos

| arquivo | conteudo |
| --- | --- |
| `structures_full.csv` | todos os registros, aprovados e rejeitados, com status e motivo |
| `structures_rejected.csv` | apenas os rejeitados, com estagio e codigo de erro |
| `run_manifest.json` | proveniencia, metricas por estagio e exclusoes |
| `parameters.yaml` | parametros e versoes usados |

## Como interpretar

Cada registro carrega o estagio em que foi rejeitado, quando aplicavel, e a lista de
transformacoes observadas. Nenhuma estrutura e descartada sem registro: a soma das
exclusoes por motivo no manifesto reconcilia com a diferenca entre entrada e saida.

Colisoes de bloco1 do InChIKey **nao** sao duplicatas - enantiomeros compartilham o
primeiro bloco. A identidade e o InChIKey completo.

## Reproducao

{verdict}

```bash
{provenance.reproduction_command(provenance.input_name)}
```

Versoes: RDKit {provenance.versions.get('rdkit')}, chembl_structure_pipeline
{provenance.versions.get('chembl_structure_pipeline')}. Mudancas de versao do RDKit
alteram percepcao de aromaticidade e regras de padronizacao - dois lotes com RDKit
diferente nao sao comparaveis mesmo com o mesmo policy_hash.
"""


def reproducibility_package(report: RunReport) -> bytes:
    """Arquivo ZIP consolidado com dados, manifesto e instruções."""
    buffer = io.BytesIO()
    root = f"structure-curation-run-{report.provenance.run_id}"
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(f"{root}/structures_full.csv", full_csv(report))
        archive.writestr(f"{root}/structures_rejected.csv", rejected_csv(report))
        archive.writestr(
            f"{root}/run_manifest.json",
            json.dumps(run_manifest(report), indent=2, ensure_ascii=False),
        )
        archive.writestr(f"{root}/parameters.yaml", parameters_yaml(report))
        archive.writestr(f"{root}/README.md", package_readme(report))
    return buffer.getvalue()


def finalize(report: RunReport, duration_seconds: float) -> RunReport:
    """Fecha a proveniência com duração e hash da saída."""
    report.provenance = report.provenance.finish(
        duration_seconds=duration_seconds, output_hash=sha256_of(full_csv(report))
    )
    return report
