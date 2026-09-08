"""Camada de adaptacao entre os relatorios da engenharia e a interface.

Nao e uma segunda fonte de verdade. Todo numero aqui vem de
:class:`~curation.reporting.RunReport`; o que esta camada acrescenta e vocabulario
de apresentacao: rotulos legiveis, texto de ajuda, arestas explicitas do grafo e o
agrupamento por no.

O motivo de existir e separar duas responsabilidades que estavam misturadas: a
interface precisava saber o nome interno dos estagios, a ordem deles e como somar
exclusoes. Isso e conhecimento sobre o pipeline, nao sobre widgets, e agora vive
fora do Streamlit, onde pode ser testado sem um navegador.

O status dos nos reutiliza :class:`~curation.reporting.StageStatus`. Nao existe
enumeracao paralela.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

from curation.models import CurationRecord
from curation.reporting import (
    DEDUP_STAGE,
    ExclusionGroup,
    RunReport,
    StageReport,
    StageStatus,
)

#: Rotulo legivel de cada estagio real do motor. E traducao de nome, nao
#: reagrupamento: cada rotulo corresponde a um estagio que o pipeline executa, para
#: que o grafo exibido case com o codigo.
NODE_LABELS: dict[str, str] = {
    "PARSE": "Leitura",
    "STANDARDIZE": "Padronizacao",
    "GET_PARENT": "Remocao de sais",
    "VALENCE_GATE": "Validacao",
    "ELIGIBILITY": "Elegibilidade",
    "CANONICALIZE": "Canonicalizacao",
    DEDUP_STAGE: "Deduplicacao",
}

NODE_HELP: dict[str, str] = {
    "PARSE": "Le o texto da estrutura e monta o grafo quimico. So recusa o que o "
    "programa nao consegue interpretar.",
    "STANDARDIZE": "Aplica as regras de padronizacao do ChEMBL: normaliza grupos "
    "funcionais, separa metais ligados de forma inadequada e ajusta cargas.",
    "GET_PARENT": "Remove contra-ions e solventes, isolando a estrutura principal. "
    "Misturas legitimas de dois principios ativos sao preservadas.",
    "VALENCE_GATE": "Verifica se a estrutura resultante e quimicamente valida. Vem "
    "depois da padronizacao de proposito, para dar ao motor a chance de corrigir.",
    "ELIGIBILITY": "Aplica os limites de escopo do estudo sobre a estrutura "
    "principal, nunca sobre o sal.",
    "CANONICALIZE": "Gera a forma canonica da estrutura e o identificador InChIKey.",
    DEDUP_STAGE: "Identifica estruturas repetidas pelo InChIKey completo. "
    "Estruturas aparentadas, que so compartilham o inicio da chave, nao sao "
    "duplicatas.",
}

#: Parametros de execucao relevantes a cada no. Exibir todos em todo lugar
#: esconderia qual deles governa o estagio que o usuario esta olhando.
NODE_PARAMETERS: dict[str, tuple[str, ...]] = {
    "PARSE": (),
    "STANDARDIZE": (),
    "GET_PARENT": (),
    "VALENCE_GATE": (),
    "ELIGIBILITY": ("max_mw", "max_ha", "require_carbon"),
    "CANONICALIZE": (),
    DEDUP_STAGE: ("deduplicate",),
}

STATUS_GLYPHS: dict[StageStatus, str] = {
    StageStatus.PENDING: "o",
    StageStatus.RUNNING: "*",
    StageStatus.SUCCESS: "v",
    StageStatus.WARNING: "!",
    StageStatus.FAILED: "x",
    StageStatus.SKIPPED: "-",
}

STATUS_TEXT: dict[StageStatus, str] = {
    StageStatus.PENDING: "Aguardando",
    StageStatus.RUNNING: "Executando",
    StageStatus.SUCCESS: "Concluido",
    StageStatus.WARNING: "Concluido com remocoes",
    StageStatus.FAILED: "Falhou",
    StageStatus.SKIPPED: "Ignorado",
}


def label_for(node_id: str) -> str:
    return NODE_LABELS.get(node_id, node_id.replace("_", " ").title())


@dataclass(frozen=True)
class NodeArtifact:
    """Conjunto de registros produzido por um no, exportavel pela interface.

    Carrega os identificadores, nao o CSV: a serializacao continua sendo de
    :mod:`curation.reporting`, para que nao existam dois geradores de arquivo.
    """

    name: str
    description: str
    record_ids: tuple[str, ...]

    @property
    def count(self) -> int:
        return len(self.record_ids)


@dataclass(frozen=True)
class NodeContract:
    """Um no do grafo, com tudo que o painel de detalhe precisa exibir."""

    id: str
    label: str
    description: str
    status: StageStatus
    input_count: int
    output_count: int
    transformed_count: int
    rejected_count: int
    duration_seconds: float
    parameters: dict[str, object] = field(default_factory=dict)
    exclusions: tuple[ExclusionGroup, ...] = ()
    artifacts: tuple[NodeArtifact, ...] = ()
    help_text: str = ""

    @property
    def glyph(self) -> str:
        return STATUS_GLYPHS[self.status]

    @property
    def status_text(self) -> str:
        """Status em texto. O estado nunca e comunicado so por cor."""
        return STATUS_TEXT[self.status]

    @property
    def has_exclusions(self) -> bool:
        return self.rejected_count > 0


@dataclass(frozen=True)
class GraphEdge:
    """Dependencia explicita entre dois nos.

    A interface nao infere a ordem pela posicao: as arestas sao dados.
    """

    source: str
    target: str


@dataclass(frozen=True)
class GraphContract:
    """O grafo completo do pipeline."""

    nodes: tuple[NodeContract, ...]
    edges: tuple[GraphEdge, ...]

    def node(self, node_id: str) -> Optional[NodeContract]:
        return next((node for node in self.nodes if node.id == node_id), None)

    @property
    def node_ids(self) -> tuple[str, ...]:
        return tuple(node.id for node in self.nodes)

    def downstream_of(self, node_id: str) -> tuple[str, ...]:
        """Nos alcancaveis a partir de ``node_id``, em ordem topologica."""
        reachable: list[str] = []
        frontier = [node_id]
        while frontier:
            current = frontier.pop(0)
            for edge in self.edges:
                if edge.source == current and edge.target not in reachable:
                    reachable.append(edge.target)
                    frontier.append(edge.target)
        return tuple(reachable)


@dataclass(frozen=True)
class RunViewModel:
    """A execucao como entidade primaria da interface."""

    id: str
    status: StageStatus
    started_at: Optional[str]
    finished_at: Optional[str]
    duration_seconds: Optional[float]
    input_name: str
    input_count: int
    output_count: int
    approved: int
    rejected: int
    transformed: int
    duplicates: int
    configuration: dict[str, object]
    pipeline_version: str
    policy_hash: str
    graph: GraphContract

    @property
    def status_text(self) -> str:
        return STATUS_TEXT[self.status]

    @property
    def glyph(self) -> str:
        return STATUS_GLYPHS[self.status]

    @property
    def approval_rate(self) -> float:
        return self.approved / self.input_count if self.input_count else 0.0


# --- Construcao a partir do relatorio -------------------------------------------


def _transformed_at(records: Sequence[CurationRecord], node_id: str) -> int:
    """Registros que sofreram transformacao observada neste estagio."""
    return sum(
        1
        for record in records
        if any(event.stage.value == node_id for event in record.transformations)
    )


def _rejected_ids_at(records: Sequence[CurationRecord], node_id: str) -> tuple[str, ...]:
    return tuple(
        record.input_id
        for record in records
        if not record.passed
        and record.rejection_stage is not None
        and record.rejection_stage.value == node_id
    )


def _artifacts_for(
    report: RunReport, stage: StageReport, rejected_ids: tuple[str, ...]
) -> tuple[NodeArtifact, ...]:
    artifacts: list[NodeArtifact] = []
    if rejected_ids:
        artifacts.append(
            NodeArtifact(
                name=f"removidas_{stage.name.lower()}.csv",
                description="Estruturas removidas neste estagio, com motivo",
                record_ids=rejected_ids,
            )
        )
    if stage.name == DEDUP_STAGE and report.index is not None:
        collisions = tuple(
            collision.incoming_id for collision in report.index.collisions
        )
        if collisions:
            artifacts.append(
                NodeArtifact(
                    name="conflicts.csv",
                    description="Colisoes de identidade, incluindo as que nao sao "
                    "duplicatas",
                    record_ids=collisions,
                )
            )
    return tuple(artifacts)


def build_node(report: RunReport, stage: StageReport) -> NodeContract:
    rejected_ids = _rejected_ids_at(report.records, stage.name)
    wanted = NODE_PARAMETERS.get(stage.name, ())
    parameters = {
        key: value
        for key, value in report.provenance.parameters.items()
        if key in wanted
    }
    return NodeContract(
        id=stage.name,
        label=label_for(stage.name),
        description=stage.description,
        status=stage.status,
        input_count=stage.n_input,
        output_count=stage.n_output,
        transformed_count=_transformed_at(report.records, stage.name),
        rejected_count=stage.n_excluded,
        duration_seconds=stage.duration_seconds,
        parameters=parameters,
        exclusions=tuple(
            group for group in report.exclusion_groups() if group.stage == stage.name
        ),
        artifacts=_artifacts_for(report, stage, rejected_ids),
        help_text=NODE_HELP.get(stage.name, stage.description),
    )


def build_graph(report: RunReport) -> GraphContract:
    """Monta o grafo na ordem real de execucao do motor.

    As arestas sao encadeadas na sequencia dos estagios; o layout e deterministico
    para uma mesma definicao de pipeline e nao depende dos dados da execucao.
    """
    stages = list(report.stages) + (
        [report.dedup] if report.dedup is not None else []
    )
    nodes = tuple(build_node(report, stage) for stage in stages)
    edges = tuple(
        GraphEdge(source=left.id, target=right.id)
        for left, right in zip(nodes, nodes[1:])
    )
    return GraphContract(nodes=nodes, edges=edges)


def run_status(graph: GraphContract) -> StageStatus:
    """Status agregado da execucao, derivado dos nos.

    Uma falha em qualquer no torna a execucao falha; qualquer remocao a torna
    concluida com avisos. Nada aqui e estimado.
    """
    statuses = {node.status for node in graph.nodes}
    if StageStatus.FAILED in statuses:
        return StageStatus.FAILED
    if StageStatus.RUNNING in statuses:
        return StageStatus.RUNNING
    if StageStatus.WARNING in statuses:
        return StageStatus.WARNING
    if statuses == {StageStatus.SKIPPED}:
        return StageStatus.SKIPPED
    return StageStatus.SUCCESS


def build_run_view(report: RunReport) -> RunViewModel:
    """Adapta um :class:`~curation.reporting.RunReport` ao modelo da interface."""
    graph = build_graph(report)
    kpis = report.kpis()
    provenance = report.provenance
    last = graph.nodes[-1] if graph.nodes else None

    return RunViewModel(
        id=provenance.run_id,
        status=run_status(graph),
        started_at=provenance.started_at,
        finished_at=provenance.finished_at,
        duration_seconds=provenance.duration_seconds,
        input_name=provenance.input_name,
        input_count=report.total,
        output_count=last.output_count if last is not None else 0,
        approved=kpis["Approved"],
        rejected=kpis["Rejected"],
        transformed=kpis["Warnings"],
        duplicates=kpis["Duplicates"],
        configuration=dict(provenance.parameters),
        pipeline_version=provenance.pipeline_version,
        policy_hash=provenance.policy_hash,
        graph=graph,
    )
