"""Testes da camada de adaptacao entre relatorios e interface.

Sao os testes que a interface nao precisa mais carregar: rotulos, arestas,
propagacao de falha e agrupamento por no ficam verificaveis sem um navegador.
"""

from __future__ import annotations

import pytest

from curation.pipeline import CurationPipeline
from curation.reporting import DEDUP_STAGE, PIPELINE_STAGES, StageStatus
from curation.viewmodel import (
    NODE_HELP,
    NODE_LABELS,
    STATUS_TEXT,
    GraphEdge,
    build_graph,
    build_run_view,
    label_for,
    run_status,
)

POLICY_HASH = "0" * 64

CORPUS = (
    "CCO\n"
    "CC(=O)O[Na]\n"
    "N[C@@H](C)C(=O)O.Cl\n"
    "C(C)(C)(C)(C)C\n"
    "C1CC\n"
    "CCO\n"
)


@pytest.fixture(scope="module")
def view():
    pipeline = CurationPipeline(POLICY_HASH)
    return build_run_view(
        pipeline.run_report(CORPUS, input_bytes=CORPUS.encode(), input_name="c.smi")
    )


# --- Grafo ------------------------------------------------------------------------


def test_graph_has_one_node_per_real_stage(view) -> None:
    """O grafo espelha os estagios que o motor executa, nao um resumo didatico."""
    expected = [name for name, _ in PIPELINE_STAGES] + [DEDUP_STAGE]
    assert list(view.graph.node_ids) == expected


def test_edges_are_explicit_and_chained(view) -> None:
    """A interface nao infere dependencia pela posicao: as arestas sao dados."""
    ids = view.graph.node_ids
    assert view.graph.edges == tuple(
        GraphEdge(source=left, target=right) for left, right in zip(ids, ids[1:])
    )


def test_every_edge_points_at_existing_nodes(view) -> None:
    for edge in view.graph.edges:
        assert view.graph.node(edge.source) is not None
        assert view.graph.node(edge.target) is not None


def test_downstream_is_topological(view) -> None:
    assert view.graph.downstream_of("ELIGIBILITY") == (
        "CANONICALIZE",
        DEDUP_STAGE,
    )
    assert view.graph.downstream_of(DEDUP_STAGE) == ()


def test_graph_layout_does_not_depend_on_the_data() -> None:
    """Topologia deterministica: dois lotes diferentes produzem o mesmo grafo."""
    pipeline = CurationPipeline(POLICY_HASH)
    left = build_graph(pipeline.run_report("CCO\n", input_bytes=b"CCO\n"))
    right = build_graph(
        pipeline.run_report("C(C)(C)(C)(C)C\n", input_bytes=b"x")
    )
    assert left.node_ids == right.node_ids
    assert left.edges == right.edges


# --- Continuidade dos numeros -------------------------------------------------------


def test_node_counts_are_continuous(view) -> None:
    nodes = view.graph.nodes
    for previous, following in zip(nodes, nodes[1:]):
        assert following.input_count == previous.output_count


def test_node_counts_come_from_the_report(view) -> None:
    """A interface nao recalcula: a soma das rejeicoes bate com o relatorio."""
    total_rejected = sum(node.rejected_count for node in view.graph.nodes)
    assert total_rejected == view.rejected + view.duplicates


def test_run_totals_match_the_first_and_last_node(view) -> None:
    assert view.input_count == view.graph.nodes[0].input_count
    assert view.output_count == view.graph.nodes[-1].output_count


# --- Status --------------------------------------------------------------------------


def test_status_is_never_conveyed_by_colour_alone(view) -> None:
    """Todo estado tem glifo e texto."""
    for node in view.graph.nodes:
        assert node.glyph
        assert node.status_text
        assert node.status_text != node.id


def test_every_status_has_a_label() -> None:
    assert set(STATUS_TEXT) == set(StageStatus)


def test_run_status_aggregates_the_nodes(view) -> None:
    assert view.status is run_status(view.graph)
    assert view.status is StageStatus.WARNING


def test_failure_propagates_to_downstream_nodes() -> None:
    """Se nenhuma estrutura sobrevive a um estagio, os seguintes ficam ignorados."""
    pipeline = CurationPipeline(POLICY_HASH)
    source = "C(C)(C)(C)(C)C\n"
    view = build_run_view(pipeline.run_report(source, input_bytes=source.encode()))

    failing = view.graph.node("STANDARDIZE")
    assert failing is not None and failing.status is StageStatus.FAILED

    for node_id in view.graph.downstream_of("STANDARDIZE"):
        node = view.graph.node(node_id)
        assert node is not None
        assert node.status is StageStatus.SKIPPED, node_id

    assert view.status is StageStatus.FAILED


def test_clean_run_reports_success() -> None:
    pipeline = CurationPipeline(POLICY_HASH)
    source = "CCO\nCCN\nCC(=O)Oc1ccccc1C(=O)O\n"
    view = build_run_view(pipeline.run_report(source, input_bytes=source.encode()))
    assert view.status is StageStatus.SUCCESS


# --- Vocabulario ------------------------------------------------------------------------


def test_every_stage_has_a_readable_label_and_help() -> None:
    for name, _ in PIPELINE_STAGES:
        assert name in NODE_LABELS
        assert name in NODE_HELP
    assert DEDUP_STAGE in NODE_LABELS and DEDUP_STAGE in NODE_HELP


def test_labels_are_not_internal_identifiers(view) -> None:
    for node in view.graph.nodes:
        assert node.label != node.id
        assert "_" not in node.label


def test_unknown_stage_falls_back_readably() -> None:
    assert label_for("NOVO_ESTAGIO") == "Novo Estagio"


# --- Parametros e artefatos ----------------------------------------------------------------


def test_only_relevant_parameters_reach_each_node(view) -> None:
    """Exibir todos os parametros em todo no esconderia qual governa o estagio."""
    eligibility = view.graph.node("ELIGIBILITY")
    assert eligibility is not None
    assert set(eligibility.parameters) == {"max_mw", "max_ha", "require_carbon"}

    parse = view.graph.node("PARSE")
    assert parse is not None and parse.parameters == {}


def test_artifacts_reference_records_not_generated_files(view) -> None:
    """O artefato aponta identificadores; a serializacao continua em reporting."""
    standardize = view.graph.node("STANDARDIZE")
    assert standardize is not None
    assert standardize.artifacts
    artifact = standardize.artifacts[0]
    assert artifact.count == standardize.rejected_count
    assert all(isinstance(identifier, str) for identifier in artifact.record_ids)


def test_nodes_without_removals_have_no_artifacts(view) -> None:
    canonicalize = view.graph.node("CANONICALIZE")
    assert canonicalize is not None
    assert canonicalize.rejected_count == 0
    assert canonicalize.artifacts == ()


def test_exclusions_are_scoped_to_their_node(view) -> None:
    for node in view.graph.nodes:
        for group in node.exclusions:
            assert group.stage == node.id


# --- Run view -----------------------------------------------------------------------------------


def test_run_view_carries_provenance_identity(view) -> None:
    assert view.id
    assert view.policy_hash == POLICY_HASH
    assert view.pipeline_version
    assert view.duration_seconds is not None and view.duration_seconds > 0


def test_configuration_is_the_stored_run_configuration(view) -> None:
    """A configuracao vem do Run, nao dos widgets da tela."""
    assert set(view.configuration) >= {"max_mw", "max_ha", "deduplicate"}


def test_approval_rate_is_derived_not_stored(view) -> None:
    assert view.approval_rate == pytest.approx(view.approved / view.input_count)


def test_empty_run_does_not_break_the_view() -> None:
    pipeline = CurationPipeline(POLICY_HASH)
    view = build_run_view(pipeline.run_report("", input_bytes=b""))
    assert view.input_count == 0
    assert view.status is StageStatus.SKIPPED
    assert view.approval_rate == 0.0
