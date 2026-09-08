"""Testes do componente de grafo interativo.

O layout e o payload sao calculados em Python de proposito: geometria
deterministica e testavel sem navegador. O que so o navegador exercita, a
interacao de zoom e arrasto, fica declarado como nao coberto no fim do arquivo.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from curation.components.dag import (
    BORDER_RADIUS,
    NODE_GAP,
    NODE_HEIGHT,
    NODE_WIDTH,
    PADDING,
    STATUS_COLORS,
    canvas_size,
    graph_payload,
    layout_nodes,
)
from curation.pipeline import CurationPipeline
from curation.reporting import StageStatus
from curation.viewmodel import build_run_view

POLICY_HASH = "0" * 64
FRONTEND = Path(__file__).resolve().parents[1] / "src/curation/components/dag/frontend"


@pytest.fixture(scope="module")
def graph():
    source = "CCO\nCC(=O)O[Na]\nC(C)(C)(C)(C)C\nCCO\n"
    pipeline = CurationPipeline(POLICY_HASH)
    return build_run_view(
        pipeline.run_report(source, input_bytes=source.encode())
    ).graph


# --- Layout -------------------------------------------------------------------


def test_layout_covers_every_node(graph) -> None:
    assert [node.id for node in layout_nodes(graph)] == list(graph.node_ids)


def test_layout_uses_the_specified_geometry(graph) -> None:
    for node in layout_nodes(graph):
        assert (node.width, node.height) == (NODE_WIDTH, NODE_HEIGHT)


def test_nodes_are_stacked_without_overlap(graph) -> None:
    positions = layout_nodes(graph)
    for previous, following in zip(positions, positions[1:]):
        assert following.y == previous.bottom_y + NODE_GAP
        assert following.y > previous.bottom_y


def test_layout_is_independent_of_the_batch_data() -> None:
    """Topologia deterministica: os dados da execucao nao movem os nos."""
    pipeline = CurationPipeline(POLICY_HASH)
    left = build_run_view(pipeline.run_report("CCO\n", input_bytes=b"a")).graph
    right = build_run_view(
        pipeline.run_report("C(C)(C)(C)(C)C\nCCO\nCCN\n", input_bytes=b"b")
    ).graph
    assert layout_nodes(left) == layout_nodes(right)


def test_canvas_encloses_every_node(graph) -> None:
    width, height = canvas_size(graph)
    for node in layout_nodes(graph):
        assert node.x + node.width + PADDING <= width
        assert node.bottom_y + PADDING <= height


def test_empty_graph_has_a_usable_canvas() -> None:
    from curation.viewmodel import GraphContract

    width, height = canvas_size(GraphContract(nodes=(), edges=()))
    assert width > 0 and height > 0
    assert layout_nodes(GraphContract(nodes=(), edges=())) == ()


# --- Payload -------------------------------------------------------------------------


def test_payload_is_json_serialisable(graph) -> None:
    """Atravessa a fronteira para o navegador, entao precisa ser JSON puro."""
    json.dumps(graph_payload(graph, None))


def test_payload_carries_no_chemical_records(graph) -> None:
    """Nenhum SMILES ou InChIKey cruza para o iframe.

    O grafo desenha contagens e status; o detalhamento quimico continua sendo
    renderizado pelo Streamlit, onde o conteudo nao sai do processo Python.
    """
    text = json.dumps(graph_payload(graph, None))
    for leaked in ("smiles", "inchikey", "curated", "raw_"):
        assert leaked not in text.lower()


def test_payload_edges_connect_declared_nodes(graph) -> None:
    payload = graph_payload(graph, None)
    ids = {node["id"] for node in payload["nodes"]}
    for edge in payload["edges"]:
        assert edge["source"] in ids and edge["target"] in ids


def test_edge_geometry_runs_downward(graph) -> None:
    for edge in graph_payload(graph, None)["edges"]:
        assert edge["y2"] > edge["y1"], "a aresta deve descer do no de origem"


def test_selection_is_passed_through(graph) -> None:
    assert graph_payload(graph, "STANDARDIZE")["selected"] == "STANDARDIZE"
    assert graph_payload(graph, None)["selected"] is None


def test_every_node_carries_glyph_and_text(graph) -> None:
    """Acessibilidade: o estado nunca depende so da cor."""
    for node in graph_payload(graph, None)["nodes"]:
        assert node["glyph"]
        assert node["statusText"]
        assert node["color"].startswith("#")


def test_every_status_has_a_colour() -> None:
    assert set(STATUS_COLORS) == set(StageStatus)


def test_payload_geometry_matches_the_module_constants(graph) -> None:
    geometry = graph_payload(graph, None)["geometry"]
    assert geometry == {
        "nodeWidth": NODE_WIDTH,
        "nodeHeight": NODE_HEIGHT,
        "radius": BORDER_RADIUS,
    }


# --- Frontend ---------------------------------------------------------------------------


def test_frontend_entry_point_exists() -> None:
    assert (FRONTEND / "index.html").is_file()


def test_frontend_declares_the_streamlit_protocol() -> None:
    """Sem estas tres mensagens o componente nunca aparece nem devolve valor."""
    source = (FRONTEND / "index.html").read_text(encoding="utf-8")
    for message in (
        "streamlit:componentReady",
        "streamlit:setComponentValue",
        "streamlit:setFrameHeight",
        "streamlit:render",
    ):
        assert message in source, message


def test_frontend_has_no_external_dependency() -> None:
    """Sem build npm e sem CDN: arquivos estaticos que o Cloud consegue servir."""
    source = (FRONTEND / "index.html").read_text(encoding="utf-8")
    assert not re.search(r'src\s*=\s*"https?://', source)
    assert not re.search(r'href\s*=\s*"https?://', source)


def test_frontend_implements_zoom_pan_and_fit() -> None:
    source = (FRONTEND / "index.html").read_text(encoding="utf-8")
    for capability in ("zoom-in", "zoom-out", "fit", "mousedown", "wheel"):
        assert capability in source, capability


def test_frontend_nodes_are_keyboard_reachable() -> None:
    """Selecao por teclado, nao so por clique."""
    source = (FRONTEND / "index.html").read_text(encoding="utf-8")
    assert 'tabindex: "0"' in source
    assert "keydown" in source
    assert "aria-label" in source


# --- NAO TESTADO -------------------------------------------------------------------------
#
# Comportamento que exige um navegador real e nao esta coberto aqui:
#   - o gesto de zoom, arrasto e enquadramento;
#   - a entrega efetiva do valor clicado ao Streamlit pelo postMessage;
#   - o reenquadramento na primeira renderizacao.
# Cobri-los exigiria um teste de navegador (Playwright ou equivalente), que
# adicionaria uma dependencia pesada e um servidor em execucao. O que da para
# verificar sem navegador esta acima: geometria, payload e contrato do protocolo.
