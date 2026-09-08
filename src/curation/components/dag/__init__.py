"""Grafo interativo do pipeline como componente Streamlit.

Implementado em JavaScript puro, sem React e sem etapa de build. A razao e
operacional: um componente com toolchain npm exigiria bundle compilado versionado
no repositorio e teria que sobreviver ao deploy no Streamlit Cloud, que instala
apenas o que esta em ``requirements.txt``. Arquivos estaticos servidos pelo proprio
Streamlit nao tem esse problema.

O componente devolve o identificador do no clicado, que a aplicacao usa como
``selected_node``. Sem esse canal de retorno o clique no no nao existiria: um
``st.components.v1.html`` comum e um iframe isolado, sem caminho de volta ao
Python.

O layout e calculado aqui, em Python, e nao no navegador: e geometria
deterministica, depende so da definicao do pipeline e nao dos dados da execucao,
entao pode ser testada sem abrir um navegador.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import streamlit.components.v1 as components

from curation.reporting import StageStatus
from curation.viewmodel import GraphContract

_FRONTEND = Path(__file__).parent / "frontend"

#: Dimensoes fixas, conforme o plano de redesenho.
NODE_WIDTH = 180
NODE_HEIGHT = 90
NODE_GAP = 46
BORDER_RADIUS = 6
PADDING = 12

#: Cor por status. Nunca e o unico portador do estado: todo no carrega glifo e
#: texto, e a interface permanece legivel em escala de cinza.
STATUS_COLORS: dict[StageStatus, str] = {
    StageStatus.PENDING: "#9aa0a6",
    StageStatus.RUNNING: "#1a73e8",
    StageStatus.SUCCESS: "#1e8e3e",
    StageStatus.WARNING: "#b06000",
    StageStatus.FAILED: "#d93025",
    StageStatus.SKIPPED: "#9aa0a6",
}


@dataclass(frozen=True)
class NodeLayout:
    """Posicao de um no no plano do grafo."""

    id: str
    x: int
    y: int
    width: int = NODE_WIDTH
    height: int = NODE_HEIGHT

    @property
    def center_x(self) -> int:
        return self.x + self.width // 2

    @property
    def bottom_y(self) -> int:
        return self.y + self.height


def layout_nodes(graph: GraphContract) -> tuple[NodeLayout, ...]:
    """Empilha os nos verticalmente, na ordem de execucao.

    Deterministico para uma mesma definicao de pipeline: a posicao depende so do
    indice do no, nunca das contagens do lote.
    """
    return tuple(
        NodeLayout(id=node.id, x=PADDING, y=PADDING + index * (NODE_HEIGHT + NODE_GAP))
        for index, node in enumerate(graph.nodes)
    )


def canvas_size(graph: GraphContract) -> tuple[int, int]:
    count = len(graph.nodes)
    if not count:
        return (NODE_WIDTH + 2 * PADDING, NODE_HEIGHT + 2 * PADDING)
    height = PADDING * 2 + count * NODE_HEIGHT + (count - 1) * NODE_GAP
    return (NODE_WIDTH + 2 * PADDING, height)


def graph_payload(graph: GraphContract, selected: Optional[str]) -> dict:
    """Serializa o grafo para o navegador.

    Envia apenas o que o desenho precisa. Nenhum registro quimico atravessa a
    fronteira: o detalhamento continua sendo renderizado pelo Streamlit.
    """
    positions = {node.id: node for node in layout_nodes(graph)}
    width, height = canvas_size(graph)

    return {
        "width": width,
        "height": height,
        "selected": selected,
        "geometry": {
            "nodeWidth": NODE_WIDTH,
            "nodeHeight": NODE_HEIGHT,
            "radius": BORDER_RADIUS,
        },
        "nodes": [
            {
                "id": node.id,
                "label": node.label,
                "glyph": node.glyph,
                "statusText": node.status_text,
                "color": STATUS_COLORS[node.status],
                "count": node.output_count,
                "removed": node.rejected_count,
                "x": positions[node.id].x,
                "y": positions[node.id].y,
            }
            for node in graph.nodes
        ],
        "edges": [
            {
                "source": edge.source,
                "target": edge.target,
                "x": positions[edge.source].center_x,
                "y1": positions[edge.source].bottom_y,
                "y2": positions[edge.target].y,
            }
            for edge in graph.edges
        ],
    }


_component = components.declare_component("curation_dag", path=str(_FRONTEND))


def render_dag(
    graph: GraphContract,
    selected: Optional[str] = None,
    height: int = 560,
    key: str = "curation_dag",
) -> Optional[str]:
    """Desenha o grafo e devolve o identificador do no clicado.

    Returns:
        O ``id`` do no selecionado pelo usuario, ou ``selected`` quando nada foi
        clicado nesta interacao. Nunca inventa uma selecao: se o valor devolvido
        pelo navegador nao corresponder a um no existente, e descartado.
    """
    clicked = _component(
        payload=graph_payload(graph, selected),
        height=height,
        key=key,
        default=selected,
    )
    if clicked in graph.node_ids:
        return clicked
    return selected
