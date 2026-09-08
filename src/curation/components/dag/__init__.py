"""Grafo interativo do pipeline como componente Streamlit usando streamlit-flow.

Refatorado para eliminar o hack de HTML/JS estático. Utiliza o pacote `streamlit-flow`
para injetar o DAG renderizado em React Flow diretamente no Streamlit sem necessidade
de build step (npm/vite) no repositório.
"""

from __future__ import annotations
from typing import Optional

from streamlit_flow import streamlit_flow, StreamlitFlowNode, StreamlitFlowEdge
from streamlit_flow.state import StreamlitFlowState

from curation.reporting import StageStatus
from curation.viewmodel import GraphContract

STATUS_BORDER_COLORS: dict[StageStatus, str] = {
    StageStatus.SUCCESS: "#10B981",   # Emerald 500
    StageStatus.WARNING: "#F59E0B",   # Amber 500
    StageStatus.FAILED: "#EF4444",    # Rose 500
    StageStatus.RUNNING: "#3B82F6",   # Blue 500
    StageStatus.PENDING: "#475569",   # Slate 630
    StageStatus.SKIPPED: "#334155",   # Slate 700
}

def render_dag(
    graph: GraphContract,
    selected: Optional[str] = None,
    height: int = 560,
    key: str = "curation_dag",
) -> Optional[str]:
    """Desenha o grafo interativo e devolve o identificador do no clicado.

    Returns:
        O ``id`` do nó selecionado pelo usuário, ou ``selected`` quando nada foi
        clicado. Preserva a seleção anterior se o clique for fora do nó.
    """
    flow_nodes = []
    
    # Layout horizontal linear com espaçamento fixo e limpo
    x_offset = 50
    y_pos = 120
    spacing = 220

    for i, node in enumerate(graph.nodes):
        border_color = STATUS_BORDER_COLORS.get(node.status, "#475569")
        is_selected = (node.id == selected)

        # Label científico enxuto
        content = (
            f"**[{node.id}]**\n\n"
            f"`In: {node.input_count} | Out: {node.output_count}`\n\n"
            f"`Rej: {node.rejected_count}`"
        )

        flow_nodes.append(
            StreamlitFlowNode(
                id=node.id,
                pos=(x_offset + (i * spacing), y_pos),
                data={"label": content},
                node_type="default",
                style={
                    "background": "#0F172A",
                    "color": "#F8FAFC",
                    "border": f"{'2px' if is_selected else '1px'} solid {border_color}",
                    "borderRadius": "6px",
                    "padding": "10px 14px",
                    "fontFamily": "ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace",
                    "fontSize": "11px",
                    "boxShadow": "0 4px 6px -1px rgba(0, 0, 0, 0.3)" if is_selected else "none",
                }
            )
        )

    flow_edges = [
        StreamlitFlowEdge(
            id=f"{edge.source}->{edge.target}",
            source=edge.source,
            target=edge.target,
            animated=True,
            style={"stroke": "#64748B", "strokeWidth": 2}
        )
        for edge in graph.edges
    ]

    state = StreamlitFlowState(nodes=flow_nodes, edges=flow_edges, selected_id=selected)

    updated_state = streamlit_flow(
        key=key,
        state=state,
        fit_view=True,
        show_controls=True,
        height=height,
        enable_pane_menu=False,
        enable_node_menu=False,
    )
    
    # Sincronização de Estado Segura: 
    # Se updated_state vier preenchido e com um nó válido, retorna ele.
    # Caso contrário (clique no vazio), preserva a visualização anterior (selected)
    new_selected = updated_state.selected_id if updated_state else None
    
    if new_selected in graph.node_ids:
        return new_selected
    return selected
