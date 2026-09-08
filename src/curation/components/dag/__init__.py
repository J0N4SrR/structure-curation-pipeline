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

#: Dimensoes fixas para simular o layout anterior
NODE_WIDTH = 250
NODE_HEIGHT = 100
NODE_GAP = 46
PADDING = 12

#: Cor por status para bordas/backgrounds
STATUS_COLORS: dict[StageStatus, str] = {
    StageStatus.PENDING: "#9aa0a6",
    StageStatus.RUNNING: "#1a73e8",
    StageStatus.SUCCESS: "#1e8e3e",
    StageStatus.WARNING: "#b06000",
    StageStatus.FAILED: "#d93025",
    StageStatus.SKIPPED: "#9aa0a6",
}

def build_flow(graph: GraphContract) -> tuple[list[StreamlitFlowNode], list[StreamlitFlowEdge]]:
    """Converte o GraphContract para nós e arestas do streamlit-flow."""
    nodes = []
    
    for index, node in enumerate(graph.nodes):
        # Layout vertical empilhado
        x = PADDING
        y = PADDING + index * (NODE_HEIGHT + NODE_GAP)
        
        # Mapeando os contadores de aprovação/rejeição no label
        counts_info = f"({node.output_count} saídas)"
        if node.rejected_count > 0:
            counts_info = f"({node.output_count} saídas | {node.rejected_count} rej.)"
            
        label = f"{node.glyph} {node.label}\n{node.status_text}\n{counts_info}"
        
        color = STATUS_COLORS.get(node.status, "#475569")
        
        nodes.append(
            StreamlitFlowNode(
                id=node.id,
                pos=(x, y),
                data={"label": label},
                node_type="default",
                style={
                    "background": "#1E293B", 
                    "color": "#F8FAFC", 
                    "border": f"2px solid {color}", 
                    "padding": "10px", 
                    "borderRadius": "8px", 
                    "width": f"{NODE_WIDTH}px",
                    "whiteSpace": "pre-wrap"
                }
            )
        )
        
    edges = [
        StreamlitFlowEdge(id=f"{e.source}-{e.target}", source=e.source, target=e.target, animated=True)
        for e in graph.edges
    ]
    
    return nodes, edges

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
    nodes, edges = build_flow(graph)
    
    state = StreamlitFlowState(nodes=nodes, edges=edges, selected_id=selected)
    
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
