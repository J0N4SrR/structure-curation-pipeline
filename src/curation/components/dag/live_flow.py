from __future__ import annotations

from typing import Optional
from streamlit_flow.elements import StreamlitFlowNode, StreamlitFlowEdge
from streamlit_flow.state import StreamlitFlowState
from curation.reporting import StageStatus

# Estilos científicos rigorosos por estado Kubeflow
STAGE_STYLES = {
    StageStatus.PENDING: {
        "background": "#0B0F19",
        "color": "#64748B",
        "border": "1px dashed #334155",
        "borderRadius": "6px",
        "padding": "10px 14px",
        "fontFamily": "ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace",
        "fontSize": "11px",
    },
    StageStatus.RUNNING: {
        "background": "#0F172A",
        "color": "#38BDF8",
        "border": "2px solid #0284C7",
        "boxShadow": "0 0 15px rgba(2, 132, 199, 0.4)",
        "borderRadius": "6px",
        "padding": "10px 14px",
        "fontFamily": "ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace",
        "fontSize": "11px",
    },
    StageStatus.SUCCESS: {
        "background": "#0F172A",
        "color": "#F8FAFC",
        "border": "1px solid #10B981",
        "borderRadius": "6px",
        "padding": "10px 14px",
        "fontFamily": "ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace",
        "fontSize": "11px",
    },
    StageStatus.WARNING: {
        "background": "#0F172A",
        "color": "#F8FAFC",
        "border": "1px solid #F59E0B",
        "borderRadius": "6px",
        "padding": "10px 14px",
        "fontFamily": "ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace",
        "fontSize": "11px",
    },
    StageStatus.FAILED: {
        "background": "#1E1215",
        "color": "#FDA4AF",
        "border": "1px solid #F43F5E",
        "borderRadius": "6px",
        "padding": "10px 14px",
        "fontFamily": "ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace",
        "fontSize": "11px",
    },
    StageStatus.SKIPPED: {
        "background": "#0F172A",
        "color": "#94A3B8",
        "border": "1px dashed #475569",
        "borderRadius": "6px",
        "padding": "10px 14px",
        "fontFamily": "ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace",
        "fontSize": "11px",
    }
}

ORDERED_STAGES = [
    ("PARSE", "Leitura SMILES"),
    ("STANDARDIZE", "Padronização"),
    ("GET_PARENT", "Estrutura-Mãe"),
    ("VALENCE_GATE", "Validação Valência"),
    ("ELIGIBILITY", "Elegibilidade"),
    ("CANONICALIZE", "Identidade InChIKey"),
    ("DEDUPLICATION", "Deduplicação")
]

def build_live_dag_state(
    stage_states: dict[str, dict], 
    selected_node: Optional[str] = None
) -> StreamlitFlowState:
    nodes = []
    edges = []
    x_offset = 40
    spacing = 220
    y_pos = 120

    for idx, (stage_id, label) in enumerate(ORDERED_STAGES):
        state_data = stage_states.get(stage_id, {
            "status": StageStatus.PENDING,
            "in": "-",
            "out": "-",
            "rej": "-"
        })
        
        status = state_data["status"]
        status_icon = "⏳" if status == StageStatus.PENDING else "🔄" if status == StageStatus.RUNNING else "✓" if status == StageStatus.SUCCESS else "!"

        card_content = (
            f"**{status_icon} [{stage_id}]**\n\n"
            f"`In: {state_data['in']} | Out: {state_data['out']}`\n\n"
            f"`Rej: {state_data['rej']}`"
        )
        
        # Apply selection highlight if this node is selected
        node_style = dict(STAGE_STYLES.get(status, STAGE_STYLES[StageStatus.PENDING]))
        if selected_node == stage_id:
            node_style["boxShadow"] = "0 4px 6px -1px rgba(0, 0, 0, 0.3)"
            if status != StageStatus.RUNNING:
                node_style["border"] = node_style["border"].replace("1px", "2px")

        nodes.append(
            StreamlitFlowNode(
                id=stage_id,
                pos=(x_offset + (idx * spacing), y_pos),
                data={"label": card_content},
                node_type="default",
                style=node_style
            )
        )

        if idx > 0:
            prev_stage = ORDERED_STAGES[idx - 1][0]
            # Aresta anima enquanto o nó destino ou fonte estiverem em execução
            is_active = status == StageStatus.RUNNING or stage_states.get(prev_stage, {}).get("status") == StageStatus.RUNNING
            edges.append(
                StreamlitFlowEdge(
                    id=f"{prev_stage}->{stage_id}",
                    source=prev_stage,
                    target=stage_id,
                    animated=is_active,
                    style={"stroke": "#0284C7" if is_active else "#334155", "strokeWidth": 2}
                )
            )

    return StreamlitFlowState(nodes=nodes, edges=edges, selected_id=selected_node)
