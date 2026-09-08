"""Interface Streamlit (Wizard UI) para o Structure Curation Pipeline.

Refatoração baseada em princípios de HCI, UX, A11y e Frontend Design:
- 3 Estados Principais: INPUT, CURATING, COMPLETE.
- Contraste de cores WCAG 2.1 AA e suporte a acessibilidade (aria-live, role="status").
- Feedback visual com pré-visualização de amostras e animações de progresso.
- Design limpo, glassmorphism sutil e apresentação guiada.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Optional

import streamlit as st
from rdkit import Chem, RDLogger

try:
    from rdkit.Chem import Draw

    DRAWING_AVAILABLE = True
    DRAWING_ERROR = ""
except ImportError as _error:
    Draw = None
    DRAWING_AVAILABLE = False
    DRAWING_ERROR = str(_error)

from curation.filters import EligibilityCriteria
from curation.io import compute_policy_hash, preview_input
from curation.pipeline import CurationPipeline
from curation.reporting import (
    DEDUP_STAGE,
    RunReport,
    StageReport,
    StageStatus,
    full_csv,
    rejected_csv,
    run_manifest,
    to_csv,
)

RDLogger.DisableLog("rdApp.*")

DEFAULT_DECISIONS = Path("docs/decisions.md")
MIN_STAGE_DISPLAY_TIME = 1.5  # segundos de apresentação por etapa na UI

# --- Definição das Etapas do Wizard ------------------------------------------------

WIZARD_STAGES = [
    {
        "id": "INPUT",
        "title": "Input",
        "description": "Recepção e triagem inicial do lote de estruturas.",
    },
    {
        "id": "STANDARDIZATION",
        "title": "Standardization",
        "description": "Normalização de grupos funcionais, cargas e aromaticidade.",
    },
    {
        "id": "PARENT",
        "title": "Parent Structure",
        "description": "Isolamento da estrutura-mãe (remoção de sais e solventes).",
    },
    {
        "id": "VALIDATION",
        "title": "Validation",
        "description": "Validação de integridade química e portão de valência.",
    },
    {
        "id": "DEDUPLICATION",
        "title": "Deduplication",
        "description": "Identificação de duplicatas por InChIKey completo.",
    },
    {
        "id": "ELIGIBILITY",
        "title": "Eligibility",
        "description": "Filtros de massa molecular e número de átomos pesados.",
    },
    {
        "id": "OUTPUT",
        "title": "Output",
        "description": "Geração do dataset canônico final e estatísticas.",
    },
]


def inject_styles() -> None:
    st.markdown(
        """
        <style>
          /* Estilos globais e acessibilidade */
          .block-container {
            max-width: 800px;
            padding-top: 2rem;
            padding-bottom: 3.5rem;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
          }
          
          /* Cabeçalho do Wizard */
          .wizard-header {
            text-align: center;
            margin-bottom: 2.2rem;
          }
          .wizard-title {
            font-size: 2.3rem;
            font-weight: 800;
            letter-spacing: -0.03em;
            margin-bottom: 0.2rem;
          }
          .wizard-sub {
            font-size: 1.05rem;
            opacity: 0.75;
          }
          
          /* Cards do Wizard com Glassmorphism sutil */
          .wizard-card {
            border: 1px solid rgba(128, 128, 128, 0.25);
            border-radius: 12px;
            padding: 1.2rem 1.5rem;
            margin-bottom: 0.8rem;
            background: rgba(255, 255, 255, 0.03);
            backdrop-filter: blur(8px);
            transition: border-color 0.2s ease, box-shadow 0.2s ease;
          }
          .wizard-card:focus-within {
            outline: 2px solid #1a73e8;
            outline-offset: 2px;
          }
          
          /* Animação pulsante acessível para etapa em execução */
          @keyframes pulse-running {
            0% { border-color: rgba(26, 115, 232, 0.4); box-shadow: 0 0 8px rgba(26, 115, 232, 0.15); }
            50% { border-color: rgba(26, 115, 232, 0.9); box-shadow: 0 0 16px rgba(26, 115, 232, 0.35); }
            100% { border-color: rgba(26, 115, 232, 0.4); box-shadow: 0 0 8px rgba(26, 115, 232, 0.15); }
          }
          
          .wizard-card-running {
            animation: pulse-running 2s infinite ease-in-out;
            background: rgba(26, 115, 232, 0.04);
          }
          .wizard-card-completed {
            border-color: rgba(30, 142, 62, 0.45);
            background: rgba(30, 142, 62, 0.02);
          }
          .wizard-card-failed {
            border-color: rgba(217, 48, 37, 0.45);
            background: rgba(217, 48, 37, 0.02);
          }
          
          /* Badges e cores de contraste alto (WCAG AA) */
          .status-badge {
            font-weight: 700;
            font-size: 1.1rem;
            margin-right: 0.7rem;
            display: inline-block;
          }
          .status-pending { color: #70757a; }
          .status-running { color: #1a73e8; }
          .status-completed { color: #1e8e3e; }
          .status-failed { color: #d93025; }
          
          /* Indicador de progresso no topo */
          .progress-tracker {
            font-size: 0.85rem;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: #1a73e8;
            margin-bottom: 0.4rem;
          }
          .summary-number {
            font-size: 2.1rem;
            font-weight: 800;
            color: #1e8e3e;
            letter-spacing: -0.02em;
          }
          .step-connector {
            text-align: center;
            font-size: 1.2rem;
            opacity: 0.35;
            margin: -0.4rem 0;
            user-select: none;
          }
        </style>
        """,
        unsafe_allow_html=True,
    )


# --- Estado da Aplicação -------------------------------------------------------------


def init_session_state() -> None:
    if "ui_state" not in st.session_state:
        st.session_state["ui_state"] = "INPUT"  # INPUT | CURATING | COMPLETE
    if "curating_step_idx" not in st.session_state:
        st.session_state["curating_step_idx"] = 0
    if "step_start_time" not in st.session_state:
        st.session_state["step_start_time"] = 0.0
    if "raw_input" not in st.session_state:
        st.session_state["raw_input"] = None
    if "input_name" not in st.session_state:
        st.session_state["input_name"] = ""
    if "config" not in st.session_state:
        st.session_state["config"] = {
            "max_mw": 1000.0,
            "max_ha": 100,
            "deduplicate": True,
            "policy_hash": compute_policy_hash(DEFAULT_DECISIONS)
            if DEFAULT_DECISIONS.is_file()
            else "UNVERSIONED_POLICY",
            "policy_path": str(DEFAULT_DECISIONS),
        }
    if "report" not in st.session_state:
        st.session_state["report"] = None


def reset_to_input() -> None:
    st.session_state["ui_state"] = "INPUT"
    st.session_state["curating_step_idx"] = 0
    st.session_state["step_start_time"] = 0.0
    st.session_state["raw_input"] = None
    st.session_state["input_name"] = ""
    st.session_state["report"] = None
    st.rerun()


# --- TELA 1: INPUT -------------------------------------------------------------------


def render_input_screen() -> None:
    st.markdown(
        """
        <div class='wizard-header' role='region' aria-label='Structure Curation Header'>
            <div class='wizard-title'>Structure Curation</div>
            <div class='wizard-sub'>Curate your molecular library</div>
            <div style='font-size:0.85rem; opacity:0.6; margin-top:0.3rem;'>Standardize, validate and deduplicate chemical structures</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    top_col1, top_col2 = st.columns([6, 1])
    with top_col2:
        with st.popover("⚙", help="Configuration Options"):
            st.markdown("### Configuration")
            cfg_mw = st.number_input(
                "Molecular weight maximum (Da)",
                min_value=50.0,
                max_value=10000.0,
                value=st.session_state["config"]["max_mw"],
                step=50.0,
                help="Maximum allowed molecular weight on parent structure.",
            )
            cfg_ha = st.number_input(
                "Heavy atoms maximum",
                min_value=5,
                max_value=1000,
                value=st.session_state["config"]["max_ha"],
                step=5,
                help="Maximum heavy atoms allowed on parent structure.",
            )
            cfg_dedup = st.checkbox(
                "Deduplicate structures",
                value=st.session_state["config"]["deduplicate"],
                help="Deduplicate exact structures based on full InChIKey.",
            )

            st.session_state["config"]["max_mw"] = float(cfg_mw)
            st.session_state["config"]["max_ha"] = int(cfg_ha)
            st.session_state["config"]["deduplicate"] = bool(cfg_dedup)

    tab_paste, tab_file = st.tabs(["Paste SMILES", "Upload file"])
    raw_bytes: Optional[bytes] = None
    input_name = ""

    with tab_paste:
        text = st.text_area(
            "SMILES Input",
            height=160,
            placeholder="CCO\nCC(=O)O[Na]\nN[C@@H](C)C(=O)O.Cl\nc1ccccc1",
            label_visibility="collapsed",
            help="Paste one SMILES string per line",
        )
        if text.strip():
            raw_bytes = text.encode("utf-8")
            input_name = "pasted_structures.smi"

    with tab_file:
        uploaded = st.file_uploader(
            "Upload file",
            type=["csv", "tsv", "smi", "smiles", "txt"],
            label_visibility="collapsed",
            help="Upload a file containing chemical structures (CSV, TSV, SMI)",
        )
        if uploaded is not None:
            raw_bytes = uploaded.getvalue()
            input_name = uploaded.name

    if raw_bytes:
        try:
            preview = preview_input(raw_bytes.decode("utf-8", errors="replace"))
            total_detected = preview.total
            valid_detected = preview.total - preview.n_invalid
        except Exception:
            total_detected = len(raw_bytes.decode("utf-8", errors="ignore").splitlines())
            valid_detected = total_detected

        if total_detected > 0:
            st.success(f"✓ **{input_name}** ({total_detected} molecules detected)")
            st.session_state["raw_input"] = raw_bytes
            st.session_state["input_name"] = input_name

            # Pré-visualização de amostras para confirmação do usuário (HCI / Error Prevention)
            with st.expander("Preview sample molecules", expanded=False):
                try:
                    lines = [line.strip() for line in raw_bytes.decode("utf-8", errors="replace").splitlines() if line.strip()]
                    sample_rows = [{"Index": i + 1, "SMILES / Structure": line} for i, line in enumerate(lines[:5])]
                    st.dataframe(sample_rows, use_container_width=True, hide_index=True)
                except Exception:
                    st.caption("Sample preview unavailable.")
        else:
            st.warning("Please provide at least one valid molecule.")
            st.session_state["raw_input"] = None
    else:
        st.session_state["raw_input"] = None

    st.markdown("<br>", unsafe_allow_html=True)
    btn_container = st.columns([1, 2, 1])
    with btn_container[1]:
        start_disabled = st.session_state["raw_input"] is None
        if st.button(
            "Start curation",
            type="primary",
            use_container_width=True,
            disabled=start_disabled,
            help="Click to start the automated curation pipeline",
        ):
            # Iniciar execução do backend e mudar estado da UI para CURATING
            run_backend_pipeline()
            st.session_state["ui_state"] = "CURATING"
            st.session_state["curating_step_idx"] = 0
            st.session_state["step_start_time"] = time.time()
            st.rerun()


# --- Execução do Backend (sem sleeps químicos) --------------------------------------


def run_backend_pipeline() -> None:
    raw = st.session_state["raw_input"]
    name = st.session_state["input_name"]
    cfg = st.session_state["config"]

    pipeline = CurationPipeline(
        policy_hash=cfg["policy_hash"],
        criteria=EligibilityCriteria(
            max_molecular_weight=cfg["max_mw"],
            max_heavy_atoms=cfg["max_ha"],
        ),
        deduplicate=cfg["deduplicate"],
    )

    report = pipeline.run_report(
        raw.decode("utf-8", errors="replace"),
        parameters={
            "max_mw": cfg["max_mw"],
            "max_ha": cfg["max_ha"],
            "deduplicate": cfg["deduplicate"],
        },
        input_bytes=raw,
        input_name=name,
        policy_path=cfg["policy_path"],
    )
    st.session_state["report"] = report


# --- TELA 2 & 3: CURATING / COMPLETE (WIZARD) ----------------------------------------


def get_step_summary(step_id: str, report: Optional[RunReport]) -> str:
    if report is None:
        return "Processing..."

    kpis = report.kpis()
    total_in = kpis["Total Processed"]

    if step_id == "INPUT":
        return f"{total_in} molecules"

    if step_id == "STANDARDIZATION":
        std_stage = next((s for s in report.stages if s.name == "STANDARDIZE"), None)
        out_n = std_stage.n_output if std_stage else total_in
        return f"{out_n} processed"

    if step_id == "PARENT":
        parent_stage = next((s for s in report.stages if s.name == "GET_PARENT"), None)
        out_n = parent_stage.n_output if parent_stage else total_in
        return f"{out_n} parent structures"

    if step_id == "VALIDATION":
        val_stage = next((s for s in report.stages if s.name == "VALENCE_GATE"), None)
        valid = val_stage.n_output if val_stage else total_in
        rejected = val_stage.n_excluded if val_stage else 0
        return f"{valid} valid · {rejected} rejected"

    if step_id == "DEDUPLICATION":
        if report.dedup:
            return f"{report.dedup.n_output} unique · {report.dedup.n_excluded} duplicates"
        return f"{total_in} unique · 0 duplicates"

    if step_id == "ELIGIBILITY":
        el_stage = next((s for s in report.stages if s.name == "ELIGIBILITY"), None)
        eligible = el_stage.n_output if el_stage else total_in
        rejected = el_stage.n_excluded if el_stage else 0
        return f"{eligible} eligible · {rejected} rejected"

    if step_id == "OUTPUT":
        return f"{kpis['Approved']} curated structures"

    return ""


def render_step_detail(step_id: str, report: RunReport) -> None:
    kpis = report.kpis()

    if step_id == "INPUT":
        st.markdown(f"**Molecules:** {kpis['Total Processed']}")
        st.markdown(f"**Source File:** `{report.provenance.input_name}`")

    elif step_id == "STANDARDIZATION":
        std_stage = next((s for s in report.stages if s.name == "STANDARDIZE"), None)
        if std_stage:
            st.markdown(f"**Processed:** {std_stage.n_input}")
            st.markdown(f"**Standardized:** {std_stage.n_output}")
            st.markdown(f"**Exclusions:** {std_stage.n_excluded}")

    elif step_id == "PARENT":
        parent_stage = next((s for s in report.stages if s.name == "GET_PARENT"), None)
        if parent_stage:
            st.markdown(f"**Structures Processed:** {parent_stage.n_input}")
            st.markdown(f"**Parent Structures:** {parent_stage.n_output}")

    elif step_id == "VALIDATION":
        val_stage = next((s for s in report.stages if s.name == "VALENCE_GATE"), None)
        if val_stage:
            st.markdown(f"**Evaluated:** {val_stage.n_input}")
            st.markdown(f"**Valid:** {val_stage.n_output}")
            st.markdown(f"**Rejected:** {val_stage.n_excluded}")

    elif step_id == "DEDUPLICATION":
        if report.dedup:
            st.markdown(f"**Evaluated:** {report.dedup.n_input}")
            st.markdown(f"**Unique:** {report.dedup.n_output}")
            st.markdown(f"**Duplicates:** {report.dedup.n_excluded}")

    elif step_id == "ELIGIBILITY":
        el_stage = next((s for s in report.stages if s.name == "ELIGIBILITY"), None)
        if el_stage:
            st.markdown(f"**Evaluated:** {el_stage.n_input}")
            st.markdown(f"**Eligible:** {el_stage.n_output}")
            st.markdown(f"**Rejected:** {el_stage.n_excluded}")

    elif step_id == "OUTPUT":
        st.markdown(f"**Curated Dataset:** {kpis['Approved']} structures")
        st.markdown(f"**Total Rejected:** {kpis['Rejected']}")
        st.markdown(f"**Total Duplicates:** {kpis['Duplicates']}")


def render_curating_and_result_screen() -> None:
    report: Optional[RunReport] = st.session_state.get("report")
    is_complete = st.session_state["ui_state"] == "COMPLETE"

    # Gerenciamento de tempo de exibição da UI para animação do Wizard
    if not is_complete:
        current_idx = st.session_state["curating_step_idx"]
        elapsed = time.time() - st.session_state["step_start_time"]
        if elapsed >= MIN_STAGE_DISPLAY_TIME:
            if current_idx < len(WIZARD_STAGES) - 1:
                st.session_state["curating_step_idx"] = current_idx + 1
                st.session_state["step_start_time"] = time.time()
                st.rerun()
            else:
                st.session_state["ui_state"] = "COMPLETE"
                st.rerun()

    # Cabeçalho do Estado com Acessibilidade (role="status", aria-live="polite")
    if is_complete and report:
        kpis = report.kpis()
        st.markdown(
            f"""
            <div class='wizard-header' role='status' aria-live='polite'>
                <div class='wizard-title' style='color:#1e8e3e;'>Curation complete</div>
                <div class='summary-number'>{kpis['Total Processed']} input → {kpis['Approved']} curated</div>
            </div>
            """,
            unsafe_allow_html=True,
        )
    else:
        cur_num = st.session_state['curating_step_idx'] + 1
        st.markdown(
            f"""
            <div class='wizard-header' role='status' aria-live='polite'>
                <div class='progress-tracker'>STEP {cur_num} OF {len(WIZARD_STAGES)}</div>
                <div class='wizard-title'>Curation</div>
                <div class='wizard-sub'>Processing molecular library...</div>
            </div>
            """,
            unsafe_allow_html=True,
        )

    # Renderização da sequência de cards do Wizard
    cur_idx = st.session_state["curating_step_idx"] if not is_complete else len(WIZARD_STAGES) - 1

    for idx, stage_info in enumerate(WIZARD_STAGES):
        if is_complete or idx < cur_idx:
            status_icon = "✓"
            status_class = "status-completed"
            card_class = "wizard-card-completed"
            status_text = "Completed"
        elif idx == cur_idx and not is_complete:
            status_icon = "◉"
            status_class = "status-running"
            card_class = "wizard-card-running"
            status_text = "Running"
        else:
            status_icon = "○"
            status_class = "status-pending"
            card_class = ""
            status_text = "Pending"

        summary = get_step_summary(stage_info["id"], report) if (idx <= cur_idx or is_complete) else ""

        with st.container():
            col1, col2 = st.columns([5, 1])
            with col1:
                st.markdown(
                    f"""
                    <div class='wizard-card {card_class}' tabindex='0' aria-label='{stage_info["title"]} stage {status_text}'>
                        <div>
                            <span class='status-badge {status_class}'>{status_icon}</span>
                            <strong>{stage_info['title']}</strong>
                            <span style='font-size:0.8rem; opacity:0.55; margin-left:0.5rem;'>{status_text}</span>
                        </div>
                        <div style='font-size:0.9rem; opacity:0.8; margin-top:0.35rem;'>
                            {summary if summary else stage_info['description']}
                        </div>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )
            with col2:
                if (idx <= cur_idx or is_complete) and report:
                    with st.popover("Details", help=f"View details for {stage_info['title']}"):
                        st.markdown(f"### {stage_info['title']}")
                        st.caption(stage_info["description"])
                        st.divider()
                        render_step_detail(stage_info["id"], report)

        if idx < len(WIZARD_STAGES) - 1:
            st.markdown("<div class='step-connector'>↓</div>", unsafe_allow_html=True)

    # Botões de Download e Nova Análise (Apenas ao finalizar)
    if is_complete and report:
        st.markdown("<br><hr><br>", unsafe_allow_html=True)
        col_d1, col_d2, col_d3 = st.columns(3)

        with col_d1:
            st.download_button(
                "Download curated dataset",
                to_csv(report.approved, ("input_id", "raw_smiles", "curated_smiles", "inchikey", "status")),
                file_name="curated_structures.csv",
                mime="text/csv",
                use_container_width=True,
            )

        with col_d2:
            st.download_button(
                "Download rejected dataset",
                rejected_csv(report),
                file_name="rejected_structures.csv",
                mime="text/csv",
                use_container_width=True,
            )

        with col_d3:
            st.download_button(
                "Download audit log",
                full_csv(report),
                file_name="audit_log.csv",
                mime="text/csv",
                use_container_width=True,
            )

        st.markdown("<br>", unsafe_allow_html=True)
        new_col = st.columns([1, 2, 1])
        with new_col[1]:
            if st.button("New analysis", use_container_width=True, help="Reset and start a new analysis"):
                reset_to_input()

    # Rerun automático para atualização suave durante o estado CURATING
    if not is_complete:
        time.sleep(0.1)
        st.rerun()


# --- Aplicação Principal -------------------------------------------------------------


def main() -> None:
    st.set_page_config(page_title="Structure Curation", layout="centered")
    inject_styles()
    init_session_state()

    ui_state = st.session_state["ui_state"]

    if ui_state == "INPUT":
        render_input_screen()
    else:
        render_curating_and_result_screen()


if __name__ == "__main__":
    main()
