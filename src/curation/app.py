"""Interface Streamlit (Wizard UI) para o Structure Curation Pipeline.

Organizada estritamente no fluxo em 8 etapas:
1. CARREGAR
2. CONFIGURAR
3. EXECUTAR
4. VER O PIPELINE
5. ENTENDER O RESULTADO
6. INVESTIGAR UMA ESTRUTURA
7. BAIXAR OS DADOS
8. REPRODUZIR A EXECUÇÃO
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
    reproducibility_package,
    run_manifest,
    to_csv,
)

RDLogger.DisableLog("rdApp.*")

DEFAULT_DECISIONS = Path("docs/decisions.md")
MIN_STAGE_DISPLAY_TIME = 1.5

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


def render_mol_image(smiles: Optional[str], legend: str = "", size: tuple[int, int] = (250, 250)):
    """Gera uma imagem 2D da molécula via RDKit caso o RDKit esteja disponível."""
    if not smiles or not DRAWING_AVAILABLE or Draw is None:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Draw.MolToImage(mol, size=size, legend=legend)
    except Exception:
        pass
    return None


def inject_styles() -> None:
    st.markdown(
        """
        <style>
          .block-container {
            max-width: 900px;
            padding-top: 1.5rem;
            padding-bottom: 3rem;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
          }
          
          /* Stepper visual de 8 passos no topo */
          .stepper-bar {
            display: flex;
            flex-wrap: wrap;
            align-items: center;
            justify-content: space-between;
            gap: 0.3rem;
            margin-bottom: 2rem;
            padding: 0.8rem 1rem;
            background: rgba(255, 255, 255, 0.03);
            border: 1px solid rgba(128, 128, 128, 0.2);
            border-radius: 12px;
          }
          .step-pill {
            font-size: 0.72rem;
            font-weight: 700;
            padding: 0.25rem 0.55rem;
            border-radius: 6px;
            color: #70757a;
            background: rgba(128, 128, 128, 0.1);
            white-space: nowrap;
          }
          .step-pill.active {
            color: #1a73e8;
            background: rgba(26, 115, 232, 0.12);
            border: 1px solid rgba(26, 115, 232, 0.35);
          }
          .step-pill.completed {
            color: #1e8e3e;
            background: rgba(30, 142, 62, 0.12);
            border: 1px solid rgba(30, 142, 62, 0.35);
          }
          .step-arrow {
            font-size: 0.7rem;
            opacity: 0.4;
            user-select: none;
          }
          
          .section-card {
            border: 1px solid rgba(128, 128, 128, 0.25);
            border-radius: 12px;
            padding: 1.5rem;
            margin-bottom: 1.5rem;
            background: rgba(255, 255, 255, 0.02);
            backdrop-filter: blur(8px);
          }
          .section-title {
            font-size: 1.25rem;
            font-weight: 800;
            letter-spacing: -0.02em;
            margin-bottom: 0.5rem;
            display: flex;
            align-items: center;
            gap: 0.5rem;
          }
          .section-sub {
            font-size: 0.88rem;
            opacity: 0.75;
            margin-bottom: 1rem;
          }
          
          /* Cards do Funil */
          .wizard-card {
            border: 1px solid rgba(128, 128, 128, 0.25);
            border-radius: 10px;
            padding: 1rem 1.2rem;
            margin-bottom: 0.6rem;
            background: rgba(255, 255, 255, 0.03);
          }
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
          .status-badge {
            font-weight: 700;
            font-size: 1rem;
            margin-right: 0.5rem;
          }
          .status-pending { color: #70757a; }
          .status-running { color: #1a73e8; }
          .status-completed { color: #1e8e3e; }
          
          .summary-number {
            font-size: 2.2rem;
            font-weight: 800;
            color: #1e8e3e;
            letter-spacing: -0.02em;
          }
          .step-connector {
            text-align: center;
            font-size: 1.1rem;
            opacity: 0.35;
            margin: -0.3rem 0;
          }
        </style>
        """,
        unsafe_allow_html=True,
    )


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


def render_stepper_header(current_step: int) -> None:
    steps = [
        "1. CARREGAR",
        "2. CONFIGURAR",
        "3. EXECUTAR",
        "4. VER O PIPELINE",
        "5. ENTENDER O RESULTADO",
        "6. INVESTIGAR UMA ESTRUTURA",
        "7. BAIXAR OS DADOS",
        "8. REPRODUZIR A EXECUÇÃO",
    ]
    pills_html = []
    for idx, name in enumerate(steps, start=1):
        if idx < current_step:
            cls = "step-pill completed"
        elif idx == current_step:
            cls = "step-pill active"
        else:
            cls = "step-pill"
        pills_html.append(f"<span class='{cls}'>{name}</span>")

    arrow = "<span class='step-arrow'>→</span>"
    bar_inner = arrow.join(pills_html)
    st.markdown(f"<div class='stepper-bar'>{bar_inner}</div>", unsafe_allow_html=True)


# --- TELA 1: INPUT (CARREGAR -> CONFIGURAR -> EXECUTAR) ------------------------------


def render_input_screen() -> None:
    render_stepper_header(current_step=1)

    st.markdown(
        """
        <div style='text-align: center; margin-bottom: 2rem;'>
            <h1 style='font-size: 2.2rem; font-weight: 800; margin-bottom: 0.3rem;'>Structure Curation Pipeline</h1>
            <p style='font-size: 1.05rem; opacity: 0.75;'>Curadoria química reprodutível, transparente e auditável</p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    # 1. CARREGAR
    st.markdown("### 1. CARREGAR")
    st.caption("Forneça o lote de estruturas químicas para triagem e curadoria.")

    tab_paste, tab_file = st.tabs(["Cole SMILES", "Upload de Arquivo"])
    raw_bytes: Optional[bytes] = None
    input_name = ""

    with tab_paste:
        text = st.text_area(
            "SMILES Input",
            height=140,
            placeholder="CCO\nCC(=O)O[Na]\nN[C@@H](C)C(=O)O.Cl\nc1ccccc1",
            label_visibility="collapsed",
            help="Cole um código SMILES por linha",
        )
        if text.strip():
            raw_bytes = text.encode("utf-8")
            input_name = "pasted_structures.smi"

    with tab_file:
        uploaded = st.file_uploader(
            "Upload file",
            type=["csv", "tsv", "smi", "smiles", "txt"],
            label_visibility="collapsed",
            help="Envie um arquivo contendo estruturas químicas (.csv, .tsv, .smi)",
        )
        if uploaded is not None:
            raw_bytes = uploaded.getvalue()
            input_name = uploaded.name

    if raw_bytes:
        try:
            preview = preview_input(raw_bytes.decode("utf-8", errors="replace"))
            total_detected = preview.total
        except Exception:
            total_detected = len(raw_bytes.decode("utf-8", errors="ignore").splitlines())

        if total_detected > 0:
            st.success(f"✓ **{input_name}** ({total_detected} moléculas detectadas)")
            st.session_state["raw_input"] = raw_bytes
            st.session_state["input_name"] = input_name

            with st.expander("Pré-visualizar amostras", expanded=False):
                try:
                    lines = [line.strip() for line in raw_bytes.decode("utf-8", errors="replace").splitlines() if line.strip()]
                    sample_rows = [{"#": i + 1, "Estrutura / SMILES": line} for i, line in enumerate(lines[:5])]
                    st.dataframe(sample_rows, width="stretch", hide_index=True)
                except Exception:
                    st.caption("Pré-visualização não disponível.")
        else:
            st.warning("Por favor, forneça pelo menos uma molécula válida.")
            st.session_state["raw_input"] = None
    else:
        st.session_state["raw_input"] = None

    st.divider()

    # 2. CONFIGURAR
    st.markdown("### 2. CONFIGURAR")
    st.caption("Ajuste os parâmetros de corte químico e deduplicação para a execução.")

    cfg_col1, cfg_col2, cfg_col3 = st.columns(3)

    with cfg_col1:
        cfg_mw = st.number_input(
            "Massa Molecular Máx. (Da)",
            min_value=50.0,
            max_value=10000.0,
            value=st.session_state["config"]["max_mw"],
            step=50.0,
            help="Massa molecular máxima permitida na estrutura-mãe.",
        )
    with cfg_col2:
        cfg_ha = st.number_input(
            "Átomos Pesados Máx.",
            min_value=5,
            max_value=1000,
            value=st.session_state["config"]["max_ha"],
            step=5,
            help="Número máximo de átomos pesados permitidos.",
        )
    with cfg_col3:
        st.markdown("<br>", unsafe_allow_html=True)
        cfg_dedup = st.checkbox(
            "Deduplicar por InChIKey",
            value=st.session_state["config"]["deduplicate"],
            help="Remove duplicatas exatas baseando-se no InChIKey completo.",
        )

    st.session_state["config"]["max_mw"] = float(cfg_mw)
    st.session_state["config"]["max_ha"] = int(cfg_ha)
    st.session_state["config"]["deduplicate"] = bool(cfg_dedup)

    st.caption(f"**Política ativa (SHA-256):** `{st.session_state['config']['policy_hash'][:16]}...` (Arquivo: `{st.session_state['config']['policy_path']}`)")

    st.divider()

    # 3. EXECUTAR
    st.markdown("### 3. EXECUTAR")
    st.caption("Inicie o pipeline automatizado de curadoria e padronização.")

    btn_container = st.columns([1, 2, 1])
    with btn_container[1]:
        start_disabled = st.session_state["raw_input"] is None
        if st.button(
            "Start curation",
            type="primary",
            width="stretch",
            disabled=start_disabled,
            help="Clique para iniciar o pipeline de curadoria",
        ):
            run_backend_pipeline()
            if st.session_state.get("fast_mode"):
                st.session_state["ui_state"] = "COMPLETE"
            else:
                st.session_state["ui_state"] = "CURATING"
                st.session_state["curating_step_idx"] = 0
                st.session_state["step_start_time"] = time.time()
            st.rerun()


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


# --- TELA 2: RESULT (VER PIPELINE -> ENTENDER -> INVESTIGAR -> BAIXAR -> REPRODUZIR) --


def get_step_summary(step_id: str, report: Optional[RunReport]) -> str:
    if report is None:
        return "Processando..."

    kpis = report.kpis()
    total_in = kpis["Total Processed"]

    if step_id == "INPUT":
        return f"{total_in} moléculas lidas da entrada"

    if step_id == "STANDARDIZATION":
        std_stage = next((s for s in report.stages if s.name == "STANDARDIZE"), None)
        out_n = std_stage.n_output if std_stage else total_in
        return f"{out_n} estruturas normalizadas"

    if step_id == "PARENT":
        parent_stage = next((s for s in report.stages if s.name == "GET_PARENT"), None)
        out_n = parent_stage.n_output if parent_stage else total_in
        return f"{out_n} estruturas-mãe isoladas"

    if step_id == "VALIDATION":
        val_stage = next((s for s in report.stages if s.name == "VALENCE_GATE"), None)
        valid = val_stage.n_output if val_stage else total_in
        rejected = val_stage.n_excluded if val_stage else 0
        return f"{valid} válidas · {rejected} rejeitadas por valência/sanitização"

    if step_id == "DEDUPLICATION":
        if report.dedup:
            return f"{report.dedup.n_output} únicas · {report.dedup.n_excluded} duplicatas exatas"
        return f"{total_in} únicas · 0 duplicatas"

    if step_id == "ELIGIBILITY":
        el_stage = next((s for s in report.stages if s.name == "ELIGIBILITY"), None)
        eligible = el_stage.n_output if el_stage else total_in
        rejected = el_stage.n_excluded if el_stage else 0
        return f"{eligible} elegíveis · {rejected} rejeitadas por cortes de escopo"

    if step_id == "OUTPUT":
        return f"{kpis['Approved']} estruturas canônicas finais no dataset curado"

    return ""


def render_results_screen() -> None:
    report: Optional[RunReport] = st.session_state.get("report")
    is_complete = st.session_state["ui_state"] == "COMPLETE"
    fast_mode = st.session_state.get("fast_mode", False)

    if not is_complete:
        if fast_mode:
            st.session_state["ui_state"] = "COMPLETE"
            st.rerun()
        else:
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

    render_stepper_header(current_step=4 if not is_complete else 8)

    top_bar_col1, top_bar_col2 = st.columns([3, 1])
    with top_bar_col1:
        st.title("Resultado da Curadoria")
    with top_bar_col2:
        st.markdown("<br>", unsafe_allow_html=True)
        if st.button("Nova análise", width="stretch", help="Reiniciar e carregar novas estruturas"):
            reset_to_input()

    # 4. VER O PIPELINE
    st.markdown("### 4. VER O PIPELINE")
    st.caption("Acompanhe o funil de execução e a passagem das moléculas em cada estágio.")

    cur_idx = st.session_state["curating_step_idx"] if not is_complete else len(WIZARD_STAGES) - 1

    for idx, stage_info in enumerate(WIZARD_STAGES):
        if is_complete or idx < cur_idx:
            status_icon = "✓"
            status_class = "status-completed"
            card_class = "wizard-card-completed"
            status_text = "Concluído"
        elif idx == cur_idx and not is_complete:
            status_icon = "◉"
            status_class = "status-running"
            card_class = "wizard-card-running"
            status_text = "Executando"
        else:
            status_icon = "○"
            status_class = "status-pending"
            card_class = ""
            status_text = "Pendente"

        summary = get_step_summary(stage_info["id"], report) if (idx <= cur_idx or is_complete) else ""

        with st.container():
            col1, col2 = st.columns([5, 1])
            with col1:
                st.markdown(
                    f"""
                    <div class='wizard-card {card_class}'>
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
                    with st.popover("Detalhes", help=f"Ver detalhes de {stage_info['title']}"):
                        st.markdown(f"#### {stage_info['title']}")
                        st.caption(stage_info["description"])
                        st.divider()
                        kpis = report.kpis()
                        if stage_info["id"] == "INPUT":
                            st.write(f"Total lido: {kpis['Total Processed']}")
                        elif stage_info["id"] == "VALIDATION":
                            val_stage = next((s for s in report.stages if s.name == "VALENCE_GATE"), None)
                            if val_stage:
                                st.write(f"Entrada: {val_stage.n_input}")
                                st.write(f"Saída: {val_stage.n_output}")
                                st.write(f"Rejeitadas: {val_stage.n_excluded}")
                        elif stage_info["id"] == "ELIGIBILITY":
                            el_stage = next((s for s in report.stages if s.name == "ELIGIBILITY"), None)
                            if el_stage:
                                st.write(f"Entrada: {el_stage.n_input}")
                                st.write(f"Saída: {el_stage.n_output}")
                                st.write(f"Rejeitadas: {el_stage.n_excluded}")

        if idx < len(WIZARD_STAGES) - 1:
            st.markdown("<div class='step-connector'>↓</div>", unsafe_allow_html=True)

    if not is_complete:
        if not fast_mode:
            time.sleep(0.1)
            st.rerun()
        return

    st.divider()

    # 5. ENTENDER O RESULTADO
    st.markdown("### 5. ENTENDER O RESULTADO")
    st.caption("Métricas consolidadas de desempenho, aprovação e causas de rejeição.")

    if report:
        kpis = report.kpis()
        approved_pct = (kpis["Approved"] / kpis["Total Processed"] * 100.0) if kpis["Total Processed"] > 0 else 0.0

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Total Processado", kpis["Total Processed"])
        m2.metric("Aprovadas / Curadas", f"{kpis['Approved']} ({approved_pct:.1f}%)")
        m3.metric("Rejeitadas", kpis["Rejected"])
        m4.metric("Duplicatas Exatas", kpis["Duplicates"])

        if kpis["Rejected"] > 0:
            st.markdown("##### Motivos de Rejeição Observados")
            rejection_rows = []
            for record in report.rejected:
                code = record.rejection_code.value if record.rejection_code else "DESCONHECIDO"
                stage = record.rejection_stage.value if record.rejection_stage else "N/A"
                detail = record.rejection_detail or ""
                rejection_rows.append({"ID": record.input_id, "Estágio": stage, "Código": code, "Detalhe": detail})
            st.dataframe(rejection_rows, width="stretch", hide_index=True)

    st.divider()

    # 6. INVESTIGAR UMA ESTRUTURA
    st.markdown("### 6. INVESTIGAR UMA ESTRUTURA")
    st.caption("Selecione uma molécula para auditar sua linhagem completa de transformações e renderização 2D.")

    if report and report.records:
        filter_opt = st.radio(
            "Filtrar moléculas",
            ["Todas", "Apenas Aprovadas", "Apenas Rejeitadas"],
            horizontal=True,
            label_visibility="collapsed",
        )

        if filter_opt == "Apenas Aprovadas":
            records_to_show = [r for r in report.records if r.passed]
        elif filter_opt == "Apenas Rejeitadas":
            records_to_show = [r for r in report.records if not r.passed]
        else:
            records_to_show = report.records

        if not records_to_show:
            st.info("Nenhuma estrutura encontrada para o filtro selecionado.")
        else:
            options_map = {
                f"{r.input_id} | Status: {r.status} | SMILES: {r.raw_smiles[:30]}": r
                for r in records_to_show
            }
            selected_key = st.selectbox(
                "Selecione uma molécula para auditar:",
                options=list(options_map.keys()),
            )
            selected_rec = options_map[selected_key]

            st.markdown(f"#### Molécula {selected_rec.input_id}")

            status_color = "#1e8e3e" if selected_rec.passed else "#d93025"
            st.markdown(
                f"**Status final:** <span style='color:{status_color}; font-weight:800;'>{selected_rec.status}</span>",
                unsafe_allow_html=True,
            )

            if not selected_rec.passed:
                st.warning(
                    f"**Motivo do descarte:** Estágio `{selected_rec.rejection_stage}` | Código `{selected_rec.rejection_code}`\n\n"
                    f"**Detalhe:** {selected_rec.rejection_detail or 'Sem detalhe adicional'}"
                )

            img_col1, img_col2 = st.columns(2)
            with img_col1:
                st.markdown("**SMILES de Entrada (Original):**")
                st.code(selected_rec.raw_smiles, language="text")
                img_in = render_mol_image(selected_rec.raw_smiles, legend="Entrada")
                if img_in:
                    st.image(img_in, width=240)

            with img_col2:
                st.markdown("**SMILES Curado (Canônico / Mãe):**")
                if selected_rec.curated_smiles:
                    st.code(selected_rec.curated_smiles, language="text")
                    img_out = render_mol_image(selected_rec.curated_smiles, legend="Curada")
                    if img_out:
                        st.image(img_out, width=240)
                else:
                    st.info("Estrutura não possui SMILES curado (foi rejeitada).")

            st.markdown("##### Linhagem de Transformações por Estágio")
            if selected_rec.transformations:
                trans_rows = [
                    {
                        "Estágio": t.stage.value,
                        "Regra": t.rule,
                        "Antes": t.before_smiles,
                        "Depois": t.after_smiles,
                        "Detalhes": t.detail,
                    }
                    for t in selected_rec.transformations
                ]
                st.dataframe(trans_rows, width="stretch", hide_index=True)
            else:
                st.caption("Nenhuma transformação alterou a conectividade desta molécula.")

            if selected_rec.inchikey:
                st.caption(f"**InChIKey:** `{selected_rec.inchikey}`")
            if selected_rec.removed_fragments:
                st.caption(f"**Fragmentos/Sais removidos:** `{selected_rec.removed_fragments}`")
            if selected_rec.delta_net_charge != 0:
                st.caption(f"**Variação de carga líquida:** `{selected_rec.delta_net_charge}`")

    st.divider()

    # 7. BAIXAR OS DADOS
    st.markdown("### 7. BAIXAR OS DADOS")
    st.caption("Exporte os datasets curados, rejeitados, logs de auditoria e manifesto completo.")

    if report:
        d_col1, d_col2, d_col3, d_col4 = st.columns(4)

        with d_col1:
            st.download_button(
                "Dataset Curado (CSV)",
                to_csv(report.approved, ("input_id", "raw_smiles", "curated_smiles", "inchikey", "status")),
                file_name="curated_structures.csv",
                mime="text/csv",
                width="stretch",
            )
        with d_col2:
            st.download_button(
                "Dataset Rejeitado (CSV)",
                rejected_csv(report),
                file_name="rejected_structures.csv",
                mime="text/csv",
                width="stretch",
            )
        with d_col3:
            st.download_button(
                "Log de Auditoria (CSV)",
                full_csv(report),
                file_name="audit_log.csv",
                mime="text/csv",
                width="stretch",
            )
        with d_col4:
            st.download_button(
                "Pacote Completo (ZIP)",
                reproducibility_package(report),
                file_name=f"curation_run_{report.provenance.run_id[:8]}.zip",
                mime="application/zip",
                width="stretch",
            )

    st.divider()

    # 8. REPRODUZIR A EXECUÇÃO
    st.markdown("### 8. REPRODUZIR A EXECUÇÃO")
    st.caption("Comando exato de terminal e proveniência para garantir reprodutibilidade determinística.")

    if report:
        prov = report.provenance
        cmd = prov.reproduction_command(prov.input_name)

        st.success("✓ **Execução 100% Reproduzível via CLI**")
        st.code(cmd, language="bash")

        with st.expander("Metadados de Proveniência Científica", expanded=False):
            prov_data = [
                {"Metadado": "Run ID", "Valor": prov.run_id},
                {"Metadado": "Política Hash (SHA-256)", "Valor": prov.policy_hash},
                {"Metadado": "Versão RDKit", "Valor": prov.versions.get("rdkit", "N/A")},
                {"Metadado": "Versão Pipeline ChEMBL", "Valor": prov.versions.get("chembl_structure_pipeline", "N/A")},
                {"Metadado": "Commit Git", "Valor": f"{prov.git.commit} (dirty: {prov.git.dirty})"},
            ]
            st.dataframe(prov_data, width="stretch", hide_index=True)


# --- APLICAÇÃO PRINCIPAL -------------------------------------------------------------


def main() -> None:
    st.set_page_config(page_title="Structure Curation Pipeline", layout="centered")
    inject_styles()
    init_session_state()

    ui_state = st.session_state["ui_state"]

    if ui_state == "INPUT":
        render_input_screen()
    else:
        render_results_screen()


if __name__ == "__main__":
    main()
