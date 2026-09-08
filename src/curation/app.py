"""Interface Streamlit - Workspace Científico de Curadoria Estrutural.

Cards reativos de pipeline com progressão sequencial em tempo real,
popovers discretos para parâmetros e downloads, e inspeção molecular detalhada.
"""

from __future__ import annotations

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
from curation.models import Stage
from curation.reporting import (
    DEDUP_STAGE,
    RunReport,
    StageStatus,
    full_csv,
    rejected_csv,
    reproducibility_package,
    to_csv,
)

RDLogger.DisableLog("rdApp.*")

DEFAULT_DECISIONS = Path("docs/decisions.md")

PIPELINE_STAGES = [
    ("PARSE", "Leitura e Validação Sintática"),
    ("STANDARDIZE", "Padronização e Neutralização (ChEMBL)"),
    ("GET_PARENT", "Estrutura-Mãe e Desagregação de Sais"),
    ("VALENCE_GATE", "Portão Estrito de Valência"),
    ("ELIGIBILITY", "Filtros de Escopo Químico (MW / Átomos Pesados)"),
    ("CANONICALIZE", "Identidade Canônica e InChIKey"),
    (DEDUP_STAGE, "Deduplicação Estrutural"),
]


def render_mol_image(smiles: Optional[str], legend: str = "", size: tuple[int, int] = (250, 250)):
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
            max-width: 1000px;
            padding-top: 4.5rem !important;
            padding-bottom: 3rem;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
          }
          
          /* Cards Científicos do Pipeline */
          .pipeline-card {
            background: #0F172A;
            border: 1px solid #334155;
            border-radius: 8px;
            padding: 14px 18px;
            margin-bottom: 10px;
            font-family: ui-monospace, monospace;
            font-size: 12px;
            display: flex;
            justify-content: space-between;
            align-items: center;
          }
          .card-running {
            border: 2px solid #0284C7 !important;
            box-shadow: 0 0 16px rgba(2, 132, 199, 0.45);
            background: #0F172A;
          }
          .card-success {
            border-left: 6px solid #10B981 !important;
          }
          .card-warning {
            border-left: 6px solid #F59E0B !important;
          }
          .badge-running {
            background: #0284C7;
            color: #FFFFFF;
            padding: 4px 8px;
            border-radius: 4px;
            font-weight: 700;
          }
          .badge-success {
            background: #065F46;
            color: #34D399;
            padding: 4px 8px;
            border-radius: 4px;
            font-weight: 700;
          }
          .badge-warning {
            background: #78350F;
            color: #FBBF24;
            padding: 4px 8px;
            border-radius: 4px;
            font-weight: 700;
          }
        </style>
        """,
        unsafe_allow_html=True,
    )


def init_session_state() -> None:
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
    if "completed_cards" not in st.session_state:
        st.session_state["completed_cards"] = []
    if "inspect_stage" not in st.session_state:
        st.session_state["inspect_stage"] = None


def render_card_html(stage_id: str, label: str, status: str, n_in: any, n_out: any, n_rej: any) -> str:
    if status == "RUNNING":
        css_class = "pipeline-card card-running"
        badge = '<span class="badge-running">RUNNING</span>'
    elif status == "WARNING":
        css_class = "pipeline-card card-warning"
        badge = '<span class="badge-warning">WARNING</span>'
    elif status == "SUCCESS":
        css_class = "pipeline-card card-success"
        badge = '<span class="badge-success">COMPLETED</span>'
    else:
        css_class = "pipeline-card"
        badge = '<span style="color: #94A3B8; font-weight: 600;">PENDING</span>'

    # Alto contraste: texto cinza claro (#E2E8F0) e valores destacados em ciano (#38BDF8)
    metrics = (
        f'<span style="color: #E2E8F0; font-size: 11px;">'
        f'IN: <b style="color: #FFFFFF;">{n_in}</b> &nbsp;|&nbsp; '
        f'OUT: <b style="color: #FFFFFF;">{n_out}</b> &nbsp;|&nbsp; '
        f'REJ: <b style="color: #F87171;">{n_rej}</b>'
        f'</span>'
    )

    return f"""
    <div class="{css_class}">
        <div>
            <span style="font-weight: 800; color: #FFFFFF; margin-right: 10px; font-size: 13px;">[{stage_id}]</span>
            <span style="color: #CBD5E1; font-weight: 500;">{label}</span>
        </div>
        <div style="display: flex; align-items: center; gap: 16px;">
            {metrics}
            {badge}
        </div>
    </div>
    """


def execute_pipeline_live(cards_placeholder):
    raw = st.session_state["raw_input"]
    cfg = st.session_state["config"]
    name = st.session_state["input_name"]

    pipeline = CurationPipeline(
        policy_hash=cfg["policy_hash"],
        criteria=EligibilityCriteria(
            max_molecular_weight=cfg["max_mw"],
            max_heavy_atoms=cfg["max_ha"]
        ),
        deduplicate=cfg["deduplicate"]
    )

    # 1. Executa a curadoria real no backend
    report = pipeline.run_report(
        raw.decode("utf-8", errors="replace"),
        parameters=cfg,
        input_bytes=raw,
        input_name=name,
        policy_path=cfg["policy_path"],
    )

    stage_map = {s.name: s for s in report.stages}
    if report.dedup:
        stage_map[DEDUP_STAGE] = report.dedup

    rendered_cards = []

    # 2. Faz os cards surgirem sequencialmente na tela
    for stage_id, label in PIPELINE_STAGES:
        # Mostra o nó atual como RUNNING
        running_html = render_card_html(stage_id, label, "RUNNING", "-", "-", "-")
        cards_placeholder.markdown("".join(rendered_cards + [running_html]), unsafe_allow_html=True)
        time.sleep(0.55)

        # Atualiza métricas reais e fixa como finalizado
        rep = stage_map.get(stage_id)
        n_in = rep.n_input if rep else "-"
        n_out = rep.n_output if rep else "-"
        n_rej = rep.n_excluded if rep else 0
        status = "WARNING" if (rep and rep.has_exclusions) else "SUCCESS"

        done_html = render_card_html(stage_id, label, status, n_in, n_out, n_rej)
        rendered_cards.append(done_html)
        cards_placeholder.markdown("".join(rendered_cards), unsafe_allow_html=True)
        time.sleep(0.2)

    st.session_state["completed_cards"] = rendered_cards
    st.session_state["report"] = report


def render_workbench():
    # --- 1. BARRA SUPERIOR DE AÇÕES E POPOVERS ---
    col_title, col_cfg, col_dl, col_run = st.columns([5, 1, 1, 2])

    with col_title:
        st.markdown("### ⚗️ Structure Curation Pipeline")
        st.caption(f"Política ativa: `{st.session_state['config']['policy_hash'][:16]}...`")

    # POPOVER: Configurações do Pipeline (ESCONDIDO)
    with col_cfg:
        with st.popover("⚙️", help="Configuração do Pipeline"):
            st.markdown("##### ⚙️ Parâmetros do Estudo")
            st.session_state["config"]["max_mw"] = st.number_input(
                "Massa Molecular Máx. (Da)", value=st.session_state["config"]["max_mw"], step=50.0
            )
            st.session_state["config"]["max_ha"] = st.number_input(
                "Átomos Pesados Máx.", value=st.session_state["config"]["max_ha"], step=5
            )
            st.session_state["config"]["deduplicate"] = st.checkbox(
                "Deduplicação InChIKey", value=st.session_state["config"]["deduplicate"]
            )

    # POPOVER: Opções de Download (ESCONDIDO)
    with col_dl:
        report: Optional[RunReport] = st.session_state.get("report")
        with st.popover("📥", help="Opções de Download e Exportação"):
            st.markdown("##### 📥 Exportar Resultados")
            
            if report and report.records:
                # 1. Escolha do Dataset Base
                target_dataset = st.radio(
                    "Base:",
                    ["Aprovadas", "Rejeitadas", "Auditoria Completa"],
                    horizontal=True,
                    label_visibility="collapsed"
                )

                # Carrega o DataFrame correspondente
                import pandas as pd
                if target_dataset == "Aprovadas":
                    base_records = report.approved
                    default_min = ["input_id", "curated_smiles", "inchikey"]
                elif target_dataset == "Rejeitadas":
                    base_records = report.rejected
                    default_min = ["input_id", "raw_smiles", "rejection_code", "rejection_stage"]
                else:
                    base_records = report.records
                    default_min = ["input_id", "curated_smiles", "status", "rejection_code"]

                # Converte para DataFrame em memória
                df_export = pd.DataFrame([r.model_dump(mode="json") for r in base_records])

                # 2. Seletor de Escopo de Colunas
                scope = st.segmented_control(
                    "Colunas:",
                    options=["Mínimo", "Completo", "Personalizado"],
                    default="Mínimo"
                )

                if scope == "Mínimo":
                    cols = [c for c in default_min if c in df_export.columns]
                elif scope == "Completo":
                    cols = df_export.columns.tolist()
                else:
                    cols = st.multiselect(
                        "Selecione as colunas:",
                        options=df_export.columns.tolist(),
                        default=[c for c in default_min if c in df_export.columns]
                    )

                df_final = df_export[cols] if cols else df_export

                # 3. Opção de Edição Rápida
                edit_mode = st.toggle("✏️ Editar dados antes de baixar", value=False)
                if edit_mode:
                    df_final = st.data_editor(
                        df_final,
                        num_rows="dynamic",
                        height=250,
                        use_container_width=True
                    )

                # 4. Download do CSV customizado
                csv_bytes = df_final.to_csv(index=False).encode("utf-8")
                st.download_button(
                    label=f"⬇️ Baixar CSV ({len(df_final)} mols)",
                    data=csv_bytes,
                    file_name=f"{target_dataset.lower()}_{scope.lower()}.csv",
                    mime="text/csv",
                    width="stretch"
                )

                st.divider()
                # Pacote completo de reprodutibilidade mantido para conformidade científica
                st.download_button(
                    "📦 Pacote de Reprodutibilidade (ZIP)",
                    reproducibility_package(report),
                    file_name=f"curation_pkg_{report.provenance.run_id[:8]}.zip",
                    mime="application/zip",
                    width="stretch"
                )
            else:
                st.caption("Execute o pipeline para habilitar as opções de exportação.")

    with col_run:
        can_run = st.session_state.get("raw_input") is not None
        start_btn = st.button("▶ Executar Pipeline", type="primary", disabled=not can_run, width="stretch")

    # --- 2. ÁREA DE ENTRADA (UPLOAD / SMILES) ---
    with st.expander("📂 Ingestão de Moléculas (SMILES)", expanded=(st.session_state.get("raw_input") is None)):
        tab_paste, tab_file = st.tabs(["Cole SMILES", "Upload de Arquivo"])
        raw_bytes: Optional[bytes] = None
        input_name = ""

        with tab_paste:
            with st.form("form_smiles_input", clear_on_submit=False):
                text = st.text_area(
                    "SMILES Input",
                    height=130,
                    placeholder="CCO\nCC(=O)O[Na]\nN[C@@H](C)C(=O)O.Cl\nc1ccccc1",
                    label_visibility="collapsed",
                    help="Insira um SMILES por linha",
                )
                btn_load = st.form_submit_button(
                    "📥 Carregar Moléculas",
                    type="secondary",
                    use_container_width=True,
                )

            if btn_load and text.strip():
                raw_bytes = text.encode("utf-8")
                input_name = "pasted_structures.smi"

        with tab_file:
            uploaded = st.file_uploader(
                "Arquivo de estruturas",
                type=["csv", "tsv", "smi", "smiles", "txt"],
                label_visibility="collapsed"
            )
            if uploaded is not None:
                raw_bytes = uploaded.getvalue()
                input_name = uploaded.name

        if raw_bytes:
            try:
                preview = preview_input(raw_bytes.decode("utf-8", errors="replace"))
                if preview.total > 0:
                    st.success(f"✓ **{input_name}** ({preview.total} moléculas detectadas)")
                    st.session_state["raw_input"] = raw_bytes
                    st.session_state["input_name"] = input_name
                else:
                    st.error("Nenhuma estrutura válida detectada.")
                    st.session_state["raw_input"] = None
            except Exception:
                st.session_state["raw_input"] = None

    st.markdown("---")

    # --- 3. CARDS DO PIPELINE (PROGRESSÃO EM TEMPO REAL) ---
    cards_slot = st.empty()

    if start_btn:
        execute_pipeline_live(cards_slot)
    elif st.session_state["completed_cards"]:
        cards_slot.markdown("".join(st.session_state["completed_cards"]), unsafe_allow_html=True)
    else:
        # Estado Inicial: Cards cinzas aguardando início
        initial_html = [
            render_card_html(s_id, s_label, "PENDING", "-", "-", "-")
            for s_id, s_label in PIPELINE_STAGES
        ]
        cards_slot.markdown("".join(initial_html), unsafe_allow_html=True)

    # --- 4. DETALHES ANALÍTICOS (APÓS A CONCLUSÃO) ---
    if st.session_state.get("report"):
        st.markdown("---")
        render_audit_details(st.session_state["report"])


def render_audit_details(report: RunReport):
    tab_overview, tab_lineage, tab_reproducibility = st.tabs([
        "Visão Geral", "Linhagem Molecular (2D)", "Manifesto de Reprodutibilidade"
    ])

    with tab_overview:
        kpis = report.kpis()
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Total Processado", kpis["Total Processed"])
        col2.metric("Aprovadas", kpis["Approved"])
        col3.metric("Rejeitadas", kpis["Rejected"])
        col4.metric("Avisos / Transformadas", kpis["Warnings"])

        st.markdown("##### Exclusões por Motivo")
        exclusions = report.exclusion_groups()
        if exclusions:
            st.dataframe(
                [{"Motivo": g.reason, "Estágio": g.stage, "Quantidade": g.count} for g in exclusions],
                width="stretch",
                hide_index=True
            )
        else:
            st.success("Nenhuma estrutura foi excluída neste lote.")

    with tab_lineage:
        records = report.records
        filter_type = st.radio(
            "Filtrar:", ["Todas", "Apenas Aprovadas", "Apenas Rejeitadas"],
            horizontal=True, label_visibility="collapsed"
        )
        if filter_type == "Apenas Aprovadas":
            records = [r for r in records if r.passed]
        elif filter_type == "Apenas Rejeitadas":
            records = [r for r in records if not r.passed]

        if not records:
            st.info("Nenhuma estrutura encontrada para o filtro.")
            return

        rec_dict = {f"{r.input_id} | {r.status} | {r.raw_smiles[:30]}": r for r in records}
        selected_key = st.selectbox("Selecione a molécula para auditar:", list(rec_dict.keys()))
        rec = rec_dict[selected_key]

        col_in, col_out = st.columns(2)
        with col_in:
            st.markdown("**Entrada Bruta:**")
            st.code(rec.raw_smiles, language="text")
            img_in = render_mol_image(rec.raw_smiles, legend="Entrada")
            if img_in:
                st.image(img_in, width=220)

        with col_out:
            st.markdown("**Estrutura-Mãe (Curada):**")
            if rec.curated_smiles:
                st.code(rec.curated_smiles, language="text")
                img_out = render_mol_image(rec.curated_smiles, legend="Curada")
                if img_out:
                    st.image(img_out, width=220)
            else:
                st.caption("Estrutura rejeitada (sem SMILES curado).")

        if rec.transformations:
            st.markdown("##### Eventos de Transformação Registrados")
            for t in rec.transformations:
                st.info(f"**[{t.stage.value}]** {t.rule} — {t.detail}")

    with tab_reproducibility:
        prov = report.provenance
        st.code(prov.reproduction_command(input_path=prov.input_name), language="bash")
        c1, c2 = st.columns(2)
        with c1:
            st.json(prov.environment)
        with c2:
            st.json(prov.versions)


def main() -> None:
    st.set_page_config(page_title="Structure Curation Pipeline", layout="centered")
    inject_styles()
    init_session_state()
    render_workbench()


if __name__ == "__main__":
    main()
