"""Interface Streamlit (Wizard UI) para o Structure Curation Pipeline.

Organizada como um Workspace contínuo, simulando a execução ao vivo do pipeline
(estilo Kubeflow / Nextflow Tower).
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

from streamlit_flow import streamlit_flow
from curation.components.dag import build_live_dag_state, ORDERED_STAGES
from curation.filters import EligibilityCriteria
from curation.io import compute_policy_hash, preview_input
from curation.pipeline import CurationPipeline
from curation.viewmodel import NodeContract, RunViewModel, build_run_view
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
    if "stage_live_states" not in st.session_state:
        st.session_state["stage_live_states"] = {
            s[0]: {"status": StageStatus.PENDING, "in": "-", "out": "-", "rej": "-"}
            for s in ORDERED_STAGES
        }
    if "selected_node" not in st.session_state:
        st.session_state["selected_node"] = None
    if "run_view" not in st.session_state:
        st.session_state["run_view"] = None

def execute_pipeline_live(dag_container):
    raw = st.session_state["raw_input"]
    cfg = st.session_state["config"]
    
    pipeline = CurationPipeline(
        policy_hash=cfg["policy_hash"],
        criteria=EligibilityCriteria(max_molecular_weight=cfg["max_mw"], max_heavy_atoms=cfg["max_ha"]),
        deduplicate=cfg["deduplicate"]
    )

    # Inicializa todos como PENDING
    states = {s[0]: {"status": StageStatus.PENDING, "in": "-", "out": "-", "rej": "-"} for s in ORDERED_STAGES}

    # Executa a curadoria real para obter as métricas do motor
    report = pipeline.run_report(
        raw.decode("utf-8", errors="replace"),
        parameters=cfg,
        input_bytes=raw,
        input_name=st.session_state["input_name"],
        policy_path=cfg["policy_path"],
    )
    
    stage_reports = {s.name: s for s in report.stages}
    if report.dedup:
        stage_reports["DEDUPLICATION"] = report.dedup

    # Animação passo a passo dos nós
    for stage_id, label in ORDERED_STAGES:
        # Marca estágio como RUNNING (Azul, pulsante, aresta animada)
        states[stage_id]["status"] = StageStatus.RUNNING
        with dag_container:
            streamlit_flow("curation_live_dag", build_live_dag_state(states), height=360, fit_view=False)
        
        # Pausa cadenciada para percepção do processamento científico
        time.sleep(0.6)

        # Atualiza métricas reais e marca como SUCCESS / WARNING
        rep = stage_reports.get(stage_id)
        if rep:
            status = StageStatus.WARNING if rep.has_exclusions else StageStatus.SUCCESS
            states[stage_id] = {
                "status": status,
                "in": rep.n_input,
                "out": rep.n_output,
                "rej": rep.n_excluded
            }
        else:
            states[stage_id]["status"] = StageStatus.SUCCESS

        with dag_container:
            streamlit_flow("curation_live_dag", build_live_dag_state(states), height=360, fit_view=False)

    st.session_state["stage_live_states"] = states
    st.session_state["report"] = report
    st.session_state["run_view"] = build_run_view(report)
    st.session_state["selected_node"] = None
    st.rerun()

def render_workbench():
    top_bar = st.container()

    # 2. Área de Entrada de Dados (Expansível / Compacta)
    with st.expander("📂 Ingestão de Moléculas (SMILES) ou Upload", expanded=(st.session_state.get("raw_input") is None)):
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
                "Enviar arquivo",
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
            except Exception as e:
                total_detected = 0

            if total_detected > 0:
                st.success(f"✓ **{input_name}** ({total_detected} moléculas detectadas)")
                st.session_state["raw_input"] = raw_bytes
                st.session_state["input_name"] = input_name
            else:
                st.error(f"**Input inválido:** Não foi possível detectar estruturas válidas em '{input_name}'. Verifique o formato do arquivo.")
                st.session_state["raw_input"] = None
        elif st.session_state.get("raw_input") and not text.strip() and not uploaded:
             st.success(f"✓ **{st.session_state['input_name']}** pronto para execução.")

    with top_bar:
        # 1. Barra de Ações Superior (Estilo Kubeflow Workspace)
        col_info, col_cfg, col_dl, col_action = st.columns([5, 1, 1, 2])
        
        with col_info:
            st.markdown("### ⚗️ Workspace de Curadoria Estrutural")
            st.caption(f"Política Ativa: `{st.session_state['config']['policy_hash'][:16]}...`")

        # Ícone de Engrenagem (Configuração do Pipeline)
        with col_cfg:
            with st.popover("⚙️", help="Configurações e Parâmetros"):
                st.markdown("##### ⚙️ Parâmetros do Estudo")
                st.session_state["config"]["max_mw"] = st.number_input(
                    "Massa Molecular Máx.", value=st.session_state["config"]["max_mw"]
                )
                st.session_state["config"]["max_ha"] = st.number_input(
                    "Átomos Pesados Máx.", value=st.session_state["config"]["max_ha"]
                )
                st.session_state["config"]["deduplicate"] = st.checkbox(
                    "Deduplicação InChIKey", value=st.session_state["config"]["deduplicate"]
                )

        # Ícone de Download (Artefatos da Execução)
        with col_dl:
            report = st.session_state.get("report")
            with st.popover("📥", help="Opções de Download"):
                st.markdown("##### 📥 Exportação de Artefatos")
                if report:
                    st.download_button(
                        "📄 Dataset Curado (CSV)",
                        to_csv(report.approved, ("input_id", "raw_smiles", "curated_smiles", "inchikey", "status")),
                        file_name="curated_structures.csv",
                        mime="text/csv",
                        width="stretch"
                    )
                    st.download_button(
                        "📄 Dataset Rejeitado (CSV)",
                        rejected_csv(report),
                        file_name="rejected_structures.csv",
                        mime="text/csv",
                        width="stretch"
                    )
                    st.download_button(
                        "📊 Log de Auditoria (CSV)",
                        full_csv(report),
                        file_name="audit_log.csv",
                        mime="text/csv",
                        width="stretch"
                    )
                    st.download_button(
                        "📦 Pacote ZIP",
                        reproducibility_package(report),
                        file_name=f"curation_pkg_{report.provenance.run_id[:8]}.zip",
                        mime="application/zip",
                        width="stretch"
                    )
                else:
                    st.caption("Nenhum dado disponível. Execute o pipeline primeiro.")

        with col_action:
            can_run = st.session_state.get("raw_input") is not None
            run_btn = st.button("▶ Executar Pipeline", type="primary", disabled=not can_run, width="stretch")

    # 3. Canvas do DAG (O Grafo Vivo)
    st.markdown("---")
    dag_placeholder = st.empty()

    # Se clicou em executar: aciona o ciclo reativo passo a passo
    if run_btn:
        execute_pipeline_live(dag_placeholder)

    # Renderização estática do frame atual / interatividade pós-execução
    flow_state = build_live_dag_state(
        st.session_state["stage_live_states"],
        selected_node=st.session_state.get("selected_node")
    )
    with dag_placeholder:
        updated = streamlit_flow("curation_live_dag", flow_state, height=360, fit_view=False)
        # Handle selection logic only if execution is finished (report exists)
        if updated and st.session_state.get("report") is not None:
            if updated.selected_id != st.session_state.get("selected_node"):
                st.session_state["selected_node"] = updated.selected_id
                st.rerun()

    # 4. Detalhes contextuais do nó clicado (Structure Lineage e Métricas)
    if st.session_state.get("report"):
        selecionado = st.session_state.get("selected_node")
        
        # Botão para limpar a seleção
        if selecionado:
            if st.button("Voltar ao resumo da execução", width="stretch"):
                st.session_state["selected_node"] = None
                st.rerun()
                
        if selecionado:
            view = st.session_state["run_view"]
            node = view.graph.node(selecionado) if selecionado else None
            if node:
                render_node_details(node, st.session_state["report"])
        else:
            render_run_details(st.session_state["run_view"], st.session_state["report"])


def render_run_details(view: RunViewModel, report: RunReport) -> None:
    """Resumo global, exibido quando nenhum no esta selecionado."""
    st.markdown(f"#### Detalhes da Execução {view.glyph} {view.status_text}")
    st.caption(
        "Nenhum estágio selecionado. Clique em um nó do grafo para ver o que "
        "aconteceu nele."
    )
    
    tab_overview, tab_lineage = st.tabs(["Visão Geral", "Linhagem de Estruturas"])
    
    with tab_overview:
        a, b, c, d = st.columns(4)
        a.metric("Entrada", view.input_count, help="Total de estruturas químicas processadas nesta execução.")
        b.metric("Aprovadas", view.approved, help="Total de estruturas que passaram por todos os testes e estão no dataset final.")
        c.metric("Rejeitadas", view.rejected, help="Total de estruturas reprovadas e enviadas para o dataset rejeitado.")
        d.metric("Taxa de Aprovação", f"{view.approval_rate:.0%}", help="Porcentagem de estruturas processadas que foram aprovadas.")

        e, f, g = st.columns(3)
        e.metric("Transformadas", view.transformed, help="Número de estruturas cujos átomos, ligações ou cargas foram modificados.")
        f.metric("Duplicatas", view.duplicates, help="Número de duplicatas estruturais exatas encontradas (baseado no InChIKey).")
        g.metric(
            "Duração",
            f"{view.duration_seconds:.2f} s" if view.duration_seconds else "-",
            help="Tempo real (wall-clock time) para executar todos os estágios do pipeline."
        )

        st.markdown("---")
        st.markdown("#### Reprodutibilidade e Manifesto")
        prov = report.provenance
        
        with st.expander("📦 Comando de Reprodução (CLI) e Ambiente", expanded=True):
            st.code(prov.reproduction_command(input_path=prov.input_name), language="bash")
            
            env_col1, env_col2 = st.columns(2)
            with env_col1:
                st.markdown("**Sistema**")
                st.json(prov.environment)
            with env_col2:
                st.markdown("**Bibliotecas**")
                st.json(prov.versions)
                
            if not prov.reproducible:
                st.warning("Esta execução possui bloqueadores de reprodução exata:")
                for blocker in prov.reproduction_blockers():
                    st.markdown(f"- {blocker}")


    with tab_lineage:
        render_structure_lineage(report)


def render_structure_lineage(report: RunReport) -> None:
    st.caption("Selecione uma molécula para auditar sua linhagem completa de transformações e renderização 2D.")

    if not report.records:
        st.info("Nenhuma estrutura disponível.")
        return

    filter_opt = st.radio(
        "Filtrar moléculas",
        ["Todas", "Apenas Aprovadas", "Apenas Rejeitadas"],
        horizontal=True,
        label_visibility="collapsed",
        help="Filtra a lista de moléculas com base no status final da execução.",
    )

    if filter_opt == "Apenas Aprovadas":
        records_to_show = [r for r in report.records if r.passed]
    elif filter_opt == "Apenas Rejeitadas":
        records_to_show = [r for r in report.records if not r.passed]
    else:
        records_to_show = report.records

    if not records_to_show:
        st.info("Nenhuma estrutura encontrada para o filtro selecionado.")
        return

    options_map = {
        f"{r.input_id} | Status: {r.status} | SMILES: {r.raw_smiles[:30]}": r
        for r in records_to_show
    }
    selected_key = st.selectbox(
        "Selecione uma molécula para auditar:",
        options=list(options_map.keys()),
        help="Pesquise ou selecione uma molécula pelo SMILES original ou identificador.",
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

    st.markdown("##### 🧬 Linhagem de Transformações (Timeline)")
    if selected_rec.transformations:
        for idx, t in enumerate(selected_rec.transformations):
            with st.expander(f"Step {idx+1}: {t.stage.value} — {t.rule}", expanded=True):
                if t.detail:
                    st.info(t.detail)
                c1, c2 = st.columns(2)
                with c1:
                    st.markdown("**Antes:**")
                    st.code(t.before_smiles, language="text")
                    img_b = render_mol_image(t.before_smiles, size=(200, 200))
                    if img_b:
                        st.image(img_b, width=200)
                with c2:
                    st.markdown("**Depois:**")
                    st.code(t.after_smiles, language="text")
                    img_a = render_mol_image(t.after_smiles, size=(200, 200))
                    if img_a:
                        st.image(img_a, width=200)
    else:
        st.caption("Nenhuma transformação alterou a conectividade desta molécula.")

    if selected_rec.inchikey:
        st.caption(f"**InChIKey:** `{selected_rec.inchikey}`")
    if selected_rec.removed_fragments:
        st.caption(f"**Fragmentos/Sais removidos:** `{selected_rec.removed_fragments}`")
    if selected_rec.delta_net_charge != 0:
        st.caption(f"**Variação de carga líquida:** `{selected_rec.delta_net_charge}`")


def render_node_details(node: NodeContract, report: RunReport) -> None:
    """Detalhe contextual do no selecionado."""
    st.markdown(f"#### {node.label}  {node.glyph} {node.status_text}")
    st.caption(node.help_text)

    a, b, c, d, e = st.columns(5)
    a.metric("Entrada", node.input_count, help="Quantidade de estruturas recebidas na entrada do estágio.")
    b.metric("Resultado", node.output_count, help="Quantidade de estruturas que avançaram para o próximo estágio.")
    c.metric("Transformadas", node.transformed_count, help="Estruturas quimicamente modificadas neste estágio.")
    d.metric("Removidas", node.rejected_count, help="Estruturas removidas e enviadas para o dataset de rejeitadas.")
    e.metric(
        "Duração",
        f"{node.duration_seconds * 1000:.0f} ms" if node.duration_seconds else "-",
        help="Tempo computacional gasto na execução exclusiva deste nó."
    )

    if node.parameters:
        with st.expander("Parâmetros que governam este estágio"):
            st.json(node.parameters)

    if node.exclusions:
        st.markdown("##### O que foi removido aqui")
        st.dataframe(
            [
                {"Motivo": group.reason, "Quantidade": group.count}
                for group in node.exclusions
            ],
            hide_index=True,
            width="stretch",
        )

    for artifact in node.artifacts:
        registros = report.records_by_id(artifact.record_ids)
        st.download_button(
            f"{artifact.name} ({artifact.count})",
            to_csv(
                registros,
                ("input_id", "raw_smiles", "rejection_code", "rejection_detail"),
            ),
            file_name=artifact.name,
            mime="text/csv",
            help=artifact.description,
            key=f"artifact_{node.id}_{artifact.name}",
        )

    if not node.exclusions and not node.artifacts:
        st.success("Nenhuma estrutura foi removida neste estagio.")


def main() -> None:
    st.set_page_config(page_title="Structure Curation Pipeline", layout="centered")
    inject_styles()
    init_session_state()
    render_workbench()


if __name__ == "__main__":
    main()
