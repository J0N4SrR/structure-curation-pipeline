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

from curation.components.dag import render_dag
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
          .status-success { color: #1e8e3e; }
          .status-warning { color: #f29900; }
          .status-failed { color: #d93025; }
          .status-skipped { color: #5f6368; font-style: italic; }
          
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
    st.session_state["selected_node"] = None
    st.session_state["ui_state"] = "INPUT"

    st.session_state["raw_input"] = None
    st.session_state["input_name"] = ""
    st.session_state["report"] = None
    st.rerun()


# --- TELA 1: INPUT ------------------------------


def render_input_screen() -> None:

    st.markdown(
        """
        <div style='text-align: center; margin-bottom: 2rem;'>
            <h1 style='font-size: 2.2rem; font-weight: 800; margin-bottom: 0.3rem;'>Structure Curation Pipeline</h1>
            <p style='font-size: 1.05rem; opacity: 0.75;'>Curadoria química reprodutível, transparente e auditável</p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("### Lote de Estruturas")
    st.caption("Forneça as estruturas químicas em formato SMILES para iniciar a execução do pipeline.")

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
        except Exception as e:
            total_detected = 0

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
            st.error(f"**Input inválido:** Não foi possível detectar estruturas válidas em '{input_name}'. Verifique o formato do arquivo.")
            st.session_state["raw_input"] = None
    else:
        st.info("Nenhum input carregado. Cole SMILES ou faça upload de um arquivo para começar.")
        st.session_state["raw_input"] = None

    st.divider()

    st.markdown("### Configuração do Pipeline")
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

    st.markdown("### Executar")
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
            st.session_state["ui_state"] = "COMPLETE"
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



def render_run_details(view: RunViewModel, report: RunReport) -> None:
    """Resumo global, exibido quando nenhum no esta selecionado."""
    st.markdown(f"#### Detalhes da execucao {view.glyph} {view.status_text}")
    st.caption(
        "Nenhum estagio selecionado. Clique em um no do grafo para ver o que "
        "aconteceu nele."
    )
    
    tab_overview, tab_lineage, tab_exports = st.tabs(["Overview & Proveniência", "Structure Lineage", "Exportar Dados"])
    
    with tab_overview:
        a, b, c, d = st.columns(4)
        a.metric("Entrada", view.input_count, help="Total de estruturas químicas processadas nesta execução.")
        b.metric("Aprovadas", view.approved, help="Total de estruturas que passaram por todos os testes e estão no dataset final.")
        c.metric("Rejeitadas", view.rejected, help="Total de estruturas reprovadas e enviadas para o dataset rejeitado.")
        d.metric("Taxa de aprovacao", f"{view.approval_rate:.0%}", help="Porcentagem de estruturas processadas que foram aprovadas.")

        e, f, g = st.columns(3)
        e.metric("Transformadas", view.transformed, help="Número de estruturas cujos átomos, ligações ou cargas foram modificados.")
        f.metric("Duplicatas", view.duplicates, help="Número de duplicatas estruturais exatas encontradas (baseado no InChIKey).")
        g.metric(
            "Duracao",
            f"{view.duration_seconds:.2f} s" if view.duration_seconds else "-",
            help="Tempo real (wall-clock time) para executar todos os estágios do pipeline."
        )
        
        st.markdown("---")
        st.markdown("#### Proveniência da Execução")
        prov = report.provenance
        st.markdown(f"**Run ID:** `{prov.run_id}`")
        st.markdown(f"**Política (SHA-256):** `{prov.policy_hash[:16]}...`")
        
        with st.expander("Ambiente e Versões"):
            st.json(prov.versions)
            st.json(prov.environment)
            if prov.git.commit != "unavailable":
                st.markdown(f"**Git Commit:** `{prov.git.commit}` (Dirty: {prov.git.dirty})")

        with st.expander("Comando de Reprodução (CLI)"):
            st.code(prov.reproduction_command(input_path=prov.input_name), language="bash")
            if not prov.reproducible:
                st.warning("Esta execução possui bloqueadores de reprodução exata:")
                for blocker in prov.reproduction_blockers():
                    st.markdown(f"- {blocker}")

    with tab_lineage:
        render_structure_lineage(report)
        
    with tab_exports:
        render_exports(report)


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


def render_exports(report: RunReport) -> None:
    st.caption("Exporte os datasets curados, rejeitados, logs de auditoria e manifesto completo.")
    d_col1, d_col2 = st.columns(2)

    with d_col1:
        st.download_button(
            "Dataset Curado (CSV)",
            to_csv(report.approved, ("input_id", "raw_smiles", "curated_smiles", "inchikey", "status")),
            file_name="curated_structures.csv",
            mime="text/csv",
            width="stretch",
            help="Arquivo contendo as moléculas finais, normalizadas e deduplicadas.",
        )
        st.download_button(
            "Log de Auditoria (CSV)",
            full_csv(report),
            file_name="audit_log.csv",
            mime="text/csv",
            width="stretch",
            help="Relatório detalhado contendo a trajetória completa e status de cada molécula.",
        )
    with d_col2:
        st.download_button(
            "Dataset Rejeitado (CSV)",
            rejected_csv(report),
            file_name="rejected_structures.csv",
            mime="text/csv",
            width="stretch",
            help="Arquivo contendo as moléculas descartadas e os motivos da rejeição.",
        )
        st.download_button(
            "Pacote Completo (ZIP)",
            reproducibility_package(report),
            file_name=f"curation_run_{report.provenance.run_id[:8]}.zip",
            mime="application/zip",
            width="stretch",
            help="Arquivo ZIP com todos os datasets, relatório consolidado em JSON e o registro das políticas aplicadas.",
        )
    
    st.markdown("---")
    st.markdown("#### Proveniência da Execução")
    prov = report.provenance
    st.markdown(f"**Run ID:** `{prov.run_id}`")
    st.markdown(f"**Política (SHA-256):** `{prov.policy_hash[:16]}...`")
    
    with st.expander("Ambiente e Versões"):
        st.json(prov.versions)
        st.json(prov.environment)
        if prov.git.commit != "unavailable":
            st.markdown(f"**Git Commit:** `{prov.git.commit}` (Dirty: {prov.git.dirty})")

    with st.expander("Comando de Reprodução (CLI)"):
        st.code(prov.reproduction_command(input_path=prov.input_name), language="bash")
        if not prov.reproducible:
            st.warning("Esta execução possui bloqueadores de reprodução exata:")
            for blocker in prov.reproduction_blockers():
                st.markdown(f"- {blocker}")


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
        "Duracao",
        f"{node.duration_seconds * 1000:.0f} ms" if node.duration_seconds else "-",
        help="Tempo computacional gasto na execução exclusiva deste nó."
    )

    if node.parameters:
        with st.expander("Parametros que governam este estagio"):
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


def render_pipeline_graph(report: RunReport) -> None:
    """Grafo interativo e painel contextual.

    O clique em um no volta do navegador pelo componente e vira
    ``selected_node``; o painel abaixo alterna entre detalhe da execucao e
    detalhe do estagio.
    """
    view = build_run_view(report)
    st.session_state["run_view"] = view

    left, right = st.columns([1, 2])
    with left:
        escolhido = render_dag(
            view.graph, selected=st.session_state.get("selected_node"), height=560
        )
        if escolhido != st.session_state.get("selected_node"):
            st.session_state["selected_node"] = escolhido
            st.rerun()
        if st.session_state.get("selected_node"):
            if st.button("Voltar ao resumo da execucao", width="stretch"):
                st.session_state["selected_node"] = None
                st.rerun()

    with right:
        selecionado = st.session_state.get("selected_node")
        node = view.graph.node(selecionado) if selecionado else None
        if node is None:
            render_run_details(view, report)
        else:
            render_node_details(node, report)


def render_results_screen() -> None:
    report: Optional[RunReport] = st.session_state.get("report")
    
    if st.session_state["ui_state"] != "COMPLETE":
        st.session_state["ui_state"] = "COMPLETE"
        st.rerun()

    top_bar_col1, top_bar_col2 = st.columns([3, 1])
    with top_bar_col1:
        st.title("Resultado da Curadoria")
    with top_bar_col2:
        st.markdown("<br>", unsafe_allow_html=True)
        if st.button("Nova análise", width="stretch", help="Reiniciar e carregar novas estruturas"):
            reset_to_input()

    if report:
        render_pipeline_graph(report)


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
