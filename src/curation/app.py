"""Interface Streamlit: cliente fino sobre a camada de engenharia.

Apresentação apenas. Nenhuma regra química, critério de aceitação, métrica ou
linhagem é calculada aqui — tudo vem de :mod:`curation.pipeline`,
:mod:`curation.reporting`, :mod:`curation.provenance` e :mod:`curation.io`.

Quando um dado não existe na engenharia, a interface diz que não existe. Não há
progresso simulado, estado estimado nem transformação inferida.

Organização em três níveis de detalhe: o essencial fica visível, a explicação abre
ao clicar numa etapa, e os metadados técnicos ficam atrás de "detalhes".

Execução::

    streamlit run streamlit_app.py
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import streamlit as st
from rdkit import Chem, RDLogger

try:
    from rdkit.Chem import Draw

    DRAWING_AVAILABLE = True
    DRAWING_ERROR = ""
except ImportError as _error:  # pragma: no cover - depende do ambiente de deploy
    # rdMolDraw2D linka contra libXrender/libX11/libXext do sistema, ausentes em
    # containers slim. A renderização é uma funcionalidade entre várias: sem ela o
    # app degrada para SMILES em texto, em vez de derrubar a auditoria inteira.
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
    record_row,
    rejected_csv,
    reproducibility_package,
    run_manifest,
    structure_lineage,
    to_csv,
)

RDLogger.DisableLog("rdApp.*")

DEFAULT_DECISIONS = Path("docs/decisions.md")


# --- Vocabulário de apresentação ------------------------------------------------
#
# Rótulos legíveis para os estágios reais do motor. É tradução de nome, não
# reagrupamento: cada rótulo corresponde a um estágio que o pipeline de fato
# executa, para que o funil exibido case com o código.

STAGE_LABELS: dict[str, str] = {
    "PARSE": "Leitura",
    "STANDARDIZE": "Padronização",
    "GET_PARENT": "Remoção de sais",
    "VALENCE_GATE": "Validação",
    "ELIGIBILITY": "Elegibilidade",
    "CANONICALIZE": "Canonicalização",
    DEDUP_STAGE: "Deduplicação",
}

STAGE_HELP: dict[str, str] = {
    "PARSE": "Lê o texto da estrutura e monta o grafo químico. Aqui só são "
    "recusadas estruturas que o programa não consegue interpretar.",
    "STANDARDIZE": "Aplica as regras de padronização do ChEMBL: normaliza grupos "
    "funcionais, separa metais ligados de forma inadequada e ajusta cargas.",
    "GET_PARENT": "Remove contra-íons e solventes, isolando a estrutura principal. "
    "Misturas legítimas de dois princípios ativos são preservadas.",
    "VALENCE_GATE": "Verifica se a estrutura resultante é quimicamente válida. Vem "
    "depois da padronização de propósito, para dar ao motor a chance de corrigir.",
    "ELIGIBILITY": "Aplica os limites de escopo do estudo (peso molecular e número "
    "de átomos) sobre a estrutura principal, nunca sobre o sal.",
    "CANONICALIZE": "Gera a forma canônica da estrutura e o identificador InChIKey.",
    DEDUP_STAGE: "Identifica estruturas repetidas pelo InChIKey completo.",
}

STATUS_LABELS: dict[StageStatus, tuple[str, str, str]] = {
    StageStatus.PENDING: ("○", "Aguardando", "#9aa0a6"),
    StageStatus.RUNNING: ("◐", "Executando", "#1a73e8"),
    StageStatus.SUCCESS: ("✓", "Concluído", "#1e8e3e"),
    StageStatus.WARNING: ("!", "Concluído com remoções", "#b06000"),
    StageStatus.FAILED: ("✗", "Falhou", "#d93025"),
    StageStatus.SKIPPED: ("–", "Ignorado", "#9aa0a6"),
}

REJECTION_LABELS: dict[str, str] = {
    "ERR_SYNTAX": "Estrutura não interpretável",
    "ERR_EMPTY": "Estrutura vazia",
    "ERR_VALENCE": "Valência inválida",
    "ERR_KEKULIZE": "Aromaticidade inconsistente",
    "ERR_SANITIZE": "Estrutura quimicamente inválida",
    "ERR_STANDARDIZE": "Falha na padronização",
    "ERR_GET_PARENT": "Falha ao isolar a estrutura principal",
    "ERR_CANONICALIZE": "Falha ao gerar a forma canônica",
    "ERR_INCHI": "Falha ao gerar o InChIKey",
    "ERR_ORGANOMETALLIC": "Organometálico descartado por política",
    "ERR_MW_LIMIT": "Peso molecular acima do limite",
    "ERR_HA_LIMIT": "Átomos demais",
    "ERR_INTERNAL": "Erro interno do processamento",
    "EXACT_DUPLICATE": "Duplicata",
    "ANNOTATION_CONFLICT": "Duplicata com anotação divergente",
    "BLOCK1_COLLISION": "Estrutura aparentada (não é duplicata)",
    "UNKNOWN": "Motivo não registrado",
}

#: Colunas de exportação com rótulo legível. Só entram colunas que a engenharia
#: realmente produz — não há "warnings" no schema, então não é oferecida.
EXPORT_LABELS: list[tuple[str, str]] = [
    ("input_id", "Identificador"),
    ("raw_smiles", "Estrutura original"),
    ("curated_smiles", "Estrutura curada"),
    ("inchikey", "InChIKey"),
    ("status", "Status"),
    ("rejection_code", "Motivo da rejeição"),
    ("rejection_detail", "Detalhe da rejeição"),
    ("inchi", "InChI"),
    ("molecular_formula", "Fórmula molecular"),
    ("parent_mw", "Peso molecular"),
    ("parent_heavy_atoms", "Átomos pesados"),
    ("removed_fragments", "Fragmentos removidos"),
    ("salt_removed", "Teve sal removido"),
    ("n_components_parent", "Componentes"),
    ("delta_stereocenters", "Variação de centros quirais"),
    ("transformations", "Transformações registradas"),
    ("policy_hash", "Hash da política"),
]

DEFAULT_EXPORT = ("input_id", "raw_smiles", "curated_smiles", "inchikey", "status")

HELP = {
    "max_mw": "Limite máximo de peso molecular do critério de elegibilidade. O "
    "valor usado nesta execução é registrado no manifesto para permitir reprodução.",
    "max_ha": "Número máximo de átomos não-hidrogênio permitido pelo critério de "
    "elegibilidade.",
    "dedup": "Identifica estruturas duplicadas usando o InChIKey completo. O "
    "primeiro bloco do InChIKey não é usado isoladamente para declarar duplicatas — "
    "enantiômeros o compartilham.",
    "policy": "Identificador da versão da política de curadoria usada nesta "
    "execução. Permite verificar se duas execuções usaram a mesma política.",
    "inchikey": "Identificador textual derivado da representação InChI da estrutura.",
    "formats": "Você pode enviar CSV, TSV ou SMI. O sistema identifica a coluna de "
    "estruturas quando aplicável e faz uma triagem sintática antes da execução.",
    "manifest": "O manifesto registra as informações necessárias para identificar e "
    "auditar esta execução: parâmetros, versões, hashes e metadados.",
}


def label_stage(name: str) -> str:
    return STAGE_LABELS.get(name, name.replace("_", " ").title())


def label_reason(code: str) -> str:
    return REJECTION_LABELS.get(code, code)


# --- Estilo ----------------------------------------------------------------------


def inject_styles() -> None:
    st.markdown(
        """
        <style>
          .block-container { padding-top: 2.2rem; max-width: 1150px; }
          .app-title { font-size: 1.5rem; font-weight: 700; margin-bottom: .1rem; }
          .app-sub { opacity: .7; font-size: .92rem; }
          .run-chip { font-size: .8rem; opacity: .75;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace; }
          .stage-box { border: 1px solid rgba(128,128,128,.3); border-radius: 6px;
            padding: .5rem .35rem; text-align: center; line-height: 1.4; }
          .stage-label { font-size: .72rem; opacity: .8; }
          .stage-mark { font-size: 1.15rem; font-weight: 700; }
          .stage-n { font-size: 1.15rem; font-weight: 600;
            font-variant-numeric: tabular-nums; }
          .stage-out { font-size: .7rem; color: #b06000; }
          .stage-ms { font-size: .66rem; opacity: .55;
            font-variant-numeric: tabular-nums; }
          .mono { font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
            font-size: .78rem; }
          .arrow { text-align:center; padding-top:2.1rem; opacity:.35; }
        </style>
        """,
        unsafe_allow_html=True,
    )


# --- Nível 1: cabeçalho e ajuda global ---------------------------------------------


def render_header(report: Optional[RunReport]) -> None:
    left, right = st.columns([3, 1])
    with left:
        st.markdown("<div class='app-title'>Structure Curation</div>", unsafe_allow_html=True)
        st.markdown(
            "<div class='app-sub'>Padronize, valide e rastreie suas estruturas "
            "químicas.</div>",
            unsafe_allow_html=True,
        )
    with right:
        if report is not None:
            st.markdown(
                f"<div class='run-chip'>Run {report.provenance.run_id}</div>"
                "<div style='color:#1e8e3e;font-weight:600'>✓ Concluído</div>",
                unsafe_allow_html=True,
            )

    with st.expander("Como funciona"):
        st.markdown(
            "O pipeline recebe estruturas químicas, executa as etapas de validação "
            "e curadoria e produz resultados rastreáveis. Cada execução registra "
            "parâmetros, versões, decisões e resultados, para facilitar auditoria e "
            "reprodução.\n\n"
            "**O caminho é sempre o mesmo:** carregar → configurar → executar → "
            "ver o pipeline → entender o resultado → investigar uma estrutura → "
            "baixar os dados → reproduzir a execução."
        )


def render_global_help() -> None:
    with st.expander("Como interpretar esta página"):
        st.markdown(
            "**Cores e símbolos das etapas** — cada estado aparece com ícone *e* "
            "texto, nunca só por cor:\n\n"
            "| | significado |\n| --- | --- |\n"
            "| ○ Aguardando | a etapa ainda não rodou |\n"
            "| ◐ Executando | a etapa está em andamento |\n"
            "| ✓ Concluído | todas as estruturas passaram |\n"
            "| ! Concluído com remoções | a etapa removeu estruturas |\n"
            "| ✗ Falhou | nenhuma estrutura passou |\n"
            "| – Ignorado | a etapa não se aplicou |\n\n"
            "**O que é uma etapa** — uma operação do pipeline. O número embaixo do "
            "ícone é quantas estruturas *saíram* dela.\n\n"
            "**Aprovada, rejeitada e duplicata são coisas diferentes.** Rejeitada é "
            "uma estrutura que não passou por algum critério, e o motivo fica "
            "registrado. Duplicata é uma estrutura válida que já havia aparecido "
            "antes — ela continua no arquivo completo, apenas não conta como "
            "identidade nova.\n\n"
            "**Estruturas aparentadas não são duplicatas.** Dois enantiômeros "
            "compartilham o começo do InChIKey mas são compostos distintos, com "
            "atividades biológicas possivelmente diferentes. O sistema os reporta "
            "como aparentados e nunca os funde.\n\n"
            "**No CSV**, cada linha é uma estrutura de entrada. A coluna `status` "
            "diz se foi aprovada, e `rejection_code` diz por que não foi.\n\n"
            "**Para reproduzir**, use a seção *Proveniência da execução*."
        )


# --- Nível 1: entrada ---------------------------------------------------------------


def render_input() -> tuple[Optional[bytes], str]:
    st.markdown("### 1. Carregar estruturas")
    upload_tab, paste_tab = st.tabs(["📁 Enviar arquivo", "✎ Colar estruturas"])

    raw: Optional[bytes] = None
    name = ""

    with upload_tab:
        st.caption("CSV · TSV · SMI", help=HELP["formats"])
        uploaded = st.file_uploader(
            "Arquivo", type=["csv", "tsv", "smi", "smiles", "txt"],
            label_visibility="collapsed",
        )
        if uploaded is not None:
            raw, name = uploaded.getvalue(), uploaded.name

    with paste_tab:
        st.caption("Uma estrutura por linha")
        text = st.text_area(
            "Estruturas", height=150, label_visibility="collapsed",
            placeholder="CC(=O)O[Na]\nN[C@@H](C)C(=O)O.Cl\nCC(=O)Oc1ccccc1C(=O)O",
        )
        if text.strip():
            raw, name = text.encode("utf-8"), "estruturas coladas"

    return raw, name


def render_input_summary(raw: bytes, name: str) -> None:
    preview = preview_input(raw.decode("utf-8", errors="replace"))
    valid = preview.total - preview.n_invalid

    st.markdown(f"**{name}** — {preview.total:,} estruturas encontradas".replace(",", "."))
    st.markdown(f"✓ {valid:,} com sintaxe válida".replace(",", "."))
    if preview.n_invalid:
        st.markdown(f"! {preview.n_invalid} precisam de atenção")
        st.caption(
            "Esta é apenas uma triagem de sintaxe. Ela não substitui a validação do "
            "pipeline: uma estrutura aprovada aqui ainda pode ser rejeitada depois."
        )

    with st.expander("Ver dados de entrada"):
        st.dataframe(
            [{"Identificador": i, "Estrutura": s} for i, s in preview.head],
            width="stretch", hide_index=True,
        )
        if preview.invalid:
            st.markdown("**Estruturas que precisam de atenção**")
            st.dataframe(
                [{"Identificador": i, "Estrutura": s} for i, s in preview.invalid],
                width="stretch", hide_index=True,
            )


# --- Nível 1: configuração ------------------------------------------------------------


def render_configuration() -> dict:
    st.markdown("### 2. Configuração da execução")
    first, second, third = st.columns(3)

    max_mw = first.number_input(
        "Peso molecular máximo (Da)", 50.0, 10000.0,
        float(EligibilityCriteria.max_molecular_weight), step=50.0,
        help=HELP["max_mw"],
    )
    max_ha = second.number_input(
        "Átomos pesados máximos", 5, 1000,
        int(EligibilityCriteria.max_heavy_atoms), step=5, help=HELP["max_ha"],
    )
    third.markdown("<div style='height:1.75rem'></div>", unsafe_allow_html=True)
    deduplicate = third.checkbox(
        "Deduplicar estruturas", value=True, help=HELP["dedup"]
    )

    with st.expander("Configuração avançada"):
        decisions = st.text_input(
            "Arquivo de política de curadoria", str(DEFAULT_DECISIONS),
            help=HELP["policy"],
        )
        path = Path(decisions)
        if path.is_file():
            policy_hash = compute_policy_hash(path)
            st.markdown(f"**Política utilizada**  \n`{decisions}`  \nVersão: identificada ✓")
            with st.expander("Detalhes"):
                st.code(f"policy_hash: {policy_hash}", language="text")
        else:
            policy_hash = "UNVERSIONED_POLICY"
            st.warning(
                "Arquivo de política não encontrado. A execução será marcada como "
                "não versionada e não poderá ser comparada a execuções versionadas."
            )

    return {
        "max_mw": max_mw,
        "max_ha": int(max_ha),
        "deduplicate": deduplicate,
        "policy_hash": policy_hash,
        "policy_path": decisions,
    }


# --- Nível 1: pipeline ------------------------------------------------------------------


def all_stages(report: RunReport) -> list[StageReport]:
    return list(report.stages) + ([report.dedup] if report.dedup else [])


def stage_box(stage: StageReport) -> str:
    mark, _, color = STATUS_LABELS[stage.status]
    removed = (
        f"<div class='stage-out'>−{stage.n_excluded}</div>"
        if stage.has_exclusions else "<div class='stage-out'>&nbsp;</div>"
    )
    duration = (
        f"<div class='stage-ms'>{stage.duration_seconds * 1000:.0f} ms</div>"
        if stage.duration_seconds else "<div class='stage-ms'>&nbsp;</div>"
    )
    return (
        f"<div class='stage-box'>"
        f"<div class='stage-mark' style='color:{color}'>{mark}</div>"
        f"<div class='stage-label'>{label_stage(stage.name)}</div>"
        f"<div class='stage-n'>{stage.n_output}</div>"
        f"{removed}{duration}</div>"
    )


def render_pipeline(report: RunReport) -> str:
    stages = all_stages(report)
    columns = st.columns(len(stages) * 2 - 1)
    for position, stage in enumerate(stages):
        with columns[position * 2]:
            st.markdown(stage_box(stage), unsafe_allow_html=True)
        if position < len(stages) - 1:
            columns[position * 2 + 1].markdown(
                "<div class='arrow'>→</div>", unsafe_allow_html=True
            )

    st.caption("Clique em uma etapa para entender o que aconteceu")
    return st.radio(
        "Etapa", [stage.name for stage in stages],
        format_func=label_stage, horizontal=True, label_visibility="collapsed",
    )


# --- Nível 1: o que aconteceu -------------------------------------------------------------


def render_summary(report: RunReport) -> None:
    kpis = report.kpis()
    st.markdown("### O que aconteceu?")

    lines = [
        f"**{kpis['Total Processed']:,} estruturas foram processadas.**".replace(",", "."),
        f"✓ {kpis['Approved']:,} aprovadas".replace(",", "."),
    ]
    if kpis["Rejected"]:
        lines.append(f"! {kpis['Rejected']:,} rejeitadas".replace(",", "."))
    if kpis["Duplicates"]:
        lines.append(f"! {kpis['Duplicates']:,} duplicatas".replace(",", "."))
    if kpis["Warnings"]:
        lines.append(
            f"· {kpis['Warnings']:,} tiveram a estrutura alterada pela curadoria".replace(",", ".")
        )
    st.markdown("  \n".join(lines))

    largest = report.largest_reduction()
    if largest is not None:
        st.caption(
            f"A maior redução ocorreu em **{label_stage(largest.name)}** "
            f"({largest.n_excluded} estruturas)."
        )


# --- Nível 2: detalhe da etapa ----------------------------------------------------------


def render_stage_detail(report: RunReport, name: str) -> None:
    stage = next(s for s in all_stages(report) if s.name == name)
    mark, status_text, color = STATUS_LABELS[stage.status]

    st.markdown(
        f"#### {label_stage(name)} "
        f"<span style='color:{color};font-size:.9rem'>{mark} {status_text}</span>",
        unsafe_allow_html=True,
    )
    st.markdown(f"*{STAGE_HELP.get(name, stage.description)}*")

    a, b, c, d = st.columns(4)
    a.metric("Entrada", stage.n_input)
    b.metric("Resultado", stage.n_output)
    c.metric("Removidas", stage.n_excluded)
    d.metric(
        "Tempo",
        f"{stage.duration_seconds * 1000:.0f} ms" if stage.duration_seconds else "—",
    )

    if name == DEDUP_STAGE:
        st.info(
            "Estruturas aparentadas — que compartilham o começo do InChIKey — não "
            "são contadas como duplicatas. Enantiômeros são compostos distintos."
        )

    removed = [
        record for record in report.records
        if not record.passed and record.rejection_stage
        and record.rejection_stage.value == name
    ]
    if removed:
        with st.expander(f"Ver as {len(removed)} estruturas removidas nesta etapa"):
            st.dataframe(
                [
                    {
                        "Identificador": r.input_id,
                        "Estrutura": r.raw_smiles,
                        "Motivo": label_reason(
                            r.rejection_code.value if r.rejection_code else "UNKNOWN"
                        ),
                        "Detalhe": r.rejection_detail or "",
                    }
                    for r in removed
                ],
                width="stretch", hide_index=True,
            )

    with st.expander("Parâmetros usados nesta execução"):
        st.json(report.provenance.parameters)


# --- Nível 2: exclusões --------------------------------------------------------------------


def render_exclusions(report: RunReport) -> None:
    groups = report.exclusion_groups()
    if not groups:
        st.success("Nenhuma estrutura foi removida.")
        return

    total = sum(group.count for group in groups)
    st.markdown(f"**{total} estruturas removidas.** Nenhuma desaparece sem registro.")
    st.dataframe(
        [
            {
                "Motivo": label_reason(group.reason),
                "Etapa": label_stage(group.stage),
                "Quantidade": group.count,
            }
            for group in groups
        ],
        width="stretch", hide_index=True,
    )

    options = [f"{label_reason(g.reason)} — {g.count}" for g in groups]
    chosen = st.selectbox("Abrir um motivo", options)
    group = groups[options.index(chosen)]
    affected = report.records_by_id(group.record_ids)

    st.dataframe(
        [
            {
                "Identificador": r.input_id,
                "Estrutura": r.raw_smiles,
                "Detalhe": r.rejection_detail or "",
            }
            for r in affected
        ],
        width="stretch", hide_index=True,
    )
    st.download_button(
        "Baixar estas estruturas (CSV)",
        to_csv(affected, ("input_id", "raw_smiles", "rejection_code", "rejection_detail")),
        file_name=f"removidas_{group.reason.lower()}.csv", mime="text/csv",
    )


# --- Nível 2: rastrear uma estrutura ----------------------------------------------------------


def draw(smiles: str, size: int = 240):
    """Imagem 2D da estrutura, ou ``None`` quando não é possível renderizar."""
    if not smiles or not DRAWING_AVAILABLE:
        return None
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None
    try:
        mol.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(mol)
        return Draw.MolToImage(mol, size=(size, size))
    except Exception:
        return None


def render_trace(report: RunReport) -> None:
    st.caption(
        "Veja como uma estrutura mudou durante o processamento e em qual etapa cada "
        "decisão foi tomada."
    )
    if not DRAWING_AVAILABLE:
        st.info(
            "Este ambiente não consegue desenhar estruturas "
            f"(`{DRAWING_ERROR}`). Os textos e a trajetória continuam disponíveis."
        )

    identifiers = [record.input_id for record in report.records]
    if not identifiers:
        return
    chosen = st.selectbox("Estrutura", identifiers)
    record = next(r for r in report.records if r.input_id == chosen)

    original, curated = st.columns(2)
    with original:
        st.markdown("**Estrutura original**")
        image = draw(record.raw_smiles)
        if image is not None:
            st.image(image)
        st.code(record.raw_smiles, language="text")
    with curated:
        st.markdown("**Estrutura final**" if record.passed else "**Rejeitada**")
        image = draw(record.curated_smiles or "")
        if image is not None:
            st.image(image)
        st.code(record.curated_smiles or "—", language="text")

    if not record.passed:
        code = record.rejection_code.value if record.rejection_code else "UNKNOWN"
        stage = record.rejection_stage.value if record.rejection_stage else ""
        st.error(
            f"**{label_reason(code)}** em *{label_stage(stage)}* — "
            f"{record.rejection_detail}"
        )

    steps = structure_lineage(record)
    st.markdown("**Trajetória**")
    if steps:
        for step in steps:
            st.markdown(
                f"↓ **{label_stage(step.stage)}** — {step.detail}  \n"
                f"<span class='mono'>{step.before_smiles} → {step.after_smiles}</span>",
                unsafe_allow_html=True,
            )
    else:
        st.caption(
            "A camada de engenharia não registrou transformações intermediárias "
            "para esta estrutura."
        )

    if record.removed_fragments:
        st.markdown(f"**Fragmentos removidos:** `{record.removed_fragments}`")
    if record.inchikey:
        st.markdown(f"**InChIKey:** `{record.inchikey}`", help=HELP["inchikey"])


# --- Nível 2: resultados ---------------------------------------------------------------------


def render_results(report: RunReport) -> None:
    scope = st.radio(
        "Filtro", ["Todos", "Aprovados", "Rejeitados", "Transformados"],
        horizontal=True, label_visibility="collapsed",
    )
    selected = {
        "Todos": report.records,
        "Aprovados": report.approved,
        "Rejeitados": report.rejected,
        "Transformados": [r for r in report.approved if r.transformations],
    }[scope]

    query = st.text_input(
        "Pesquisar estrutura, identificador ou InChIKey", "",
        label_visibility="collapsed",
        placeholder="Pesquisar estrutura, identificador ou InChIKey",
    )
    if query:
        needle = query.lower()
        selected = [
            record for record in selected
            if needle in record.input_id.lower()
            or needle in record.raw_smiles.lower()
            or needle in (record.curated_smiles or "").lower()
            or needle in (record.inchikey or "").lower()
        ]

    st.caption(f"{len(selected)} estruturas")
    st.dataframe(
        [
            {
                "Identificador": r.input_id,
                "Original": r.raw_smiles,
                "Curada": r.curated_smiles or "",
                "InChIKey": r.inchikey or "",
                "Status": "Aprovada" if r.passed else "Rejeitada",
                "Motivo": label_reason(r.rejection_code.value)
                if r.rejection_code else "",
            }
            for r in selected
        ],
        width="stretch", hide_index=True,
    )


# --- Nível 2: exportação -------------------------------------------------------------------------


def render_downloads(report: RunReport) -> None:
    st.caption(
        "Escolha entre baixar tudo ou selecionar somente as informações necessárias."
    )

    st.markdown("**📦 Todos os dados**")
    st.caption("Inclui estruturas aprovadas, rejeitadas, status e motivos.")
    st.download_button(
        "Baixar CSV completo", full_csv(report),
        file_name="structures_full.csv", mime="text/csv",
    )

    st.divider()
    st.markdown("**Seleção personalizada**")
    labels = {label: column for column, label in EXPORT_LABELS}
    chosen = st.multiselect(
        "Colunas", list(labels),
        [label for column, label in EXPORT_LABELS if column in DEFAULT_EXPORT],
        label_visibility="collapsed",
    )
    if chosen:
        st.download_button(
            "Baixar CSV selecionado",
            to_csv(report.records, [labels[label] for label in chosen]),
            file_name="structures_custom.csv", mime="text/csv",
        )

    st.divider()
    st.markdown("**Estruturas rejeitadas**")
    st.download_button(
        "Baixar rejeitadas", rejected_csv(report),
        file_name="rejected_structures.csv", mime="text/csv",
    )

    st.divider()
    st.markdown("**Reprodutibilidade**")
    manifest = run_manifest(report)
    first, second = st.columns(2)
    with first:
        st.caption("📋 Manifesto da execução", help=HELP["manifest"])
        st.download_button(
            "Baixar manifest.json",
            json.dumps(manifest, indent=2, ensure_ascii=False),
            file_name="run_manifest.json", mime="application/json",
        )
    with second:
        st.caption("📦 Pacote completo — dados, manifesto e instruções")
        st.download_button(
            "Baixar pacote de reprodução", reproducibility_package(report),
            file_name=f"structure-curation-run-{report.provenance.run_id}.zip",
            mime="application/zip",
        )


# --- Nível 3: proveniência ----------------------------------------------------------------------------


def render_provenance(report: RunReport) -> None:
    provenance = report.provenance
    blockers = provenance.reproduction_blockers()

    if blockers:
        st.markdown("**! Reprodução exata não garantida**")
        st.markdown("\n".join(f"- {item}" for item in blockers))
    else:
        st.markdown("**✓ Esta execução pode ser reproduzida**")

    a, b, c, d = st.columns(4)
    a.metric("Pipeline", provenance.pipeline_version)
    b.metric("Git commit", provenance.git.commit[:7])
    c.metric(
        "Política",
        "versionada" if provenance.policy_hash != "UNVERSIONED_POLICY" else "ausente",
    )
    d.metric("RDKit", provenance.versions.get("rdkit", "—"))

    st.markdown("**Comando de reprodução**")
    st.code(provenance.reproduction_command(provenance.input_name), language="bash")

    with st.expander("Ver detalhes técnicos"):
        st.markdown("**Execução**")
        st.json(
            {
                "run_id": provenance.run_id,
                "started_at": provenance.started_at,
                "finished_at": provenance.finished_at,
                "duration_seconds": provenance.duration_seconds,
                "input_name": provenance.input_name,
                "input_sha256": provenance.input_hash,
                "output_sha256": provenance.output_hash,
            }
        )
        st.markdown("**Código e política**")
        st.json(
            {
                "pipeline_version": provenance.pipeline_version,
                "git_commit": provenance.git.commit,
                "git_branch": provenance.git.branch,
                "git_dirty": provenance.git.dirty,
                "policy_hash": provenance.policy_hash,
                "policy_path": provenance.policy_path,
            }
        )
        st.markdown("**Parâmetros**")
        st.json(provenance.parameters)
        st.markdown("**Versões e ambiente**")
        st.json({**provenance.versions, **provenance.environment})
        st.caption(
            "Mudanças de versão do RDKit alteram percepção de aromaticidade e regras "
            "de padronização: duas execuções com RDKit diferente não são comparáveis "
            "mesmo com a mesma política."
        )


# --- Aplicação --------------------------------------------------------------------------------------------


def execute(raw: bytes, name: str, configuration: dict) -> RunReport:
    pipeline = CurationPipeline(
        policy_hash=configuration["policy_hash"],
        criteria=EligibilityCriteria(
            max_molecular_weight=configuration["max_mw"],
            max_heavy_atoms=configuration["max_ha"],
        ),
        deduplicate=configuration["deduplicate"],
    )

    placeholder = st.empty()
    bar = st.progress(0.0)

    def report_progress(position: int, total: int, _identifier: str) -> None:
        # Progresso real: posição do registro no lote, informada pelo pipeline.
        # Não há percentual por etapa porque a engenharia não fornece um.
        placeholder.markdown(f"Processando estrutura {position} de {total}…")
        bar.progress(position / total if total else 1.0)

    report = pipeline.run_report(
        raw.decode("utf-8", errors="replace"),
        parameters={
            "max_mw": configuration["max_mw"],
            "max_ha": configuration["max_ha"],
            "deduplicate": configuration["deduplicate"],
        },
        input_bytes=raw,
        input_name=name,
        policy_path=configuration["policy_path"],
        progress=report_progress,
    )
    placeholder.empty()
    bar.empty()
    return report


def main() -> None:
    st.set_page_config(page_title="Structure Curation", layout="wide")
    inject_styles()

    report: Optional[RunReport] = st.session_state.get("report")
    render_header(report)
    st.divider()

    raw, name = render_input()
    if raw:
        render_input_summary(raw, name)

    st.divider()
    configuration = render_configuration()

    st.divider()
    if st.button(
        "▶ Executar curadoria", type="primary", disabled=raw is None,
        width="content",
    ):
        st.session_state["report"] = execute(raw, name, configuration)
        st.rerun()

    if report is None:
        st.info("Carregue um arquivo ou cole estruturas para começar.")
        render_global_help()
        return

    st.divider()
    st.markdown("### 3. Pipeline")
    selected = render_pipeline(report)
    st.divider()
    render_stage_detail(report, selected)

    st.divider()
    render_summary(report)

    st.divider()
    exclusions, trace, results, downloads, provenance = st.tabs(
        [
            "O que foi removido",
            "Rastrear uma estrutura",
            "Resultados",
            "Baixar resultados",
            "Proveniência da execução",
        ]
    )
    with exclusions:
        render_exclusions(report)
    with trace:
        render_trace(report)
    with results:
        render_results(report)
    with downloads:
        render_downloads(report)
    with provenance:
        render_provenance(report)

    st.divider()
    render_global_help()


if __name__ == "__main__":
    main()
