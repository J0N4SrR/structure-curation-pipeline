"""Interface Streamlit: cliente fino sobre a camada de engenharia.

Este módulo é **exclusivamente apresentação**. Nenhuma regra química, nenhum
critério de aceitação e nenhuma métrica são calculados aqui: tudo vem de
:mod:`curation.pipeline`, :mod:`curation.reporting` e :mod:`curation.provenance`.

Quando uma informação não existe na camada de engenharia, a interface diz que não
existe. Não há estado estimado, progresso simulado nem linhagem inferida.

Execução::

    streamlit run src/curation/app.py
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
    # rdMolDraw2D linka contra libXrender/libX11/libXext do sistema, que faltam em
    # containers slim. A renderizacao e uma funcionalidade entre varias: sem ela o
    # app degrada para SMILES em texto, em vez de derrubar a auditoria inteira.
    Draw = None
    DRAWING_AVAILABLE = False
    DRAWING_ERROR = str(_error)

from curation.filters import EligibilityCriteria
from curation.io import compute_policy_hash, preview_input
from curation.pipeline import CurationPipeline
from curation.reporting import (
    DEDUP_STAGE,
    EXPORTABLE_COLUMNS,
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

STATUS_STYLE: dict[StageStatus, tuple[str, str]] = {
    StageStatus.PENDING: ("○", "#9aa0a6"),
    StageStatus.RUNNING: ("◐", "#1a73e8"),
    StageStatus.SUCCESS: ("✓", "#1e8e3e"),
    StageStatus.WARNING: ("!", "#f9ab00"),
    StageStatus.FAILED: ("✗", "#d93025"),
    StageStatus.SKIPPED: ("–", "#9aa0a6"),
}

DEFAULT_EXPORT_COLUMNS = (
    "input_id",
    "raw_smiles",
    "curated_smiles",
    "inchikey",
    "status",
)


# --- Estilo ---------------------------------------------------------------------


def inject_styles() -> None:
    st.markdown(
        """
        <style>
          .block-container { padding-top: 2rem; max-width: 1200px; }
          .stage-card {
            border: 1px solid rgba(128,128,128,.28); border-radius: 6px;
            padding: .55rem .5rem; text-align: center; line-height: 1.35;
          }
          .stage-name { font-size: .68rem; letter-spacing: .06em;
            text-transform: uppercase; opacity: .75; }
          .stage-mark { font-size: 1.25rem; font-weight: 700; }
          .stage-count { font-size: 1.1rem; font-weight: 600; font-variant-numeric: tabular-nums; }
          .stage-drop { font-size: .7rem; color: #d93025; }
          .stage-time { font-size: .66rem; opacity: .6; font-variant-numeric: tabular-nums; }
          .run-banner { display:flex; justify-content:space-between; align-items:baseline;
            border-bottom:1px solid rgba(128,128,128,.28); padding-bottom:.5rem; margin-bottom:1rem; }
          .run-title { font-size:1.05rem; font-weight:700; letter-spacing:.04em; }
          .mono { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size:.78rem; }
        </style>
        """,
        unsafe_allow_html=True,
    )


def stage_card(stage: StageReport) -> str:
    mark, color = STATUS_STYLE[stage.status]
    drop = (
        f"<div class='stage-drop'>−{stage.n_excluded}</div>"
        if stage.has_exclusions
        else "<div class='stage-drop'>&nbsp;</div>"
    )
    duration = (
        f"<div class='stage-time'>{stage.duration_seconds * 1000:.0f} ms</div>"
        if stage.duration_seconds
        else "<div class='stage-time'>&nbsp;</div>"
    )
    return (
        f"<div class='stage-card'>"
        f"<div class='stage-name'>{stage.name.replace('_', ' ')}</div>"
        f"<div class='stage-mark' style='color:{color}'>{mark}</div>"
        f"<div class='stage-count'>{stage.n_output}</div>"
        f"{drop}{duration}</div>"
    )


# --- Barra lateral ---------------------------------------------------------------


def sidebar() -> dict:
    """Coleta parâmetros. Não valida regras — só encaminha à engenharia."""
    st.sidebar.markdown("### Parâmetros de curadoria")
    st.sidebar.caption(
        "Os cortes incidem sobre a estrutura-mãe isolada, nunca sobre o sal (D-10)."
    )
    max_mw = st.sidebar.number_input(
        "Peso molecular máximo (Da)", 50.0, 10000.0,
        float(EligibilityCriteria.max_molecular_weight), step=50.0,
    )
    max_ha = st.sidebar.number_input(
        "Átomos pesados máximos", 5, 1000,
        int(EligibilityCriteria.max_heavy_atoms), step=5,
    )
    deduplicate = st.sidebar.checkbox(
        "Deduplicar por InChIKey", value=True,
        help="Identidade é o InChIKey completo; o bloco1 apenas agrupa (D-07).",
    )
    decisions = st.sidebar.text_input(
        "Arquivo de política (ADRs)", str(DEFAULT_DECISIONS)
    )

    path = Path(decisions)
    if path.is_file():
        policy_hash = compute_policy_hash(path)
        st.sidebar.success(f"policy_hash `{policy_hash[:16]}…`")
    else:
        policy_hash = "UNVERSIONED_POLICY"
        st.sidebar.warning(
            "Política não encontrada. O lote será marcado como não versionado e "
            "não será comparável a lotes versionados."
        )

    return {
        "max_mw": max_mw,
        "max_ha": max_ha,
        "deduplicate": deduplicate,
        "policy_hash": policy_hash,
        "policy_path": decisions,
    }


# --- Entrada ---------------------------------------------------------------------


def input_area() -> tuple[Optional[bytes], str]:
    upload_tab, paste_tab = st.tabs(["Upload de arquivo", "Colar SMILES"])

    with upload_tab:
        uploaded = st.file_uploader(
            "CSV, TSV ou SMI", type=["csv", "tsv", "smi", "smiles", "txt"]
        )
        if uploaded is not None:
            return uploaded.getvalue(), uploaded.name

    with paste_tab:
        text = st.text_area(
            "Uma estrutura por linha", height=170,
            placeholder="CC(=O)O[Na]\nN[C@@H](C)C(=O)O.Cl\nCC(=O)Oc1ccccc1C(=O)O",
        )
        if text.strip():
            return text.encode("utf-8"), "colado.smi"

    return None, ""


def render_preview(raw: bytes, name: str) -> None:
    preview = preview_input(raw.decode("utf-8", errors="replace"))

    left, middle, right = st.columns(3)
    left.metric("Registros", preview.total)
    middle.metric("Inválidos (sintaxe)", preview.n_invalid)
    right.metric("Fonte", name)

    st.caption(f"Coluna de estrutura: {preview.smiles_column}")
    if preview.head:
        st.dataframe(
            [{"input_id": i, "raw_smiles": s} for i, s in preview.head],
            use_container_width=True, hide_index=True,
        )
    if preview.invalid:
        with st.expander(f"{preview.n_invalid} inválidos na triagem preliminar"):
            st.caption(
                "Triagem apenas sintática. Não substitui a validação do pipeline: "
                "uma estrutura aprovada aqui ainda pode ser rejeitada por valência."
            )
            st.dataframe(
                [{"input_id": i, "raw_smiles": s} for i, s in preview.invalid],
                use_container_width=True, hide_index=True,
            )


# --- Pipeline --------------------------------------------------------------------


def render_funnel(report: RunReport) -> str:
    stages = list(report.stages) + ([report.dedup] if report.dedup else [])
    columns = st.columns(len(stages) * 2 - 1)

    for position, stage in enumerate(stages):
        with columns[position * 2]:
            st.markdown(stage_card(stage), unsafe_allow_html=True)
        if position < len(stages) - 1:
            columns[position * 2 + 1].markdown(
                "<div style='text-align:center;padding-top:2.4rem;opacity:.4'>→</div>",
                unsafe_allow_html=True,
            )

    return st.radio(
        "Etapa em detalhe",
        [stage.name for stage in stages],
        horizontal=True,
        label_visibility="collapsed",
    )


def render_stage_detail(report: RunReport, name: str) -> None:
    stages = {s.name: s for s in list(report.stages) + ([report.dedup] if report.dedup else [])}
    stage = stages[name]
    mark, color = STATUS_STYLE[stage.status]

    st.markdown(
        f"#### {name.replace('_', ' ')} "
        f"<span style='color:{color}'>{mark} {stage.status.value}</span>",
        unsafe_allow_html=True,
    )
    st.caption(stage.description)

    summary, records_tab, parameters_tab, exclusions_tab, provenance_tab = st.tabs(
        ["Summary", "Registros", "Parameters", "Exclusões", "Provenance"]
    )

    with summary:
        a, b, c, d = st.columns(4)
        a.metric("Entrada", stage.n_input)
        b.metric("Saída", stage.n_output)
        c.metric("Excluídos", stage.n_excluded)
        d.metric(
            "Duração",
            f"{stage.duration_seconds * 1000:.0f} ms" if stage.duration_seconds else "—",
        )
        if stage.name == DEDUP_STAGE:
            st.info(
                "Colisões de bloco1 **não** são duplicatas: enantiômeros "
                "compartilham o primeiro bloco do InChIKey (D-07)."
            )

    with records_tab:
        affected = [
            record for record in report.records
            if not record.passed
            and record.rejection_stage
            and record.rejection_stage.value == name
        ]
        if affected:
            st.dataframe(
                [record_row(r, ("input_id", "raw_smiles", "rejection_code", "rejection_detail")) for r in affected],
                use_container_width=True, hide_index=True,
            )
        else:
            st.caption("Nenhum registro foi rejeitado nesta etapa.")

    with parameters_tab:
        st.json(report.provenance.parameters)

    with exclusions_tab:
        groups = [g for g in report.exclusion_groups() if g.stage == name]
        if not groups:
            st.caption("Nenhuma exclusão nesta etapa.")
        for group in groups:
            with st.expander(f"{group.count}  {group.reason}"):
                st.dataframe(
                    [record_row(r, ("input_id", "raw_smiles", "rejection_detail"))
                     for r in report.records_by_id(group.record_ids)],
                    use_container_width=True, hide_index=True,
                )

    with provenance_tab:
        st.json(
            {
                "rdkit": report.provenance.versions.get("rdkit"),
                "chembl_structure_pipeline": report.provenance.versions.get(
                    "chembl_structure_pipeline"
                ),
                "policy_hash": report.provenance.policy_hash,
                "pipeline_version": report.provenance.pipeline_version,
            }
        )


def render_exclusions(report: RunReport) -> None:
    groups = report.exclusion_groups()
    if not groups:
        st.success("Nenhuma estrutura foi excluída.")
        return

    total = sum(group.count for group in groups)
    st.markdown(f"**{total} excluídas** — nenhuma desaparece sem registro.")
    st.dataframe(
        [{"motivo": g.reason, "etapa": g.stage, "n": g.count} for g in groups],
        use_container_width=True, hide_index=True,
    )
    chosen = st.selectbox(
        "Abrir registros de uma categoria",
        [f"{g.reason} @ {g.stage} ({g.count})" for g in groups],
    )
    group = groups[[f"{g.reason} @ {g.stage} ({g.count})" for g in groups].index(chosen)]
    st.dataframe(
        [record_row(r, ("input_id", "raw_smiles", "rejection_stage", "rejection_detail"))
         for r in report.records_by_id(group.record_ids)],
        use_container_width=True, hide_index=True,
    )


# --- Rastreabilidade e diff ---------------------------------------------------------


def draw(smiles: str, size: int = 260):
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


def render_lineage(report: RunReport) -> None:
    if not DRAWING_AVAILABLE:
        st.warning(
            "Renderização de estruturas indisponível neste ambiente: "
            f"`{DRAWING_ERROR}`. Os SMILES e a trajetória continuam abaixo — "
            "faltam apenas as imagens."
        )

    identifiers = [record.input_id for record in report.records]
    if not identifiers:
        return
    chosen = st.selectbox("Estrutura", identifiers)
    record = next(r for r in report.records if r.input_id == chosen)

    original, curated = st.columns(2)
    with original:
        st.caption("Original")
        image = draw(record.raw_smiles)
        if image is not None:
            st.image(image)
        st.code(record.raw_smiles, language="text")
    with curated:
        st.caption("Curada" if record.passed else "Rejeitada")
        image = draw(record.curated_smiles or "")
        if image is not None:
            st.image(image)
        st.code(record.curated_smiles or "—", language="text")

    if not record.passed:
        st.error(
            f"**{record.rejection_code.value}** em `{record.rejection_stage.value}` — "
            f"{record.rejection_detail}"
        )

    steps = structure_lineage(record)
    if steps:
        st.markdown("**Trajetória observada**")
        for step in steps:
            st.markdown(
                f"`{step.stage}` · **{step.rule}** — {step.detail}  \n"
                f"<span class='mono'>{step.before_smiles} → {step.after_smiles}</span>",
                unsafe_allow_html=True,
            )
    else:
        st.caption(
            "Nenhuma transformação observada. O motor é uma caixa-preta de oito "
            "operações internas: etapas sem efeito detectado pelas sondas não "
            "produzem passos de linhagem, e a interface não os inventa."
        )

    if record.removed_fragments:
        st.markdown(f"**Fragmentos removidos:** `{record.removed_fragments}`")
        image = draw(record.removed_fragments, size=180)
        if image is not None:
            st.image(image)

    if record.inchikey:
        st.markdown(
            f"**InChIKey** `{record.inchikey}` · bloco1 `{record.inchikey_block1}`"
        )


# --- Proveniência e exportação --------------------------------------------------------


def render_provenance(report: RunReport) -> None:
    provenance = report.provenance
    blockers = provenance.reproduction_blockers()

    if blockers:
        st.warning(
            "**Esta execução não é exatamente reproduzível.**\n\n"
            + "\n".join(f"- {item}" for item in blockers)
        )
    else:
        st.success("Execução reproduzível a partir dos metadados registrados.")

    left, right = st.columns(2)
    with left:
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
    with right:
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

    st.markdown("**Versões e ambiente**")
    st.json({**provenance.versions, **provenance.environment})

    st.markdown("**Reproduce this run**")
    st.code(provenance.reproduction_command(provenance.input_name), language="bash")
    st.caption(
        "Mudanças de versão do RDKit alteram percepção de aromaticidade e regras de "
        "padronização: dois lotes com RDKit diferente não são comparáveis mesmo com "
        "o mesmo policy_hash."
    )


def render_exports(report: RunReport) -> None:
    run_id = report.provenance.run_id
    complete, custom, rejected_tab, manifest_tab, package_tab = st.tabs(
        ["Completo", "Personalizado", "Rejeitados", "Manifesto", "Pacote"]
    )

    with complete:
        st.caption("Todos os registros, aprovados e rejeitados, com status e motivo.")
        st.download_button(
            "structures_full.csv", full_csv(report),
            file_name="structures_full.csv", mime="text/csv",
        )

    with custom:
        columns = st.multiselect(
            "Colunas", list(EXPORTABLE_COLUMNS), list(DEFAULT_EXPORT_COLUMNS)
        )
        st.caption(
            "`inchi` e `molecular_formula` são calculadas sob demanda. Colunas que a "
            "camada de engenharia não produz não aparecem nesta lista."
        )
        if columns:
            st.download_button(
                "structures_custom.csv", to_csv(report.records, columns),
                file_name="structures_custom.csv", mime="text/csv",
            )

    with rejected_tab:
        st.caption("Identificador, SMILES, etapa de falha, código e detalhe.")
        st.download_button(
            "rejected_structures.csv", rejected_csv(report),
            file_name="rejected_structures.csv", mime="text/csv",
        )

    with manifest_tab:
        manifest = run_manifest(report)
        st.json(manifest, expanded=False)
        st.download_button(
            "run_manifest.json", json.dumps(manifest, indent=2, ensure_ascii=False),
            file_name="run_manifest.json", mime="application/json",
        )

    with package_tab:
        st.caption(
            "CSV completo, rejeitados, manifesto, parâmetros e um README explicando "
            "como interpretar e reproduzir."
        )
        st.download_button(
            f"structure-curation-run-{run_id}.zip",
            reproducibility_package(report),
            file_name=f"structure-curation-run-{run_id}.zip",
            mime="application/zip",
        )


def render_results_table(report: RunReport) -> None:
    scope = st.radio(
        "Escopo", ["Todos", "Aprovados", "Rejeitados", "Com transformações"],
        horizontal=True, label_visibility="collapsed",
    )
    selected = {
        "Todos": report.records,
        "Aprovados": report.approved,
        "Rejeitados": report.rejected,
        "Com transformações": [r for r in report.approved if r.transformations],
    }[scope]

    query = st.text_input("Buscar (id, SMILES ou InChIKey)", "")
    if query:
        needle = query.lower()
        selected = [
            record for record in selected
            if needle in record.input_id.lower()
            or needle in record.raw_smiles.lower()
            or needle in (record.curated_smiles or "").lower()
            or needle in (record.inchikey or "").lower()
        ]

    st.caption(f"{len(selected)} registros")
    st.dataframe(
        [record_row(record, (
            "input_id", "raw_smiles", "curated_smiles", "inchikey", "status",
            "rejection_code", "removed_fragments", "parent_mw", "n_components_parent",
        )) for record in selected],
        use_container_width=True, hide_index=True,
    )


# --- Aplicação -----------------------------------------------------------------------


def main() -> None:
    st.set_page_config(page_title="Structure Curation Pipeline", layout="wide")
    inject_styles()

    parameters = sidebar()
    report: Optional[RunReport] = st.session_state.get("report")

    status = (
        f"RUN {report.provenance.run_id} · ✓ SUCCESS" if report else "nenhuma execução"
    )
    st.markdown(
        f"<div class='run-banner'><span class='run-title'>STRUCTURE CURATION "
        f"PIPELINE</span><span class='mono'>{status}</span></div>",
        unsafe_allow_html=True,
    )

    raw, name = input_area()
    if raw:
        render_preview(raw, name)
        if st.button("Executar pipeline", type="primary"):
            pipeline = CurationPipeline(
                policy_hash=parameters["policy_hash"],
                criteria=EligibilityCriteria(
                    max_molecular_weight=parameters["max_mw"],
                    max_heavy_atoms=int(parameters["max_ha"]),
                ),
                deduplicate=parameters["deduplicate"],
            )
            with st.spinner("Executando…"):
                st.session_state["report"] = pipeline.run_report(
                    raw.decode("utf-8", errors="replace"),
                    parameters={
                        "max_mw": parameters["max_mw"],
                        "max_ha": int(parameters["max_ha"]),
                        "deduplicate": parameters["deduplicate"],
                    },
                    input_bytes=raw,
                    input_name=name,
                    policy_path=parameters["policy_path"],
                )
            st.rerun()

    if report is None:
        st.info("Carregue um arquivo ou cole SMILES para executar o pipeline.")
        return

    st.divider()
    kpis = report.kpis()
    for column, (label, value) in zip(st.columns(len(kpis)), kpis.items()):
        column.metric(label, value)

    st.divider()
    selected_stage = render_funnel(report)
    st.divider()
    render_stage_detail(report, selected_stage)

    st.divider()
    exclusions, lineage, results, provenance_section, exports = st.tabs(
        ["Exclusões", "Rastreabilidade", "Resultados", "Proveniência", "Exportação"]
    )
    with exclusions:
        render_exclusions(report)
    with lineage:
        render_lineage(report)
    with results:
        render_results_table(report)
    with provenance_section:
        render_provenance(report)
    with exports:
        render_exports(report)


if __name__ == "__main__":
    main()
