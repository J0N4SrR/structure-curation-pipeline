"""Testes da camada de relatório e proveniência.

Cobrem o que o app exibe. O app em si é apresentação e não é testado aqui —
Streamlit não está instalado neste ambiente.
"""

from __future__ import annotations

import io
import json
import zipfile

import pytest

from curation import CurationPipeline
from curation.provenance import GitState, RunProvenance, sha256_of
from curation.reporting import (
    DEDUP_STAGE,
    EXPORTABLE_COLUMNS,
    StageStatus,
    build_stage_reports,
    full_csv,
    package_readme,
    parameters_yaml,
    record_row,
    rejected_csv,
    reproducibility_package,
    run_manifest,
    structure_lineage,
    to_csv,
)

POLICY_HASH = "0" * 64
CORPUS = (
    "CC(=O)O[Na]\nN[C@@H](C)C(=O)O.Cl\nCCO\nOCC\nN[C@H](C)C(=O)O\n"
    "C(C)(C)(C)(C)C\nC1CC\n[Pt](Cl)(Cl)(N)N\n"
)


@pytest.fixture(scope="module")
def report():
    pipeline = CurationPipeline(POLICY_HASH)
    return pipeline.run_report(
        CORPUS, input_bytes=CORPUS.encode(), input_name="corpus.smi"
    )


# --- Funil ---------------------------------------------------------------------------


def test_funnel_is_continuous(report) -> None:
    """A saída de um estágio é a entrada do seguinte — sem buracos."""
    stages = report.stages
    assert stages[0].n_input == report.total
    for previous, following in zip(stages, stages[1:]):
        assert following.n_input == previous.n_output


def test_every_drop_is_explained(report) -> None:
    """Nenhuma estrutura desaparece silenciosamente."""
    dropped = sum(stage.n_excluded for stage in report.stages)
    assert dropped == len(report.rejected)

    accounted = sum(
        group.count
        for group in report.exclusion_groups()
        if group.stage != DEDUP_STAGE
    )
    assert accounted == dropped


def test_stage_status_reflects_exclusions(report) -> None:
    for stage in report.stages:
        expected = StageStatus.WARNING if stage.n_excluded else StageStatus.SUCCESS
        assert stage.status is expected, stage.name


def test_durations_are_measured_not_estimated(report) -> None:
    assert sum(stage.duration_seconds for stage in report.stages) > 0
    assert all(stage.duration_seconds >= 0 for stage in report.stages)


def test_empty_run_produces_skipped_stages() -> None:
    for stage in build_stage_reports([]):
        assert stage.status is StageStatus.SKIPPED


def test_dedup_stage_excludes_only_real_duplicates(report) -> None:
    """Colisões de bloco1 não reduzem a contagem: não são duplicatas (D-07)."""
    assert report.dedup is not None
    assert report.dedup.n_excluded == report.n_duplicates
    assert "BLOCK1_COLLISION" in report.dedup.rejection_counts


# --- KPIs -----------------------------------------------------------------------------


def test_kpis_are_internally_consistent(report) -> None:
    kpis = report.kpis()
    assert kpis["Approved"] + kpis["Rejected"] == kpis["Total Processed"]
    assert kpis["Warnings"] <= kpis["Approved"]


# --- Linhagem ---------------------------------------------------------------------------


def test_lineage_comes_from_observed_events(report) -> None:
    record = next(r for r in report.records if r.raw_smiles == "CC(=O)O[Na]")
    steps = structure_lineage(record)

    assert steps
    assert {step.rule for step in steps} <= {e.rule for e in record.transformations}
    assert all(step.before_smiles and step.after_smiles for step in steps)


def test_lineage_is_empty_when_nothing_was_observed(report) -> None:
    """Sem evento observado não há passo inventado."""
    record = next(r for r in report.records if r.raw_smiles == "CCO")
    assert structure_lineage(record) == []


# --- Exportação ---------------------------------------------------------------------------


def test_full_csv_includes_rejected_records(report) -> None:
    text = full_csv(report)
    assert "C(C)(C)(C)(C)C" in text
    assert "ERR_VALENCE" in text
    assert len(text.splitlines()) == report.total + 1


def test_custom_csv_honours_column_selection(report) -> None:
    text = to_csv(report.records, ("input_id", "inchikey"))
    assert text.splitlines()[0] == "input_id,inchikey"


def test_derived_columns_are_computed(report) -> None:
    record = next(r for r in report.approved if r.curated_smiles)
    row = record_row(record, ("inchi", "molecular_formula"))
    assert row["molecular_formula"]
    assert row["inchi"].startswith("InChI=")


def test_exportable_columns_only_expose_real_data(report) -> None:
    """Colunas que a engenharia não produz não podem ser oferecidas."""
    assert "tautomer_status" not in EXPORTABLE_COLUMNS
    row = record_row(report.records[0], EXPORTABLE_COLUMNS)
    assert set(row) == set(EXPORTABLE_COLUMNS)


def test_rejected_csv_carries_stage_and_reason(report) -> None:
    text = rejected_csv(report)
    assert "rejection_stage" in text and "rejection_code" in text
    assert len(text.splitlines()) == len(report.rejected) + 1


# --- Manifesto e pacote ---------------------------------------------------------------------


def test_manifest_reconciles_with_records(report) -> None:
    manifest = run_manifest(report)
    assert manifest["metrics"]["Total Processed"] == report.total
    assert manifest["run_id"] == report.provenance.run_id
    assert len(manifest["stages"]) == len(report.stages) + 1
    assert manifest["provenance"]["input_hash"] == sha256_of(CORPUS.encode())


def test_package_contains_all_five_artifacts(report) -> None:
    archive = zipfile.ZipFile(io.BytesIO(reproducibility_package(report)))
    names = {name.split("/")[-1] for name in archive.namelist()}
    assert names == {
        "structures_full.csv",
        "structures_rejected.csv",
        "run_manifest.json",
        "parameters.yaml",
        "README.md",
    }


def test_package_readme_states_reproducibility_verdict(report) -> None:
    text = package_readme(report)
    assert report.provenance.run_id in text
    assert "Reproducao" in text
    assert ("NAO e exatamente reproduzivel" in text) or ("e reproduzivel" in text)


def test_parameters_yaml_records_versions(report) -> None:
    text = parameters_yaml(report)
    assert "policy_hash:" in text and "rdkit:" in text


# --- Proveniência ---------------------------------------------------------------------------


def test_output_hash_is_set_after_finalize(report) -> None:
    assert report.provenance.output_hash
    assert report.provenance.duration_seconds > 0


def test_dirty_tree_blocks_reproducibility_claim() -> None:
    dirty = RunProvenance(
        run_id="r", started_at="t", policy_hash="a" * 64,
        git=GitState(commit="abc123", branch="main", dirty=True),
    )
    assert not dirty.reproducible
    assert any("nao commitadas" in b or "não commitadas" in b
               for b in dirty.reproduction_blockers())


def test_unversioned_policy_blocks_reproducibility_claim() -> None:
    provenance = RunProvenance(
        run_id="r", started_at="t", policy_hash="UNVERSIONED_POLICY",
        git=GitState(commit="abc", branch="main", dirty=False),
    )
    assert not provenance.reproducible
    assert provenance.reproduction_blockers()


def test_clean_state_is_reproducible() -> None:
    provenance = RunProvenance(
        run_id="r", started_at="t", policy_hash="a" * 64,
        git=GitState(commit="abc", branch="main", dirty=False),
    )
    assert provenance.reproducible
    assert provenance.reproduction_blockers() == []


def test_reproduction_command_includes_parameters() -> None:
    provenance = RunProvenance(
        run_id="r", started_at="t", policy_hash="a" * 64,
        parameters={"max_mw": 500.0, "max_ha": 50, "deduplicate": False},
        git=GitState(commit="abc123", branch="main", dirty=False),
        policy_path="docs/decisions.md",
    )
    command = provenance.reproduction_command("data/in.smi")
    for fragment in ("git checkout abc123", "--max-mw 500.0", "--max-ha 50",
                     "--no-dedup", "--decisions docs/decisions.md"):
        assert fragment in command


def test_run_id_is_deterministic_for_same_input() -> None:
    """O sufixo do run_id vem do hash da entrada, não de aleatoriedade."""
    data = b"CCO\n"
    a = RunProvenance.start("h", {}, data)
    b = RunProvenance.start("h", {}, data)
    assert a.run_id.split("-")[1] == b.run_id.split("-")[1]


# --- Regressões de renderização e alinhamento ----------------------------------------


def test_numeric_columns_stay_numeric_for_arrow(report) -> None:
    """``None`` não pode virar string: Arrow recusa colunas mistas.

    Um registro rejeitado não tem peso molecular. Preencher com ``""`` produzia
    uma coluna de ``float`` e ``str`` e quebrava a renderização da tabela inteira.
    """
    import pandas as pd
    import pyarrow as pa

    frame = pd.DataFrame(
        [
            record_row(record, ("input_id", "parent_mw", "parent_heavy_atoms"))
            for record in report.records
        ]
    )
    assert frame["parent_mw"].dtype.kind == "f"
    pa.Table.from_pandas(frame)


def test_csv_still_writes_empty_cells_for_missing_values(report) -> None:
    """Preservar ``None`` não pode mudar a saída em CSV."""
    text = to_csv(report.rejected, ("input_id", "curated_smiles", "parent_mw"))
    body = text.splitlines()[1]
    assert body.endswith(",,"), body


# --- Cobertura das funções auxiliares ---------------------------------------------


def test_largest_reduction_points_at_the_worst_stage(report) -> None:
    stage = report.largest_reduction()
    assert stage is not None
    worst = max(
        (s.n_excluded for s in report.stages + ([report.dedup] if report.dedup else [])),
    )
    assert stage.n_excluded == worst


def test_largest_reduction_is_none_when_nothing_was_removed() -> None:
    """A frase de resumo não pode aparecer quando não houve redução."""
    pipeline = CurationPipeline(POLICY_HASH)
    clean = "CCO\nCCN\nCC(=O)Oc1ccccc1C(=O)O\n"
    clean_report = pipeline.run_report(clean, input_bytes=clean.encode())
    assert clean_report.largest_reduction() is None


def test_progress_callback_reports_real_positions() -> None:
    """O progresso vem do pipeline, não é estimado pela interface."""
    pipeline = CurationPipeline(POLICY_HASH)
    seen: list[tuple[int, int, str]] = []
    source = "CCO\nCCN\nC(C)(C)(C)(C)C\n"

    pipeline.run_report(
        source, input_bytes=source.encode(),
        progress=lambda position, total, identifier: seen.append(
            (position, total, identifier)
        ),
    )

    assert [position for position, _, _ in seen] == [1, 2, 3]
    assert {total for _, total, _ in seen} == {3}


def test_progress_is_optional() -> None:
    pipeline = CurationPipeline(POLICY_HASH)
    assert pipeline.run_report("CCO\n", input_bytes=b"CCO\n").total == 1


def test_stage_timings_reset_between_runs() -> None:
    """Tempos acumulados de um lote não podem vazar para o seguinte."""
    pipeline = CurationPipeline(POLICY_HASH)
    long_source = "\n".join(["CC(=O)Oc1ccccc1C(=O)O"] * 40)
    short_source = "CCO"

    pipeline.run_report(long_source, input_bytes=long_source.encode())
    second = pipeline.run_report(short_source, input_bytes=short_source.encode())

    total = sum(stage.duration_seconds for stage in second.stages)
    assert total > 0
    assert second.total == 1
