"""Testes do harness de validação científica (Fase 5)."""

from __future__ import annotations

from pathlib import Path

import pytest

from curation import CurationPipeline
from curation.validation.ablation import (
    canonicalize_tautomer,
    default_variants,
    run_ablation,
    standardize_without_tartrate_flattening,
    to_markdown,
    write_report,
)
from curation.validation.benchmarking import (
    ConfusionMatrix,
    GoldStandardEntry,
    GoldStandardError,
    agreement_report,
    benchmark_report,
    cohens_kappa,
    evaluate_pipeline,
    load_gold_standard,
    write_annotation_template,
)

POLICY_HASH = "0" * 64


# --- Kappa de Cohen ------------------------------------------------------------------


def test_kappa_perfect_agreement() -> None:
    assert cohens_kappa(["ACCEPT", "REJECT"] * 5, ["ACCEPT", "REJECT"] * 5) == 1.0


def test_kappa_total_disagreement_is_negative() -> None:
    assert cohens_kappa(["ACCEPT", "REJECT"] * 4, ["REJECT", "ACCEPT"] * 4) < 0


def test_kappa_single_category_returns_one() -> None:
    """Pe = 1 torna κ indefinido; a convenção adotada é 1.0."""
    assert cohens_kappa(["ACCEPT"] * 5, ["ACCEPT"] * 5) == 1.0


def test_kappa_discounts_chance_agreement() -> None:
    """Concordância alta com classes desbalanceadas produz κ baixo."""
    a = ["ACCEPT"] * 18 + ["REJECT"] * 2
    b = ["ACCEPT"] * 19 + ["REJECT"]
    kappa = cohens_kappa(a, b)
    observed = sum(x == y for x, y in zip(a, b)) / len(a)
    assert observed > 0.9
    assert kappa < observed


def test_kappa_rejects_mismatched_lengths() -> None:
    with pytest.raises(ValueError):
        cohens_kappa(["ACCEPT"], ["ACCEPT", "REJECT"])


# --- Carregamento do padrão-ouro -------------------------------------------------------


def write_gold(path: Path, rows: list[str]) -> Path:
    header = "compound_id,raw_smiles,annotator_a,annotator_b,consensus,reason"
    path.write_text("\n".join([header, *rows]) + "\n", encoding="utf-8")
    return path


def test_load_gold_standard(tmp_path: Path) -> None:
    source = write_gold(
        tmp_path / "gold.csv",
        [
            "G1,CCO,ACCEPT,ACCEPT,,limpo",
            "G2,C(C)(C)(C)(C)C,REJECT,REJECT,,valencia",
            "G3,CC(=O)O[Na],ACCEPT,REJECT,ACCEPT,resolvido em revisao",
        ],
    )
    entries = load_gold_standard(source)

    assert len(entries) == 3
    assert entries[0].annotators_agree
    assert entries[2].resolved_label == "ACCEPT"


def test_unresolved_disagreement_has_no_truth(tmp_path: Path) -> None:
    source = write_gold(tmp_path / "g.csv", ["G1,CCO,ACCEPT,REJECT,,em aberto"])
    assert load_gold_standard(source)[0].resolved_label is None


def test_invalid_label_is_rejected(tmp_path: Path) -> None:
    source = write_gold(tmp_path / "g.csv", ["G1,CCO,TALVEZ,ACCEPT,,"])
    with pytest.raises(GoldStandardError, match="ACCEPT ou REJECT"):
        load_gold_standard(source)


def test_missing_annotator_is_rejected(tmp_path: Path) -> None:
    source = write_gold(tmp_path / "g.csv", ["G1,CCO,ACCEPT,,,"])
    with pytest.raises(GoldStandardError, match="ambos os anotadores"):
        load_gold_standard(source)


def test_missing_columns_are_reported(tmp_path: Path) -> None:
    (tmp_path / "g.csv").write_text("compound_id,raw_smiles\nG1,CCO\n", encoding="utf-8")
    with pytest.raises(GoldStandardError, match="colunas ausentes"):
        load_gold_standard(tmp_path / "g.csv")


def test_annotation_template_round_trips(tmp_path: Path) -> None:
    target = tmp_path / "template.csv"
    assert write_annotation_template([("A1", "CCO"), ("A2", "CCN")], target) == 2

    text = target.read_text(encoding="utf-8")
    assert "annotator_a" in text and "A1,CCO" in text


# --- Matriz de confusão -----------------------------------------------------------------


def test_confusion_matrix_metrics() -> None:
    matrix = ConfusionMatrix(
        true_positive=8, false_positive=2, true_negative=88, false_negative=2
    )
    assert matrix.precision == pytest.approx(0.8)
    assert matrix.recall == pytest.approx(0.8)
    assert matrix.f1 == pytest.approx(0.8)
    assert matrix.accuracy == pytest.approx(0.96)


def test_false_rejection_rate_denominator_is_expert_accepted() -> None:
    """FRR = FP / (FP + VN): fração dos aceitáveis que o pipeline descartou."""
    matrix = ConfusionMatrix(
        true_positive=10, false_positive=5, true_negative=95, false_negative=0
    )
    assert matrix.false_rejection_rate == pytest.approx(5 / 100)


def test_metrics_are_zero_safe() -> None:
    empty = ConfusionMatrix(0, 0, 0, 0)
    assert (empty.precision, empty.recall, empty.f1, empty.false_rejection_rate) == (
        0.0, 0.0, 0.0, 0.0,
    )


def test_evaluate_pipeline_against_gold_standard() -> None:
    pipeline = CurationPipeline(POLICY_HASH)
    entries = [
        GoldStandardEntry("G1", "CCO", "ACCEPT", "ACCEPT"),
        GoldStandardEntry("G2", "C(C)(C)(C)(C)C", "REJECT", "REJECT"),
        GoldStandardEntry("G3", "C1CC", "REJECT", "REJECT"),
        GoldStandardEntry("G4", "CC(=O)Oc1ccccc1C(=O)O", "ACCEPT", "ACCEPT"),
    ]
    matrix = evaluate_pipeline(entries, pipeline)

    assert (matrix.true_positive, matrix.true_negative) == (2, 2)
    assert (matrix.false_positive, matrix.false_negative) == (0, 0)
    assert matrix.false_rejection_rate == 0.0


def test_unresolved_entries_are_excluded_from_the_matrix() -> None:
    """Sem verdade acordada não há contra o que comparar."""
    pipeline = CurationPipeline(POLICY_HASH)
    entries = [
        GoldStandardEntry("G1", "CCO", "ACCEPT", "ACCEPT"),
        GoldStandardEntry("G2", "CCN", "ACCEPT", "REJECT"),
    ]
    matrix = evaluate_pipeline(entries, pipeline)

    assert matrix.total == 1
    assert matrix.excluded_unresolved == 1


def test_agreement_report_counts_disagreements() -> None:
    entries = [
        GoldStandardEntry("G1", "CCO", "ACCEPT", "ACCEPT"),
        GoldStandardEntry("G2", "CCN", "ACCEPT", "REJECT"),
        GoldStandardEntry("G3", "CCC", "REJECT", "ACCEPT", consensus="REJECT"),
    ]
    report = agreement_report(entries)

    assert (report.n, report.disagreements, report.unresolved) == (3, 2, 1)
    assert report.interpretation


def test_benchmark_report_warns_on_low_kappa() -> None:
    pipeline = CurationPipeline(POLICY_HASH)
    entries = [
        GoldStandardEntry(f"G{i}", "CCO", "ACCEPT", "REJECT" if i % 2 else "ACCEPT")
        for i in range(6)
    ]
    report = benchmark_report(entries, pipeline)

    assert "κ de Cohen" in report
    assert "concordância inter-anotador" in report.lower()
    assert "Atenção" in report


# --- Ablação ------------------------------------------------------------------------------

CORPUS = [
    ("C1", "OC(=O)[C@H](O)[C@@H](O)C(=O)O"),
    ("C2", "N[C@@H](C)C(=O)O"),
    ("C3", "[Pt](Cl)(Cl)(N)N"),
    ("C4", "CC(=O)Oc1ccccc1C(=O)O"),
    ("C5", "CC(=O)O[Na]"),
]


def test_default_variants_cover_the_open_decisions() -> None:
    adrs = {variant.adr for variant in default_variants()}
    assert {"D-02", "D-03", "D-04"} <= adrs


def test_tartrate_variant_preserves_stereochemistry() -> None:
    """Ablação da D-02: o efeito precisa ser mensurável, não apenas declarado."""
    baseline, no_flatten = run_ablation(
        CORPUS, [v for v in default_variants() if v.name in ("baseline", "no_tartrate_flattening")]
    )
    assert no_flatten.defined_stereocenters > baseline.defined_stereocenters


def test_alternative_standardizer_matches_reference_without_tartrate() -> None:
    """A sequência recomposta só pode divergir no passo omitido.

    Se a referência mudar a ordem interna, esta asserção quebra - que é o objetivo.
    """
    import chembl_structure_pipeline as csp
    from rdkit import Chem

    for smiles in ("CC(=O)Oc1ccccc1C(=O)O", "N[C@@H](C)C(=O)O", "CC(=O)O[Na]"):
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        reference = Chem.MolToSmiles(csp.standardize_mol(Chem.Mol(mol)))
        recomposed = Chem.MolToSmiles(
            standardize_without_tartrate_flattening(Chem.Mol(mol))
        )
        assert reference == recomposed, smiles


def test_tautomer_canonicalization_removes_sp3_stereo() -> None:
    """``RemoveSp3Stereo`` é True por padrão - a base empírica da D-04."""
    from rdkit import Chem

    canonical = canonicalize_tautomer(Chem.MolFromSmiles("N[C@@H](C)C(=O)O"))
    assert "@" not in Chem.MolToSmiles(canonical)


def test_organometallic_variant_discards_excluded() -> None:
    baseline, discard = run_ablation(
        CORPUS,
        [v for v in default_variants() if v.name in ("baseline", "discard_organometallics")],
    )
    assert discard.n_passed < baseline.n_passed
    assert "ERR_ORGANOMETALLIC" in discard.rejection_counts


def test_markdown_report_contrasts_against_baseline() -> None:
    report = to_markdown(run_ablation(CORPUS))

    for heading in ("Aprovação e deduplicação", "Distribuição físico-química",
                    "Conservação de centros quirais", "Leitura"):
        assert heading in report
    assert "baseline" in report


def test_write_report_emits_markdown_and_figures(tmp_path: Path) -> None:
    path = write_report(run_ablation(CORPUS), tmp_path)

    assert path.is_file() and path.name == "ablation.md"
    assert list(tmp_path.glob("*.png")), "figuras devem ser geradas"
    assert "Figuras" in path.read_text(encoding="utf-8")


def test_variant_metrics_are_self_consistent() -> None:
    for result in run_ablation(CORPUS):
        assert result.n_passed + result.n_rejected == result.n_input
        assert result.n_unique <= result.n_passed
        assert 0.0 <= result.pass_rate <= 1.0
        assert 0.0 <= result.dedup_rate <= 1.0
