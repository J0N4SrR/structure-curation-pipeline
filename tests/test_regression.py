"""Regressão dos modos de falha empíricos, através do orquestrador.

Cada caso aqui trava um comportamento medido contra a implementação de referência.
Diferente de ``test_engine.py``, que exercita o wrapper, esta suíte entra pela API
pública do pipeline — é o contrato que a CLI e a interface vão consumir.

Todos os valores esperados foram obtidos executando o ``chembl_structure_pipeline``.
Nenhum foi escrito de memória: o teste original de benzoato de sódio, redigido a
partir da intuição de que o sódio sairia, estava errado.
"""

from __future__ import annotations

import pytest

from curation import CurationPipeline
from curation.models import RejectionCode, Stage

POLICY_HASH = "0" * 64

COTRIMOXAZOLE = (
    "COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC.Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1"
)


@pytest.fixture(scope="module")
def pipeline() -> CurationPipeline:
    return CurationPipeline(policy_hash=POLICY_HASH)


# --- 1. Precedência iônica (D-01) -------------------------------------------------


def test_covalent_sodium_acetate_keeps_ionic_precedence(
    pipeline: CurationPipeline,
) -> None:
    """``normalize`` antes de ``desalt``.

    ``GetMolFrags`` vê ``CC(=O)O[Na]`` como fragmento único. Com a ordem invertida o
    sódio permaneceria covalente na saída, aprovado e sem sinal de erro.
    """
    record = pipeline.process_single("CC(=O)O[Na]", "REG_001")

    assert record.passed
    assert record.curated_smiles == "CC(=O)[O-].[Na+]"
    assert record.n_components_parent == 2
    assert record.has_rule("alkali_metal_ionized")
    assert not record.has_rule("alkali_metal_residual")
    assert record.removed_fragments is None


# --- 2 e 3. Portão de valência posterior (D-09) -----------------------------------


def test_neutral_quaternary_ammonium_is_repaired(pipeline: CurationPipeline) -> None:
    """Sanitização estrita antecipada levantaria ``AtomValenceException``."""
    record = pipeline.process_single("C[N](C)(C)C", "REG_002")

    assert record.passed
    assert record.curated_smiles == "C[N+](C)(C)C"
    assert record.delta_net_charge == 1


def test_neutral_diazonium_is_repaired(pipeline: CurationPipeline) -> None:
    record = pipeline.process_single("c1ccccc1[N]#N", "REG_003")

    assert record.passed
    assert record.curated_smiles == "N#[N+]c1ccccc1"


# --- 4. Multicomponentes (D-06) ----------------------------------------------------


def test_cotrimoxazole_keeps_both_active_components(
    pipeline: CurationPipeline,
) -> None:
    """``LargestFragmentChooser`` devolveria só a trimetoprima."""
    record = pipeline.process_single(COTRIMOXAZOLE, "REG_004")

    assert record.passed
    assert record.n_components_parent == 2
    assert record.removed_fragments is None
    assert "S(=O)(=O)" in record.curated_smiles


# --- 5. Guard "tudo é sal" (D-06) ---------------------------------------------------


def test_sodium_benzoate_is_left_intact(pipeline: CurationPipeline) -> None:
    """Benzoato e sódio constam ambos de ``salts.smi``.

    Quando todos os componentes casam com a tabela, a referência não remove nada.
    O valor esperado é o medido, não o intuitivo.
    """
    record = pipeline.process_single("[Na+].[O-]C(=O)c1ccccc1", "REG_005")

    assert record.passed
    assert record.curated_smiles == "O=C([O-])c1ccccc1.[Na+]"
    assert record.n_components_parent == 2
    assert record.removed_fragments is None
    assert not record.has_rule("fragments_removed")


# --- 6. Dessalificação com neutralização ---------------------------------------------


def test_alanine_hydrochloride_yields_neutral_amino_acid(
    pipeline: CurationPipeline,
) -> None:
    """Sal removido, estrutura-mãe neutra e quiralidade preservada."""
    record = pipeline.process_single("N[C@@H](C)C(=O)O.Cl", "REG_006")

    assert record.passed
    assert record.curated_smiles == "C[C@H](N)C(=O)O"
    assert record.removed_fragments == "Cl"
    assert record.has_rule("fragments_removed"), "equivalente a flag_salt_removed=True"
    assert record.n_components_parent == 1
    assert record.delta_stereocenters == 0
    assert record.inchikey == "QNAYBMKLOCPYGJ-REOHCLBHSA-N"


# --- 7. Rejeição legítima --------------------------------------------------------


def test_pentavalent_carbon_is_rejected(pipeline: CurationPipeline) -> None:
    record = pipeline.process_single("C(C)(C)(C)(C)C", "REG_007")

    assert not record.passed
    assert record.rejection_code is RejectionCode.ERR_VALENCE
    assert record.rejection_stage is not Stage.PARSE, "rejeicao so depois do motor"
    assert record.rejection_detail


# --- Portão de valência x exclude_flag (D-03) --------------------------------------


def test_valence_gate_does_not_override_exclusion_policy(
    pipeline: CurationPipeline,
) -> None:
    """Carborano de 8 boros: o portão não pode revogar a D-03.

    A referência pula a sanitização de compostos com ``exclude_flag`` ativo. Um
    portão incondicional descartaria esta estrutura, presente no corpus da própria
    referência sob a anotação ">7 Boron atoms".
    """
    carborane = (
        "CC1=CN([C@H]2C[C@@H](O)[C@@H](CO)O2)C(=O)N(CC234C56B78B29%10B32%11C43B54%12"
        "B765B896B%1027B3%114C5%1267)C1=O"
    )
    record = pipeline.process_single(carborane, "REG_008")

    assert record.excluded_flag, "estrutura deve acionar o exclude_flag"
    assert record.passed, "composto excluido nao pode cair no portao de valencia"
    assert record.rejection_code is None


def test_cisplatin_is_preserved(pipeline: CurationPipeline) -> None:
    record = pipeline.process_single("[Pt](Cl)(Cl)(N)N", "REG_009")

    assert record.passed
    assert record.excluded_flag
    assert record.has_rule("standardization_skipped_excluded")


# --- Barreira global -------------------------------------------------------------


@pytest.mark.parametrize(
    "raw_smiles",
    ["", "   ", "C1CC", "[[[", "nao e um smiles", "C" * 5000, "\x00", "C(", "]["],
)
def test_no_input_can_break_the_batch(
    pipeline: CurationPipeline, raw_smiles: str
) -> None:
    record = pipeline.process_single(raw_smiles, "REG_HOSTILE")
    assert record.status in ("PASSED", "REJECTED")
    if not record.passed:
        assert record.rejection_code is not None


def test_engine_failure_is_contained_by_the_barrier() -> None:
    """Um motor com defeito degrada um registro, nunca derruba o lote."""

    class BrokenEngine:
        _rdkit_version = "x"
        _csp_version = "y"

        def curate(self, input_id: str, raw_smiles: str):
            raise RuntimeError("motor com defeito")

    pipeline = CurationPipeline(POLICY_HASH, engine=BrokenEngine())
    record = pipeline.process_single("CCO", "REG_010")

    assert not record.passed
    assert record.rejection_code is RejectionCode.ERR_INTERNAL
    assert "motor com defeito" in record.rejection_detail


# --- Idempotência estrita ----------------------------------------------------------

CORPUS = (
    "CC(=O)O[Na]",
    "CCO[Na]",
    "C[N](C)(C)C",
    "c1ccccc1[N]#N",
    COTRIMOXAZOLE,
    "[Na+].[O-]C(=O)c1ccccc1",
    "N[C@@H](C)C(=O)O.Cl",
    "OC(=O)[C@H](O)[C@@H](O)C(=O)O",
    "CC(=O)Oc1ccccc1C(=O)O",
    "CS(=O)C",
    "C[N+](C)(C)CCC[O-]",
    "CC(=O)O.CC(=O)O",
    "[Pt](Cl)(Cl)(N)N",
    "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
)


@pytest.mark.parametrize("raw_smiles", CORPUS)
def test_strict_idempotency(pipeline: CurationPipeline, raw_smiles: str) -> None:
    """f(x) == f(f(x)), em SMILES e em InChIKey."""
    once = pipeline.process_single(raw_smiles, "IDEM")
    assert once.passed, raw_smiles

    twice = pipeline.process_single(once.curated_smiles, "IDEM")
    assert twice.passed
    assert twice.curated_smiles == once.curated_smiles
    assert twice.inchikey == once.inchikey


def test_idempotency_holds_for_transformation_events(
    pipeline: CurationPipeline,
) -> None:
    """A segunda passagem não pode disparar transformações novas."""
    once = pipeline.process_single("CC(=O)O[Na]", "IDEM")
    twice = pipeline.process_single(once.curated_smiles, "IDEM")
    assert not twice.has_rule("alkali_metal_ionized")


# --- Lote de ponta a ponta ---------------------------------------------------------


def test_run_produces_atomic_output_and_summary(tmp_path) -> None:
    pipeline = CurationPipeline(policy_hash=POLICY_HASH)
    source = "\n".join(["CC(=O)O[Na]", "N[C@@H](C)C(=O)O.Cl", "C(C)(C)(C)(C)C"])

    summary = pipeline.run(source, tmp_path)

    assert summary.total == 3
    assert summary.passed == 2
    assert summary.rejected == 1
    assert summary.rejection_counts == {"ERR_VALENCE": 1}
    assert 0.0 < summary.pass_rate < 1.0
    assert "aprovados" in summary.format_report()

    from curation.io import read_manifest

    manifest = read_manifest(tmp_path)
    assert manifest is not None
    assert manifest["counts"]["total"] == 3
    assert not list(tmp_path.glob("*.tmp"))


def test_process_many_is_lazy(pipeline: CurationPipeline) -> None:
    consumed: list[str] = []

    def source():
        for index in range(500):
            consumed.append(str(index))
            yield (f"CMPD_{index}", "CCO")

    stream = pipeline.process_many(source())
    next(stream)
    assert len(consumed) < 500


def test_summarize_matches_run_counts(pipeline: CurationPipeline) -> None:
    records = [
        pipeline.process_single(smiles, f"S{index}")
        for index, smiles in enumerate(["CCO", "C(C)(C)(C)(C)C"])
    ]
    summary = CurationPipeline.summarize(records)
    assert (summary.total, summary.passed, summary.rejected) == (2, 1, 1)
