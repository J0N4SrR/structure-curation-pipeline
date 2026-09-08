"""Regressão do wrapper do motor.

Cada teste trava um modo de falha estabelecido empiricamente contra a implementação
de referência. Os valores esperados foram obtidos executando o
``chembl_structure_pipeline``, nunca escritos de memória - o atalho oposto foi o que
produziu o teste incorreto de benzoato de sódio no roadmap original.
"""

from __future__ import annotations

import pytest
from rdkit import Chem

from curation import EngineWrapper, Stage
from curation.models import RejectionCode

POLICY_HASH = "0" * 64


@pytest.fixture(scope="module")
def engine() -> EngineWrapper:
    return EngineWrapper(policy_hash=POLICY_HASH)


# --- D-01: precedência normalize -> desalt --------------------------------------


def test_covalent_alkali_metal_is_ionized(engine: EngineWrapper) -> None:
    """Trava a inversão de ordem que corromperia dados silenciosamente.

    ``GetMolFrags`` enxerga ``CC(=O)O[Na]`` como fragmento único. Se o stripping
    rodasse antes da normalização, o sódio atravessaria ligado covalentemente, com
    ``status=PASSED`` e nenhum sinal de erro.
    """
    record = engine.curate("T", "CC(=O)O[Na]")

    assert record.passed
    assert record.curated_smiles == "CC(=O)[O-].[Na+]"
    assert record.has_rule("alkali_metal_ionized")
    assert record.n_components_parent == 2


def test_no_residual_covalent_alkali_metal(engine: EngineWrapper) -> None:
    """A sonda residual é o alarme direto para a regressão da D-01."""
    record = engine.curate("T", "CC(=O)O[Na]")
    assert not record.has_rule("alkali_metal_residual")


def test_ionization_is_not_reported_as_fragment_removal(
    engine: EngineWrapper,
) -> None:
    """Regressão: transformação não pode ser contabilizada como remoção.

    A linha de base de ``removed_fragments`` é a molécula padronizada, não a entrada
    crua - a normalização reescreve os fragmentos retidos.
    """
    record = engine.curate("T", "CC(=O)O[Na]")
    assert record.removed_fragments is None


# --- D-09: nenhuma rejeição antes do reparo -------------------------------------


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("C[N](C)(C)C", "C[N+](C)(C)C"),
        ("C[N](C)(C)CCO", "C[N+](C)(C)CCO"),
        ("c1ccccc1[N]#N", "N#[N+]c1ccccc1"),
    ],
)
def test_engine_repairs_what_strict_sanitization_would_reject(
    engine: EngineWrapper, smiles: str, expected: str
) -> None:
    """Amônios quaternários neutros e diazônios são reparados, não descartados.

    Uma sanitização estrita antes do motor levantaria ``AtomValenceException`` nos
    três casos, descartando classes inteiras de compostos válidos.
    """
    record = engine.curate("T", smiles)
    assert record.passed
    assert record.curated_smiles == expected


def test_valence_rejection_occurs_after_the_engine(engine: EngineWrapper) -> None:
    """Valência inválida é rejeitada, mas só depois da tentativa de reparo."""
    record = engine.curate("T", "C(C)(C)(C)(C)C")

    assert not record.passed
    assert record.rejection_code is RejectionCode.ERR_VALENCE
    assert record.rejection_stage is Stage.STANDARDIZE


def test_parse_stage_never_rejects_on_valence(engine: EngineWrapper) -> None:
    chemical_errors = {
        RejectionCode.ERR_VALENCE,
        RejectionCode.ERR_KEKULIZE,
        RejectionCode.ERR_SANITIZE,
    }
    for smiles in ("C(C)(C)(C)(C)C", "C[N](C)(C)C", "c1ccccc1[N]#N", "O[N](=O)O"):
        record = engine.curate("T", smiles)
        if not record.passed and record.rejection_stage is Stage.PARSE:
            assert record.rejection_code not in chemical_errors, smiles


# --- D-06: entidades multicomponentes -------------------------------------------


def test_combination_drug_is_preserved(engine: EngineWrapper) -> None:
    """Cotrimoxazol mantém os dois princípios ativos.

    ``LargestFragmentChooser`` devolveria apenas a trimetoprima.
    """
    cotrimoxazole = (
        "COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC."
        "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1"
    )
    record = engine.curate("T", cotrimoxazole)

    assert record.passed
    assert record.n_components_parent == 2
    assert record.removed_fragments is None


def test_all_components_matching_salt_list_are_kept(engine: EngineWrapper) -> None:
    """Benzoato de sódio: guard "tudo é sal, então mantém tudo".

    Tanto ``Sodium`` quanto ``Benzoate`` constam de ``salts.smi``. O valor esperado
    aqui é o da referência, não a intuição de que o sódio sairia.
    """
    record = engine.curate("T", "[Na+].[O-]C(=O)c1ccccc1")

    assert record.passed
    assert record.curated_smiles == "O=C([O-])c1ccccc1.[Na+]"
    assert record.n_components_parent == 2


def test_salt_is_stripped(engine: EngineWrapper) -> None:
    record = engine.curate("T", "N[C@@H](C)C(=O)O.Cl")

    assert record.passed
    assert record.curated_smiles == "C[C@H](N)C(=O)O"
    assert record.removed_fragments == "Cl"
    assert record.has_rule("fragments_removed")


def test_identical_components_are_deduplicated(engine: EngineWrapper) -> None:
    record = engine.curate("T", "CC(=O)O.CC(=O)O")

    assert record.passed
    assert record.curated_smiles == "CC(=O)O"
    assert record.n_components_parent == 1


# --- D-02: achatamento de tartaratos ---------------------------------------------


def test_tartrate_flattening_is_recorded(engine: EngineWrapper) -> None:
    """A perda de quiralidade é aceita, mas nunca silenciosa."""
    record = engine.curate("T", "OC(=O)[C@H](O)[C@@H](O)C(=O)O")

    assert record.passed
    assert record.has_rule("tartrate_flattened")
    assert record.delta_stereocenters == -2
    assert "@" not in record.curated_smiles


# --- D-03: exclude_flag ----------------------------------------------------------


def test_excluded_compound_is_kept_not_discarded(engine: EngineWrapper) -> None:
    """Cisplatina atravessa marcada, sem padronização."""
    record = engine.curate("T", "[Pt](Cl)(Cl)(N)N")

    assert record.passed
    assert record.excluded_flag
    assert record.has_rule("standardization_skipped_excluded")


# --- D-05 e D-07: estereoquímica e identidade ------------------------------------


def test_chirality_is_preserved(engine: EngineWrapper) -> None:
    record = engine.curate("T", "N[C@@H](C)C(=O)O")

    assert record.curated_smiles == "C[C@H](N)C(=O)O"
    assert record.inchikey == "QNAYBMKLOCPYGJ-REOHCLBHSA-N"
    assert record.delta_stereocenters == 0


def test_undefined_stereocenters_are_counted(engine: EngineWrapper) -> None:
    record = engine.curate("T", "CC(O)C(N)C(=O)O")
    assert record.n_undefined_stereocenters == 2
    assert record.n_stereocenters_total == 2


def test_enantiomers_get_distinct_identity_but_share_block1(
    engine: EngineWrapper,
) -> None:
    """Justifica a D-07: o bloco1 não distingue enantiômeros."""
    left = engine.curate("T", "N[C@@H](C)C(=O)O")
    right = engine.curate("T", "N[C@H](C)C(=O)O")

    assert left.inchikey != right.inchikey
    assert left.inchikey_block1 == right.inchikey_block1


# --- Rejeições de entrada ---------------------------------------------------------


@pytest.mark.parametrize(
    ("smiles", "code"),
    [
        ("C1CC", RejectionCode.ERR_SYNTAX),
        ("", RejectionCode.ERR_EMPTY),
        ("   ", RejectionCode.ERR_EMPTY),
    ],
)
def test_input_rejections(
    engine: EngineWrapper, smiles: str, code: RejectionCode
) -> None:
    record = engine.curate("T", smiles)

    assert not record.passed
    assert record.rejection_code is code
    assert record.rejection_stage is Stage.PARSE
    assert record.rejection_detail


def test_rejected_records_still_carry_provenance(engine: EngineWrapper) -> None:
    record = engine.curate("T", "C1CC")
    assert record.policy_hash == POLICY_HASH
    assert record.rdkit_version and record.csp_version and record.pipeline_version


# --- Propriedades formais ---------------------------------------------------------

CORPUS = (
    "CC(=O)O[Na]",
    "CCO[Na]",
    "C[N](C)(C)C",
    "c1ccccc1[N]#N",
    "N[C@@H](C)C(=O)O.Cl",
    "OC(=O)[C@H](O)[C@@H](O)C(=O)O",
    "[Na+].[O-]C(=O)c1ccccc1",
    "CC(=O)Oc1ccccc1C(=O)O",
    "CS(=O)C",
    "C[N+](C)(C)CCC[O-]",
    "CC(=O)O.CC(=O)O",
    "[Pt](Cl)(Cl)(N)N",
)


@pytest.mark.parametrize("smiles", CORPUS)
def test_idempotency(engine: EngineWrapper, smiles: str) -> None:
    """f(x) == f(f(x)): reprocessar a saída não pode alterá-la."""
    once = engine.curate("T", smiles)
    assert once.passed

    twice = engine.curate("T", once.curated_smiles)
    assert twice.passed
    assert twice.curated_smiles == once.curated_smiles
    assert twice.inchikey == once.inchikey


def test_caller_molecule_is_not_mutated(engine: EngineWrapper) -> None:
    """O wrapper passa cópias defensivas para a biblioteca."""
    mol = Chem.MolFromSmiles("CC(=O)O[Na]", sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    before = Chem.MolToSmiles(mol)

    engine.curate("T", "CC(=O)O[Na]")

    assert Chem.MolToSmiles(mol) == before


def test_determinism(engine: EngineWrapper) -> None:
    def signature() -> tuple:
        return tuple(
            (record.curated_smiles, tuple(record.rule_names()))
            for record in (engine.curate("T", smiles) for smiles in CORPUS)
        )

    assert len({signature() for _ in range(5)}) == 1


def test_curate_never_raises_on_hostile_input(engine: EngineWrapper) -> None:
    """Nenhuma entrada pode derrubar o processamento do lote."""
    for smiles in ("", "   ", "not a smiles", "C1CC", "[[[", "C" * 5000, "\x00"):
        record = engine.curate("T", smiles)
        assert record.status in ("PASSED", "REJECTED")


def test_record_serializes_round_trip(engine: EngineWrapper) -> None:
    record = engine.curate("CMPD_0000001", "CC(=O)O[Na]")
    restored = type(record).model_validate_json(record.model_dump_json())
    assert restored == record
