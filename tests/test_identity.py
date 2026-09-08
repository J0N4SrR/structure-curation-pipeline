"""Preservação de identidade estrutural em fármacos-controle.

A asserção é sobre **InChIKey**, nunca sobre igualdade literal de SMILES: a
representação textual muda sem que o composto mude. ``CC(=O)O`` e ``OC(C)=O`` são
o mesmo ácido acético, e o pipeline devolve a forma canônica do RDKit — comparar
strings de entrada com strings de saída falharia por motivo errado.

A distinção que organiza este arquivo:

* **Identidade** — o composto entra e sai sendo o mesmo. O InChIKey se conserva.
* **Normalização** — o composto entra e sai *diferente por decisão de política*
  (sal removido, carga neutralizada, tautômero preservado). Aqui o InChIKey **não**
  se conserva, e comparar contra a entrada seria o teste errado; compara-se contra
  a estrutura-mãe esperada.

Todos os SMILES de referência foram verificados por fórmula molecular contra os
valores de literatura antes de entrarem aqui.
"""

from __future__ import annotations

import pytest
from rdkit import Chem

from curation import CurationPipeline

POLICY_HASH = "0" * 64


@pytest.fixture(scope="module")
def pipeline() -> CurationPipeline:
    return CurationPipeline(policy_hash=POLICY_HASH)


def inchikey(smiles: str) -> str:
    """InChIKey de um SMILES, para comparar identidade química."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {smiles!r}")
    return Chem.MolToInchiKey(mol)


#: (nome, ChEMBL ID, SMILES, fórmula molecular de literatura).
#: A fórmula é verificada em teste próprio: um controle cuja estrutura esteja
#: errada validaria o pipeline contra a molécula errada, sem que nada acusasse.
CONTROL_DRUGS: tuple[tuple[str, str, str, str], ...] = (
    ("Paracetamol", "CHEMBL112", "CC(=O)NC1=CC=C(O)C=C1", "C8H9NO2"),
    ("Aspirina", "CHEMBL25", "CC(=O)OC1=CC=CC=C1C(=O)O", "C9H8O4"),
    ("Ibuprofeno", "CHEMBL521", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "C13H18O2"),
    ("Naproxeno", "CHEMBL154", "COc1ccc2cc(C(C)C(=O)O)ccc2c1", "C14H14O3"),
    ("Diclofenaco", "CHEMBL139", "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl", "C14H11Cl2NO2"),
    ("Metformina", "CHEMBL1431", "CN(C)C(=N)NC(=N)N", "C4H11N5"),
    ("Metoprolol", "CHEMBL13", "CC(C)NCC(COC1=CC=C(C=C1)CCOC)O", "C15H25NO3"),
    ("Propranolol", "CHEMBL27", "CC(C)NCC(COC1=CC=CC2=CC=CC=C21)O", "C16H21NO2"),
    ("Lidocaina", "CHEMBL79", "CCN(CC)CC(=O)Nc1c(C)cccc1C", "C14H22N2O"),
    ("Captopril", "CHEMBL1560", "CC(CS)C(=O)N1CCCC1C(=O)O", "C9H15NO3S"),
    (
        "Losartana", "CHEMBL191",
        "CCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1", "C22H23ClN6O",
    ),
    (
        "Omeprazol", "CHEMBL1503",
        "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1", "C17H19N3O3S",
    ),
    ("Fluconazol", "CHEMBL106", "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F", "C13H12F2N6O"),
    ("Warfarina", "CHEMBL1464", "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O", "C19H16O4"),
    (
        "Atorvastatina", "CHEMBL1487",
        "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)"
        "C[C@@H](O)CC(=O)O", "C33H35FN2O5",
    ),
    (
        "Sildenafil", "CHEMBL192",
        "CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12",
        "C22H30N6O4S",
    ),
    ("Cafeina", "CHEMBL113", "CN1C(=O)N(C)c2ncn(C)c2C1=O", "C8H10N4O2"),
    (
        "Amoxicilina", "CHEMBL1082",
        "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)C3=CC=C(O)C=C3)C(=O)N2[C@H]1C(=O)O",
        "C16H19N3O5S",
    ),
    ("Dopamina", "CHEMBL59", "C1=CC(=C(C=C1CCN)O)O", "C8H11NO2"),
    ("Clorpromazina", "CHEMBL71", "CN(C)CCCN1c2ccccc2Sc2ccc(Cl)cc21", "C17H19ClN2S"),
)

IDS = [name for name, _, _, _ in CONTROL_DRUGS]


# --- Integridade do próprio conjunto de controle -------------------------------


@pytest.mark.parametrize(
    ("name", "chembl_id", "smiles", "formula"), CONTROL_DRUGS, ids=IDS
)
def test_control_structure_matches_its_declared_formula(
    name: str, chembl_id: str, smiles: str, formula: str
) -> None:
    """O controle precisa ser a molécula que diz ser.

    Sem esta verificação, um SMILES errado rotulado com o nome de um fármaco
    passaria em todos os testes de identidade — validando o pipeline contra a
    molécula errada e sem nada acusar. Metade dos SMILES propostos originalmente
    para este conjunto falhava aqui.
    """
    from rdkit.Chem import rdMolDescriptors

    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"{name} ({chembl_id}) não parseia"
    assert rdMolDescriptors.CalcMolFormula(mol) == formula


# --- Happy path: identidade preservada ------------------------------------------


@pytest.mark.parametrize(
    ("name", "chembl_id", "smiles", "formula"), CONTROL_DRUGS, ids=IDS
)
def test_control_drug_keeps_its_identity(
    pipeline: CurationPipeline, name: str, chembl_id: str, smiles: str, formula: str
) -> None:
    """Um fármaco neutro e sem sal atravessa o pipeline sendo o mesmo composto."""
    record = pipeline.process_single(smiles, chembl_id)

    assert record.passed, f"{name}: {record.rejection_detail}"
    assert record.inchikey == inchikey(smiles), f"{name} mudou de identidade"


@pytest.mark.parametrize(
    ("name", "chembl_id", "smiles", "formula"), CONTROL_DRUGS, ids=IDS
)
def test_control_drug_is_idempotent(
    pipeline: CurationPipeline, name: str, chembl_id: str, smiles: str, formula: str
) -> None:
    once = pipeline.process_single(smiles, chembl_id)
    twice = pipeline.process_single(once.curated_smiles, chembl_id)
    assert twice.inchikey == once.inchikey, name


@pytest.mark.parametrize(
    ("name", "chembl_id", "smiles", "formula"), CONTROL_DRUGS, ids=IDS
)
def test_defined_stereocenters_survive(
    pipeline: CurationPipeline, name: str, chembl_id: str, smiles: str, formula: str
) -> None:
    """Nenhum controle pode perder estereoquímica.

    Nenhum deles é tartarato, então a exceção da D-02 não se aplica: qualquer
    perda aqui é regressão.
    """
    record = pipeline.process_single(smiles, chembl_id)
    assert record.delta_stereocenters == 0, name


def test_representation_changes_but_identity_does_not(
    pipeline: CurationPipeline,
) -> None:
    """A mesma molécula escrita de formas diferentes converge.

    É a razão de a asserção ser sobre InChIKey e não sobre a string do SMILES.
    """
    variants = [
        "CC(=O)Oc1ccccc1C(=O)O",
        "O=C(O)c1ccccc1OC(C)=O",
        "CC(=O)OC1=CC=CC=C1C(O)=O",
    ]
    keys = {pipeline.process_single(v, "ASA").inchikey for v in variants}
    assert len(keys) == 1


# --- Normalização: identidade muda por decisão de política ------------------------


@pytest.mark.parametrize(
    ("name", "salt_form", "expected_parent"),
    [
        ("cloridrato de propranolol", "CC(C)NCC(COC1=CC=CC2=CC=CC=C21)O.Cl",
         "CC(C)NCC(COC1=CC=CC2=CC=CC=C21)O"),
        ("cloridrato de metformina", "CN(C)C(=N)NC(=N)N.Cl", "CN(C)C(=N)NC(=N)N"),
        ("cloridrato de dopamina", "C1=CC(=C(C=C1CCN)O)O.Cl", "C1=CC(=C(C=C1CCN)O)O"),
        ("cloridrato de lidocaína", "CCN(CC)CC(=O)Nc1c(C)cccc1C.Cl",
         "CCN(CC)CC(=O)Nc1c(C)cccc1C"),
    ],
)
def test_salt_form_converges_to_the_free_base(
    pipeline: CurationPipeline, name: str, salt_form: str, expected_parent: str
) -> None:
    """A comparação é contra a estrutura-mãe esperada, não contra a entrada.

    Aqui a identidade **deve** mudar: é remoção de sal, não preservação. Exigir
    que o InChIKey da entrada se conservasse seria o teste errado.
    """
    record = pipeline.process_single(salt_form, "SALT")

    assert record.passed, name
    assert record.inchikey == inchikey(expected_parent), name
    assert record.removed_fragments == "Cl", name


def test_free_base_and_its_salt_share_the_final_identity(
    pipeline: CurationPipeline,
) -> None:
    """Depois da dessalificação, sal e base livre são a mesma identidade."""
    base = pipeline.process_single("CN(C)C(=N)NC(=N)N", "BASE")
    salt = pipeline.process_single("CN(C)C(=N)NC(=N)N.Cl", "SALT")
    assert base.inchikey == salt.inchikey


# --- Negative testing: entradas que devem falhar ------------------------------------


@pytest.mark.parametrize(
    ("label", "smiles", "expected_code"),
    [
        ("anel não fechado", "C1CC", "ERR_SYNTAX"),
        ("parêntese aberto", "CC(", "ERR_SYNTAX"),
        ("colchete solto", "][", "ERR_SYNTAX"),
        ("texto livre", "não é um smiles", "ERR_SYNTAX"),
        ("string vazia", "", "ERR_EMPTY"),
        ("apenas espaços", "     ", "ERR_EMPTY"),
        ("carbono pentavalente", "C(C)(C)(C)(C)C", "ERR_VALENCE"),
        ("oxigênio trivalente", "CO(C)C", "ERR_VALENCE"),
    ],
)
def test_invalid_input_is_rejected_with_the_right_code(
    pipeline: CurationPipeline, label: str, smiles: str, expected_code: str
) -> None:
    record = pipeline.process_single(smiles, "NEG")

    assert not record.passed, label
    assert record.rejection_code is not None
    assert record.rejection_code.value == expected_code, label
    assert record.rejection_detail, "toda rejeição precisa de motivo legível"


def test_rejected_records_carry_no_identity(pipeline: CurationPipeline) -> None:
    """Um registro rejeitado não pode sair com InChIKey preenchido."""
    record = pipeline.process_single("C(C)(C)(C)(C)C", "NEG")
    assert record.inchikey is None
    assert record.curated_smiles is None


def test_valence_rejection_is_not_premature(pipeline: CurationPipeline) -> None:
    """Estruturas reparáveis não podem cair no filtro de valência.

    Amônio quaternário neutro e diazônio neutro falham numa sanitização estrita,
    mas o normalizador da referência os corrige.
    """
    for smiles in ("C[N](C)(C)C", "c1ccccc1[N]#N"):
        assert pipeline.process_single(smiles, "NEG").passed, smiles


# --- Edge cases: limites -------------------------------------------------------------


def test_single_atom_is_valid(pipeline: CurationPipeline) -> None:
    record = pipeline.process_single("C", "EDGE")
    assert record.passed
    assert record.parent_heavy_atoms == 1


def test_cutoff_boundary_is_inclusive(pipeline: CurationPipeline) -> None:
    """No limite exato o composto passa; um passo além, não."""
    from curation.filters import EligibilityCriteria

    reference = pipeline.process_single("CC(=O)Oc1ccccc1C(=O)O", "EDGE")
    weight = reference.parent_mw
    assert weight is not None

    at_limit = CurationPipeline(
        POLICY_HASH, criteria=EligibilityCriteria(max_molecular_weight=weight)
    )
    just_below = CurationPipeline(
        POLICY_HASH, criteria=EligibilityCriteria(max_molecular_weight=weight - 0.01)
    )
    assert at_limit.process_single("CC(=O)Oc1ccccc1C(=O)O", "EDGE").passed
    assert not just_below.process_single("CC(=O)Oc1ccccc1C(=O)O", "EDGE").passed


def test_very_long_input_does_not_break(pipeline: CurationPipeline) -> None:
    record = pipeline.process_single("C" * 400, "EDGE")
    assert record.status in ("PASSED", "REJECTED")


def test_whitespace_around_structure_is_tolerated(
    pipeline: CurationPipeline,
) -> None:
    padded = pipeline.process_single("  CC(=O)Oc1ccccc1C(=O)O  ", "EDGE")
    clean = pipeline.process_single("CC(=O)Oc1ccccc1C(=O)O", "EDGE")
    assert padded.inchikey == clean.inchikey


def test_all_fragments_are_salts(pipeline: CurationPipeline) -> None:
    """Quando tudo casa com a tabela de sais, a referência não remove nada."""
    record = pipeline.process_single("[Na+].[O-]C(=O)c1ccccc1", "EDGE")
    assert record.passed
    assert record.n_components_parent == 2


def test_repeated_identical_fragments_collapse(pipeline: CurationPipeline) -> None:
    record = pipeline.process_single("CC(=O)O.CC(=O)O.CC(=O)O", "EDGE")
    assert record.inchikey == inchikey("CC(=O)O")


def test_isotope_is_stripped(pipeline: CurationPipeline) -> None:
    assert pipeline.process_single("[13CH4]", "EDGE").inchikey == inchikey("C")


def test_excluded_compound_survives_the_valence_gate(
    pipeline: CurationPipeline,
) -> None:
    """Compostos com exclude_flag não passam por sanitização estrita (D-03)."""
    record = pipeline.process_single("[Pt](Cl)(Cl)(N)N", "EDGE")
    assert record.passed and record.excluded_flag


# --- Error handling: nada derruba o lote -----------------------------------------------


@pytest.mark.parametrize(
    "hostile",
    ["", "   ", "\x00", "\n\t", "[[[", "C(((", "%%%%", "C" * 5000, "🧪", "NaN", "None"],
)
def test_hostile_input_yields_a_record_instead_of_an_exception(
    pipeline: CurationPipeline, hostile: str
) -> None:
    record = pipeline.process_single(hostile, "HOSTILE")
    assert record.status in ("PASSED", "REJECTED")
    if not record.passed:
        assert record.rejection_code is not None
        assert record.rejection_stage is not None


def test_a_broken_record_does_not_stop_the_batch(pipeline: CurationPipeline) -> None:
    """Uma entrada ruim no meio não pode impedir as seguintes de processar."""
    corpus = [
        ("A", "CC(=O)Oc1ccccc1C(=O)O"),
        ("B", "C(C)(C)(C)(C)C"),
        ("C", "não é um smiles"),
        ("D", "CN1C(=O)N(C)c2ncn(C)c2C1=O"),
    ]
    records = list(pipeline.process_many(corpus))

    assert len(records) == 4
    assert [r.passed for r in records] == [True, False, False, True]


def test_engine_failure_is_contained(pipeline: CurationPipeline) -> None:
    class Exploding:
        _rdkit_version = "x"
        _csp_version = "y"

        def curate(self, input_id: str, raw_smiles: str):
            raise MemoryError("falha imprevista")

    contained = CurationPipeline(POLICY_HASH, engine=Exploding())
    record = contained.process_single("CCO", "ERR")

    assert not record.passed
    assert record.rejection_code.value == "ERR_INTERNAL"
    assert "falha imprevista" in record.rejection_detail


def test_every_record_carries_provenance_even_when_rejected(
    pipeline: CurationPipeline,
) -> None:
    for smiles in ("CC(=O)Oc1ccccc1C(=O)O", "C(C)(C)(C)(C)C", ""):
        record = pipeline.process_single(smiles, "PROV")
        assert record.policy_hash == POLICY_HASH
        assert record.rdkit_version and record.csp_version


def test_batch_of_control_drugs_is_fully_accounted_for(
    pipeline: CurationPipeline,
) -> None:
    """Nenhum controle desaparece: entradas = aprovados + rejeitados."""
    corpus = "\n".join(smiles for _, _, smiles, _ in CONTROL_DRUGS)
    report = pipeline.run_report(corpus, input_bytes=corpus.encode())

    assert report.total == len(CONTROL_DRUGS)
    assert len(report.approved) == len(CONTROL_DRUGS)
    assert report.n_unique == len(CONTROL_DRUGS), "controles são todos distintos"
