"""Contratos das decisões D-11 (isótopos) e D-12 (entradas sem carbono).

Os dois comportamentos já existiam no código; o que faltava era a decisão
registrada e o teste que a trava. Uma política implementada mas não documentada
não é auditável, e uma política não travada muda sem que ninguém perceba.

Inclui também testes de propriedade determinísticos — sem Hypothesis, ver a
justificativa em ``test_canonicalisation_converges_under_atom_reordering``.
"""

from __future__ import annotations

import random

import pytest
from rdkit import Chem

from curation import CurationPipeline, EligibilityCriteria
from curation.dedup import CollisionType, DedupIndex
from curation.models import RejectionCode, Stage

POLICY_HASH = "0" * 64


@pytest.fixture(scope="module")
def pipeline() -> CurationPipeline:
    return CurationPipeline(policy_hash=POLICY_HASH)


def key_of(smiles: str) -> str:
    return Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))


# --- D-11: política isotópica é REMOVE -------------------------------------------


@pytest.mark.parametrize(
    ("label", "labelled", "unlabelled"),
    [
        ("clorofórmio-d", "[2H]C(Cl)(Cl)Cl", "ClC(Cl)Cl"),
        ("etano 12C/13C", "[12CH3][13CH3]", "CC"),
        ("benzeno-d1", "[2H]c1ccccc1", "c1ccccc1"),
        ("etano-d6", "[2H]C([2H])([2H])C([2H])([2H])[2H]", "CC"),
        ("metano-13C", "[13CH4]", "C"),
    ],
)
def test_isotope_labels_are_removed(
    pipeline: CurationPipeline, label: str, labelled: str, unlabelled: str
) -> None:
    """D-11: a marcação isotópica é descartada na estrutura-mãe."""
    record = pipeline.process_single(labelled, "ISO")

    assert record.passed, label
    assert record.inchikey == key_of(unlabelled), label


def test_isotope_removal_preserves_connectivity(
    pipeline: CurationPipeline,
) -> None:
    """Remover isótopo não pode alterar o esqueleto."""
    labelled = pipeline.process_single("[2H]C([2H])([2H])C([2H])([2H])[2H]", "ISO")
    plain = pipeline.process_single("CC", "ISO")
    assert labelled.curated_smiles == plain.curated_smiles


def test_isotope_removal_is_idempotent(pipeline: CurationPipeline) -> None:
    once = pipeline.process_single("[2H]C(Cl)(Cl)Cl", "ISO")
    twice = pipeline.process_single(once.curated_smiles, "ISO")
    assert twice.inchikey == once.inchikey


def test_isotopologues_collapse_into_one_identity(
    pipeline: CurationPipeline,
) -> None:
    """Consequência declarada da D-11, e a razão de ela precisar ser registrada.

    Um estudo com padrões internos deuterados precisa saber que eles serão
    contabilizados como duplicatas do análogo comum.
    """
    index = DedupIndex()
    index.add(pipeline.process_single("CC", "PLAIN"))
    collision = index.add(
        pipeline.process_single("[2H]C([2H])([2H])C([2H])([2H])[2H]", "D6")
    )

    assert index.unique_count == 1
    assert collision is not None
    assert collision.collision_type is CollisionType.EXACT_DUPLICATE


# --- D-12: entradas sem carbono ------------------------------------------------------


@pytest.mark.parametrize(
    "inorganic", ["[Na+].[Cl-]", "[Cu+2]", "[OH2]", "O=S(=O)(O)O"]
)
def test_carbon_free_input_passes_by_default(
    pipeline: CurationPipeline, inorganic: str
) -> None:
    """D-12: por padrão o motor não rejeita estruturas sem carbono.

    Padrão desligado de propósito: ligá-lo por omissão mudaria silenciosamente o
    resultado de lotes já processados.
    """
    assert pipeline.process_single(inorganic, "INORG").passed, inorganic


@pytest.mark.parametrize("inorganic", ["[Na+].[Cl-]", "[Cu+2]", "[OH2]"])
def test_require_carbon_rejects_inorganic_input(inorganic: str) -> None:
    strict = CurationPipeline(
        POLICY_HASH, criteria=EligibilityCriteria(require_carbon=True)
    )
    record = strict.process_single(inorganic, "INORG")

    assert not record.passed, inorganic
    assert record.rejection_code is RejectionCode.ERR_NO_CARBON
    assert record.rejection_stage is Stage.ELIGIBILITY


@pytest.mark.parametrize(
    "organic",
    ["CCO", "CC(=O)Oc1ccccc1C(=O)O", "CC(=O)[O-].[Cu+2].CC(=O)[O-]"],
)
def test_require_carbon_keeps_anything_containing_carbon(organic: str) -> None:
    """Inclui o acetato de cobre: tem carbono, logo não é alcançado pelo critério."""
    strict = CurationPipeline(
        POLICY_HASH, criteria=EligibilityCriteria(require_carbon=True)
    )
    assert strict.process_single(organic, "ORG").passed, organic


def test_sodium_chloride_survives_because_every_fragment_is_a_salt(
    pipeline: CurationPipeline,
) -> None:
    """A origem do comportamento é o guard da D-06, não um descuido.

    Sodium e Chloride constam ambos de ``salts.smi``; quando todos os componentes
    casam, a referência não remove nada.
    """
    record = pipeline.process_single("[Na+].[Cl-]", "NACL")
    assert record.passed
    assert record.n_components_parent == 2
    assert record.removed_fragments is None


def test_require_carbon_is_recorded_in_the_parameters() -> None:
    """O critério aplicado precisa aparecer no manifesto para o lote ser reproduzível."""
    strict = CurationPipeline(
        POLICY_HASH, criteria=EligibilityCriteria(require_carbon=True)
    )
    assert strict.default_parameters()["require_carbon"] is True
    assert CurationPipeline(POLICY_HASH).default_parameters()["require_carbon"] is False


def test_criteria_description_mentions_the_active_criterion() -> None:
    assert "carbono" in EligibilityCriteria(require_carbon=True).describe()
    assert "carbono" not in EligibilityCriteria().describe()


# --- Propriedades ------------------------------------------------------------------------

#: Corpus pequeno e quimicamente escolhido, cobrindo as classes cujas políticas
#: divergem entre si: simples, sal, carga, estereocentro, ligação dupla estéreo,
#: isótopo, multicomponente e excluído.
PROPERTY_CORPUS = (
    "CC(=O)Oc1ccccc1C(=O)O",
    "N[C@@H](C)C(=O)O.Cl",
    "C[N+](C)(C)CCC[O-]",
    "N[C@@H](C)C(=O)O",
    "OC(=O)/C=C/C(=O)O",
    "[2H]C(Cl)(Cl)Cl",
    "COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC.Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1",
    "[Pt](Cl)(Cl)(N)N",
    "CC(=O)O[Na]",
)


@pytest.mark.parametrize("smiles", PROPERTY_CORPUS)
def test_normalisation_is_idempotent(pipeline: CurationPipeline, smiles: str) -> None:
    """N(N(x)) == N(x) sobre as classes cujas políticas diferem."""
    once = pipeline.process_single(smiles, "PROP")
    assert once.passed, smiles

    twice = pipeline.process_single(once.curated_smiles, "PROP")
    assert twice.curated_smiles == once.curated_smiles
    assert twice.inchikey == once.inchikey


@pytest.mark.parametrize("smiles", PROPERTY_CORPUS)
def test_pipeline_is_deterministic(pipeline: CurationPipeline, smiles: str) -> None:
    """P(x) == P(x) em identidade, status, parent e exclusão."""
    runs = [pipeline.process_single(smiles, "PROP") for _ in range(4)]
    signatures = {
        (
            record.status,
            record.curated_smiles,
            record.inchikey,
            record.removed_fragments,
            record.excluded_flag,
            tuple(record.rule_names()),
        )
        for record in runs
    }
    assert len(signatures) == 1, smiles


@pytest.mark.parametrize("smiles", PROPERTY_CORPUS)
def test_canonicalisation_converges_under_atom_reordering(
    pipeline: CurationPipeline, smiles: str
) -> None:
    """A identidade não pode depender da ordem dos átomos na entrada.

    Este é o único gerador de propriedades que vale a pena aqui, e é por isso que
    Hypothesis não foi introduzido: um gerador genérico de strings produziria quase
    só SMILES inválidos, exercitando o caminho de erro que já está coberto, ao passo
    que ``RenumberAtoms`` gera entradas **quimicamente idênticas e sintaticamente
    distintas** — exatamente a invariante que importa. Semente fixa, para que uma
    falha seja reproduzível.
    """
    reference = pipeline.process_single(smiles, "PERM")
    assert reference.passed, smiles

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(mol)

    generator = random.Random(20260907)
    for _ in range(5):
        order = list(range(mol.GetNumAtoms()))
        generator.shuffle(order)
        shuffled = Chem.MolToSmiles(
            Chem.RenumberAtoms(mol, order), canonical=False
        )
        record = pipeline.process_single(shuffled, "PERM")
        assert record.passed, shuffled
        assert record.inchikey == reference.inchikey, (
            f"a identidade mudou com a ordem dos átomos: {shuffled}"
        )


def test_duplicate_copies_collapse_to_one_identity(
    pipeline: CurationPipeline,
) -> None:
    """N cópias da mesma identidade produzem uma identidade."""
    index = DedupIndex()
    for position in range(6):
        index.add(pipeline.process_single("CC(=O)Oc1ccccc1C(=O)O", f"C{position}"))
    assert index.unique_count == 1
    assert len(index.collisions) == 5
