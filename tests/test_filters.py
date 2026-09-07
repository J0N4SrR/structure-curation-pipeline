"""Testes dos critérios de elegibilidade (D-10)."""

from __future__ import annotations

import pytest
from rdkit import Chem

from curation import CurationPipeline, EligibilityCriteria
from curation.filters import evaluate, molecular_properties
from curation.models import RejectionCode, Stage

POLICY_HASH = "0" * 64


def mol(smiles: str) -> Chem.Mol:
    parsed = Chem.MolFromSmiles(smiles, sanitize=False)
    parsed.UpdatePropertyCache(strict=False)
    return parsed


def test_defaults_match_the_adr() -> None:
    criteria = EligibilityCriteria()
    assert criteria.max_molecular_weight == 1000.0
    assert criteria.max_heavy_atoms == 100


def test_small_molecule_is_eligible() -> None:
    verdict = evaluate(mol("CC(=O)Oc1ccccc1C(=O)O"), EligibilityCriteria())
    assert verdict.eligible
    assert verdict.rejection_code is None
    assert 179 < verdict.molecular_weight < 181
    assert verdict.heavy_atoms == 13


def test_molecular_weight_cutoff() -> None:
    verdict = evaluate(mol("CC(=O)Oc1ccccc1C(=O)O"), EligibilityCriteria(max_molecular_weight=100))
    assert not verdict.eligible
    assert verdict.rejection_code is RejectionCode.ERR_MW_LIMIT
    assert "excede o limite" in verdict.detail


def test_heavy_atom_cutoff() -> None:
    verdict = evaluate(mol("CC(=O)Oc1ccccc1C(=O)O"), EligibilityCriteria(max_heavy_atoms=5))
    assert not verdict.eligible
    assert verdict.rejection_code is RejectionCode.ERR_HA_LIMIT
    assert "atomos pesados" in verdict.detail or "átomos pesados" in verdict.detail


def test_boundary_is_inclusive() -> None:
    weight, heavy = molecular_properties(mol("CCO"))
    assert evaluate(mol("CCO"), EligibilityCriteria(max_molecular_weight=weight)).eligible
    assert evaluate(mol("CCO"), EligibilityCriteria(max_heavy_atoms=heavy)).eligible


def test_properties_are_reported_even_on_rejection() -> None:
    """Rejeitar a 1002 Da é diferente de rejeitar a 4000 Da."""
    verdict = evaluate(mol("CC(=O)Oc1ccccc1C(=O)O"), EligibilityCriteria(max_molecular_weight=50))
    assert verdict.molecular_weight > 50
    assert verdict.heavy_atoms > 0


def test_works_on_unsanitizable_molecule() -> None:
    """Compostos com exclude_flag nunca são sanitizados (D-03), mas são filtrados."""
    verdict = evaluate(mol("[Pt](Cl)(Cl)(N)N"), EligibilityCriteria())
    assert verdict.eligible
    assert verdict.molecular_weight > 290


def test_multicomponent_context_appears_in_detail() -> None:
    """Interação D-10 x D-06: o corte incide sobre a soma dos componentes."""
    verdict = evaluate(
        mol("CC(=O)O.CC(=O)O"), EligibilityCriteria(max_molecular_weight=50), n_components=2
    )
    assert not verdict.eligible
    assert "componentes" in verdict.detail


# --- Integração com o pipeline -------------------------------------------------------


def test_pipeline_rejects_at_eligibility_stage() -> None:
    pipeline = CurationPipeline(
        POLICY_HASH, criteria=EligibilityCriteria(max_molecular_weight=100)
    )
    record = pipeline.process_single("CC(=O)Oc1ccccc1C(=O)O", "F1")

    assert not record.passed
    assert record.rejection_code is RejectionCode.ERR_MW_LIMIT
    assert record.rejection_stage is Stage.ELIGIBILITY
    assert record.parent_mw is not None


def test_cutoff_applies_to_parent_not_to_the_salt() -> None:
    """D-10: medir a massa do contra-íon rejeitaria compostos válidos.

    Cloridrato de alanina pesa ~125 Da; a estrutura-mãe, ~89 Da. Com o corte em
    100 Da, só passa se o filtro incidir sobre a estrutura-mãe.
    """
    pipeline = CurationPipeline(
        POLICY_HASH, criteria=EligibilityCriteria(max_molecular_weight=100)
    )
    record = pipeline.process_single("N[C@@H](C)C(=O)O.Cl", "F2")

    assert record.passed, "o corte deve incidir sobre o parent, nao sobre o sal"
    assert record.curated_smiles == "C[C@H](N)C(=O)O"
    assert record.parent_mw is not None and record.parent_mw < 100


def test_properties_recorded_on_passing_records() -> None:
    pipeline = CurationPipeline(POLICY_HASH)
    record = pipeline.process_single("CC(=O)Oc1ccccc1C(=O)O", "F3")
    assert record.parent_mw is not None and record.parent_heavy_atoms == 13


def test_excluded_compound_is_still_filtered() -> None:
    """D-03 preserva o composto, mas não o isenta dos critérios de escopo."""
    pipeline = CurationPipeline(
        POLICY_HASH, criteria=EligibilityCriteria(max_heavy_atoms=2)
    )
    record = pipeline.process_single("[Pt](Cl)(Cl)(N)N", "F4")
    assert not record.passed
    assert record.rejection_code is RejectionCode.ERR_HA_LIMIT
