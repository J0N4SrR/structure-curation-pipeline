"""Testes unitários das sondas post-hoc.

As sondas observam fronteiras de estágio. Aqui os estados antes/depois são
construídos à mão, isolando a lógica de detecção do comportamento do motor.
"""

from __future__ import annotations

from rdkit import Chem

from curation.models import Stage
from curation.probes import (
    DEFAULT_PROBES,
    CovalentAlkaliMetalProbe,
    FragmentCountProbe,
    NetChargeProbe,
    TartrateFlattenProbe,
    count_defined_stereocenters,
    fragment_smiles,
    run_probes,
    safe_smiles,
)


def mol(smiles: str) -> Chem.Mol:
    parsed = Chem.MolFromSmiles(smiles, sanitize=False)
    parsed.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(parsed)
    return parsed


# --- Metal alcalino ---------------------------------------------------------------


def test_alkali_probe_detects_ionization() -> None:
    events = CovalentAlkaliMetalProbe().inspect(
        mol("CC(=O)O[Na]"), mol("CC(=O)[O-].[Na+]"), Stage.STANDARDIZE
    )
    assert [event.rule for event in events] == ["alkali_metal_ionized"]


def test_alkali_probe_silent_when_no_metal() -> None:
    assert (
        CovalentAlkaliMetalProbe().inspect(
            mol("CCO"), mol("CCO"), Stage.STANDARDIZE
        )
        == []
    )


def test_alkali_probe_flags_residual_at_final_stage() -> None:
    """Alarme para a regressão de precedência da D-01."""
    events = CovalentAlkaliMetalProbe().inspect(
        mol("CC(=O)O[Na]"), mol("CC(=O)O[Na]"), Stage.GET_PARENT
    )
    assert [event.rule for event in events] == ["alkali_metal_residual"]


def test_alkali_probe_does_not_flag_residual_before_final_stage() -> None:
    events = CovalentAlkaliMetalProbe().inspect(
        mol("CC(=O)O[Na]"), mol("CC(=O)O[Na]"), Stage.STANDARDIZE
    )
    assert events == []


# --- Tartarato ---------------------------------------------------------------------


def test_tartrate_probe_detects_flattening() -> None:
    before = mol("OC(=O)[C@H](O)[C@@H](O)C(=O)O")
    after = mol("OC(=O)C(O)C(O)C(=O)O")
    events = TartrateFlattenProbe().inspect(before, after, Stage.STANDARDIZE)

    assert [event.rule for event in events] == ["tartrate_flattened"]
    assert "2 -> 0" in events[0].detail


def test_tartrate_probe_ignores_stereo_loss_elsewhere() -> None:
    """Exige as duas condições: padrão de tartarato e perda de centros.

    Sem isso, qualquer perda de estereoquímica seria atribuída ao achatamento.
    """
    before = mol("N[C@@H](C)C(=O)O")
    after = mol("NC(C)C(=O)O")
    assert TartrateFlattenProbe().inspect(before, after, Stage.STANDARDIZE) == []


def test_tartrate_probe_silent_when_stereo_preserved() -> None:
    structure = mol("OC(=O)[C@H](O)[C@@H](O)C(=O)O")
    assert TartrateFlattenProbe().inspect(structure, structure, Stage.STANDARDIZE) == []


# --- Carga líquida -----------------------------------------------------------------


def test_charge_probe_detects_neutralization() -> None:
    events = NetChargeProbe().inspect(
        mol("CC(=O)[O-]"), mol("CC(=O)O"), Stage.GET_PARENT
    )
    assert [event.rule for event in events] == ["net_charge_changed"]
    assert "-1 -> +0" in events[0].detail


def test_charge_probe_silent_on_balanced_salt() -> None:
    """O ``Uncharger`` não neutraliza enquanto houver contra-íon."""
    structure = mol("CC(=O)[O-].[Na+]")
    assert NetChargeProbe().inspect(structure, structure, Stage.STANDARDIZE) == []


# --- Contagem de fragmentos --------------------------------------------------------


def test_fragment_probe_reports_removal_with_smiles() -> None:
    events = FragmentCountProbe().inspect(
        mol("N[C@@H](C)C(=O)O.Cl"), mol("N[C@@H](C)C(=O)O"), Stage.GET_PARENT
    )
    assert [event.rule for event in events] == ["fragments_removed"]
    assert "Cl" in events[0].detail


def test_fragment_probe_reports_split() -> None:
    """Aumento de componentes é auditoria legítima, não anomalia."""
    events = FragmentCountProbe().inspect(
        mol("CC(=O)O[Na]"), mol("CC(=O)[O-].[Na+]"), Stage.STANDARDIZE
    )
    assert [event.rule for event in events] == ["fragments_split"]


def test_fragment_probe_silent_when_count_unchanged() -> None:
    structure = mol("CCO")
    assert FragmentCountProbe().inspect(structure, structure, Stage.PARSE) == []


# --- Utilitários e isolamento -------------------------------------------------------


def test_safe_smiles_tolerates_none_and_unsanitized() -> None:
    assert safe_smiles(None) == ""
    assert safe_smiles(mol("C[N](C)(C)C"))


def test_fragment_smiles_is_order_stable() -> None:
    left = fragment_smiles(mol("CCO.Cl"))
    right = fragment_smiles(mol("Cl.CCO"))
    assert left == right


def test_count_defined_stereocenters() -> None:
    assert count_defined_stereocenters(mol("N[C@@H](C)C(=O)O")) == 1
    assert count_defined_stereocenters(mol("NC(C)C(=O)O")) == 0
    assert count_defined_stereocenters(None) == 0


def test_a_failing_probe_does_not_break_the_record() -> None:
    """Sonda defeituosa degrada a auditoria, nunca derruba a curadoria."""

    class Exploding:
        rule = "exploding"

        def inspect(self, before, after, stage):
            raise ValueError("sonda com defeito")

    probes = (Exploding(), NetChargeProbe())
    events = run_probes(probes, mol("CC(=O)[O-]"), mol("CC(=O)O"), Stage.GET_PARENT)

    assert [event.rule for event in events] == ["net_charge_changed"]


def test_default_probe_order_is_deterministic() -> None:
    assert [type(probe).__name__ for probe in DEFAULT_PROBES] == [
        "CovalentAlkaliMetalProbe",
        "TartrateFlattenProbe",
        "NetChargeProbe",
        "FragmentCountProbe",
    ]
