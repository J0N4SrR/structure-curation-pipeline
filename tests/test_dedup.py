"""Testes da deduplicação e do relatório de conflitos (D-07)."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from curation import CurationPipeline
from curation.dedup import CollisionType, DedupIndex

POLICY_HASH = "0" * 64


@pytest.fixture(scope="module")
def pipeline() -> CurationPipeline:
    return CurationPipeline(policy_hash=POLICY_HASH)


def curate(pipeline: CurationPipeline, smiles: str, identifier: str):
    record = pipeline.process_single(smiles, identifier)
    assert record.passed, smiles
    return record


# --- Identidade -----------------------------------------------------------------------


def test_first_occurrence_is_indexed(pipeline: CurationPipeline) -> None:
    index = DedupIndex()
    record = curate(pipeline, "CCO", "A1")

    assert index.add(record) is None
    assert index.unique_count == 1
    assert index.first_occurrence(record.inchikey) == "A1"


def test_exact_duplicate_is_detected(pipeline: CurationPipeline) -> None:
    index = DedupIndex()
    index.add(curate(pipeline, "CCO", "A1"))
    collision = index.add(curate(pipeline, "OCC", "A2"))

    assert collision is not None
    assert collision.collision_type is CollisionType.EXACT_DUPLICATE
    assert (collision.existing_id, collision.incoming_id) == ("A1", "A2")
    assert index.unique_count == 1, "duplicata nao cria identidade nova"


def test_bare_anion_converges_to_the_free_acid(pipeline: CurationPipeline) -> None:
    """O ``Uncharger`` neutraliza o ânion isolado, tornando-o o mesmo composto.

    Sem contra-íon para balancear a carga, ``[O-]C(=O)c1ccccc1`` vira ácido
    benzoico - mesma identidade, portanto duplicata exata.
    """
    index = DedupIndex()
    index.add(curate(pipeline, "O=C(O)c1ccccc1", "A1"))
    collision = index.add(curate(pipeline, "[O-]C(=O)c1ccccc1", "A2"))

    assert collision is not None
    assert collision.collision_type is CollisionType.EXACT_DUPLICATE
    assert index.unique_count == 1


def test_sodium_salt_is_a_different_compound_entirely(
    pipeline: CurationPipeline,
) -> None:
    """Interação D-06 x D-07: o sódio é retido, então a molécula é outra.

    Como o guard "tudo é sal" preserva o ``[Na+]``, o benzoato de sódio não é uma
    forma protonada do ácido benzoico - é outro composto, com outro esqueleto de
    conectividade. Não há colisão a reportar.
    """
    acid = curate(pipeline, "O=C(O)c1ccccc1", "A1")
    salt = curate(pipeline, "[Na+].[O-]C(=O)c1ccccc1", "A2")

    assert acid.inchikey_block1 != salt.inchikey_block1

    index = DedupIndex()
    index.add(acid)
    assert index.add(salt) is None
    assert index.unique_count == 2


# --- Bloco 1 --------------------------------------------------------------------------


def test_enantiomers_collide_on_block1_but_are_not_duplicates(
    pipeline: CurationPipeline,
) -> None:
    """A razão de ser da D-07: bloco1 ignora estereoquímica."""
    index = DedupIndex()
    left = curate(pipeline, "N[C@@H](C)C(=O)O", "L")
    right = curate(pipeline, "N[C@H](C)C(=O)O", "D")

    index.add(left)
    collision = index.add(right)

    assert left.inchikey != right.inchikey
    assert left.inchikey_block1 == right.inchikey_block1
    assert collision is not None
    assert collision.collision_type is CollisionType.BLOCK1_COLLISION
    assert "nao e duplicata" in collision.detail or "não é duplicata" in collision.detail
    assert index.unique_count == 2, "enantiomeros nao podem ser fundidos"


def test_block1_collision_is_subclassified(pipeline: CurationPipeline) -> None:
    index = DedupIndex()
    index.add(curate(pipeline, "N[C@@H](C)C(=O)O", "L"))
    collision = index.add(curate(pipeline, "N[C@H](C)C(=O)O", "D"))

    assert collision is not None
    assert "estereoquímica" in collision.detail or "estereoq" in collision.detail


def test_group_returns_related_identities(pipeline: CurationPipeline) -> None:
    index = DedupIndex()
    left = curate(pipeline, "N[C@@H](C)C(=O)O", "L")
    right = curate(pipeline, "N[C@H](C)C(=O)O", "D")
    index.add(left)
    index.add(right)

    assert index.group(left.inchikey_block1) == {left.inchikey, right.inchikey}


# --- Conflito de anotação --------------------------------------------------------------


def test_annotation_conflict_is_distinguished_from_plain_duplicate(
    pipeline: CurationPipeline,
) -> None:
    """Mesma estrutura com anotações divergentes exige decisão humana."""
    index = DedupIndex()
    index.add(curate(pipeline, "CCO", "A1"), annotation="ativo")
    collision = index.add(curate(pipeline, "OCC", "A2"), annotation="inativo")

    assert collision is not None
    assert collision.collision_type is CollisionType.ANNOTATION_CONFLICT
    assert "ativo" in collision.detail and "inativo" in collision.detail


def test_matching_annotations_are_plain_duplicates(
    pipeline: CurationPipeline,
) -> None:
    index = DedupIndex()
    index.add(curate(pipeline, "CCO", "A1"), annotation="ativo")
    collision = index.add(curate(pipeline, "OCC", "A2"), annotation="ativo")

    assert collision is not None
    assert collision.collision_type is CollisionType.EXACT_DUPLICATE


# --- Registros ignorados ---------------------------------------------------------------


def test_rejected_records_are_not_indexed(pipeline: CurationPipeline) -> None:
    index = DedupIndex()
    assert index.add(pipeline.process_single("C(C)(C)(C)(C)C", "R1")) is None
    assert index.unique_count == 0


# --- Relatório --------------------------------------------------------------------------


def test_conflicts_file_is_written_even_when_empty(tmp_path: Path) -> None:
    """Ausência de conflitos é um resultado, distinto de relatório não gerado."""
    target = tmp_path / "conflicts.csv"
    assert DedupIndex().write_conflicts(target) == 0
    assert target.is_file()

    with target.open(encoding="utf-8") as handle:
        assert next(csv.reader(handle))[0] == "collision_type"


def test_conflicts_file_records_each_collision(
    pipeline: CurationPipeline, tmp_path: Path
) -> None:
    index = DedupIndex()
    index.add(curate(pipeline, "CCO", "A1"))
    index.add(curate(pipeline, "OCC", "A2"))
    index.add(curate(pipeline, "N[C@@H](C)C(=O)O", "L"))
    index.add(curate(pipeline, "N[C@H](C)C(=O)O", "D"))

    target = tmp_path / "conflicts.csv"
    assert index.write_conflicts(target) == 2

    with target.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    kinds = {row["collision_type"] for row in rows}
    assert kinds == {"EXACT_DUPLICATE", "BLOCK1_COLLISION"}
    assert index.collision_counts() == {"EXACT_DUPLICATE": 1, "BLOCK1_COLLISION": 1}


# --- Integração com o lote ---------------------------------------------------------------


def test_run_emits_conflicts_under_the_atomicity_guarantee(tmp_path: Path) -> None:
    pipeline = CurationPipeline(policy_hash=POLICY_HASH)
    source = "\n".join(["CCO", "OCC", "N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"])

    summary = pipeline.run(source, tmp_path)

    assert (tmp_path / "conflicts.csv").is_file()
    assert not list(tmp_path.glob("*.tmp"))
    assert summary.unique == 3
    assert summary.conflicts == 2

    from curation.io import read_manifest

    manifest = read_manifest(tmp_path)
    assert manifest["deduplication"]["unique"] == 3
    assert manifest["deduplication"]["conflicts"] == 2
    assert manifest["outputs"]["conflicts"] == "conflicts.csv"


def test_dedup_can_be_disabled(tmp_path: Path) -> None:
    pipeline = CurationPipeline(policy_hash=POLICY_HASH, deduplicate=False)
    summary = pipeline.run("CCO\nOCC\n", tmp_path)
    assert summary.conflicts == 0
    assert summary.unique == 0
