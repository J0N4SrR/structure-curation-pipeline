"""Testes do subsistema de I/O.

Dois eixos: ingestão não posicional (a coluna de estrutura é descoberta, não
assumida) e atomicidade de escrita (um lote interrompido não deixa saída
consumível).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from curation.io import (
    MANIFEST_NAME,
    BatchWriter,
    compute_policy_hash,
    read_input,
    read_manifest,
)
from curation.models import CurationRecord, RejectionCode, Stage, TransformationEvent

POLICY_HASH = "0" * 64


# --- Fábricas -----------------------------------------------------------------------


def make_record(
    input_id: str = "CMPD_0000001",
    status: str = "PASSED",
    policy_hash: str = POLICY_HASH,
) -> CurationRecord:
    common = dict(
        input_id=input_id,
        raw_smiles="CC(=O)O[Na]",
        rdkit_version="2026.03.5",
        csp_version="1.2.4",
        pipeline_version="0.1.0",
        policy_hash=policy_hash,
    )
    if status == "PASSED":
        return CurationRecord(
            status="PASSED",
            curated_smiles="CC(=O)[O-].[Na+]",
            inchikey="VMHLLURERBWHNL-UHFFFAOYSA-M",
            inchikey_block1="VMHLLURERBWHNL",
            n_components_parent=2,
            transformations=[
                TransformationEvent(
                    stage=Stage.STANDARDIZE,
                    rule="alkali_metal_ionized",
                    before_smiles="CC(=O)O[Na]",
                    after_smiles="CC(=O)[O-].[Na+]",
                    detail="ligacoes covalentes: 1 -> 0",
                )
            ],
            **common,
        )
    return CurationRecord(
        status="REJECTED",
        rejection_code=RejectionCode.ERR_VALENCE,
        rejection_stage=Stage.STANDARDIZE,
        rejection_detail="valencia explicita excedida",
        **common,
    )


# --- Ingestão -----------------------------------------------------------------------


def test_reads_csv_with_header_and_id_column(tmp_path: Path) -> None:
    source = tmp_path / "in.csv"
    source.write_text("compound_id,smiles\nA1,CCO\nA2,c1ccccc1\n", encoding="utf-8")
    assert list(read_input(source)) == [("A1", "CCO"), ("A2", "c1ccccc1")]


def test_reads_tsv(tmp_path: Path) -> None:
    source = tmp_path / "in.tsv"
    source.write_text("id\tsmiles\nX\tCCO\nY\tCCN\n", encoding="utf-8")
    assert list(read_input(source)) == [("X", "CCO"), ("Y", "CCN")]


def test_detects_smiles_column_when_name_comes_first(tmp_path: Path) -> None:
    """Formato ``NOME<TAB>SMILES`` — o usado pelos ``.smi`` do ChEMBL.

    É o inverso da convenção Daylight. Uma leitura posicional interpretaria
    ``Ethanolamine`` como estrutura.
    """
    source = tmp_path / "salts.smi"
    source.write_text(
        "Ethanolamine\tNCCO\nDeanol\tCN(C)CCO\nArginine\tNC(CCCNC(=N)N)C(=O)O\n",
        encoding="utf-8",
    )
    assert list(read_input(source)) == [
        ("Ethanolamine", "NCCO"),
        ("Deanol", "CN(C)CCO"),
        ("Arginine", "NC(CCCNC(=N)N)C(=O)O"),
    ]


def test_detects_smiles_column_when_smiles_comes_first(tmp_path: Path) -> None:
    """Formato Daylight ``SMILES<espaço>id``, o layout oposto do teste anterior."""
    source = tmp_path / "daylight.smi"
    source.write_text("NCCO ethanolamine\nCN(C)CCO deanol\n", encoding="utf-8")
    assert list(read_input(source)) == [
        ("ethanolamine", "NCCO"),
        ("deanol", "CN(C)CCO"),
    ]


def test_generates_deterministic_ids_when_id_absent(tmp_path: Path) -> None:
    source = tmp_path / "only_smiles.smi"
    source.write_text("CCO\nCCN\nc1ccccc1\n", encoding="utf-8")
    first = list(read_input(source))
    assert first == [
        ("CMPD_0000001", "CCO"),
        ("CMPD_0000002", "CCN"),
        ("CMPD_0000003", "c1ccccc1"),
    ]
    assert list(read_input(source)) == first, "leituras repetidas devem coincidir"


def test_skips_blank_lines_and_comments(tmp_path: Path) -> None:
    source = tmp_path / "messy.smi"
    source.write_text(
        "# comentario\n\nCCO\n   \n# outro\nCCN\n\n", encoding="utf-8"
    )
    assert [smiles for _, smiles in read_input(source)] == ["CCO", "CCN"]


def test_strips_surrounding_whitespace(tmp_path: Path) -> None:
    source = tmp_path / "spaced.csv"
    source.write_text("id,smiles\n  A1  ,  CCO  \n", encoding="utf-8")
    assert list(read_input(source)) == [("A1", "CCO")]


def test_reads_pasted_text() -> None:
    assert list(read_input("CCO\nCCN\n")) == [
        ("CMPD_0000001", "CCO"),
        ("CMPD_0000002", "CCN"),
    ]


def test_reads_file_handle(tmp_path: Path) -> None:
    source = tmp_path / "in.smi"
    source.write_text("CCO\nCCN\n", encoding="utf-8")
    with source.open(encoding="utf-8") as handle:
        assert [smiles for _, smiles in read_input(handle)] == ["CCO", "CCN"]


def test_handles_utf8_bom(tmp_path: Path) -> None:
    source = tmp_path / "bom.csv"
    source.write_bytes("id,smiles\nA1,CCO\n".encode("utf-8-sig"))
    assert list(read_input(source)) == [("A1", "CCO")]


def test_falls_back_gracefully_on_non_utf8(tmp_path: Path) -> None:
    """Bytes não-UTF-8 num campo de texto não podem derrubar o lote."""
    source = tmp_path / "latin1.csv"
    source.write_bytes("id,smiles\nMalv\xe3o,CCO\n".encode("latin-1"))
    rows = list(read_input(source))
    assert len(rows) == 1
    assert rows[0][1] == "CCO"


def test_empty_source_yields_nothing(tmp_path: Path) -> None:
    source = tmp_path / "empty.csv"
    source.write_text("", encoding="utf-8")
    assert list(read_input(source)) == []


def test_reading_is_lazy() -> None:
    """A leitura não pode materializar a fonte inteira em memória."""
    consumed: list[int] = []

    def lines():
        for index in range(1000):
            consumed.append(index)
            yield "CCO\n"

    stream = read_input(lines())
    next(stream)
    assert len(consumed) < 1000


# --- Escrita e atomicidade -----------------------------------------------------------


def test_commit_promotes_files_and_writes_manifest(tmp_path: Path) -> None:
    with BatchWriter(tmp_path, POLICY_HASH) as writer:
        writer.write(make_record("A1", "PASSED"))
        writer.write(make_record("A2", "REJECTED"))

    for name in ("curated.csv", "rejected.csv", "audit.csv", MANIFEST_NAME):
        assert (tmp_path / name).is_file(), f"{name} ausente"
    assert not list(tmp_path.glob("*.tmp")), "nenhum temporario deve sobreviver"

    manifest = read_manifest(tmp_path)
    assert manifest is not None
    assert manifest["counts"] == {"total": 2, "passed": 1, "rejected": 1}
    assert manifest["rejection_counts"] == {"ERR_VALENCE": 1}
    assert manifest["policy_hash"] == POLICY_HASH
    assert manifest["versions"]["rdkit"]
    assert manifest["versions"]["chembl_structure_pipeline"]
    assert manifest["started_at"] and manifest["finished_at"]


def test_exception_leaves_no_consumable_output(tmp_path: Path) -> None:
    """Um lote interrompido não pode deixar arquivo final nem manifesto."""
    with pytest.raises(RuntimeError):
        with BatchWriter(tmp_path, POLICY_HASH) as writer:
            writer.write(make_record("A1", "PASSED"))
            raise RuntimeError("falha simulada no meio do lote")

    assert not (tmp_path / "curated.csv").exists()
    assert not (tmp_path / "rejected.csv").exists()
    assert not (tmp_path / "audit.csv").exists()
    assert not (tmp_path / MANIFEST_NAME).exists()
    assert not list(tmp_path.glob("*.tmp")), "temporarios devem ser removidos"
    assert read_manifest(tmp_path) is None


def test_manifest_absence_marks_incomplete_run(tmp_path: Path) -> None:
    """Sem manifesto, a saída é inválida mesmo que os CSVs existam."""
    (tmp_path / "curated.csv").write_text("input_id\nA1\n", encoding="utf-8")
    assert read_manifest(tmp_path) is None


def test_stale_temporaries_are_cleaned_on_enter(tmp_path: Path) -> None:
    """Limpeza pós-crash: ``__exit__`` não roda quando o SO cai."""
    stale = tmp_path / "curated.csv.tmp"
    stale.write_text("lixo de execucao anterior", encoding="utf-8")
    (tmp_path / f"{MANIFEST_NAME}.tmp").write_text("{", encoding="utf-8")

    with BatchWriter(tmp_path, POLICY_HASH) as writer:
        writer.write(make_record("A1", "PASSED"))

    assert "lixo" not in (tmp_path / "curated.csv").read_text(encoding="utf-8")
    assert not list(tmp_path.glob("*.tmp"))


def test_io_failure_during_commit_aborts_cleanly(tmp_path: Path, monkeypatch) -> None:
    """Falha de disco na promoção não pode deixar saída pela metade."""
    import curation.io as io_module

    def exploding_replace(src, dst):
        raise OSError("disco cheio (simulado)")

    with pytest.raises(OSError):
        with BatchWriter(tmp_path, POLICY_HASH) as writer:
            writer.write(make_record("A1", "PASSED"))
            monkeypatch.setattr(io_module.os, "replace", exploding_replace)

    assert not (tmp_path / "curated.csv").exists()
    assert not (tmp_path / MANIFEST_NAME).exists()
    assert not list(tmp_path.glob("*.tmp"))


def test_output_dir_is_created(tmp_path: Path) -> None:
    target = tmp_path / "nested" / "run01"
    with BatchWriter(target, POLICY_HASH) as writer:
        writer.write(make_record("A1", "PASSED"))
    assert (target / MANIFEST_NAME).is_file()


def test_segregation_between_curated_and_rejected(tmp_path: Path) -> None:
    with BatchWriter(tmp_path, POLICY_HASH) as writer:
        writer.write(make_record("A1", "PASSED"))
        writer.write(make_record("A2", "REJECTED"))

    curated = (tmp_path / "curated.csv").read_text(encoding="utf-8")
    rejected = (tmp_path / "rejected.csv").read_text(encoding="utf-8")
    audit = (tmp_path / "audit.csv").read_text(encoding="utf-8")

    assert "A1" in curated and "A2" not in curated
    assert "A2" in rejected and "A1" not in rejected
    assert "A1" in audit and "A2" in audit


def test_audit_serializes_transformations_as_json(tmp_path: Path) -> None:
    import csv as csv_module

    with BatchWriter(tmp_path, POLICY_HASH) as writer:
        writer.write(make_record("A1", "PASSED"))

    with (tmp_path / "audit.csv").open(encoding="utf-8") as handle:
        row = next(iter(csv_module.DictReader(handle)))

    events = json.loads(row["transformations"])
    assert events[0]["rule"] == "alkali_metal_ionized"
    assert events[0]["stage"] == "STANDARDIZE"


def test_policy_hash_mismatch_is_counted(tmp_path: Path) -> None:
    """Misturar lotes de políticas distintas invalidaria a comparação."""
    with BatchWriter(tmp_path, POLICY_HASH) as writer:
        writer.write(make_record("A1", "PASSED"))
        writer.write(make_record("A2", "PASSED", policy_hash="f" * 64))

    manifest = read_manifest(tmp_path)
    assert manifest is not None
    assert manifest["policy_hash_mismatches"] == 1


def test_write_after_close_is_rejected(tmp_path: Path) -> None:
    with BatchWriter(tmp_path, POLICY_HASH) as writer:
        writer.write(make_record("A1", "PASSED"))
    with pytest.raises(RuntimeError):
        writer.write(make_record("A2", "PASSED"))


def test_empty_batch_still_produces_manifest(tmp_path: Path) -> None:
    with BatchWriter(tmp_path, POLICY_HASH):
        pass
    manifest = read_manifest(tmp_path)
    assert manifest is not None
    assert manifest["counts"] == {"total": 0, "passed": 0, "rejected": 0}


# --- Hash de política ----------------------------------------------------------------


def test_policy_hash_is_stable_and_content_sensitive(tmp_path: Path) -> None:
    document = tmp_path / "decisions.md"
    document.write_text("# ADR\nD-01\n", encoding="utf-8")
    original = compute_policy_hash(document)

    assert original == compute_policy_hash(document)
    assert len(original) == 64

    document.write_text("# ADR\nD-01 revisada\n", encoding="utf-8")
    assert compute_policy_hash(document) != original


def test_row_with_empty_smiles_cell_is_preserved(tmp_path: Path) -> None:
    """Alinhamento linha a linha: uma célula vazia é um registro, não um descarte.

    Descartá-la na leitura fazia um dado sumir sem decisão registrada, violando a
    regra de que nenhuma estrutura desaparece silenciosamente.
    """
    source = tmp_path / "com_vazio.csv"
    source.write_text(
        "SMILES,Tipo\nCCO,ok\n,vazio\nCCN,ok\n", encoding="utf-8"
    )
    rows = list(read_input(source))

    assert len(rows) == 3, "as tres linhas de dados devem virar registros"
    assert rows[1][1] == ""


def test_fully_blank_lines_are_still_skipped(tmp_path: Path) -> None:
    """Quebra de linha final não é dado."""
    source = tmp_path / "trailing.smi"
    source.write_text("CCO\n\nCCN\n\n\n", encoding="utf-8")
    assert [s for _, s in read_input(source)] == ["CCO", "CCN"]


def test_uppercase_smiles_header_is_detected(tmp_path: Path) -> None:
    """Cabeçalhos são comparados em minúsculas, então ``SMILES`` casa.

    E a coluna ``Tipo`` não vira identificador: só cabeçalhos reconhecidos como
    identificador assumem esse papel, o resto recebe id gerado.
    """
    source = tmp_path / "maiusculas.csv"
    source.write_text("SMILES,Tipo\nCCO,solvente\nCCN,amina\n", encoding="utf-8")
    assert list(read_input(source)) == [
        ("CMPD_0000001", "CCO"),
        ("CMPD_0000002", "CCN"),
    ]
