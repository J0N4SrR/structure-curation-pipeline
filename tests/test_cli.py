"""Testes da interface de linha de comando."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from curation.cli import build_parser, run
from curation.pipeline import BatchSummary

CORPUS = "\n".join(
    [
        "CC(=O)O[Na]",
        "N[C@@H](C)C(=O)O.Cl",
        "CC(=O)Oc1ccccc1C(=O)O",
        "C(C)(C)(C)(C)C",
        "C1CC",
        "CCO",
        "OCC",
    ]
)


@pytest.fixture
def corpus(tmp_path: Path) -> Path:
    source = tmp_path / "input.smi"
    source.write_text(CORPUS + "\n", encoding="utf-8")
    return source


# --- Parser ---------------------------------------------------------------------


def test_parser_defaults_match_the_adr() -> None:
    args = build_parser().parse_args(["-i", "x.smi", "-o", "out"])
    assert args.max_mw == 1000.0
    assert args.max_ha == 100


def test_input_and_out_dir_are_required() -> None:
    with pytest.raises(SystemExit):
        build_parser().parse_args(["-i", "x.smi"])
    with pytest.raises(SystemExit):
        build_parser().parse_args(["-o", "out"])


# --- Execução --------------------------------------------------------------------


def test_run_produces_all_outputs(corpus: Path, tmp_path: Path) -> None:
    out_dir = tmp_path / "out"
    assert run(["-i", str(corpus), "-o", str(out_dir), "--quiet"]) == 0

    for name in ("curated.csv", "rejected.csv", "audit.csv", "conflicts.csv", "manifest.json"):
        assert (out_dir / name).is_file(), name
    assert not list(out_dir.glob("*.tmp"))


def test_cutoffs_are_honoured(corpus: Path, tmp_path: Path) -> None:
    out_dir = tmp_path / "out"
    run(["-i", str(corpus), "-o", str(out_dir), "--max-mw", "50", "--quiet"])

    manifest = json.loads((out_dir / "manifest.json").read_text(encoding="utf-8"))
    assert "ERR_MW_LIMIT" in manifest["rejection_counts"]


def test_heavy_atom_cutoff_is_honoured(corpus: Path, tmp_path: Path) -> None:
    out_dir = tmp_path / "out"
    run(["-i", str(corpus), "-o", str(out_dir), "--max-ha", "3", "--quiet"])

    manifest = json.loads((out_dir / "manifest.json").read_text(encoding="utf-8"))
    assert "ERR_HA_LIMIT" in manifest["rejection_counts"]


def test_no_dedup_flag_disables_conflict_detection(
    corpus: Path, tmp_path: Path
) -> None:
    out_dir = tmp_path / "out"
    run(["-i", str(corpus), "-o", str(out_dir), "--no-dedup", "--quiet"])

    manifest = json.loads((out_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["deduplication"]["conflicts"] == 0


def test_missing_input_returns_error_code(tmp_path: Path) -> None:
    assert run(["-i", str(tmp_path / "ausente.smi"), "-o", str(tmp_path / "o")]) == 2


def test_missing_decisions_marks_batch_as_unversioned(
    corpus: Path, tmp_path: Path, capsys
) -> None:
    """Sem o documento de política o lote é marcado, nunca recebe hash falso."""
    out_dir = tmp_path / "out"
    run(
        [
            "-i", str(corpus), "-o", str(out_dir),
            "--decisions", str(tmp_path / "ausente.md"), "--quiet",
        ]
    )
    manifest = json.loads((out_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["policy_hash"] == "UNVERSIONED_POLICY"
    assert "nao sera comparavel" in capsys.readouterr().err


def test_real_decisions_file_yields_a_hash(corpus: Path, tmp_path: Path) -> None:
    decisions = tmp_path / "decisions.md"
    decisions.write_text("# ADR\nD-01\n", encoding="utf-8")
    out_dir = tmp_path / "out"
    run(["-i", str(corpus), "-o", str(out_dir), "--decisions", str(decisions), "--quiet"])

    manifest = json.loads((out_dir / "manifest.json").read_text(encoding="utf-8"))
    assert len(manifest["policy_hash"]) == 64


# --- Sumário ----------------------------------------------------------------------


def test_summary_is_printed(corpus: Path, tmp_path: Path, capsys) -> None:
    run(["-i", str(corpus), "-o", str(tmp_path / "out")])
    out = capsys.readouterr().out

    assert "SUMARIO DA EXECUCAO" in out
    assert "entradas processadas" in out
    assert "top motivos de rejeicao" in out
    assert "criterios:" in out


def test_quiet_suppresses_summary(corpus: Path, tmp_path: Path, capsys) -> None:
    run(["-i", str(corpus), "-o", str(tmp_path / "out"), "--quiet"])
    assert capsys.readouterr().out == ""


def test_top_rejections_is_capped_at_five() -> None:
    summary = BatchSummary(
        total=60,
        passed=0,
        rejected=60,
        rejection_counts={f"ERR_{i}": 10 - i for i in range(8)},
    )
    top = summary.top_rejections()
    assert len(top) == 5
    assert [count for _, count in top] == sorted(
        [count for _, count in top], reverse=True
    )


def test_report_includes_pass_rate_and_dedup_block() -> None:
    summary = BatchSummary(
        total=10, passed=8, rejected=2,
        rejection_counts={"ERR_VALENCE": 2},
        unique=7, conflicts=1, collision_counts={"EXACT_DUPLICATE": 1},
    )
    report = summary.format_report()

    assert "80.0%" in report
    assert "deduplicacao" in report
    assert "EXACT_DUPLICATE" in report


def test_stdin_is_accepted(tmp_path: Path, monkeypatch) -> None:
    import io as std_io

    monkeypatch.setattr("sys.stdin", std_io.StringIO("CCO\nCCN\n"))
    out_dir = tmp_path / "out"
    assert run(["-i", "-", "-o", str(out_dir), "--quiet"]) == 0

    manifest = json.loads((out_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["counts"]["total"] == 2


def test_empty_input_is_a_failure_not_a_result(tmp_path: Path) -> None:
    source = tmp_path / "vazio.smi"
    source.write_text("", encoding="utf-8")
    assert run(["-i", str(source), "-o", str(tmp_path / "out"), "--quiet"]) == 1
