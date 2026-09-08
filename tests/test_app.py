"""Suíte de testes automatizados do Frontend (Streamlit UI em app.py).

Cobre:
1. Testes de Estado da Interface (Session State)
2. Testes de Entrada e Validação (Workbench)
3. Testes do Fluxo do Orquestrador (Live DAG)
4. Testes de Configuração (Popover)
5. Testes de Download e Resultados
6. Teste de Integração Nativa (Streamlit AppTest)
"""

from __future__ import annotations

import io
import time
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from curation.app import (
    init_session_state,
    execute_pipeline_live,
    PIPELINE_STAGES,
)
import streamlit as st

APP_PATH = Path(__file__).parent.parent / "src" / "curation" / "app.py"


# --- 1. Testes de Estado da Interface (Session State) -------------------------------


def test_init_session_state(monkeypatch) -> None:
    """Garante que todas as chaves requeridas são criadas."""
    fake_state: dict[str, object] = {}
    monkeypatch.setattr(st, "session_state", fake_state)

    init_session_state()

    assert fake_state["raw_input"] is None
    assert fake_state["report"] is None
    assert "config" in fake_state
    assert "completed_cards" in fake_state


# --- 2. Testes de Entrada e Validação (Workbench) ---------------------------------


def test_valid_smiles_enables_start_button() -> None:
    """Verifica se SMILES válidos ativam o botão de início no AppTest."""
    at = AppTest.from_file(str(APP_PATH))
    at.run()

    # Expand the expander "📂 Ingestão de Moléculas (SMILES)" which is the first expander
    # Actually just write to the text area
    at.text_area[0].input("CCO\nCC(=O)O[Na]").run()
    at.run() # Rerun to update button disabled state since it is rendered above the input area

    assert at.session_state["raw_input"] is not None
    assert "2 moléculas detectadas" in at.success[0].value
    assert not at.button[0].disabled


def test_empty_input_disables_start_button() -> None:
    """Garante que entrada vazia mantém o botão de início desativado."""
    at = AppTest.from_file(str(APP_PATH))
    at.run()

    assert at.session_state["raw_input"] is None
    assert at.button[0].disabled


def test_file_upload_handling() -> None:
    """Garante o suporte a uploads de arquivos (.smi, .csv, .tsv, .txt)."""
    at = AppTest.from_file(str(APP_PATH))
    at.run()

    # Simular upload de arquivo .csv
    csv_bytes = b"SMILES\nCCO\nCC(=O)O\n"
    at.file_uploader[0].upload("molecules.csv", csv_bytes).run()
    at.run() # Rerun to update button disabled state since it is rendered above the input area

    assert at.session_state["raw_input"] == csv_bytes
    assert at.session_state["input_name"] == "molecules.csv"
    assert not at.button[0].disabled


# --- 3. Testes do Fluxo do Orquestrador (Live DAG) ---------------------------


def test_ordered_stages_sequence() -> None:
    """Verifica se os 7 estágios do DAG estão definidos na ordem correta."""
    stage_ids = [stage[0] for stage in PIPELINE_STAGES]
    expected_order = [
        "PARSE",
        "STANDARDIZE",
        "GET_PARENT",
        "VALENCE_GATE",
        "ELIGIBILITY",
        "CANONICALIZE",
        "DEDUPLICATION",
    ]
    assert stage_ids == expected_order
    assert len(PIPELINE_STAGES) == 7


def test_execute_pipeline_live_sets_report(monkeypatch) -> None:
    """Garante que a execução gera os report states."""
    fake_state = {
        "raw_input": b"CCO\nCC(=O)O",
        "input_name": "test.smi",
        "config": {
            "max_mw": 500.0,
            "max_ha": 50,
            "deduplicate": False,
            "policy_hash": "0" * 64,
            "policy_path": "docs/decisions.md",
        },
        "report": None,
    }
    monkeypatch.setattr("streamlit.session_state", fake_state)
    monkeypatch.setattr("streamlit.rerun", lambda: None)
    # Patch time.sleep to run fast
    monkeypatch.setattr(time, "sleep", lambda x: None)
    
    class DummyDagContainer:
        def __enter__(self): return self
        def __exit__(self, *args): pass
        def markdown(self, *args, **kwargs): pass
    
    execute_pipeline_live(DummyDagContainer())

    assert fake_state["report"] is not None
    assert "completed_cards" in fake_state


# --- 4. Testes de Configuração (Configuration Popover) ------------------------------


def test_config_updates(monkeypatch) -> None:
    """Verifica se alterações de configuração atualizam os parâmetros do pipeline."""
    fake_state = {
        "raw_input": b"CCO",
        "input_name": "test.smi",
        "config": {
            "max_mw": 500.0,
            "max_ha": 50,
            "deduplicate": False,
            "policy_hash": "0" * 64,
            "policy_path": "docs/decisions.md",
        },
        "report": None,
        "input_name": "manual_input.smi"
    }
    monkeypatch.setattr("streamlit.session_state", fake_state)
    monkeypatch.setattr("streamlit.rerun", lambda: None)
    monkeypatch.setattr(time, "sleep", lambda x: None)

    class DummyDagContainer:
        def __enter__(self): return self
        def __exit__(self, *args): pass
        def markdown(self, *args, **kwargs): pass

    execute_pipeline_live(DummyDagContainer())

    report = fake_state["report"]
    assert report is not None
    assert report.provenance.parameters["max_mw"] == 500.0
    assert report.provenance.parameters["max_ha"] == 50
    assert report.provenance.parameters["deduplicate"] is False


# --- 5. Testes de Download e Resultados (Results & Downloads) -----------------------


def test_summary_kpi_matching(monkeypatch) -> None:
    """Garante que a contagem dos KPIs calculados pelo backend corresponde ao total."""
    fake_state = {
        "raw_input": b"CCO\nCC(=O)O[Na]\nC(C)(C)(C)(C)C",
        "input_name": "test.smi",
        "config": {
            "max_mw": 1000.0,
            "max_ha": 100,
            "deduplicate": True,
            "policy_hash": "0" * 64,
            "policy_path": "docs/decisions.md",
        },
        "report": None,
        "input_name": "manual_input.smi"
    }
    monkeypatch.setattr("streamlit.session_state", fake_state)
    monkeypatch.setattr("streamlit.rerun", lambda: None)
    monkeypatch.setattr(time, "sleep", lambda x: None)

    class DummyDagContainer:
        def __enter__(self): return self
        def __exit__(self, *args): pass
        def markdown(self, *args, **kwargs): pass

    execute_pipeline_live(DummyDagContainer())

    report = fake_state["report"]
    kpis = report.kpis()

    assert kpis["Total Processed"] == 3
    assert kpis["Approved"] == 2
    assert kpis["Rejected"] == 1


# --- 6. Teste de Integração Nativa (Streamlit AppTest) -----------------------------

@pytest.mark.skip(reason="time.sleep inside AppTest causes unpredictable timeouts, logic tested in isolated unit test")
def test_app_full_lifecycle_without_exceptions() -> None:
    """Simula uma execução completa do aplicativo do início ao fim usando AppTest."""
    at = AppTest.from_file(str(APP_PATH))
    at.run()
    
    assert not at.exception


def test_selected_node_starts_empty() -> None:
    """Sem selecao, o painel mostra o resumo da execucao."""
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_file(str(APP_PATH), default_timeout=90).run()
    assert not app.exception
    assert "inspect_stage" not in app.session_state or (
        app.session_state["inspect_stage"] is None
    )

def test_node_details_panel_uses_the_view_model() -> None:
    """O painel consome NodeContract; nao recalcula metrica nenhuma."""
    from curation.pipeline import CurationPipeline
    from curation.viewmodel import build_run_view

    source = "CCO\nCC(=O)O[Na]\nC(C)(C)(C)(C)C\n"
    report = CurationPipeline("0" * 64).run_report(
        source, input_bytes=source.encode()
    )
    view = build_run_view(report)

    node = view.graph.node("STANDARDIZE")
    assert node is not None
    assert node.input_count == report.stages[1].n_input
    assert node.output_count == report.stages[1].n_output
    assert node.rejected_count == report.stages[1].n_excluded
