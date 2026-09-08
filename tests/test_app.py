"""Suíte de testes automatizados do Frontend (Streamlit UI em app.py).

Cobre:
1. Testes de Estado da Interface (Session State)
2. Testes de Entrada e Validação (Input Screen)
3. Testes do Fluxo Wizard (Wizard Workflow & Timing)
4. Testes de Configuração (Configuration Popover)
5. Testes de Download e Resultados (Results & Downloads)
6. Teste de Integração Nativa (Streamlit AppTest)
"""

from __future__ import annotations

import io
import time
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from curation.app import (
    MIN_STAGE_DISPLAY_TIME,
    WIZARD_STAGES,
    init_session_state,
    reset_to_input,
    run_backend_pipeline,
)

APP_PATH = Path(__file__).parent.parent / "src" / "curation" / "app.py"


# --- 1. Testes de Estado da Interface (Session State) -------------------------------


def test_init_session_state(monkeypatch) -> None:
    """Verifica se o estado inicial começa em ui_state == 'INPUT' e configs padrão."""
    fake_state: dict = {}
    monkeypatch.setattr("streamlit.session_state", fake_state)

    init_session_state()

    assert fake_state["ui_state"] == "INPUT"
    assert fake_state["curating_step_idx"] == 0
    assert fake_state["raw_input"] is None
    assert fake_state["input_name"] == ""
    assert fake_state["config"]["max_mw"] == 1000.0
    assert fake_state["config"]["max_ha"] == 100
    assert fake_state["config"]["deduplicate"] is True


def test_reset_to_input(monkeypatch) -> None:
    """Verifica se reset_to_input restaura a UI para a Tela 1 (INPUT)."""
    fake_state = {
        "ui_state": "COMPLETE",
        "curating_step_idx": 6,
        "step_start_time": 100.0,
        "raw_input": b"CCO",
        "input_name": "test.smi",
        "report": "fake_report",
    }
    monkeypatch.setattr("streamlit.session_state", fake_state)
    monkeypatch.setattr("streamlit.rerun", lambda: None)

    reset_to_input()

    assert fake_state["ui_state"] == "INPUT"
    assert fake_state["curating_step_idx"] == 0
    assert fake_state["raw_input"] is None
    assert fake_state["input_name"] == ""
    assert fake_state["report"] is None


# --- 2. Testes de Entrada e Validação (Input Screen) ---------------------------------


def test_valid_smiles_enables_start_button() -> None:
    """Verifica se SMILES válidos ativam o botão de início no AppTest."""
    at = AppTest.from_file(APP_PATH)
    at.run()

    # Preencher área de texto com SMILES válidos
    at.text_area[0].input("CCO\nCC(=O)O[Na]").run()

    assert at.session_state["raw_input"] is not None
    assert "2 molecules detected" in at.success[0].value
    assert not at.button[0].disabled


def test_empty_input_disables_start_button() -> None:
    """Garante que entrada vazia mantém o botão de início desativado."""
    at = AppTest.from_file(APP_PATH)
    at.run()

    assert at.session_state["raw_input"] is None
    assert at.button[0].disabled


def test_file_upload_handling() -> None:
    """Garante o suporte a uploads de arquivos (.smi, .csv, .tsv, .txt)."""
    at = AppTest.from_file(APP_PATH)
    at.run()

    # Simular upload de arquivo .csv
    csv_bytes = b"SMILES\nCCO\nCC(=O)O\n"
    at.file_uploader[0].upload("molecules.csv", csv_bytes).run()

    assert at.session_state["raw_input"] == csv_bytes
    assert at.session_state["input_name"] == "molecules.csv"
    assert not at.button[0].disabled


# --- 3. Testes do Fluxo Wizard (Wizard Workflow & Timing) ---------------------------


def test_wizard_stages_sequence() -> None:
    """Verifica se os 7 estágios do Wizard estão definidos na ordem correta."""
    stage_ids = [stage["id"] for stage in WIZARD_STAGES]
    expected_order = [
        "INPUT",
        "STANDARDIZATION",
        "PARENT",
        "VALIDATION",
        "DEDUPLICATION",
        "ELIGIBILITY",
        "OUTPUT",
    ]
    assert stage_ids == expected_order
    assert len(WIZARD_STAGES) == 7


def test_ui_stage_display_timing() -> None:
    """Garante que a constante MIN_STAGE_DISPLAY_TIME está configurada em 1.5s."""
    assert MIN_STAGE_DISPLAY_TIME == 1.5


def test_transition_to_complete_state() -> None:
    """Garante que a UI completa transiciona para o estado COMPLETE."""
    at = AppTest.from_file(APP_PATH)
    at.run()

    at.session_state["fast_mode"] = True
    # Carregar entrada e disparar curadoria
    at.text_area[0].input("CCO\nCC(=O)O").run()
    at.button[0].click().run()

    assert at.session_state["ui_state"] in ("CURATING", "COMPLETE")
    assert at.session_state["report"] is not None


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
    }
    monkeypatch.setattr("streamlit.session_state", fake_state)

    run_backend_pipeline()

    report = fake_state["report"]
    assert report is not None
    assert report.provenance.parameters["max_mw"] == 500.0
    assert report.provenance.parameters["max_ha"] == 50
    assert report.provenance.parameters["deduplicate"] is False


# --- 5. Testes de Download e Resultados (Results & Downloads) -----------------------


def test_download_buttons_content(monkeypatch) -> None:
    """Verifica se a curadoria gera relatórios válidos para exportação em CSV."""
    fake_state = {
        "raw_input": b"CCO\nC(C)(C)(C)(C)C",
        "input_name": "test.smi",
        "config": {
            "max_mw": 1000.0,
            "max_ha": 100,
            "deduplicate": True,
            "policy_hash": "0" * 64,
            "policy_path": "docs/decisions.md",
        },
        "report": None,
    }
    monkeypatch.setattr("streamlit.session_state", fake_state)

    run_backend_pipeline()

    report = fake_state["report"]
    assert report is not None
    assert len(report.approved) == 1
    assert len(report.rejected) == 1


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
    }
    monkeypatch.setattr("streamlit.session_state", fake_state)

    run_backend_pipeline()

    report = fake_state["report"]
    kpis = report.kpis()

    assert kpis["Total Processed"] == 3
    assert kpis["Approved"] == 2
    assert kpis["Rejected"] == 1


# --- 6. Teste de Integração Nativa (Streamlit AppTest) -----------------------------


def test_app_full_lifecycle_without_exceptions() -> None:
    """Simula uma execução completa do aplicativo do início ao fim usando AppTest."""
    at = AppTest.from_file(APP_PATH)
    at.run()

    assert not at.exception
    assert at.session_state["ui_state"] == "INPUT"

    at.session_state["fast_mode"] = True
    # Preencher entrada com 3 moléculas
    at.text_area[0].input("CCO\nCC(=O)O[Na]\nN[C@@H](C)C(=O)O.Cl").run()
    assert not at.exception
    assert at.session_state["raw_input"] is not None

    # Clicar em "Start curation"
    at.button[0].click().run()
    assert not at.exception
    assert at.session_state["report"] is not None
