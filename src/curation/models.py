"""Schema de proveniência do pipeline de curadoria.

Os modelos aqui definem o contrato de saída. Toda decisão de política registrada em
``docs/decisions.md`` que produza efeito observável precisa ter um campo ou um evento
correspondente neste módulo — caso contrário a decisão não é auditável.
"""

from __future__ import annotations

from enum import Enum
from typing import Literal, Optional

from pydantic import BaseModel, ConfigDict, Field


class Stage(str, Enum):
    """Fronteiras do pipeline nas quais as sondas são executadas.

    A ordem dos membros é a ordem de execução e é significativa: ``rejection_stage``
    é comparado contra ela para verificar a restrição de ordem da D-09 (instrumentar
    antes do motor, rejeitar depois).
    """

    PARSE = "PARSE"
    STANDARDIZE = "STANDARDIZE"
    GET_PARENT = "GET_PARENT"
    VALENCE_GATE = "VALENCE_GATE"
    CANONICALIZE = "CANONICALIZE"


class RejectionCode(str, Enum):
    """Causas de descarte.

    Conjunto fechado: um código novo é uma mudança de contrato e exige atualização
    do relatório de lote, que agrega por esta enumeração.
    """

    ERR_SYNTAX = "ERR_SYNTAX"
    ERR_VALENCE = "ERR_VALENCE"
    ERR_KEKULIZE = "ERR_KEKULIZE"
    ERR_SANITIZE = "ERR_SANITIZE"
    ERR_STANDARDIZE = "ERR_STANDARDIZE"
    ERR_GET_PARENT = "ERR_GET_PARENT"
    ERR_CANONICALIZE = "ERR_CANONICALIZE"
    ERR_INCHI = "ERR_INCHI"
    ERR_EMPTY = "ERR_EMPTY"
    ERR_INTERNAL = "ERR_INTERNAL"


class TransformationEvent(BaseModel):
    """Uma transformação observada na fronteira entre dois estágios.

    Registra apenas o que foi *observado* comparando os estados antes e depois. O
    motor é uma caixa-preta de oito passos internos (D-01), então nenhum evento aqui
    afirma qual regra interna disparou — apenas que o efeito é compatível com ela.
    """

    model_config = ConfigDict(frozen=True)

    stage: Stage
    rule: str
    before_smiles: str
    after_smiles: str
    detail: str = ""


class CurationRecord(BaseModel):
    """Resultado da curadoria de um único composto.

    Um registro é emitido para *toda* entrada, aprovada ou rejeitada. A segregação
    em ``curated.csv`` / ``rejected.csv`` é uma projeção desta estrutura, nunca a
    fonte primária.
    """

    model_config = ConfigDict(frozen=True)

    # --- Identidade ---
    input_id: str
    raw_smiles: str
    curated_smiles: Optional[str] = None
    inchikey: Optional[str] = None
    inchikey_block1: Optional[str] = None

    # --- Decisão ---
    status: Literal["PASSED", "REJECTED"]
    rejection_code: Optional[RejectionCode] = None
    rejection_stage: Optional[Stage] = None
    rejection_detail: Optional[str] = None

    # --- Rastreabilidade ---
    transformations: list[TransformationEvent] = Field(default_factory=list)

    # --- Diferenciais (entrada -> parent) ---
    removed_fragments: Optional[str] = None
    delta_net_charge: int = 0
    delta_formula: str = ""
    delta_stereocenters: int = 0
    n_undefined_stereocenters: int = 0
    n_stereocenters_total: int = 0
    n_components_parent: int = 0

    # --- Metadados de proveniência ---
    excluded_flag: bool = False
    rdkit_version: str
    csp_version: str
    pipeline_version: str
    policy_hash: str

    @property
    def passed(self) -> bool:
        return self.status == "PASSED"

    def rule_names(self) -> list[str]:
        """Regras observadas, na ordem em que ocorreram."""
        return [event.rule for event in self.transformations]

    def has_rule(self, rule: str) -> bool:
        return any(event.rule == rule for event in self.transformations)
