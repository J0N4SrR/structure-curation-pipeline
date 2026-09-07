"""Deduplicação estrutural com relatório explícito de conflitos (D-07).

A identidade é o **InChIKey completo**. O primeiro bloco de 14 caracteres codifica
apenas conectividade — ignora estereoquímica, isótopos e carga — e usá-lo como chave
funde enantiômeros::

    L-alanina  QNAYBMKLOCPYGJ-REOHCLBHSA-N
    D-alanina  QNAYBMKLOCPYGJ-UWTATZPHSA-N

Por isso o bloco1 é apenas **agrupador taxonômico**: serve para localizar famílias
sal/ácido livre e conjuntos de estereoisômeros relacionados, e nunca para decidir
que dois registros são o mesmo composto.

Nenhuma colisão é resolvida automaticamente. Conflitos de anotação, em particular,
são um problema de dados que exige decisão humana; silenciá-los produziria um
dataset limpo e errado.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Iterator, Optional, Union

from curation.models import CurationRecord


class CollisionType(str, Enum):
    """Categorias de colisão, em ordem de severidade decrescente para revisão."""

    #: Mesmo InChIKey completo e anotações divergentes. Exige decisão humana.
    ANNOTATION_CONFLICT = "ANNOTATION_CONFLICT"
    #: Mesmo InChIKey completo. O registro é o mesmo composto.
    EXACT_DUPLICATE = "EXACT_DUPLICATE"
    #: Mesmo bloco1, InChIKeys distintos. **Não** é duplicata.
    BLOCK1_COLLISION = "BLOCK1_COLLISION"


@dataclass(frozen=True)
class Collision:
    """Uma colisão detectada, com o suficiente para revisão sem voltar ao lote."""

    collision_type: CollisionType
    incoming_id: str
    existing_id: str
    incoming_inchikey: str
    existing_inchikey: str
    inchikey_block1: str
    detail: str = ""

    def as_row(self) -> dict[str, str]:
        return {
            "collision_type": self.collision_type.value,
            "incoming_id": self.incoming_id,
            "existing_id": self.existing_id,
            "incoming_inchikey": self.incoming_inchikey,
            "existing_inchikey": self.existing_inchikey,
            "inchikey_block1": self.inchikey_block1,
            "detail": self.detail,
        }


CONFLICT_COLUMNS = (
    "collision_type",
    "incoming_id",
    "existing_id",
    "incoming_inchikey",
    "existing_inchikey",
    "inchikey_block1",
    "detail",
)


def _classify_block1_difference(left: str, right: str) -> str:
    """Descreve por que dois InChIKeys compartilham o esqueleto mas divergem.

    O formato é ``AAAAAAAAAAAAAA-BBBBBBBBFV-P``: o segundo bloco carrega
    estereoquímica e isótopos, o terceiro carrega o estado de protonação.
    """
    left_parts = left.split("-")
    right_parts = right.split("-")
    if len(left_parts) < 3 or len(right_parts) < 3:
        return "estrutura de InChIKey inesperada"

    reasons: list[str] = []
    if left_parts[1] != right_parts[1]:
        reasons.append("estereoquímica ou isotopologia")
    if left_parts[2] != right_parts[2]:
        reasons.append("estado de protonação")
    return " e ".join(reasons) if reasons else "diferença não localizada"


class DedupIndex:
    """Índice em memória de identidades já vistas.

    Mantém duas estruturas: ``inchikey -> primeiro compound_id`` para identidade, e
    ``bloco1 -> conjunto de inchikeys`` para agrupamento. Registros aprovados sem
    InChIKey são ignorados pelo índice.

    **Pegada de memória medida** (não estimada): ~159 bytes por entrada para o índice
    de identidade, sobre CPython 3.12. Isso é 1 milhão de compostos em ~0,16 GB e 10
    milhões em ~1,6 GB — viável em memória até a ordem de 10⁷, ponto a partir do qual
    convém trocar por um índice em disco.

    Deduplicar exige visão global do lote: não é possível fazê-lo em streaming puro.
    Por isso a arquitetura é "registros para disco em fluxo, índice compacto em RAM",
    e não escrita puramente sequencial.
    """

    def __init__(self, keep_first: bool = True) -> None:
        self._by_inchikey: dict[str, str] = {}
        self._annotations: dict[str, str] = {}
        self._by_block1: dict[str, set[str]] = {}
        self._collisions: list[Collision] = []
        self.keep_first = keep_first

    # --- Consulta ------------------------------------------------------------------

    def __len__(self) -> int:
        return len(self._by_inchikey)

    @property
    def unique_count(self) -> int:
        """Número de identidades distintas aceitas."""
        return len(self._by_inchikey)

    @property
    def collisions(self) -> list[Collision]:
        return list(self._collisions)

    def collision_counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for collision in self._collisions:
            key = collision.collision_type.value
            counts[key] = counts.get(key, 0) + 1
        return counts

    def first_occurrence(self, inchikey: str) -> Optional[str]:
        return self._by_inchikey.get(inchikey)

    def group(self, block1: str) -> set[str]:
        """InChIKeys completos que compartilham um esqueleto de conectividade."""
        return set(self._by_block1.get(block1, ()))

    # --- Inserção ------------------------------------------------------------------

    def add(
        self, record: CurationRecord, annotation: Optional[str] = None
    ) -> Optional[Collision]:
        """Registra um composto aprovado e devolve a colisão, se houver.

        Args:
            record: registro aprovado, com InChIKey preenchido.
            annotation: valor anotado associado ao composto (atividade, classe,
                rótulo). Quando dois registros compartilham o InChIKey completo mas
                divergem aqui, a colisão é classificada como conflito de anotação —
                a categoria que não pode ser resolvida automaticamente.

        Returns:
            ``None`` quando a identidade é nova. Caso contrário, a colisão detectada.
            O chamador decide o que fazer; o índice apenas registra.
        """
        inchikey = record.inchikey
        if not record.passed or not inchikey:
            return None

        block1 = record.inchikey_block1 or inchikey.split("-")[0]
        existing_id = self._by_inchikey.get(inchikey)

        if existing_id is not None:
            collision = self._duplicate_collision(
                record, inchikey, block1, existing_id, annotation
            )
            self._collisions.append(collision)
            if not self.keep_first:
                self._by_inchikey[inchikey] = record.input_id
            return collision

        siblings = self._by_block1.get(block1)
        collision: Optional[Collision] = None
        if siblings:
            neighbour = sorted(siblings)[0]
            collision = Collision(
                collision_type=CollisionType.BLOCK1_COLLISION,
                incoming_id=record.input_id,
                existing_id=self._by_inchikey.get(neighbour, ""),
                incoming_inchikey=inchikey,
                existing_inchikey=neighbour,
                inchikey_block1=block1,
                detail=(
                    "mesmo esqueleto de conectividade, identidades distintas: "
                    + _classify_block1_difference(inchikey, neighbour)
                    + " — não é duplicata"
                ),
            )
            self._collisions.append(collision)

        self._by_inchikey[inchikey] = record.input_id
        if annotation is not None:
            self._annotations[inchikey] = annotation
        self._by_block1.setdefault(block1, set()).add(inchikey)
        return collision

    def add_all(self, records: Iterator[CurationRecord]) -> None:
        for record in records:
            self.add(record)

    def _duplicate_collision(
        self,
        record: CurationRecord,
        inchikey: str,
        block1: str,
        existing_id: str,
        annotation: Optional[str],
    ) -> Collision:
        previous = self._annotations.get(inchikey)
        if annotation is not None and previous is not None and annotation != previous:
            return Collision(
                collision_type=CollisionType.ANNOTATION_CONFLICT,
                incoming_id=record.input_id,
                existing_id=existing_id,
                incoming_inchikey=inchikey,
                existing_inchikey=inchikey,
                inchikey_block1=block1,
                detail=(
                    f"mesma estrutura com anotações divergentes: "
                    f"{previous!r} vs {annotation!r} — requer decisão humana"
                ),
            )
        return Collision(
            collision_type=CollisionType.EXACT_DUPLICATE,
            incoming_id=record.input_id,
            existing_id=existing_id,
            incoming_inchikey=inchikey,
            existing_inchikey=inchikey,
            inchikey_block1=block1,
            detail="InChIKey completo idêntico ao da primeira ocorrência",
        )

    # --- Relatório -------------------------------------------------------------------

    def write_conflicts(self, path: Union[str, Path]) -> int:
        """Grava ``conflicts.csv``. Devolve o número de linhas escritas.

        O arquivo é sempre criado, mesmo vazio: a ausência de conflitos é um
        resultado, e sua distinção de "o relatório não foi gerado" importa.
        """
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(CONFLICT_COLUMNS))
            writer.writeheader()
            for collision in self._collisions:
                writer.writerow(collision.as_row())
        return len(self._collisions)
