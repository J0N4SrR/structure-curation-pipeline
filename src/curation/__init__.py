"""Pipeline de curadoria e padronização estrutural.

As decisões de política que governam o comportamento químico estão em
``docs/decisions.md`` e são identificadas pelo ``policy_hash`` carimbado em cada
registro produzido.
"""

from curation.engine import PIPELINE_VERSION, EngineWrapper
from curation.dedup import Collision, CollisionType, DedupIndex
from curation.filters import EligibilityCriteria
from curation.pipeline import BatchSummary, CurationPipeline
from curation.models import (
    CurationRecord,
    RejectionCode,
    Stage,
    TransformationEvent,
)

__all__ = [
    "BatchSummary",
    "Collision",
    "CollisionType",
    "DedupIndex",
    "EligibilityCriteria",
    "CurationPipeline",
    "CurationRecord",
    "EngineWrapper",
    "PIPELINE_VERSION",
    "RejectionCode",
    "Stage",
    "TransformationEvent",
]
