"""Harness do estudo de validação científica (Fase 5).

Dois módulos independentes:

* :mod:`~curation.validation.benchmarking` — desempenho contra um padrão-ouro
  curado por dois anotadores humanos, com κ de Cohen e taxa de falsos rejeitados.
* :mod:`~curation.validation.ablation` — impacto a jusante de cada decisão de
  política, medido em deduplicação, propriedades físico-químicas e conservação de
  estereoquímica.

Nenhum dos dois compara a saída do pipeline com a da biblioteca de referência: como
o motor é o mesmo, essa comparação é tautológica e mede o wrapper, não a curadoria.
"""

from curation.validation.ablation import (
    PolicyVariant,
    VariantResult,
    default_variants,
    run_ablation,
    write_report,
)
from curation.validation.benchmarking import (
    AgreementReport,
    ConfusionMatrix,
    GoldStandardEntry,
    agreement_report,
    benchmark_report,
    cohens_kappa,
    evaluate_pipeline,
    load_gold_standard,
    write_annotation_template,
)

__all__ = [
    "AgreementReport",
    "ConfusionMatrix",
    "GoldStandardEntry",
    "PolicyVariant",
    "VariantResult",
    "agreement_report",
    "benchmark_report",
    "cohens_kappa",
    "default_variants",
    "evaluate_pipeline",
    "load_gold_standard",
    "run_ablation",
    "write_report",
]
