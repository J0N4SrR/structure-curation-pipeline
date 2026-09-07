"""Benchmark contra padrão-ouro curado manualmente.

Este módulo mede **curadoria**, não encanamento. Comparar a saída do pipeline com a
saída do ``chembl_structure_pipeline`` seria tautológico: o motor é o mesmo, então
100% de coincidência é o resultado esperado de chamar a mesma função duas vezes, e
qualquer divergência seria bug do wrapper. O único referencial com conteúdo é o
julgamento humano.

Convenção da matriz de confusão — declarada explicitamente porque a métrica
principal depende dela. A **classe positiva é REJECT**: a tarefa do pipeline é
detectar estruturas problemáticas.

===================  ==========================  ==========================
                     especialistas: REJECT       especialistas: ACCEPT
===================  ==========================  ==========================
pipeline: REJECT     verdadeiro positivo         **falso positivo**
pipeline: ACCEPT     falso negativo              verdadeiro negativo
===================  ==========================  ==========================

A *taxa de falsos rejeitados* (FRR) é ``FP / (FP + TN)``: a fração dos compostos
que os especialistas consideraram aceitáveis e que o pipeline descartou. É a
métrica mais consequente do trabalho, porque um composto descartado desaparece do
dataset sem deixar rastro no modelo treinado depois.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal, Optional, Sequence, Union

Label = Literal["ACCEPT", "REJECT"]
VALID_LABELS = ("ACCEPT", "REJECT")

GOLD_COLUMNS = (
    "compound_id",
    "raw_smiles",
    "annotator_a",
    "annotator_b",
    "consensus",
    "reason",
)


class GoldStandardError(ValueError):
    """Problema estrutural no arquivo de padrão-ouro."""


@dataclass(frozen=True)
class GoldStandardEntry:
    """Um composto julgado independentemente por dois anotadores.

    Attributes:
        consensus: rótulo acordado. Fica vazio quando os anotadores divergiram e a
            divergência ainda não foi resolvida — esses registros entram no cálculo
            de κ, mas não na matriz de confusão, porque não existe verdade contra a
            qual comparar.
    """

    compound_id: str
    raw_smiles: str
    annotator_a: Label
    annotator_b: Label
    consensus: Optional[Label] = None
    reason: str = ""

    @property
    def annotators_agree(self) -> bool:
        return self.annotator_a == self.annotator_b

    @property
    def resolved_label(self) -> Optional[Label]:
        """Verdade utilizável: o consenso, ou a concordância espontânea."""
        if self.consensus:
            return self.consensus
        return self.annotator_a if self.annotators_agree else None


def _normalize_label(value: str, field: str, row: int) -> Optional[Label]:
    text = (value or "").strip().upper()
    if not text:
        return None
    if text not in VALID_LABELS:
        raise GoldStandardError(
            f"linha {row}: {field}={value!r} invalido; use ACCEPT ou REJECT"
        )
    return text  # type: ignore[return-value]


def load_gold_standard(path: Union[str, Path]) -> list[GoldStandardEntry]:
    """Carrega o conjunto anotado.

    O arquivo é um CSV com as colunas de :data:`GOLD_COLUMNS`. ``consensus`` e
    ``reason`` podem ficar em branco; os dois rótulos de anotador são obrigatórios,
    porque sem eles não há concordância inter-anotador a medir.
    """
    source = Path(path)
    entries: list[GoldStandardEntry] = []

    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = {"compound_id", "raw_smiles", "annotator_a", "annotator_b"} - set(
            reader.fieldnames or []
        )
        if missing:
            raise GoldStandardError(f"colunas ausentes: {sorted(missing)}")

        for number, row in enumerate(reader, start=2):
            label_a = _normalize_label(row.get("annotator_a", ""), "annotator_a", number)
            label_b = _normalize_label(row.get("annotator_b", ""), "annotator_b", number)
            if label_a is None or label_b is None:
                raise GoldStandardError(
                    f"linha {number}: ambos os anotadores precisam ter rotulo"
                )
            entries.append(
                GoldStandardEntry(
                    compound_id=(row.get("compound_id") or "").strip(),
                    raw_smiles=(row.get("raw_smiles") or "").strip(),
                    annotator_a=label_a,
                    annotator_b=label_b,
                    consensus=_normalize_label(
                        row.get("consensus", ""), "consensus", number
                    ),
                    reason=(row.get("reason") or "").strip(),
                )
            )

    if not entries:
        raise GoldStandardError("conjunto vazio")
    return entries


def write_annotation_template(
    records: Iterable[tuple[str, str]], path: Union[str, Path]
) -> int:
    """Gera o CSV em branco para os anotadores preencherem.

    Existe para que a anotação possa começar **em paralelo** ao desenvolvimento. É a
    única etapa do projeto com dependência humana de prazo longo — recrutar, treinar
    e calibrar um segundo anotador leva semanas — e adiá-la para depois do código
    pronto é o que comprime a validação até ela perder valor.
    """
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(GOLD_COLUMNS))
        writer.writeheader()
        for compound_id, smiles in records:
            writer.writerow(
                {
                    "compound_id": compound_id,
                    "raw_smiles": smiles,
                    "annotator_a": "",
                    "annotator_b": "",
                    "consensus": "",
                    "reason": "",
                }
            )
            count += 1
    return count


# --- Concordância inter-anotador ----------------------------------------------------


@dataclass(frozen=True)
class AgreementReport:
    """Concordância entre os dois anotadores."""

    n: int
    observed_agreement: float
    expected_agreement: float
    kappa: float
    disagreements: int
    unresolved: int

    @property
    def interpretation(self) -> str:
        """Faixas de Landis & Koch (1977), a convenção usual na literatura."""
        thresholds = (
            (0.81, "quase perfeita"),
            (0.61, "substancial"),
            (0.41, "moderada"),
            (0.21, "razoavel"),
            (0.00, "leve"),
        )
        for cutoff, label in thresholds:
            if self.kappa >= cutoff:
                return label
        return "pobre (pior que o acaso)"

    def to_markdown(self) -> str:
        return "\n".join(
            [
                "| métrica | valor |",
                "| --- | ---: |",
                f"| compostos anotados | {self.n} |",
                f"| concordância observada | {self.observed_agreement:.3f} |",
                f"| concordância esperada por acaso | {self.expected_agreement:.3f} |",
                f"| **κ de Cohen** | **{self.kappa:.3f}** ({self.interpretation}) |",
                f"| divergências | {self.disagreements} |",
                f"| divergências não resolvidas | {self.unresolved} |",
            ]
        )


def cohens_kappa(labels_a: Sequence[str], labels_b: Sequence[str]) -> float:
    """κ de Cohen para dois anotadores sobre categorias nominais.

    ``κ = (Po - Pe) / (1 - Pe)``, onde ``Po`` é a concordância observada e ``Pe`` a
    esperada por acaso. Vale 1 para concordância perfeita, 0 para concordância no
    nível do acaso e negativo para pior que o acaso.

    Quando os dois anotadores usam uma única categoria e concordam em tudo, ``Pe``
    vale 1 e κ é indefinido; a convenção adotada é devolver 1.0, porque a
    concordância é de fato total.
    """
    if len(labels_a) != len(labels_b):
        raise ValueError("as duas sequencias precisam ter o mesmo tamanho")
    total = len(labels_a)
    if total == 0:
        raise ValueError("nenhum par de rotulos")

    observed = sum(a == b for a, b in zip(labels_a, labels_b)) / total

    counts_a = Counter(labels_a)
    counts_b = Counter(labels_b)
    expected = sum(
        (counts_a[category] / total) * (counts_b[category] / total)
        for category in set(counts_a) | set(counts_b)
    )

    if expected >= 1.0:
        return 1.0
    return (observed - expected) / (1.0 - expected)


def agreement_report(entries: Sequence[GoldStandardEntry]) -> AgreementReport:
    labels_a = [entry.annotator_a for entry in entries]
    labels_b = [entry.annotator_b for entry in entries]
    disagreements = sum(1 for entry in entries if not entry.annotators_agree)
    unresolved = sum(
        1 for entry in entries if not entry.annotators_agree and not entry.consensus
    )

    total = len(entries)
    observed = sum(a == b for a, b in zip(labels_a, labels_b)) / total
    counts_a, counts_b = Counter(labels_a), Counter(labels_b)
    expected = sum(
        (counts_a[c] / total) * (counts_b[c] / total)
        for c in set(counts_a) | set(counts_b)
    )

    return AgreementReport(
        n=total,
        observed_agreement=observed,
        expected_agreement=expected,
        kappa=cohens_kappa(labels_a, labels_b),
        disagreements=disagreements,
        unresolved=unresolved,
    )


# --- Desempenho do pipeline ----------------------------------------------------------


@dataclass(frozen=True)
class ConfusionMatrix:
    """Desempenho do pipeline contra o consenso dos especialistas.

    Classe positiva: **REJECT**. Ver a convenção no docstring do módulo.
    """

    true_positive: int
    false_positive: int
    true_negative: int
    false_negative: int
    excluded_unresolved: int = 0

    @property
    def total(self) -> int:
        return (
            self.true_positive
            + self.false_positive
            + self.true_negative
            + self.false_negative
        )

    @staticmethod
    def _ratio(numerator: int, denominator: int) -> float:
        return numerator / denominator if denominator else 0.0

    @property
    def precision(self) -> float:
        """Das rejeições do pipeline, quantas os especialistas confirmam."""
        return self._ratio(
            self.true_positive, self.true_positive + self.false_positive
        )

    @property
    def recall(self) -> float:
        """Dos problemas reais, quantos o pipeline detecta."""
        return self._ratio(
            self.true_positive, self.true_positive + self.false_negative
        )

    @property
    def f1(self) -> float:
        denominator = self.precision + self.recall
        return 2 * self.precision * self.recall / denominator if denominator else 0.0

    @property
    def false_rejection_rate(self) -> float:
        """Fração dos compostos aceitáveis que o pipeline descartou.

        A métrica mais consequente: cada unidade aqui é um composto válido que
        sumiu do dataset sem deixar rastro no modelo treinado depois.
        """
        return self._ratio(
            self.false_positive, self.false_positive + self.true_negative
        )

    @property
    def accuracy(self) -> float:
        return self._ratio(self.true_positive + self.true_negative, self.total)

    def to_markdown(self) -> str:
        lines = [
            "| | especialistas: REJECT | especialistas: ACCEPT |",
            "| --- | ---: | ---: |",
            f"| **pipeline: REJECT** | {self.true_positive} (VP) | "
            f"{self.false_positive} (FP) |",
            f"| **pipeline: ACCEPT** | {self.false_negative} (FN) | "
            f"{self.true_negative} (VN) |",
            "",
            "| métrica | valor |",
            "| --- | ---: |",
            f"| precisão | {self.precision:.3f} |",
            f"| recall | {self.recall:.3f} |",
            f"| F1 | {self.f1:.3f} |",
            f"| acurácia | {self.accuracy:.3f} |",
            f"| **taxa de falsos rejeitados (FRR)** | **{self.false_rejection_rate:.3f}** |",
        ]
        if self.excluded_unresolved:
            lines.append(
                f"| divergências não resolvidas (fora da matriz) | "
                f"{self.excluded_unresolved} |"
            )
        return "\n".join(lines)


def evaluate_pipeline(
    entries: Sequence[GoldStandardEntry], pipeline
) -> ConfusionMatrix:
    """Roda o pipeline sobre o conjunto anotado e monta a matriz de confusão.

    Registros cuja divergência entre anotadores não foi resolvida são **excluídos**
    da matriz e contados à parte: sem verdade acordada não há contra o que comparar,
    e escolher um dos anotadores arbitrariamente inflaria a métrica.
    """
    counts = {"tp": 0, "fp": 0, "tn": 0, "fn": 0}
    unresolved = 0

    for entry in entries:
        truth = entry.resolved_label
        if truth is None:
            unresolved += 1
            continue

        record = pipeline.process_single(entry.raw_smiles, entry.compound_id)
        predicted: Label = "ACCEPT" if record.passed else "REJECT"

        if truth == "REJECT":
            counts["tp" if predicted == "REJECT" else "fn"] += 1
        else:
            counts["fp" if predicted == "REJECT" else "tn"] += 1

    return ConfusionMatrix(
        true_positive=counts["tp"],
        false_positive=counts["fp"],
        true_negative=counts["tn"],
        false_negative=counts["fn"],
        excluded_unresolved=unresolved,
    )


def benchmark_report(
    entries: Sequence[GoldStandardEntry], pipeline, title: str = "Benchmark"
) -> str:
    """Relatório completo em Markdown, pronto para a documentação da tese."""
    agreement = agreement_report(entries)
    matrix = evaluate_pipeline(entries, pipeline)

    sections = [
        f"# {title}",
        "",
        "## Concordância inter-anotador",
        "",
        agreement.to_markdown(),
        "",
        "## Desempenho do pipeline contra o consenso",
        "",
        "Classe positiva: `REJECT`.",
        "",
        matrix.to_markdown(),
    ]

    if agreement.kappa < 0.61:
        sections += [
            "",
            "> **Atenção:** κ abaixo de 0,61 indica concordância inter-anotador "
            "insuficiente. As métricas de desempenho acima carregam essa "
            "incerteza — calibrar os anotadores antes de reportá-las.",
        ]
    return "\n".join(sections)
