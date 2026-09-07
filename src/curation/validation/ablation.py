"""Estudo de ablação: qual decisão de curadoria importa, e quanto.

Cada variante altera **uma** política de ``docs/decisions.md`` e mede o efeito a
jusante. É a pergunta publicável do trabalho, e é diferente de "a curadoria mudou os
descritores?" — que é aritmética, já que remover um contra-íon obviamente reduz o
peso molecular. O que tem conteúdo é comparar *escolhas de política entre si*, com
tudo o mais mantido constante.

Métricas por variante:

(a) **Deduplicação** — quantas identidades distintas restam. Sensível a tautomeria e
    a estereoquímica, portanto é onde D-02 e D-04 aparecem.
(b) **Distribuição físico-química** — MW, LogP e TPSA da estrutura-mãe.
(c) **Conservação de centros quirais** — quantos centros definidos sobrevivem.

Sobre a variante de tartaratos: ``flatten_tartrate_mol`` é chamada
incondicionalmente por ``standardize_mol`` e a API pública não permite desativá-la.
Em vez de monkeypatch, a variante **recompõe a sequência** a partir das próprias
funções do módulo de referência, omitindo um passo. É um desvio controlado, isolado
no harness experimental, que nunca toca o caminho de produção — que é justamente o
que a D-01 protege.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from statistics import mean, median, pstdev
from typing import Callable, Iterable, Optional, Sequence, Union

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

import chembl_structure_pipeline as csp
from chembl_structure_pipeline import standardizer as _ref
from chembl_structure_pipeline.exclude_flag import exclude_flag

from curation.dedup import DedupIndex
from curation.engine import EngineWrapper
from curation.filters import EligibilityCriteria
from curation.models import CurationRecord
from curation.pipeline import CurationPipeline

POLICY_HASH_PLACEHOLDER = "ABLATION"


# --- Padronizadores alternativos -------------------------------------------------


def standardize_without_tartrate_flattening(mol: Chem.Mol) -> Chem.Mol:
    """``standardize_mol`` sem o passo ``flatten_tartrate_mol`` (ablação da D-02).

    Reproduz fielmente a sequência da referência, omitindo um único passo. Se a
    biblioteca mudar a ordem interna, esta função diverge silenciosamente — daí o
    teste que compara as duas saídas em compostos sem tartarato, onde elas devem
    ser idênticas.
    """
    if exclude_flag(mol, includeRDKitSanitization=False):
        return mol

    mol = _ref.update_mol_valences(mol)
    mol = _ref.remove_sgroups_from_mol(mol)
    mol = _ref.kekulize_mol(mol)
    mol = _ref.remove_hs_from_mol(mol)
    mol = _ref.normalize_mol(mol)
    mol = _ref.uncharge_mol(mol)
    # flatten_tartrate_mol omitido deliberadamente
    mol = _ref.cleanup_drawing_mol(mol)
    Chem.SanitizeMol(mol)
    return mol


#: Enumerador reutilizado entre chamadas: instanciar por molécula desperdiça o
#: mesmo custo de compilação que torna ``get_parent_mol`` o gargalo da referência.
_TAUTOMER_ENUMERATOR = rdMolStandardize.TautomerEnumerator()


def canonicalize_tautomer(mol: Chem.Mol) -> Chem.Mol:
    """Forma tautomérica canônica (ablação da D-04).

    Nota de API: ``rdMolStandardize.TautomerCanonicalizer`` **não existe** — é o nome
    do MolVS, que está obsoleto. A interface atual é
    ``TautomerEnumerator().Canonicalize()``.
    """
    return _TAUTOMER_ENUMERATOR.Canonicalize(mol)


# --- Definição das variantes --------------------------------------------------------


@dataclass(frozen=True)
class PolicyVariant:
    """Uma política a comparar contra a linha de base."""

    name: str
    adr: str
    description: str
    build: Callable[[], CurationPipeline]

    def pipeline(self) -> CurationPipeline:
        return self.build()


def _pipeline(**engine_kwargs) -> CurationPipeline:
    engine = EngineWrapper(
        policy_hash=POLICY_HASH_PLACEHOLDER,
        criteria=EligibilityCriteria(),
        **engine_kwargs,
    )
    return CurationPipeline(POLICY_HASH_PLACEHOLDER, engine=engine)


def default_variants() -> list[PolicyVariant]:
    """As variantes derivadas das decisões em aberto de ``docs/decisions.md``."""
    return [
        PolicyVariant(
            name="baseline",
            adr="D-01..D-10",
            description="políticas vigentes, biblioteca de referência intacta",
            build=lambda: _pipeline(),
        ),
        PolicyVariant(
            name="no_tartrate_flattening",
            adr="D-02",
            description="preserva a estereoquímica de tartaratos",
            build=lambda: _pipeline(
                standardize_fn=standardize_without_tartrate_flattening
            ),
        ),
        PolicyVariant(
            name="canonical_tautomer",
            adr="D-04",
            description="canonicaliza o tautômero da estrutura-mãe",
            build=lambda: _pipeline(parent_transform=canonicalize_tautomer),
        ),
        PolicyVariant(
            name="discard_organometallics",
            adr="D-03",
            description="descarta compostos com exclude_flag em vez de preservá-los",
            build=lambda: _pipeline(discard_excluded=True),
        ),
    ]


# --- Métricas ------------------------------------------------------------------------


@dataclass(frozen=True)
class Distribution:
    """Resumo de uma propriedade contínua."""

    n: int
    mean: float
    median: float
    stdev: float
    minimum: float
    maximum: float

    @classmethod
    def from_values(cls, values: Sequence[float]) -> "Distribution":
        if not values:
            return cls(0, 0.0, 0.0, 0.0, 0.0, 0.0)
        return cls(
            n=len(values),
            mean=mean(values),
            median=median(values),
            stdev=pstdev(values) if len(values) > 1 else 0.0,
            minimum=min(values),
            maximum=max(values),
        )


@dataclass(frozen=True)
class VariantResult:
    """Resultado de uma variante sobre o corpus."""

    variant: PolicyVariant
    n_input: int
    n_passed: int
    n_rejected: int
    n_unique: int
    n_collisions: int
    collision_counts: dict[str, int]
    molecular_weight: Distribution
    logp: Distribution
    tpsa: Distribution
    defined_stereocenters: int
    molecules_with_stereo: int
    rejection_counts: dict[str, int] = field(default_factory=dict)

    @property
    def pass_rate(self) -> float:
        return self.n_passed / self.n_input if self.n_input else 0.0

    @property
    def dedup_rate(self) -> float:
        """Fração dos aprovados que foi absorvida como duplicata."""
        if not self.n_passed:
            return 0.0
        return 1.0 - (self.n_unique / self.n_passed)


def _parent_mol(record: CurationRecord) -> Optional[Chem.Mol]:
    """Reconstrói a estrutura-mãe a partir do SMILES curado.

    Usa ``sanitize=False``: compostos com ``exclude_flag`` nunca passam pela
    sanitização da referência e ``MolFromSmiles`` devolveria ``None`` para eles.
    """
    if not record.curated_smiles:
        return None
    mol = Chem.MolFromSmiles(record.curated_smiles, sanitize=False)
    if mol is None:
        return None
    try:
        mol.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(mol)
    except Exception:
        return None
    return mol


def _descriptors(mol: Chem.Mol) -> Optional[tuple[float, float, float]]:
    try:
        return (
            Descriptors.MolWt(mol),
            Crippen.MolLogP(mol),
            rdMolDescriptors.CalcTPSA(mol),
        )
    except Exception:
        return None


def _count_defined_stereocenters(mol: Chem.Mol) -> int:
    try:
        total = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
        unspecified = rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol)
        return total - unspecified
    except Exception:
        return 0


def run_variant(
    variant: PolicyVariant, corpus: Sequence[tuple[str, str]]
) -> VariantResult:
    """Executa uma variante sobre o corpus e coleta as três famílias de métricas."""
    pipeline = variant.pipeline()
    index = DedupIndex()

    weights: list[float] = []
    logps: list[float] = []
    tpsas: list[float] = []
    stereo_total = 0
    with_stereo = 0
    passed = 0
    rejection_counts: dict[str, int] = {}

    for compound_id, smiles in corpus:
        record = pipeline.process_single(smiles, compound_id)
        if not record.passed:
            code = record.rejection_code.value if record.rejection_code else "UNKNOWN"
            rejection_counts[code] = rejection_counts.get(code, 0) + 1
            continue

        passed += 1
        index.add(record)

        mol = _parent_mol(record)
        if mol is None:
            continue
        values = _descriptors(mol)
        if values is not None:
            weights.append(values[0])
            logps.append(values[1])
            tpsas.append(values[2])

        centers = _count_defined_stereocenters(mol)
        stereo_total += centers
        if centers:
            with_stereo += 1

    return VariantResult(
        variant=variant,
        n_input=len(corpus),
        n_passed=passed,
        n_rejected=len(corpus) - passed,
        n_unique=index.unique_count,
        n_collisions=len(index.collisions),
        collision_counts=index.collision_counts(),
        molecular_weight=Distribution.from_values(weights),
        logp=Distribution.from_values(logps),
        tpsa=Distribution.from_values(tpsas),
        defined_stereocenters=stereo_total,
        molecules_with_stereo=with_stereo,
        rejection_counts=rejection_counts,
    )


def run_ablation(
    corpus: Sequence[tuple[str, str]],
    variants: Optional[Sequence[PolicyVariant]] = None,
) -> list[VariantResult]:
    return [run_variant(v, corpus) for v in (variants or default_variants())]


# --- Relatório -------------------------------------------------------------------------


def _delta(value: float, baseline: float, digits: int = 1) -> str:
    difference = value - baseline
    if abs(difference) < 10**-digits:
        return "—"
    return f"{difference:+.{digits}f}"


def to_markdown(results: Sequence[VariantResult]) -> str:
    """Tabelas consolidadas, com cada variante contrastada à linha de base."""
    if not results:
        return "_sem resultados_"
    base = results[0]

    lines = [
        "# Estudo de ablação de políticas de curadoria",
        "",
        f"Corpus: {base.n_input} estruturas. Cada variante altera **uma** decisão "
        "de `docs/decisions.md`; tudo o mais permanece constante.",
        "",
        "## Variantes",
        "",
        "| variante | ADR | descrição |",
        "| --- | --- | --- |",
    ]
    for result in results:
        lines.append(
            f"| `{result.variant.name}` | {result.variant.adr} | "
            f"{result.variant.description} |"
        )

    lines += [
        "",
        "## (a) Aprovação e deduplicação",
        "",
        "| variante | aprovados | taxa | identidades únicas | Δ vs base | colisões |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for result in results:
        delta = "—" if result is base else f"{result.n_unique - base.n_unique:+d}"
        lines.append(
            f"| `{result.variant.name}` | {result.n_passed} | "
            f"{result.pass_rate:.1%} | {result.n_unique} | {delta} | "
            f"{result.n_collisions} |"
        )

    lines += [
        "",
        "## (b) Distribuição físico-química da estrutura-mãe",
        "",
        "| variante | MW médio | Δ | LogP médio | Δ | TPSA média | Δ |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for result in results:
        is_base = result is base
        lines.append(
            f"| `{result.variant.name}` | {result.molecular_weight.mean:.1f} | "
            f"{'—' if is_base else _delta(result.molecular_weight.mean, base.molecular_weight.mean)} | "
            f"{result.logp.mean:.2f} | "
            f"{'—' if is_base else _delta(result.logp.mean, base.logp.mean, 2)} | "
            f"{result.tpsa.mean:.1f} | "
            f"{'—' if is_base else _delta(result.tpsa.mean, base.tpsa.mean)} |"
        )

    lines += [
        "",
        "## (c) Conservação de centros quirais definidos",
        "",
        "| variante | centros definidos | Δ vs base | moléculas com estereoquímica |",
        "| --- | ---: | ---: | ---: |",
    ]
    for result in results:
        delta = (
            "—"
            if result is base
            else f"{result.defined_stereocenters - base.defined_stereocenters:+d}"
        )
        lines.append(
            f"| `{result.variant.name}` | {result.defined_stereocenters} | "
            f"{delta} | {result.molecules_with_stereo} |"
        )

    interpretation = _interpret(results)
    if interpretation:
        lines += ["", "## Leitura", "", *interpretation]

    return "\n".join(lines)


def _interpret(results: Sequence[VariantResult]) -> list[str]:
    """Aponta as diferenças materiais, sem concluir por elas.

    O harness mede; a interpretação causal é do pesquisador. Estas linhas apenas
    sinalizam onde olhar.
    """
    base = results[0]
    notes: list[str] = []
    for result in results[1:]:
        deltas = []
        if result.defined_stereocenters != base.defined_stereocenters:
            deltas.append(
                f"{result.defined_stereocenters - base.defined_stereocenters:+d} "
                "centros quirais definidos"
            )
        if result.n_unique != base.n_unique:
            deltas.append(f"{result.n_unique - base.n_unique:+d} identidades únicas")
        if result.n_passed != base.n_passed:
            deltas.append(f"{result.n_passed - base.n_passed:+d} compostos aprovados")
        if deltas:
            notes.append(
                f"- `{result.variant.name}` ({result.variant.adr}): "
                + ", ".join(deltas)
                + "."
            )
    if not notes:
        notes.append(
            "- Nenhuma variante produziu diferença mensurável neste corpus. "
            "Corpus pequeno ou pobre nas classes sensíveis às políticas testadas "
            "produz esse resultado — verificar a composição antes de concluir que "
            "as decisões são indiferentes."
        )
    return notes


def write_plots(results: Sequence[VariantResult], out_dir: Union[str, Path]) -> list[Path]:
    """Gera as figuras comparativas. Devolve os caminhos escritos.

    Se ``matplotlib`` não estiver disponível, devolve lista vazia em vez de falhar:
    as tabelas em Markdown são o entregável primário, as figuras são acessórias.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []

    target = Path(out_dir)
    target.mkdir(parents=True, exist_ok=True)
    names = [result.variant.name for result in results]
    written: list[Path] = []

    panels = (
        ("identidades_unicas", "Identidades únicas", [r.n_unique for r in results]),
        (
            "centros_quirais",
            "Centros quirais definidos",
            [r.defined_stereocenters for r in results],
        ),
        ("aprovados", "Compostos aprovados", [r.n_passed for r in results]),
    )
    for filename, title, values in panels:
        figure, axes = plt.subplots(figsize=(7, 4))
        axes.bar(names, values, color="#4C72B0")
        axes.set_title(title)
        axes.set_ylabel(title)
        axes.tick_params(axis="x", rotation=20)
        for index, value in enumerate(values):
            axes.text(index, value, str(value), ha="center", va="bottom", fontsize=9)
        figure.tight_layout()
        path = target / f"ablation_{filename}.png"
        figure.savefig(path, dpi=150)
        plt.close(figure)
        written.append(path)

    figure, axes = plt.subplots(figsize=(7, 4))
    axes.bar(names, [r.molecular_weight.mean for r in results], color="#DD8452")
    axes.set_title("Peso molecular médio da estrutura-mãe")
    axes.set_ylabel("MW (Da)")
    axes.tick_params(axis="x", rotation=20)
    figure.tight_layout()
    path = target / "ablation_peso_molecular.png"
    figure.savefig(path, dpi=150)
    plt.close(figure)
    written.append(path)

    return written


def write_report(
    results: Sequence[VariantResult], out_dir: Union[str, Path]
) -> Path:
    """Escreve ``ablation.md`` e as figuras no diretório indicado."""
    target = Path(out_dir)
    target.mkdir(parents=True, exist_ok=True)
    figures = write_plots(results, target)

    body = to_markdown(results)
    if figures:
        body += "\n\n## Figuras\n\n" + "\n".join(
            f"![{path.stem}]({path.name})" for path in figures
        )

    report = target / "ablation.md"
    report.write_text(body + "\n", encoding="utf-8")
    return report
