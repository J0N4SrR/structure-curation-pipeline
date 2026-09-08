"""Interface de linha de comando.

A CLI é a via de execução em lote de alto rendimento; o Streamlit, quando existir,
será um cliente fino sobre esta mesma camada. Nenhuma lógica química vive aqui.

``--out-dir`` é obrigatório por decisão de projeto: o pipeline nunca escreve no
diretório corrente por padrão, para que saídas de execuções distintas não se
misturem e o ``manifest.json`` identifique sempre um lote inteiro.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

from curation.engine import PIPELINE_VERSION
from curation.filters import EligibilityCriteria
from curation.io import compute_policy_hash
from curation.pipeline import BatchSummary, CurationPipeline

DEFAULT_DECISIONS = Path("docs/decisions.md")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="curation",
        description=(
            "Curadoria, validação e padronização estrutural de moléculas pequenas, "
            "com proveniência auditável."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "A saída só é válida quando manifest.json existe: ele é escrito por "
            "último e marca a execução como completa."
        ),
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="arquivo de entrada (.csv, .tsv, .smi) ou '-' para ler da entrada padrão",
    )
    parser.add_argument(
        "--out-dir",
        "-o",
        required=True,
        type=Path,
        help="diretório de saída (criado se não existir)",
    )
    parser.add_argument(
        "--max-mw",
        type=float,
        default=EligibilityCriteria.max_molecular_weight,
        help="peso molecular máximo da estrutura-mãe, em Da",
    )
    parser.add_argument(
        "--max-ha",
        type=int,
        default=EligibilityCriteria.max_heavy_atoms,
        help="número máximo de átomos pesados na estrutura-mãe",
    )
    parser.add_argument(
        "--require-carbon",
        action="store_true",
        help=(
            "rejeita estruturas sem carbono (sais inorgânicos, íons metálicos "
            "isolados). Desligado por padrão — ver D-12"
        ),
    )
    parser.add_argument(
        "--decisions",
        type=Path,
        default=DEFAULT_DECISIONS,
        help="arquivo de ADRs cujo SHA-256 identifica a política do lote",
    )
    parser.add_argument(
        "--no-dedup",
        action="store_true",
        help="desliga a deduplicação por InChIKey e o relatório de conflitos",
    )
    parser.add_argument(
        "--rdkit-logs",
        action="store_true",
        help=(
            "mostra os avisos do RDKit no stderr; por padrão são suprimidos, pois "
            "a causa de cada rejeição já é registrada em rejected.csv"
        ),
    )
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="suprime o sumário no terminal",
    )
    parser.add_argument(
        "--version", action="version", version=f"curation {PIPELINE_VERSION}"
    )
    return parser


def _resolve_policy_hash(path: Path) -> str:
    """SHA-256 do arquivo de decisões, ou um marcador explícito de ausência.

    Executar sem o documento de política é possível, mas o lote fica marcado como
    tal — jamais com um hash falso que sugira reprodutibilidade inexistente.
    """
    try:
        return compute_policy_hash(path)
    except OSError:
        print(
            f"aviso: {path} nao encontrado; o lote sera marcado como "
            "'UNVERSIONED_POLICY' e nao sera comparavel a lotes versionados",
            file=sys.stderr,
        )
        return "UNVERSIONED_POLICY"


def run(argv: Optional[Sequence[str]] = None) -> int:
    """Ponto de entrada. Devolve o código de saída do processo."""
    args = build_parser().parse_args(argv)

    if not args.rdkit_logs:
        # Silenciar aqui, e não na biblioteca: um CLI pode calar o logger global,
        # uma biblioteca importada por terceiros não deve.
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.*")

    criteria = EligibilityCriteria(
        max_molecular_weight=args.max_mw,
        max_heavy_atoms=args.max_ha,
        require_carbon=args.require_carbon,
    )
    pipeline = CurationPipeline(
        policy_hash=_resolve_policy_hash(args.decisions),
        criteria=criteria,
        deduplicate=not args.no_dedup,
    )

    source = sys.stdin if args.input == "-" else Path(args.input)
    if isinstance(source, Path) and not source.is_file():
        print(f"erro: entrada nao encontrada: {source}", file=sys.stderr)
        return 2

    try:
        summary = pipeline.run(source, args.out_dir)
    except KeyboardInterrupt:
        print(
            "\ninterrompido: nenhum arquivo final foi gravado e nao ha manifesto",
            file=sys.stderr,
        )
        return 130
    except OSError as error:
        print(f"erro de E/S: {error}", file=sys.stderr)
        return 1

    if not args.quiet:
        print(f"criterios: {criteria.describe()}")
        print(summary.format_report())

    return 0 if _batch_is_usable(summary) else 1


def _batch_is_usable(summary: BatchSummary) -> bool:
    """Um lote sem nenhuma entrada legível é falha de execução, não resultado."""
    return summary.total > 0


def main() -> None:
    sys.exit(run())


if __name__ == "__main__":
    main()
