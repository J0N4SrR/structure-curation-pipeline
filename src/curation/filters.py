"""Critérios de elegibilidade físico-química (D-10).

Nota terminológica: estes são **critérios de elegibilidade / escopo químico**, não
"domínio de aplicabilidade". Domínio de aplicabilidade é conceito dependente de
modelo — cobertura do espaço de descritores, *leverage*, distância ao modelo — e
pertence à etapa de modelagem, não à curadoria.

Os cortes incidem sobre a **estrutura-mãe isolada**, nunca sobre a forma salificada:
medir a massa do contra-íon rejeitaria compostos válidos. Os padrões são ponto de
partida documentado, não constante universal — produtos naturais, macrociclos e
peptídeos exigem revisão explícita.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from rdkit import Chem
from rdkit.Chem import Descriptors

from curation.models import RejectionCode


@dataclass(frozen=True)
class EligibilityCriteria:
    """Limites de escopo químico do estudo.

    Attributes:
        max_molecular_weight: peso molecular médio máximo, em Da.
        max_heavy_atoms: número máximo de átomos pesados (não-hidrogênio).
    """

    max_molecular_weight: float = 1000.0
    max_heavy_atoms: int = 100
    require_carbon: bool = False

    def describe(self) -> str:
        text = (
            f"MW <= {self.max_molecular_weight:g} Da, "
            f"átomos pesados <= {self.max_heavy_atoms}"
        )
        if self.require_carbon:
            text += ", apenas estruturas contendo carbono"
        return text


@dataclass(frozen=True)
class EligibilityVerdict:
    """Resultado da checagem, com as propriedades medidas sempre presentes.

    As propriedades são reportadas mesmo quando o composto é reprovado: saber que
    uma rejeição ocorreu a 1002 Da é diferente de saber que ocorreu a 4000 Da, e a
    diferença orienta a revisão dos cortes.
    """

    eligible: bool
    molecular_weight: float
    heavy_atoms: int
    carbon_atoms: int = 0
    rejection_code: Optional[RejectionCode] = None
    detail: str = ""


def molecular_properties(mol: Chem.Mol) -> tuple[float, int]:
    """Peso molecular médio e contagem de átomos pesados.

    Tolera moléculas não sanitizadas: compostos com ``exclude_flag`` ativo (D-03)
    nunca passam pela sanitização da referência, e ainda assim precisam ser
    filtrados. Reparsear o SMILES curado não é alternativa — para essas estruturas
    ``MolFromSmiles`` devolve ``None``.
    """
    try:
        weight = Descriptors.MolWt(mol)
    except Exception:
        copy = Chem.Mol(mol)
        copy.UpdatePropertyCache(strict=False)
        weight = Descriptors.MolWt(copy)
    return float(weight), int(mol.GetNumHeavyAtoms())


def count_carbon(mol: Chem.Mol) -> int:
    """Átomos de carbono na estrutura. Base do critério opcional da D-12."""
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)


def evaluate(
    mol: Chem.Mol,
    criteria: EligibilityCriteria,
    n_components: int = 1,
) -> EligibilityVerdict:
    """Avalia a estrutura-mãe contra os critérios.

    Args:
        mol: a estrutura-mãe, depois de ``get_parent_mol``.
        criteria: limites vigentes.
        n_components: número de componentes retidos na estrutura-mãe.

    O ``n_components`` entra no detalhe da rejeição por causa da interação entre
    D-10 e D-06: quando a estrutura-mãe é multicomponente, o corte incide sobre a
    **soma** dos componentes retidos. Um medicamento combinado de dois fármacos
    pesados pode ser reprovado como "molécula grande" sendo, na verdade, duas
    moléculas médias. Sem esse registro o falso positivo é indistinguível na revisão.
    """
    weight, heavy_atoms = molecular_properties(mol)
    carbon_atoms = count_carbon(mol)
    context = f" (componentes na estrutura-mãe: {n_components})" if n_components > 1 else ""

    if criteria.require_carbon and not carbon_atoms:
        return EligibilityVerdict(
            eligible=False,
            molecular_weight=weight,
            heavy_atoms=heavy_atoms,
            carbon_atoms=carbon_atoms,
            rejection_code=RejectionCode.ERR_NO_CARBON,
            detail=(
                "estrutura sem átomos de carbono: fora do escopo de moléculas "
                f"pequenas orgânicas{context}"
            ),
        )

    if weight > criteria.max_molecular_weight:
        return EligibilityVerdict(
            eligible=False,
            molecular_weight=weight,
            heavy_atoms=heavy_atoms,
            carbon_atoms=carbon_atoms,
            rejection_code=RejectionCode.ERR_MW_LIMIT,
            detail=(
                f"peso molecular {weight:.2f} Da excede o limite de "
                f"{criteria.max_molecular_weight:g} Da{context}"
            ),
        )

    if heavy_atoms > criteria.max_heavy_atoms:
        return EligibilityVerdict(
            eligible=False,
            molecular_weight=weight,
            heavy_atoms=heavy_atoms,
            carbon_atoms=carbon_atoms,
            rejection_code=RejectionCode.ERR_HA_LIMIT,
            detail=(
                f"{heavy_atoms} átomos pesados excedem o limite de "
                f"{criteria.max_heavy_atoms}{context}"
            ),
        )

    return EligibilityVerdict(
        eligible=True,
        molecular_weight=weight,
        heavy_atoms=heavy_atoms,
        carbon_atoms=carbon_atoms,
    )
