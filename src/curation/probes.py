"""Sondas post-hoc: instrumentação do motor sem violar seu isolamento.

``standardize_mol`` é uma caixa-preta de oito passos internos e a D-01 proíbe fork ou
monkeypatch. Atribuição por regra interna é, portanto, inobservável de fora. A
estratégia adotada é diff de fronteira mais sondas dirigidas: cada sonda observa os
estados antes e depois de um estágio e reporta um efeito *compatível* com uma regra
conhecida, sem alegar acesso ao interior do motor.

Cada sonda existe para tornar mensurável uma decisão de ``docs/decisions.md``. Sondas
não são adicionadas por conveniência: o critério é epistêmico, não de desempenho -
medições mostram que o custo das sondas é de 2 a 4% do tempo total, enquanto o motor
responde por ~85%.

Nota de implementação: todo SMARTS é compilado **uma vez** no nível do módulo. A
referência recompila 172 SMARTS por molécula dentro do corpo da função, o que responde
por ~73% do tempo do pipeline; esse erro não é reproduzido aqui.
"""

from __future__ import annotations

from typing import Optional, Protocol, Sequence, runtime_checkable

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from curation.models import Stage, TransformationEvent

# --- Padrões compilados uma única vez -------------------------------------------

#: Espelha ``_alkoxide_pattern`` de ``chembl_structure_pipeline.standardizer``.
#: Deliberadamente idêntico ao da referência: a sonda deve detectar exatamente o que
#: o motor corrige, nem mais nem menos.
_ALKALI_COVALENT = Chem.MolFromSmarts("[Li,Na,K;+0]-[#7,#8;+0]")


def _build_tartrate_query() -> Chem.Mol:
    """Reproduz a consulta de ``flatten_tartrate_mol``.

    A referência restringe o casamento a fragmentos livres de tartarato via
    ``adjustDegree``; sem esse ajuste a sonda marcaria ésteres e amidas de tartarato
    que o motor não achata.
    """
    query = Chem.MolFromSmarts("OC(=O)C(O)C(O)C(=O)O")
    params = Chem.AdjustQueryParameters.NoAdjustments()
    params.adjustDegree = True
    params.adjustDegreeFlags = Chem.AdjustQueryWhichFlags.ADJUST_IGNORENONE
    return Chem.AdjustQueryProperties(query, params)


_TARTRATE = _build_tartrate_query()


# --- Utilitários tolerantes a moléculas não sanitizadas --------------------------


def safe_smiles(mol: Optional[Chem.Mol]) -> str:
    """SMILES de uma molécula possivelmente não sanitizada.

    Estados intermediários do pipeline não são necessariamente sanitizáveis. Uma
    falha de escrita nunca pode derrubar a auditoria, então degrada para string vazia.
    """
    if mol is None:
        return ""
    try:
        return Chem.MolToSmiles(mol)
    except Exception:
        try:
            copy = Chem.Mol(mol)
            copy.UpdatePropertyCache(strict=False)
            Chem.FastFindRings(copy)
            return Chem.MolToSmiles(copy)
        except Exception:
            return ""


def count_matches(mol: Optional[Chem.Mol], query: Chem.Mol) -> int:
    if mol is None:
        return 0
    try:
        return len(mol.GetSubstructMatches(query))
    except Exception:
        return 0


def net_charge(mol: Optional[Chem.Mol]) -> int:
    if mol is None:
        return 0
    try:
        return Chem.GetFormalCharge(mol)
    except Exception:
        return sum(atom.GetFormalCharge() for atom in mol.GetAtoms())


def fragment_smiles(mol: Optional[Chem.Mol]) -> list[str]:
    """SMILES de cada componente desconectado, em ordem canônica estável."""
    if mol is None:
        return []
    try:
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    except Exception:
        return []
    return sorted(safe_smiles(frag) for frag in frags)


def count_fragments(mol: Optional[Chem.Mol]) -> int:
    if mol is None:
        return 0
    try:
        return len(Chem.GetMolFrags(mol))
    except Exception:
        return 0


def count_defined_stereocenters(mol: Optional[Chem.Mol]) -> int:
    """Centros estereoquímicos com configuração atribuída."""
    if mol is None:
        return 0
    try:
        total = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
        unspecified = rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol)
        return total - unspecified
    except Exception:
        return sum(
            1
            for atom in mol.GetAtoms()
            if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
        )


def molecular_formula(mol: Optional[Chem.Mol]) -> str:
    if mol is None:
        return ""
    try:
        return rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        return ""


# --- Protocolo das sondas --------------------------------------------------------


@runtime_checkable
class Probe(Protocol):
    """Observa uma fronteira de estágio e reporta os efeitos detectados.

    Implementações devem ser puras e não podem modificar as moléculas recebidas.
    """

    rule: str

    def inspect(
        self, before: Optional[Chem.Mol], after: Optional[Chem.Mol], stage: Stage
    ) -> list[TransformationEvent]:
        ...


class CovalentAlkaliMetalProbe:
    """Ionização de metal alcalino covalentemente ligado - D-01, canário da ordem.

    O motor quebra ligações Li/Na/K–N/O dentro de ``normalize_mol``, *antes* do
    stripping de fragmentos. Se a ordem for invertida, ``GetMolFrags`` enxerga um
    único fragmento e o metal atravessa o pipeline ligado covalentemente, sem erro
    de execução.

    Além de registrar a ionização, esta sonda emite ``alkali_metal_residual`` quando
    resta ligação covalente no estado final - é o alarme direto para essa regressão.
    """

    rule = "alkali_metal_ionized"
    residual_rule = "alkali_metal_residual"

    def inspect(
        self, before: Optional[Chem.Mol], after: Optional[Chem.Mol], stage: Stage
    ) -> list[TransformationEvent]:
        n_before = count_matches(before, _ALKALI_COVALENT)
        n_after = count_matches(after, _ALKALI_COVALENT)
        events: list[TransformationEvent] = []

        if n_before > n_after:
            events.append(
                TransformationEvent(
                    stage=stage,
                    rule=self.rule,
                    before_smiles=safe_smiles(before),
                    after_smiles=safe_smiles(after),
                    detail=f"ligacoes covalentes de metal alcalino: {n_before} -> {n_after}",
                )
            )

        if stage is Stage.GET_PARENT and n_after > 0:
            events.append(
                TransformationEvent(
                    stage=stage,
                    rule=self.residual_rule,
                    before_smiles=safe_smiles(before),
                    after_smiles=safe_smiles(after),
                    detail=(
                        f"{n_after} ligacao(oes) covalente(s) de metal alcalino "
                        "remanescente(s) na estrutura-mae"
                    ),
                )
            )

        return events


class TartrateFlattenProbe:
    """Perda de quiralidade em tartaratos - D-02.

    A D-02 aceita ``flatten_tartrate_mol`` de forma incondicional, porque a API
    pública não permite desativá-la. A contrapartida é que toda ocorrência seja
    registrada, para que a perda vire um número no relatório de lote em vez de um
    efeito invisível.

    Exige as duas condições - casamento do padrão de tartarato livre *e* redução de
    centros definidos - para não atribuir ao achatamento uma perda de estereoquímica
    de outra origem.
    """

    rule = "tartrate_flattened"

    def inspect(
        self, before: Optional[Chem.Mol], after: Optional[Chem.Mol], stage: Stage
    ) -> list[TransformationEvent]:
        if not count_matches(before, _TARTRATE):
            return []

        defined_before = count_defined_stereocenters(before)
        defined_after = count_defined_stereocenters(after)
        if defined_after >= defined_before:
            return []

        return [
            TransformationEvent(
                stage=stage,
                rule=self.rule,
                before_smiles=safe_smiles(before),
                after_smiles=safe_smiles(after),
                detail=(
                    "centros estereoquimicos definidos: "
                    f"{defined_before} -> {defined_after}"
                ),
            )
        ]


class NetChargeProbe:
    """Variação de carga formal líquida - D-01 e D-06.

    O ``Uncharger`` da referência é ciente de balanço de carga: não neutraliza um
    ânion enquanto houver contra-íon presente. Por isso o motor desprotona duas
    vezes, uma antes e outra depois do stripping. Esta sonda torna as duas
    observáveis separadamente.
    """

    rule = "net_charge_changed"

    def inspect(
        self, before: Optional[Chem.Mol], after: Optional[Chem.Mol], stage: Stage
    ) -> list[TransformationEvent]:
        charge_before = net_charge(before)
        charge_after = net_charge(after)
        if charge_before == charge_after:
            return []

        return [
            TransformationEvent(
                stage=stage,
                rule=self.rule,
                before_smiles=safe_smiles(before),
                after_smiles=safe_smiles(after),
                detail=f"carga liquida: {charge_before:+d} -> {charge_after:+d}",
            )
        ]


class FragmentCountProbe:
    """Variação no número de componentes desconectados - D-06.

    Reporta os dois sentidos. Redução é remoção de sal ou solvente, com os SMILES
    removidos no detalhe. **Aumento** ocorre quando a normalização ioniza uma ligação
    covalente e é informação de auditoria legítima, não anomalia.

    A sonda nunca escolhe fragmento: a heurística de maior fragmento é expressamente
    proibida pela D-06, e a preservação de misturas legítimas é responsabilidade do
    motor.
    """

    rule = "fragments_removed"
    split_rule = "fragments_split"

    def inspect(
        self, before: Optional[Chem.Mol], after: Optional[Chem.Mol], stage: Stage
    ) -> list[TransformationEvent]:
        n_before = count_fragments(before)
        n_after = count_fragments(after)
        if n_before == n_after:
            return []

        smiles_before = fragment_smiles(before)
        smiles_after = fragment_smiles(after)

        if n_after < n_before:
            remaining = list(smiles_after)
            removed: list[str] = []
            for smiles in smiles_before:
                if smiles in remaining:
                    remaining.remove(smiles)
                else:
                    removed.append(smiles)
            detail = f"componentes: {n_before} -> {n_after}"
            if removed:
                detail += f"; removidos: {'.'.join(removed)}"
            return [
                TransformationEvent(
                    stage=stage,
                    rule=self.rule,
                    before_smiles=safe_smiles(before),
                    after_smiles=safe_smiles(after),
                    detail=detail,
                )
            ]

        return [
            TransformationEvent(
                stage=stage,
                rule=self.split_rule,
                before_smiles=safe_smiles(before),
                after_smiles=safe_smiles(after),
                detail=f"componentes: {n_before} -> {n_after} (ionizacao de ligacao covalente)",
            )
        ]


#: Conjunto padrão. Ordem fixa para que a lista de eventos seja determinística.
DEFAULT_PROBES: tuple[Probe, ...] = (
    CovalentAlkaliMetalProbe(),
    TartrateFlattenProbe(),
    NetChargeProbe(),
    FragmentCountProbe(),
)


def run_probes(
    probes: Sequence[Probe],
    before: Optional[Chem.Mol],
    after: Optional[Chem.Mol],
    stage: Stage,
) -> list[TransformationEvent]:
    """Executa as sondas numa fronteira, isolando falhas individuais.

    Uma sonda defeituosa degrada a auditoria daquele registro; jamais derruba a
    curadoria dele.
    """
    events: list[TransformationEvent] = []
    for probe in probes:
        try:
            events.extend(probe.inspect(before, after, stage))
        except Exception:
            continue
    return events
