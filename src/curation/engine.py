"""Wrapper de auditoria sobre o motor químico do ChEMBL.

Responsabilidade única: chamar ``standardize_mol`` e ``get_parent_mol`` da biblioteca
oficial e observar o que aconteceu. Nenhuma química é implementada aqui — D-01.

Restrição de ordem (D-09): o estágio de parsing **instrumenta e registra, nunca
rejeita por valência**. O normalizador da referência repara classes inteiras que uma
sanitização estrita descartaria — amônios quaternários neutros e sais de diazônio,
entre outras. A rejeição por valência só é legítima depois do motor, e é por isso que
ela emerge de ``standardize_mol``, que sanitiza ao final.
"""

from __future__ import annotations

import re
import time
from collections import defaultdict
from dataclasses import dataclass
from typing import Callable, Optional, Sequence

import rdkit
from rdkit import Chem

import chembl_structure_pipeline as csp

from curation.filters import EligibilityCriteria, EligibilityVerdict, evaluate
from curation.models import (
    CurationRecord,
    RejectionCode,
    Stage,
    TransformationEvent,
)
from curation.probes import (
    DEFAULT_PROBES,
    Probe,
    count_defined_stereocenters,
    count_fragments,
    fragment_smiles,
    molecular_formula,
    net_charge,
    run_probes,
    safe_smiles,
)

PIPELINE_VERSION = "0.1.0"

#: Assinatura de um padronizador substituível.
StandardizeFn = Callable[[Chem.Mol], Chem.Mol]

#: Transformação opcional aplicada à estrutura-mãe, depois de ``get_parent_mol``.
ParentTransform = Callable[[Chem.Mol], Chem.Mol]

_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


@dataclass(frozen=True)
class _Failure:
    """Falha terminal de um estágio."""

    code: RejectionCode
    stage: Stage
    detail: str


def _classify(error: Exception) -> RejectionCode:
    """Mapeia exceções do RDKit para códigos de rejeição.

    ``AtomValenceException`` e ``KekulizeException`` são subclasses de
    ``MolSanitizeException``, então a ordem dos testes importa.
    """
    if isinstance(error, Chem.AtomValenceException):
        return RejectionCode.ERR_VALENCE
    if isinstance(error, Chem.KekulizeException):
        return RejectionCode.ERR_KEKULIZE
    if isinstance(error, Chem.MolSanitizeException):
        return RejectionCode.ERR_SANITIZE
    return RejectionCode.ERR_STANDARDIZE


def _parse_formula(formula: str) -> dict[str, int]:
    """Decompõe uma fórmula molecular em contagens por elemento.

    Sufixos de carga (``+``, ``-``, ``2+``) são ignorados: a variação de carga é
    reportada separadamente em ``delta_net_charge``.
    """
    counts: dict[str, int] = {}
    base = formula.split("+")[0].split("-")[0]
    for element, digits in _FORMULA_TOKEN.findall(base):
        counts[element] = counts.get(element, 0) + (int(digits) if digits else 1)
    return counts


def _formula_delta(before: str, after: str) -> str:
    """Diferença elemento a elemento, com sinal, em ordem estável.

    Ex.: ``C7H5NaO2 -> C7H6O2`` produz ``"+H1 -Na1"``. String vazia quando não há
    mudança, o que torna o campo diretamente filtrável.
    """
    if not before or not after:
        return ""
    counts_before = _parse_formula(before)
    counts_after = _parse_formula(after)
    parts: list[str] = []
    for element in sorted(set(counts_before) | set(counts_after)):
        delta = counts_after.get(element, 0) - counts_before.get(element, 0)
        if delta:
            sign = "+" if delta > 0 else "-"
            parts.append(f"{sign}{element}{abs(delta)}")
    return " ".join(parts)


class EngineWrapper:
    """Executa o motor da referência e produz um ``CurationRecord`` auditado.

    O wrapper é sem estado entre chamadas e seguro para reuso: as moléculas recebidas
    pela biblioteca são sempre cópias defensivas, de modo que nenhum estado do
    chamador é mutado.

    Args:
        policy_hash: SHA-256 de ``docs/decisions.md``. É o que torna dois lotes
            comparáveis; sem ele a proveniência não identifica sob qual política o
            resultado foi produzido.
        probes: sondas a executar em cada fronteira. O padrão cobre as decisões
            atualmente em vigor.
    """

    def __init__(
        self,
        policy_hash: str,
        probes: Sequence[Probe] = DEFAULT_PROBES,
        pipeline_version: str = PIPELINE_VERSION,
        criteria: Optional[EligibilityCriteria] = None,
        standardize_fn: Optional[StandardizeFn] = None,
        parent_transform: Optional[ParentTransform] = None,
        discard_excluded: bool = False,
    ) -> None:
        self._policy_hash = policy_hash
        self._probes = tuple(probes)
        self._criteria = criteria or EligibilityCriteria()
        self._standardize_fn = standardize_fn or csp.standardize_mol
        self._parent_transform = parent_transform
        self._discard_excluded = discard_excluded
        self._pipeline_version = pipeline_version
        self._rdkit_version = rdkit.__version__
        self._csp_version = getattr(csp, "__version__", "unknown")
        self.stage_seconds: dict[str, float] = defaultdict(float)

    def reset_timings(self) -> None:
        """Zera os acumuladores de tempo. Chamado no inicio de cada lote."""
        self.stage_seconds = defaultdict(float)

    def _timed(self, stage: Stage, action):
        """Executa ``action`` medindo o tempo real gasto no estagio.

        Os tempos sao acumulados por estagio ao longo do lote. Nao ha estimativa:
        o que a interface exibe e tempo de parede medido.
        """
        started = time.perf_counter()
        try:
            return action()
        finally:
            self.stage_seconds[stage.value] += time.perf_counter() - started

    # --- API pública -------------------------------------------------------------

    def curate(self, input_id: str, raw_smiles: str) -> CurationRecord:
        """Cura um composto, sempre devolvendo um registro.

        Nunca levanta exceção por conteúdo da entrada: qualquer falha vira um
        registro ``REJECTED`` com estágio e causa. A barreira de exceção do lote
        (Fase 3) é uma segunda linha de defesa, não a primeira.
        """
        events: list[TransformationEvent] = []

        parsed = self._timed(Stage.PARSE, lambda: self._parse(raw_smiles))
        if isinstance(parsed, _Failure):
            return self._rejected(input_id, raw_smiles, parsed, events)
        original = parsed

        standardized = self._timed(
            Stage.STANDARDIZE, lambda: self._standardize(original)
        )
        if isinstance(standardized, _Failure):
            return self._rejected(input_id, raw_smiles, standardized, events)
        events.extend(
            run_probes(self._probes, original, standardized, Stage.STANDARDIZE)
        )

        parent_result = self._timed(
            Stage.GET_PARENT, lambda: self._get_parent(standardized)
        )
        if isinstance(parent_result, _Failure):
            return self._rejected(input_id, raw_smiles, parent_result, events)
        parent, excluded = parent_result
        events.extend(run_probes(self._probes, standardized, parent, Stage.GET_PARENT))

        if excluded and safe_smiles(original) == safe_smiles(standardized):
            events.insert(
                0,
                TransformationEvent(
                    stage=Stage.STANDARDIZE,
                    rule="standardization_skipped_excluded",
                    before_smiles=safe_smiles(original),
                    after_smiles=safe_smiles(standardized),
                    detail="exclude_flag ativo: motor devolveu a estrutura inalterada (D-03)",
                ),
            )

        if excluded and self._discard_excluded:
            return self._rejected(
                input_id,
                raw_smiles,
                _Failure(
                    RejectionCode.ERR_ORGANOMETALLIC,
                    Stage.GET_PARENT,
                    "composto com exclude_flag descartado por politica de ablacao",
                ),
                events,
            )

        if self._parent_transform is not None:
            try:
                transformed = self._parent_transform(Chem.Mol(parent))
                if transformed is not None and transformed.GetNumAtoms():
                    parent = transformed
            except Exception:
                pass

        gate_failure = self._timed(
            Stage.VALENCE_GATE, lambda: self._valence_gate(parent, excluded)
        )
        if gate_failure is not None:
            return self._rejected(input_id, raw_smiles, gate_failure, events)

        verdict = self._timed(
            Stage.ELIGIBILITY,
            lambda: evaluate(parent, self._criteria, count_fragments(parent)),
        )
        if not verdict.eligible:
            return self._rejected(
                input_id,
                raw_smiles,
                _Failure(
                    verdict.rejection_code or RejectionCode.ERR_MW_LIMIT,
                    Stage.ELIGIBILITY,
                    verdict.detail,
                ),
                events,
                verdict=verdict,
            )

        identity = self._timed(
            Stage.CANONICALIZE, lambda: self._canonicalize(parent)
        )
        if isinstance(identity, _Failure):
            return self._rejected(input_id, raw_smiles, identity, events)
        curated_smiles, inchikey = identity

        return self._passed(
            input_id=input_id,
            raw_smiles=raw_smiles,
            original=original,
            standardized=standardized,
            parent=parent,
            curated_smiles=curated_smiles,
            inchikey=inchikey,
            excluded=excluded,
            events=events,
            verdict=verdict,
        )

    # --- Estágios ----------------------------------------------------------------

    def _parse(self, raw_smiles: str) -> Chem.Mol | _Failure:
        """Parsing não-restritivo.

        ``sanitize=False`` deixa a molécula em estado inconsistente — sem percepção de
        anéis nem valências calculadas — daí o ``UpdatePropertyCache(strict=False)``
        obrigatório na sequência. ``strict=False`` é o ponto central: valências
        anômalas são toleradas aqui para que o motor tenha a chance de repará-las.
        """
        text = raw_smiles.strip()
        if not text:
            return _Failure(RejectionCode.ERR_EMPTY, Stage.PARSE, "SMILES vazio")

        try:
            mol = Chem.MolFromSmiles(text, sanitize=False)
        except Exception as error:
            return _Failure(RejectionCode.ERR_SYNTAX, Stage.PARSE, str(error))

        if mol is None:
            return _Failure(
                RejectionCode.ERR_SYNTAX, Stage.PARSE, "SMILES nao interpretavel"
            )
        if mol.GetNumAtoms() == 0:
            return _Failure(
                RejectionCode.ERR_EMPTY, Stage.PARSE, "estrutura sem atomos"
            )

        try:
            mol.UpdatePropertyCache(strict=False)
        except Exception as error:
            return _Failure(RejectionCode.ERR_SYNTAX, Stage.PARSE, str(error))

        return mol

    def _standardize(self, mol: Chem.Mol) -> Chem.Mol | _Failure:
        """Padronização pela referência, sem intervenção.

        É aqui que a rejeição por valência legitimamente aparece: ``standardize_mol``
        sanitiza ao final, depois de ter tentado normalizar. Uma exceção neste ponto
        significa que o motor não conseguiu reparar a estrutura — diferente de
        rejeitá-la antes de tentar.
        """
        try:
            return self._standardize_fn(Chem.Mol(mol))
        except Exception as error:
            return _Failure(_classify(error), Stage.STANDARDIZE, str(error))

    def _get_parent(self, mol: Chem.Mol) -> tuple[Chem.Mol, bool] | _Failure:
        """Isolamento da estrutura-mãe pela referência.

        A flag de exclusão devolvida aqui é a única exposta pela API pública, e é a
        usada em ``excluded_flag`` (D-03).
        """
        try:
            parent, excluded = csp.get_parent_mol(Chem.Mol(mol))
        except Exception as error:
            return _Failure(_classify(error), Stage.GET_PARENT, str(error))

        if parent is None or parent.GetNumAtoms() == 0:
            return _Failure(
                RejectionCode.ERR_EMPTY,
                Stage.GET_PARENT,
                "estrutura-mae vazia apos remocao de fragmentos",
            )
        return parent, bool(excluded)

    def _valence_gate(
        self, parent: Chem.Mol, excluded: bool
    ) -> Optional[_Failure]:
        """Sanitização estrita, aplicada **depois** do motor (D-09).

        Aqui é o único lugar legítimo para reprovar por valência: o normalizador já
        teve a chance de reparar a estrutura. Um portão anterior descartaria amônios
        quaternários neutros e diazônios, que a referência conserta.

        Compostos com ``exclude_flag`` ativo são dispensados do portão, porque o
        próprio motor pula a sanitização deles (D-03) — aplicá-la aqui revogaria a
        decisão de preservá-los. O caso concreto é um carborano de 8 boros presente
        no corpus da referência, anotado lá como ">7 Boron atoms": ele atravessa o
        motor e seria descartado por um portão incondicional.
        """
        if excluded:
            return None
        try:
            Chem.SanitizeMol(Chem.Mol(parent))
        except Exception as error:
            return _Failure(_classify(error), Stage.VALENCE_GATE, str(error))
        return None

    def _canonicalize(self, parent: Chem.Mol) -> tuple[str, str] | _Failure:
        """SMILES canônico isomérico e InChIKey.

        A estereoquímica é preservada como recebida (D-05); o InChIKey completo é a
        chave de identidade (D-07).
        """
        try:
            smiles = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
        except Exception as error:
            return _Failure(
                RejectionCode.ERR_CANONICALIZE, Stage.CANONICALIZE, str(error)
            )
        if not smiles:
            return _Failure(
                RejectionCode.ERR_CANONICALIZE,
                Stage.CANONICALIZE,
                "SMILES canonico vazio",
            )

        try:
            inchikey = Chem.MolToInchiKey(parent)
        except Exception as error:
            return _Failure(RejectionCode.ERR_INCHI, Stage.CANONICALIZE, str(error))
        if not inchikey:
            return _Failure(
                RejectionCode.ERR_INCHI, Stage.CANONICALIZE, "InChIKey vazio"
            )

        return smiles, inchikey

    # --- Construção dos registros ------------------------------------------------

    def _rejected(
        self,
        input_id: str,
        raw_smiles: str,
        failure: _Failure,
        events: list[TransformationEvent],
        verdict: Optional[EligibilityVerdict] = None,
    ) -> CurationRecord:
        return CurationRecord(
            input_id=input_id,
            raw_smiles=raw_smiles,
            status="REJECTED",
            rejection_code=failure.code,
            rejection_stage=failure.stage,
            rejection_detail=failure.detail,
            transformations=events,
            parent_mw=verdict.molecular_weight if verdict else None,
            parent_heavy_atoms=verdict.heavy_atoms if verdict else None,
            **self._provenance(),
        )

    def _passed(
        self,
        input_id: str,
        raw_smiles: str,
        original: Chem.Mol,
        standardized: Chem.Mol,
        parent: Chem.Mol,
        curated_smiles: str,
        inchikey: str,
        excluded: bool,
        events: list[TransformationEvent],
        verdict: EligibilityVerdict,
    ) -> CurationRecord:
        removed = self._removed_fragments(standardized, parent)
        defined_before = count_defined_stereocenters(original)
        defined_after = count_defined_stereocenters(parent)

        return CurationRecord(
            input_id=input_id,
            raw_smiles=raw_smiles,
            curated_smiles=curated_smiles,
            inchikey=inchikey,
            inchikey_block1=inchikey.split("-")[0],
            status="PASSED",
            transformations=events,
            removed_fragments=removed,
            delta_net_charge=net_charge(parent) - net_charge(original),
            delta_formula=_formula_delta(
                molecular_formula(original), molecular_formula(parent)
            ),
            delta_stereocenters=defined_after - defined_before,
            n_undefined_stereocenters=self._undefined_stereocenters(parent),
            n_stereocenters_total=self._total_stereocenters(parent),
            n_components_parent=count_fragments(parent),
            parent_mw=verdict.molecular_weight,
            parent_heavy_atoms=verdict.heavy_atoms,
            excluded_flag=excluded,
            **self._provenance(),
        )

    def _provenance(self) -> dict[str, str]:
        return {
            "rdkit_version": self._rdkit_version,
            "csp_version": self._csp_version,
            "pipeline_version": self._pipeline_version,
            "policy_hash": self._policy_hash,
        }

    @staticmethod
    def _removed_fragments(
        standardized: Chem.Mol, parent: Chem.Mol
    ) -> Optional[str]:
        """Componentes descartados por ``get_parent_mol``.

        A linha de base é a molécula **padronizada**, não a entrada crua. A remoção de
        fragmentos ocorre exclusivamente em ``get_parent_mol``; comparar contra a
        entrada crua confundiria transformação com remoção, porque a normalização
        reescreve os fragmentos retidos — ``CC(=O)O[Na]`` vira ``CC(=O)[O-].[Na+]``
        sem que nada tenha sido removido.

        Comparação por multiconjunto de SMILES: a referência deduplica componentes
        idênticos, então um dímero de ácido acético colapsa em um componente e a cópia
        excedente aparece aqui corretamente como removida.
        """
        before = fragment_smiles(standardized)
        after = fragment_smiles(parent)
        remaining = list(after)
        removed: list[str] = []
        for smiles in before:
            if smiles in remaining:
                remaining.remove(smiles)
            else:
                removed.append(smiles)
        return ".".join(removed) if removed else None

    @staticmethod
    def _undefined_stereocenters(mol: Chem.Mol) -> int:
        from rdkit.Chem import rdMolDescriptors

        try:
            return rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol)
        except Exception:
            return 0

    @staticmethod
    def _total_stereocenters(mol: Chem.Mol) -> int:
        from rdkit.Chem import rdMolDescriptors

        try:
            return rdMolDescriptors.CalcNumAtomStereoCenters(mol)
        except Exception:
            return 0
