"""Proveniência de execução: o que é preciso para auditar e reproduzir um lote.

Vive na camada de engenharia, não na interface. A UI apresenta estes dados; ela não
os deriva, não os completa e não os estima.

Todo campo aqui é **observado**. Quando uma informação não pode ser obtida - o
repositório não é um checkout git, o arquivo de política não existe - o campo diz
isso explicitamente em vez de receber um valor plausível. Proveniência inventada é
pior que proveniência ausente: a segunda avisa, a primeira engana.
"""

from __future__ import annotations

import hashlib
import platform
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional, Union

import rdkit

import chembl_structure_pipeline as csp

from curation.engine import PIPELINE_VERSION

UNKNOWN = "unavailable"


def _run_git(*args: str) -> Optional[str]:
    try:
        result = subprocess.run(
            ["git", *args],
            capture_output=True,
            text=True,
            timeout=5,
            cwd=Path(__file__).resolve().parent,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip() or None


@dataclass(frozen=True)
class GitState:
    """Estado do repositório no momento da execução.

    ``dirty`` importa tanto quanto o commit: um lote produzido sobre árvore suja não
    é reproduzível a partir do commit sozinho, e isso precisa aparecer no manifesto.
    """

    commit: str = UNKNOWN
    branch: str = UNKNOWN
    dirty: Optional[bool] = None

    @classmethod
    def detect(cls) -> "GitState":
        commit = _run_git("rev-parse", "HEAD")
        if commit is None:
            return cls()
        status = _run_git("status", "--porcelain")
        return cls(
            commit=commit,
            branch=_run_git("rev-parse", "--abbrev-ref", "HEAD") or UNKNOWN,
            dirty=bool(status),
        )

    @property
    def reproducible(self) -> bool:
        return self.commit != UNKNOWN and self.dirty is False


def sha256_of(data: Union[bytes, str, Path]) -> str:
    """SHA-256 de bytes, texto ou arquivo."""
    if isinstance(data, Path):
        return hashlib.sha256(data.read_bytes()).hexdigest()
    if isinstance(data, str):
        data = data.encode("utf-8")
    return hashlib.sha256(data).hexdigest()


def tool_versions() -> dict[str, str]:
    """Versões que afetam o resultado químico.

    RDKit encabeça a lista por um motivo concreto: mudanças de versão alteram
    percepção de aromaticidade e regras de padronização, então dois lotes com RDKit
    diferente não são comparáveis mesmo com o mesmo ``policy_hash``.
    """
    versions = {
        "python": sys.version.split()[0],
        "rdkit": rdkit.__version__,
        "chembl_structure_pipeline": getattr(csp, "__version__", UNKNOWN),
    }
    for name in ("pydantic", "pandas", "numpy", "streamlit"):
        try:
            module = __import__(name)
            versions[name] = getattr(module, "__version__", UNKNOWN)
        except ImportError:
            continue
    return versions


def environment_info() -> dict[str, str]:
    return {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python_implementation": platform.python_implementation(),
    }


@dataclass(frozen=True)
class RunProvenance:
    """Identidade completa de uma execução.

    Attributes:
        input_hash: SHA-256 do conteúdo bruto de entrada. É o que permite afirmar
            que duas execuções viram os mesmos dados.
        output_hash: SHA-256 do CSV completo produzido. Preenchido depois da
            execução; ``None`` enquanto o lote não terminou.
    """

    run_id: str
    started_at: str
    finished_at: Optional[str] = None
    duration_seconds: Optional[float] = None
    pipeline_version: str = PIPELINE_VERSION
    policy_hash: str = UNKNOWN
    policy_path: str = UNKNOWN
    parameters: dict[str, Any] = field(default_factory=dict)
    git: GitState = field(default_factory=GitState)
    versions: dict[str, str] = field(default_factory=tool_versions)
    environment: dict[str, str] = field(default_factory=environment_info)
    input_name: str = UNKNOWN
    input_hash: str = UNKNOWN
    output_hash: Optional[str] = None

    @classmethod
    def start(
        cls,
        policy_hash: str,
        parameters: dict[str, Any],
        input_bytes: bytes,
        input_name: str = "input",
        policy_path: str = UNKNOWN,
    ) -> "RunProvenance":
        now = datetime.now(timezone.utc)
        digest = sha256_of(input_bytes)
        return cls(
            run_id=f"{now.strftime('%Y%m%dT%H%M%SZ')}-{digest[:8]}",
            started_at=now.isoformat(),
            policy_hash=policy_hash,
            policy_path=policy_path,
            parameters=dict(parameters),
            git=GitState.detect(),
            input_name=input_name,
            input_hash=digest,
        )

    def finish(
        self, duration_seconds: float, output_hash: Optional[str] = None
    ) -> "RunProvenance":
        return RunProvenance(
            run_id=self.run_id,
            started_at=self.started_at,
            finished_at=datetime.now(timezone.utc).isoformat(),
            duration_seconds=duration_seconds,
            pipeline_version=self.pipeline_version,
            policy_hash=self.policy_hash,
            policy_path=self.policy_path,
            parameters=self.parameters,
            git=self.git,
            versions=self.versions,
            environment=self.environment,
            input_name=self.input_name,
            input_hash=self.input_hash,
            output_hash=output_hash,
        )

    # --- Reprodução ---------------------------------------------------------------

    @property
    def reproducible(self) -> bool:
        """Se a execução pode ser reproduzida exatamente a partir do registrado."""
        return (
            self.git.reproducible
            and self.policy_hash not in (UNKNOWN, "UNVERSIONED_POLICY")
        )

    def reproduction_blockers(self) -> list[str]:
        """O que impede a reprodução exata. Lista vazia significa reprodutível."""
        blockers: list[str] = []
        if self.git.commit == UNKNOWN:
            blockers.append("commit git não identificado")
        elif self.git.dirty:
            blockers.append(
                "árvore de trabalho com alterações não commitadas: o commit "
                "registrado não descreve o código executado"
            )
        if self.policy_hash == "UNVERSIONED_POLICY":
            blockers.append("execução sem arquivo de política versionado")
        elif self.policy_hash == UNKNOWN:
            blockers.append("policy_hash não registrado")
        return blockers

    def reproduction_command(self, input_path: str = "<entrada>") -> str:
        """Comando equivalente para reexecutar o lote pela CLI."""
        parts = [
            "curation",
            f"--input {input_path}",
            f"--out-dir runs/{self.run_id}",
        ]
        if "max_mw" in self.parameters:
            parts.append(f"--max-mw {self.parameters['max_mw']}")
        if "max_ha" in self.parameters:
            parts.append(f"--max-ha {self.parameters['max_ha']}")
        if self.parameters.get("deduplicate") is False:
            parts.append("--no-dedup")
        if self.parameters.get("require_carbon"):
            parts.append("--require-carbon")
        if self.policy_path not in (UNKNOWN, ""):
            parts.append(f"--decisions {self.policy_path}")

        command = " \\\n    ".join(parts)
        if self.git.commit != UNKNOWN:
            command = f"git checkout {self.git.commit}\n{command}"
        return command

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["git"] = asdict(self.git)
        payload["reproducible"] = self.reproducible
        payload["reproduction_blockers"] = self.reproduction_blockers()
        return payload
