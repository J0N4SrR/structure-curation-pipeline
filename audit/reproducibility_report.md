# Relatório de Reprodutibilidade — Verificação Byte a Byte

**Data do Relatório:** 2026-09-07  
**Agente Responsável:** AGENT 6 — REPRODUCIBILITY AUDITOR  
**Status de Reprodutibilidade:** REPRODUCIBLE  

---

## 1. Verificação de Reprodutibilidade de Execução

Execuções repetidas sobre o mesmo dataset geram resultados **idênticos byte a byte** em todos os artefatos de saída:

- `curated.csv` -> Hash de saída determinístico.
- `rejected.csv` -> Razões e códigos de rejeição estáveis.
- `manifest.json` -> Registra versões exatas das bibliotecas e o `policy_hash`.

---

## 2. Metadados de Proveniência Registrados

Cada lote de saída carimba obrigatoriamente:
- **`policy_hash`**: SHA-256 do documento [docs/decisions.md](file:///home/jrr/Documentos/Doutorado/structure-curation-pipeline/docs/decisions.md).
- **`pipeline_version`**: `0.1.0`.
- **`rdkit_version`**: `2026.3.6`.
- **`csp_version`**: `1.2.4`.
