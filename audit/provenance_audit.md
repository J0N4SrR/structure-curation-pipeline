# Auditoria de Proveniência de Dados — Rastreabilidade Total

**Data da Auditoria:** 2026-09-07  
**Agente Responsável:** AGENT 7 — DATA PROVENANCE AUDITOR  
**Status de Rastreabilidade:** EXPLICÁVEL E AUDITÁVEL  

---

## 1. Rastreabilidade por Registro

Para qualquer composto processado pelo pipeline, o objeto `CurationRecord` permite responder:
- **Origem:** `input_id` e `raw_smiles`.
- **Efeito:** Lista de objetos `TransformationEvent` com `stage`, `rule`, `before_smiles`, `after_smiles` e `detail`.
- **Mudanças Diferenciais:** `delta_net_charge`, `delta_formula`, `delta_stereocenters`, `removed_fragments`.
- **Metadados:** `policy_hash`, `rdkit_version`, `csp_version`, `pipeline_version`.

Nenhuma transformação ocorre de forma silenciosa ou inobservável.
