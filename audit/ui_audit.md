# Auditoria da Interface de Usuário — Transparência

**Data da Auditoria:** 2026-09-07  
**Agente Responsável:** AGENT 13 — UI / TRANSPARENCY AUDITOR  
**Aplicação Inspecionada:** [src/curation/app.py](file:///home/jrr/Documentos/Doutorado/structure-curation-pipeline/src/curation/app.py) (Streamlit)  

---

## 1. Princípios de Transparência Avaliados

- **Exposição de Transformações:** A interface permite visualizar a molécula original vs. curada, a lista de transformações e o motivo exato da decisão.
- **Rastreabilidade de Política:** Exibe o `policy_hash` e as versões do RDKit, ChEMBL Pipeline e Curation Pipeline em cada execução.
- **Exportação de Artefatos:** Permite ao usuário baixar os relatórios completos, o manifesto e o log de proveniência em CSV/JSON.
