# Auditoria de Deduplicação — Identidade Estrutural

**Data da Auditoria:** 2026-09-07  
**Agente Responsável:** AGENT 8 — DEDUPLICATION AUDITOR  
**Status de Deduplicação:** COMPLIANT  

---

## 1. Classificação de Duplicatas

O módulo [src/curation/dedup.py](file:///home/jrr/Documentos/Doutorado/structure-curation-pipeline/src/curation/dedup.py) classifica duplicatas em três níveis:
- **`EXACT_DUPLICATE`:** Mesmo InChIKey completo e mesmas anotações.
- **`BLOCK1_COLLISION`:** Mesmo Bloco 1 do InChIKey (mesmo esqueleto de conectividade), mas InChIKey completo diferente (diferente estereoquímica ou grau de ionização).
- **`ANNOTATION_CONFLICT`:** Mesmo composto estrutural, mas com metadados ou anotações conflitantes.

Todos os conflitos são registrados e expostos em relatórios sem descarte silencioso de dados.
