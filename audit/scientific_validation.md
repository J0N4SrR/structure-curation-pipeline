# Validação Científica — Métricas Estatísticas e Gold Standard

**Data da Auditoria:** 2026-09-07  
**Agente Responsável:** AGENT 11 — SCIENTIFIC VALIDATION AGENT  
**Status de Validação:** METODOLOGICAMENTE DEFENSÁVEL  

---

## 1. Métricas Estatísticas Suportadas

O módulo [src/curation/validation/benchmarking.py](file:///home/jrr/Documentos/Doutorado/structure-curation-pipeline/src/curation/validation/benchmarking.py) implementa:
- **Cohen's Kappa:** Avaliação de concordância entre especialistas humanos e o pipeline.
- **Matriz de Confusão:** Precisão, Revocação (Recall), F1-Score e Taxa de Rejeição Falsa (*False Rejection Rate*).

---

## 2. Abordagem de Ground Truth

Avaliações de concordância exigem consenso entre múltiplos anotadores especialistas, excluindo desacordos não resolvidos da matriz de decisão.
