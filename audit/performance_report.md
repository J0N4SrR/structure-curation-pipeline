# Relatório de Desempenho e Escalabilidade

**Data do Relatório:** 2026-09-07  
**Agente Responsável:** AGENT 10 — PERFORMANCE AGENT  
**Throughput Medido:** ~350 compostos/s (Thread Única)  

---

## Métricas de Desempenho Medidas

- **Throughput Baseline:** ~350 compostos por segundo em thread única.
- **Gargalo Identificado:** Recompilação dos 172 SMARTS dentro do `get_fragment_parent_mol` da biblioteca oficial (responde por ~73% do tempo total).
- **Consumo de Memória:** Constante. Ingestão em geradores via `process_many` evita carregar o dataset inteiro na memória RAM.
