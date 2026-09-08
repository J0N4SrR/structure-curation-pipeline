# Relatório de Estudos de Ablação — Impacto Metodológico

**Data do Relatório:** 2026-09-07  
**Agente Responsável:** AGENT 12 — ABLATION AGENT  

---

## Variantes Metodológicas Avaliadas

O módulo [src/curation/validation/ablation.py](file:///home/jrr/Documentos/Doutorado/structure-curation-pipeline/src/curation/validation/ablation.py) avalia empiricamente variações de política contra a linha de base (*baseline*):
1. `flatten_tartrates_off`: Mede a preservação da estereoquímica de tartaratos.
2. `discard_organometallics`: Mede o impacto do descarte compulsório de compostos organometálicos.
3. `tautomer_canonicalization_on`: Mede a perda de estereoquímica $sp^3$ ao aplicar canonicalização tautomérica.

Gráficos e tabelas comparativas são gerados em Markdown e PNG para o relatório final.
