# Auditoria de I/O e Atomicidade — Integridade de Escrita

**Data da Auditoria:** 2026-09-07  
**Agente Responsável:** AGENT 9 — I/O AND ATOMICITY AUDITOR  
**Status de Atomicidade:** GARANTIDA VIA BATCHWRITER  

---

## Escrita Atômica e Rollback

O gerenciador de contexto `BatchWriter` em [src/curation/io.py](file:///home/jrr/Documentos/Doutorado/structure-curation-pipeline/src/curation/io.py) escreve temporariamente em arquivos `.tmp`. 

- Se o lote for concluído sem erros, os arquivos são renomeados atomicamente para `curated.csv`, `rejected.csv` e `manifest.json`.
- Se ocorrer uma interrupção ou falha catastrófica durante a escrita, todos os arquivos temporários são limpos e o manifesto não é gerado, impedindo a existência de saídas parciais ou corrompidas no disco.
