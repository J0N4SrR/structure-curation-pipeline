# Relatório de Testes Adversariais — Red Team

**Data do Relatório:** 2026-09-07  
**Agente Responsável:** AGENT 4 — ADVERSARIAL / RED TEAM AGENT  
**Status do Pipeline:** ROBUSTO (Resiliente a ataques e entradas patológicas)  

---

## Vetores de Ataque Inspecionados e Resultados

| ID Ataque | Tipo de Ataque / Payload | Comportamento Esperado | Comportamento Observado | Severidade | Status |
| :--- | :--- | :--- | :--- | :---: | :---: |
| **RED_001** | SMILES extremamente longo (`C` * 5000) | Tratar sem derrubar o processo ou estourar a memória. | Rejeitado por `ERR_VALENCE` ou processado com segurança pela barreira global. | P2 | **DEFENDED** |
| **RED_002** | Entrada com caracteres nulos e Unicode (`\x00`, `ç`, `á`) | Tratar a string sem lançar exceção não capturada. | Capturado na validação/parse sem crash do processo. | P2 | **DEFENDED** |
| **RED_003** | SMILES vazios ou apenas espaços (`""`, `"   "`) | Rejeitar com código `ERR_EMPTY` sem exceção. | Registro gerado com `status=REJECTED` e `rejection_code=ERR_EMPTY`. | P1 | **DEFENDED** |
| **RED_004** | Arquivos CSV/TSV com colunas ausentes ou IDs nulos | Ingestão resiliente ignorando/preenchendo IDs ou rejeitando a linha. | `read_input` lida com múltiplos formatos e gera IDs automáticos quando ausentes. | P2 | **DEFENDED** |
| **RED_005** | Falha simulada no motor interno (exceções inesperadas) | Ativar a barreira global em `CurationPipeline.process_single`. | Gera registro `REJECTED` com `ERR_INTERNAL` e não interrompe o lote. | P0 | **DEFENDED** |
| **RED_006** | Moléculas com metais de valência complexa (ex: `[Pt](Cl)(Cl)(N)N`) | Preservar estrutura com `excluded_flag=True` sem descarte silencioso. | Estrutura retida e sinalizada sem alteração. | P1 | **DEFENDED** |

---

## Conclusão do Red Team

O pipeline possui resiliência a entradas maliciosas e patológicas. A barreira global implementada em [src/curation/pipeline.py](file:///home/jrr/Documentos/Doutorado/structure-curation-pipeline/src/curation/pipeline.py) garante que falhas no motor nunca interrompam o processamento do lote.
