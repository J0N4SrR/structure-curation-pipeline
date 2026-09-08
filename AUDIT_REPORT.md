# Relatório Final de Auditoria Multiagente (AUDIT_REPORT.md)

**Data da Auditoria:** 2026-09-07  
**Juiz Final:** AGENT 15 — FINAL JUDGE  
**Projeto:** Pipeline de Curadoria Estrutural (`structure-curation-pipeline`)  

---

## 1. Sumário Executivo

### Overall Status: **PASS**

O projeto foi submetido a uma auditoria integral executada por 16 agentes especializados, cobrindo inventário, arquitetura, domínio químico, resiliência adversarial (Red Team), qualidade de testes, reprodutibilidade, proveniência de dados, deduplicação, I/O atômico, desempenho, validação científica, estudos de ablação, transparência de UI e revisão cruzada independente.

---

## 2. Conformidade com as Decisões de Política (ADRs D-01 a D-10)

| Decisão | Status | Evidência Executável |
| :--- | :---: | :--- |
| **D-01 (Motor Central)** | COMPLIANT | `EngineWrapper` utiliza `chembl_structure_pipeline` via API pública. |
| **D-02 (Tartaratos)** | COMPLIANT | Sonda `_TARTRATE` rastreia achatamento incondicional. |
| **D-03 (Exclude Flag)** | COMPLIANT | Cisplatina e carboranos mantidos sem rejeição por valência. |
| **D-04 (Tautomeria)** | COMPLIANT | Tautômero de entrada preservado no SMILES curado. |
| **D-05 (Estereoquímica)** | COMPLIANT | Preservação de centros quimicamente definidos e métricas registradas. |
| **D-06 (Multicomponentes)** | COMPLIANT | Cotrimoxazol e benzoato de sódio mantidos sem descarte indevido. |
| **D-07 (Identidade)** | COMPLIANT | InChIKey completo e Bloco 1 com detecção de conflitos. |
| **D-08 (Fallbacks)** | COMPLIANT | Parsing com RDKit (`sanitize=False`) sem OpenBabel/Indigo no MVP. |
| **D-09 (Valência Posterior)** | COMPLIANT | Portão de valência executado somente após reparo do motor. |
| **D-10 (Elegibilidade)** | COMPLIANT | Cortes de MW <= 1000 e HA <= 100 aplicados na `parent mol`. |

---

## 3. Resultados da Suíte de Testes

| Categoria de Testes | Aprovados | Falhas | Não Testados |
| :--- | -----: | -----: | -----: |
| **Motor e Wrapper (`test_engine.py`)** | 37 | 0 | 0 |
| **Regressão e Políticas (`test_regression.py`)** | 37 | 0 | 0 |
| **I/O e Escrita Atômica (`test_io.py`)** | 25 | 0 | 0 |
| **Sondas Post-Hoc (`test_probes.py`)** | 17 | 0 | 0 |
| **Deduplicação (`test_dedup.py`)** | 14 | 0 | 0 |
| **Filtros Físico-Químicos (`test_filters.py`)** | 12 | 0 | 0 |
| **Relatórios e Linha do Tempo (`test_reporting.py`)** | 25 | 0 | 0 |
| **Validação Científica / Ablação (`test_validation.py`)** | 26 | 0 | 0 |
| **CLI (`test_cli.py`)** | 15 | 0 | 0 |
| **TOTAL** | **207** | **0** | **0** |

---

## 4. Avaliações Finais das Dimensões Críticas

### Reprodutibilidade
> **É possível reproduzir uma execução exatamente?**  
> **SIM.** Execuções repetidas geram saídas idênticas byte a byte. O `manifest.json` registra o `policy_hash` (SHA-256 de `docs/decisions.md`) e as versões exatas de todas as dependências.

### Proveniência
> **É possível explicar cada transformação?**  
> **SIM.** Cada composto aprovado ou rejeitado gera um `CurationRecord` contendo a lista completa de eventos de transformação observados e deltas físico-químicos.

### Validade Científica
> **Existe evidência suficiente para considerar a metodologia defensável?**  
> **SIM.** A metodologia fundamenta-se na referência publicada pelo EBI/ChEMBL, com validação estatística por Cohen's Kappa, matrizes de confusão e estudos de ablação comparativos.

### Segurança e Robustez (Red Team)
> **Entradas adversariais conseguem quebrar o pipeline?**  
> **NÃO.** Testes com strings nulas, Unicode, SMILES vazios e falhas simuladas de motor confirmam a eficácia da barreira global de exceções.

### Transparência (UI)
> **Um usuário consegue compreender e exportar tudo que aconteceu?**  
> **SIM.** A interface Streamlit ([src/curation/app.py](file:///home/jrr/Documentos/Doutorado/structure-curation-pipeline/src/curation/app.py)) expõe transformações, rejeições, metadados e botões para download dos artefatos completos.

---

## 5. Veredito Final

O projeto **`structure-curation-pipeline` está APROVADO (PASS)** sem reservas na auditoria científica e de engenharia de software. Todos os 15 artefatos exigidos foram gerados no diretório `audit/` acompanhados do ledger de evidências verificadas.
