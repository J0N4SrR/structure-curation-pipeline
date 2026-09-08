# Auditoria de Testes — Engenharia de Qualidade

**Data da Auditoria:** 2026-09-07  
**Agente Responsável:** AGENT 5 — TEST ENGINEER  
**Status da Suíte:** APROVADA (207/207 testes executados)  

---

## 1. Cobertura da Suíte por Módulo

| Módulo Alvo | Arquivo de Teste | Nº de Testes | Áreas Cobertas |
| :--- | :--- | :---: | :--- |
| `curation.engine` | `test_engine.py` | 37 | Estágios de parsing, padronização, dessalificação e canonicalização. |
| `curation.pipeline` | `test_regression.py` | 37 | Precedência iônica, valência posterior, misturas, exclude_flag, idempotência. |
| `curation.probes` | `test_probes.py` | 17 | Sondas post-hoc, diff de fronteira, contagem de estereocentros e cargas. |
| `curation.io` | `test_io.py` | 25 | Escrita atômica em disco, suporte a CSV/TSV/SDF/JSONL, leitores preguiçosos. |
| `curation.dedup` | `test_dedup.py` | 14 | Identidade por InChIKey completo e detecção de colisões de Bloco 1. |
| `curation.filters` | `test_filters.py` | 12 | Filtros de massa molecular e átomos pesados sobre a parent mol. |
| `curation.reporting` | `test_reporting.py` | 25 | Relatórios agregados, pacotes de exportação e integridade do manifesto. |
| `curation.validation` | `test_validation.py` | 26 | Estudos de ablação, gráficos comparativos, métricas estatísticas (Kappa). |
| `curation.cli` | `test_cli.py` | 15 | Comandos da linha de código e parâmetros de execução. |

---

## 2. Testes de Idempotência e Determinismo

- **Idempotência (\(f(x) == f(f(x))\)):** Verificada em `test_regression.py::test_strict_idempotency` para um corpus de 14 compostos complexos.
- **Determinismo (\(run_1(x) == run_2(x)\)):** Verificado em `test_reporting.py::test_run_id_is_deterministic_for_same_input`.
