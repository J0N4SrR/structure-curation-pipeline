# Architecture Decision Records — Pipeline de Curadoria Estrutural

Registro formal das decisões de política química e de arquitetura que governam o
pipeline. Cada ADR fixa uma escolha, sua evidência empírica, o impacto no schema de
dados e a métrica que torna a consequência auditável.

Estas decisões são **anteriores ao código**: sete delas alteram o schema de saída e
três alteram o resultado químico. Alterar qualquer ADR invalida comparações entre
execuções, motivo pelo qual o hash deste arquivo é carimbado em cada lote
processado (ver [Rastreabilidade](#rastreabilidade)).

## Ambiente de referência

| Componente | Versão |
| --- | --- |
| Python | ≥ 3.10 |
| RDKit | 2026.03.5 |
| `chembl_structure_pipeline` | 1.2.4 |

Toda evidência empírica citada neste documento foi produzida contra este ambiente.
Mudanças de versão do RDKit alteram percepção de aromaticidade e regras de
padronização; a revalidação das evidências é obrigatória em qualquer upgrade.

## Legenda de status

- **Aceita** — vigente, implementada ou pendente de implementação.
- **Provisória** — vigente no MVP, com reavaliação agendada.
- **Substituída** — revogada por ADR posterior (referência cruzada obrigatória).

## Índice

| ID | Decisão | Status | Altera schema | Altera química |
| --- | --- | --- | :---: | :---: |
| [D-01](#d-01--motor-químico-central) | `chembl_structure_pipeline` intacto | Aceita | não | — |
| [D-02](#d-02--achatamento-de-tartaratos) | Aceitar `flatten_tartrate_mol` | Aceita | sim | sim |
| [D-03](#d-03--semântica-do-exclude_flag) | Preservar excluídos, não descartar | Aceita | sim | sim |
| [D-04](#d-04--política-tautomérica) | Preservar tautômero de entrada | Aceita | sim | não |
| [D-05](#d-05--política-estereoquímica) | Reter estereoquímica de entrada | Aceita | sim | não |
| [D-06](#d-06--entidades-multicomponentes) | Preservar misturas legítimas | Aceita | sim | sim |
| [D-07](#d-07--chave-de-identidade-estrutural) | InChIKey completo como identidade | Aceita | sim | não |
| [D-08](#d-08--fallback-de-parsing) | Sem OpenBabel/Indigo no MVP | Provisória | não | — |
| [D-09](#d-09--formato-de-entrada-e-escopo-do-checker) | SMILES primário | Provisória | sim | — |
| [D-10](#d-10--critérios-de-elegibilidade) | Cortes aplicados ao *parent* | Aceita | sim | — |
| [D-11](#d-11--política-isotópica) | Marcação isotópica removida | Aceita | sim | sim |
| [D-12](#d-12--entradas-sem-carbono) | Inorgânicos aceitos por padrão | Provisória | sim | não |

---

## D-01 — Motor químico central

**Status:** Aceita · 2026-09-07

### Contexto

A lógica de padronização e isolamento de estrutura-mãe do ChEMBL está publicada,
revisada por pares (Bento et al., *J Cheminform* 2020) e implementada em
`chembl_structure_pipeline` sob licença MIT. Reimplementá-la significa reproduzir
~1.040 linhas de química depurada ao longo de anos pelo EBI.

Três tentativas de reimplementação foram avaliadas e todas introduziram regressões
silenciosas: inversão da ordem `normalize → desalt` (perde ionização de metal
alcalino covalente), substituição do stripping por `LargestFragmentChooser` (destrói
misturas legítimas) e pré-normalização própria de nitro (redundante com o parser do
RDKit).

### Decisão

Consumir `chembl_structure_pipeline` como dependência direta, exclusivamente pela
API pública a nível de `Mol`:

```python
std       = standardize_mol(mol)
parent, _ = get_parent_mol(std)
```

Sem fork, sem monkeypatch e sem reimplementação de etapas internas no MVP. O código
próprio limita-se a parsing, instrumentação, proveniência, I/O e validação.

A sequência interna invocada é, na ordem (`standardizer.py`):

```
update_mol_valences → remove_sgroups → kekulize → remove_hs → normalize
→ uncharge → flatten_tartrate → cleanup_drawing → sanitize
→ get_isotope_parent → get_fragment_parent (solventes → sais → uncharge)
```

### Consequência conhecida e aceita

`get_fragment_parent_mol` abre `salts.smi` e `solvents.smi` e recompila 172 SMARTS
**a cada molécula** (custo medido: 1,62 ms/chamada). Isso responde por ~73% do tempo
total do pipeline.

| Configuração | Throughput |
| --- | --- |
| Biblioteca intacta | ~350 compostos/s |
| Com memoização dos SMARTS | ~1.084 compostos/s (2,6×) |

A 350 compostos/s, 1 milhão de estruturas processam em ~48 min em thread única, e o
problema é trivialmente paralelizável por registro. **O ganho de 2,6× não justifica
violar esta ADR no MVP.** O caminho correto é submeter PR ao projeto upstream (MIT).
Se a memoização vier a ser adotada localmente, torna-se um desvio formal, exigindo
nova ADR e teste de regressão provando saída idêntica à da biblioteca original.

### Impacto no schema

Nenhum campo novo. Exige o carimbo de `csp_version` e `rdkit_version` por registro.

### Consequência mensurável

Divergência estrutural entre a saída do pipeline e a chamada direta da biblioteca
sobre o mesmo corpus: **deve ser exatamente zero**. Verificada por teste de regressão
a cada build.

### Alternativas rejeitadas

- **Reimplementar a química em código próprio** — rejeitada: sem contribuição
  científica e com regressões demonstradas.
- **Fork com customizações** — rejeitada no MVP: custo de manutenção e perda do
  alinhamento com a referência, que é o baseline da validação (D-08, Fase 5).

---

## D-02 — Achatamento de tartaratos

**Status:** Aceita · 2026-09-07

### Contexto

`standardize_mol` chama `flatten_tartrate_mol` de forma **incondicional**. A assinatura
pública é `standardize_mol(m, check_exclusion=True, sanitize=True)` — não existe
parâmetro para desativá-la. Efeito medido:

```
entrada : OC(=O)[C@H](O)[C@@H](O)C(=O)O
saída   : O=C(O)C(O)C(O)C(=O)O          # dois centros quirais removidos
```

Trata-se de convenção interna de registro do ChEMBL, não de química universal, e
contraria o objetivo declarado de preservar centros quirais. Desativá-la exigiria
fork ou monkeypatch, violando [D-01](#d-01--motor-químico-central).

### Decisão

Aceitar o achatamento de tartaratos, com **registro obrigatório** no log de
proveniência. A perda é tornada visível e mensurável, não evitada.

O evento é detectado por sonda post-hoc comparando a contagem de centros
estereoquímicos antes e depois de `standardize_mol`, com confirmação por padrão
estrutural de tartarato.

### Impacto no schema

- `transformations[]` recebe evento `stage="standardize"`, `rule="tartrate_flattened"`,
  com SMILES antes/depois.
- `delta_stereocenters` registra a variação numérica.

### Consequência mensurável

Número absoluto e fração do dataset com estereoquímica perdida por esta regra,
reportados no sumário de cada lote. Se a fração for material para o endpoint em
estudo, esta ADR deve ser reaberta — e a reabertura é, ela própria, um resultado da
Fase 5.

### Alternativas rejeitadas

- **Monkeypatch para remover a chamada** — rejeitada: viola D-01 e cria divergência
  não rastreável da referência.
- **Ignorar silenciosamente** — rejeitada: a perda seria invisível na auditoria,
  exatamente o modo de falha que o pipeline existe para eliminar.

---

## D-03 — Semântica do `exclude_flag`

**Status:** Aceita · 2026-09-07

### Contexto

`exclude_flag.py` marca compostos contendo metais de transição (`METAL_LIST`, ~50
elementos; Zn é deliberadamente omitido pelo ChEMBL). O comportamento da referência é
**marcar e não padronizar**, retornando a molécula inalterada:

```
entrada     : [NH2][Pt]([NH2])([Cl])[Cl]      exclude_flag = True
standardize : [NH2][Pt]([NH2])([Cl])[Cl]      (idêntica)
```

Descartar esses registros produziria um dataset diferente do da referência e
eliminaria classes farmacologicamente relevantes (complexos de platina, por exemplo).

### Decisão

Preservar o registro com `status="PASSED"` e a flag de exclusão marcada. A
padronização é pulada pela própria biblioteca; o pipeline não acrescenta descarte.

Compostos excluídos **permanecem sujeitos** aos filtros de elegibilidade
([D-10](#d-10--critérios-de-elegibilidade)) e à deduplicação
([D-07](#d-07--chave-de-identidade-estrutural)).

### Impacto no schema

- `excluded_flag: bool` — campo obrigatório.
- `transformations[]` recebe evento `rule="standardization_skipped_excluded"`.

### Consequência mensurável

Fração do dataset com `excluded_flag=True`, ou seja, a fração que atravessou o
pipeline **sem padronização alguma**. Este número precisa constar de qualquer
publicação que use o dataset curado, sob risco de sobrestimar a homogeneidade do
conjunto.

### Alternativas rejeitadas

- **Rejeitar organometálicos** — rejeitada: divergiria da referência e removeria
  classes válidas sem justificativa no escopo declarado.

---

## D-04 — Política tautomérica

**Status:** Aceita · 2026-09-07

### Contexto

A canonicalização tautomérica é instável entre versões, altera a forma como
depositada pelo autor e pode reorganizar centros estereoquímicos. O ChEMBL não a
aplica.

Contudo, a premissa comum de que "sem canonicalização tautomérica não se detectam
duplicatas tautoméricas" é **falsa em parte**: a camada de hidrogênio móvel do InChI
padrão já normaliza tautomeria prototrópica em heteroátomos. Medição:

| Par tautomérico | InChIKey coincide? | `tautomer_key` coincide? |
| --- | :---: | :---: |
| 4-hidroxipiridina / 4-piridona | **sim** | sim |
| 2-hidroxipiridina / 2-piridona | **sim** | sim |
| amida / imidol | **sim** | sim |
| acetilacetona ceto / enol | **não** | **sim** |

O ganho real da coluna auxiliar é, portanto, estreito e bem definido: **tautomeria
com deslocamento em carbono** (1,3-dicarbonílicos e afins), que o InChI não trata
como tautomérica.

**Segunda razão, medida na ablação:** `TautomerEnumerator` tem `RemoveSp3Stereo`
ligado por padrão e **apaga centros quirais** que participem do sistema tautomérico:

```
N[C@@H](C)C(=O)O                 -> CC(N)C(=O)O            (quiralidade perdida)
OC(=O)[C@H](O)[C@@H](O)C(=O)O    -> O=C(O)C(O)C(O)C(=O)O   (quiralidade perdida)
C[C@H](O)CC(=O)C                 -> CC(=O)C[C@H](C)O       (preservada: centro isolado)
```

No corpus de ablação, a variante `canonical_tautomer` reduziu os centros quirais
definidos de 3 para 0. Canonicalizar tautômeros como política de produção seria,
portanto, incompatível com o objetivo declarado de preservar quiralidade — a menos
que `SetRemoveSp3Stereo(False)` seja aplicado explicitamente, o que altera o
resultado da canonicalização e exigiria ADR própria.

### Decisão

Preservar o tautômero de entrada na estrutura curada. Adicionar coluna auxiliar
`tautomer_key`, calculada via `rdMolStandardize.TautomerEnumerator().Canonicalize()`,
usada **exclusivamente para análise comparativa** na Fase 5 — nunca como chave de
deduplicação em produção.

### Impacto no schema

- `tautomer_key: Optional[str]` — coluna auxiliar, não participa da identidade.

### Consequência mensurável

Número de agrupamentos de duplicatas detectáveis por `tautomer_key` **e não** por
InChIKey completo. Dada a tabela acima, espera-se um número pequeno; se for
desprezível no corpus real, a coluna deve ser removida em ADR posterior por não
pagar seu custo computacional.

### Riscos operacionais

`TautomerEnumerator` tem custo combinatório em moléculas com muitos centros
tautoméricos. Aplicar limite explícito de tautômeros e tratar o estouro como valor
nulo com evento registrado, jamais como falha do registro.

---

## D-05 — Política estereoquímica

**Status:** Aceita · 2026-09-07

### Contexto

Curadoria estrutural não tem informação para atribuir configuração a um centro
indefinido. Qualquer atribuição automática inventa dado. Ao mesmo tempo, centros
indefinidos são material para modelagem: um dataset com alta fração de estereoquímica
parcial suporta conclusões mais fracas.

### Decisão

Reter a representação estereoquímica exatamente como recebida. Não atribuir, não
enumerar e não remover centros. Quantificar explicitamente o que está indefinido:

```python
rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol)   # indefinidos
rdMolDescriptors.CalcNumAtomStereoCenters(mol)              # total
```

Saída canônica em *isomeric SMILES*, preservando a informação existente.

Esta ADR interage com [D-02](#d-02--achatamento-de-tartaratos): tartaratos perdem
estereoquímica por regra do motor, e a perda aparece em `delta_stereocenters`.

### Impacto no schema

- `n_stereocenters_total: int`
- `n_stereocenters_undefined: int`
- `delta_stereocenters: int` — variação entrada → parent.

### Consequência mensurável

Fração de registros com pelo menos um centro indefinido, e fração com
`delta_stereocenters ≠ 0` (perda induzida pelo pipeline, não pela entrada). Os dois
números são distintos e ambos obrigatórios no relatório de lote.

---

## D-06 — Entidades multicomponentes

**Status:** Aceita · 2026-09-07

### Contexto

O ChEMBL rejeita explicitamente o removedor de sais nativo do RDKit — comentário em
`standardizer.py`: *"there are a number of special cases for the ChEMBL salt
stripping, so we can't use the salt remover that's built into the RDKit
standardizer."* Quando restam fragmentos distintos, ele os **combina** via
`CombineMols`, não escolhe o maior.

Medição com cotrimoxazol (trimetoprima + sulfametoxazol, medicamento combinado real):

```
LargestFragmentChooser : COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC          # perdeu o sulfametoxazol
ChEMBL get_parent_mol  : COc1cc(...)c1OC.Cc1cc(NS(=O)(=O)...)no1  # preserva ambos
```

Duas regras adicionais da referência devem ser conhecidas: fragmentos idênticos são
deduplicados por SMILES (2 ácidos acéticos → 1), e existe um guard de segurança —
se **todos** os componentes casarem com a tabela de sais, nada é removido. Este guard
faz com que benzoato de sódio (`[Na+].[O-]C(=O)c1ccccc1`) atravesse o pipeline
intacto, pois tanto `Sodium` quanto `Benzoate` constam de `salts.smi`.

### Decisão

Preservar o comportamento da referência integralmente. **Proibido** o uso de
`LargestFragmentChooser` ou de qualquer heurística de "maior fragmento" em qualquer
ponto do pipeline.

### Impacto no schema

- `removed_fragments: Optional[str]` — SMILES dos componentes removidos.
- `n_components_parent: int` — número de componentes na estrutura-mãe.

### Consequência mensurável

Número de entidades com `n_components_parent > 1` preservadas. Um pipeline com
`LargestFragmentChooser` reportaria zero, o que serve como teste de regressão direto.

### Interação crítica

Ver [D-10](#d-10--critérios-de-elegibilidade): quando o *parent* é multicomponente, os
cortes físico-químicos incidem sobre a soma dos componentes retidos.

---

## D-07 — Chave de identidade estrutural

**Status:** Aceita · 2026-09-07

### Contexto

O InChIKey tem estrutura em blocos: os 14 primeiros caracteres codificam apenas
conectividade, **ignorando estereoquímica, isótopos e carga**. Usá-los como chave de
deduplicação funde enantiômeros:

```
L-alanina            QNAYBMKLOCPYGJ-REOHCLBHSA-N   bloco1 = QNAYBMKLOCPYGJ
D-alanina            QNAYBMKLOCPYGJ-UWTATZPHSA-N   bloco1 = QNAYBMKLOCPYGJ
alanina sem estereo  QNAYBMKLOCPYGJ-UHFFFAOYSA-N   bloco1 = QNAYBMKLOCPYGJ
```

Fundir enantiômeros é inaceitável em um pipeline cujo objetivo declarado é preservar
quiralidade, e é o modo de falha clássico em QSAR quando os isômeros têm atividades
distintas.

O último bloco codifica protonação: ácido benzoico termina em `-N` e benzoato em `-M`,
com bloco1 idêntico.

### Decisão

- **Identidade / deduplicação:** InChIKey **completo**, obrigatório.
- **Agrupamento taxonômico:** primeiro bloco de 14 caracteres, apenas como
  agrupador para inspeção — permite localizar famílias sal/ácido livre e conjuntos de
  estereoisômeros relacionados.

Colisões de bloco1 com InChIKeys completos distintos **não** são duplicatas; são
relatadas para revisão, com o tipo de colisão classificado (estereoisomeria,
protonação, isotopologia).

### Impacto no schema

- `inchikey: Optional[str]` — chave de identidade.
- `inchikey_block1: Optional[str]` — agrupador.

### Consequência mensurável

Contagem de colisões de bloco1 por tipo. Cada tipo tem política de resolução própria
e nenhuma é resolvida automaticamente.

---

## D-08 — Fallback de parsing

**Status:** Provisória · 2026-09-07 · reavaliação na Fase 5

### Contexto

OpenBabel e Indigo foram cogitados como recuperadores de SMILES que o RDKit rejeita.
Três problemas: (i) OpenBabel é GPL-2.0, com implicações de distribuição que RDKit
(BSD) e Indigo (Apache-2.0) não têm; (ii) OpenBabel corrige valências
permissivamente, podendo devolver molécula quimicamente distinta da entrada; (iii)
não há evidência de que a taxa de recuperação justifique a dependência.

Trata-se de uma **pergunta empírica não respondida**, tratada até aqui como premissa.

### Decisão

Excluir ambos do MVP. A questão é realocada como estudo comparativo na Fase 5:
qual a taxa de recuperação, e que fração das recuperações é quimicamente correta sob
curadoria manual?

Se reintroduzido, o composto recuperado carrega flag permanente de proveniência e
revalidação estrita pelo RDKit é condição de aceitação.

### Impacto no schema

Nenhum no MVP. O campo `parser_source` é reservado para reintrodução futura.

### Consequência mensurável

Contagem de `ERR_SYNTAX` no corpus real — é o teto do ganho possível. Se for
desprezível, a dependência não se justifica e o resultado negativo é publicável.

---

## D-09 — Formato de entrada e escopo do checker

**Status:** Provisória · 2026-09-07 · reavaliação após o MVP

### Contexto

O `checker.py` do ChEMBL opera sobre molblocks V2000 e atribui penalidades de 2 a 7
(7 = *Illegal input*; exclusão do ChEMBL a partir de 6). Dos 17 checkers, apenas 3
podem disparar a partir de SMILES — *InChI warnings*, `num_atoms_equals_zero` e
`disallowed_radical`. Os outros 14 dependem de informação inexistente em SMILES:
coordenadas 3D, flag 3D, tipos e estéreo-flags de ligação, átomos sobrepostos,
coordenadas zeradas, ligações cruzadas em anel, estéreo-ligação em anel, informação
de polímero e blocos V3000.

Invocar a escala 2–7 sobre entrada SMILES seria indefensável: ~82% do instrumento
estaria inativo.

### Decisão

SMILES como formato primário de entrada. A validação do MVP é **topológica e de
valência**, nomeada como tal, sem qualquer referência à escala de penalidades do
ChEMBL.

Suporte a SDF/molblock e ativação do checker geométrico completo ficam adiados para
ADR posterior. A arquitetura reserva o ponto de extensão, mas o MVP não o promete.

**Restrição de ordem, derivada de evidência.** A validação de valência **não pode
rejeitar antes** de `standardize_mol`, porque o normalizador da referência repara
classes inteiras que a sanitização estrita rejeitaria:

| Entrada | Sanitização estrita | ChEMBL `standardize_mol` |
| --- | --- | --- |
| `C[N](C)(C)C` (N quaternário neutro) | `AtomValenceException` | `C[N+](C)(C)C` |
| `C[N](C)(C)CCO` | `AtomValenceException` | `C[N+](C)(C)CCO` |
| `c1ccccc1[N]#N` (diazônio neutro) | `AtomValenceException` | `N#[N+]c1ccccc1` |

Um portão anterior ao motor descartaria sais de amônio quaternário e de diazônio
válidos. A fase de parsing **instrumenta e registra**; a rejeição por valência ocorre
somente após o motor.

### Impacto no schema

- `rejection_stage: Optional[str]` — distingue rejeição no parsing (sintaxe) de
  rejeição pós-motor (valência), tornando a restrição de ordem auditável.
- `checker_penalties` fica reservado, inativo enquanto D-09 estiver vigente.

### Consequência mensurável

Contagem de rejeições por estágio. Rejeição por valência **antes** do motor é bug de
implementação, não resultado — o teste de regressão com os três casos acima trava a
ordem correta.

---

## D-10 — Critérios de elegibilidade

**Status:** Aceita · 2026-09-07

### Contexto

Cortes físico-químicos delimitam o escopo químico do estudo. Dois cuidados: aplicá-los
ao sal, e não à estrutura-mãe, mede massa de contra-íon e rejeita compostos válidos;
e cortes fixos em código impedem que o escopo seja um parâmetro do experimento.

Nota terminológica: estes são **critérios de elegibilidade / escopo químico**, não
"domínio de aplicabilidade". Domínio de aplicabilidade é conceito dependente de modelo
(cobertura do espaço de descritores, *leverage*, distância ao modelo) e pertence à
etapa de modelagem.

### Decisão

Aplicar os cortes **estritamente sobre a estrutura-mãe isolada**, após
`get_parent_mol`, nunca sobre a forma salificada. Valores parametrizáveis, com
padrões declarados:

| Critério | Padrão | Cálculo |
| --- | --- | --- |
| Peso molecular | ≤ 1000 Da | `Descriptors.MolWt` (massa média) |
| Átomos pesados | ≤ 100 | `mol.GetNumHeavyAtoms()` |

Os padrões são ponto de partida documentado, não constante universal. Domínios como
produtos naturais, macrociclos e peptídeos exigem revisão explícita.

### Interação crítica com D-06

Quando o *parent* é multicomponente, os cortes incidem sobre a **soma dos componentes
retidos**. Exemplo medido — cotrimoxazol:

```
parent: MW = 543,6 Da   heavy = 38   componentes = 2
```

Uma combinação de dois fármacos pesados pode exceder o corte e ser rejeitada como
"molécula grande", quando na verdade são duas moléculas médias. O `rejection_detail`
deve registrar `n_components_parent` sempre que a rejeição for por peso, para que o
caso seja distinguível na revisão.

### Impacto no schema

- `parent_mw: Optional[float]`
- `parent_heavy_atoms: Optional[int]`
- `rejection_code` admite `ERR_MW_CUTOFF` e `ERR_HEAVY_ATOMS_CUTOFF`.

### Consequência mensurável

Número de rejeições por critério, e — separadamente — número de rejeições por peso em
que `n_components_parent > 1`. O segundo número identifica falsos positivos
introduzidos pela interação com D-06.

---

## D-11 — Política isotópica

**Status:** Aceita · 2026-09-07

### Contexto

`get_parent_mol` chama `get_isotope_parent_mol` antes do stripping de fragmentos,
zerando o número de massa de todo átomo marcado. A política estava **implementada
mas não registrada**, e uma política não registrada não é auditável: quem lê o
dataset não tem como saber se um composto deuterado foi preservado ou colapsado no
seu análogo comum.

Comportamento medido, consistente em todos os casos testados:

| entrada | saída | igual ao análogo não marcado? |
| --- | --- | :---: |
| `[2H]C(Cl)(Cl)Cl` | `ClC(Cl)Cl` | sim |
| `[12CH3][13CH3]` | `CC` | sim |
| `[2H]c1ccccc1` | `c1ccccc1` | sim |
| `[2H]C([2H])([2H])C([2H])([2H])[2H]` | `CC` | sim |
| `[13CH4]` | `C` | sim |

### Decisão

**REMOVE.** A marcação isotópica é descartada na obtenção da estrutura-mãe,
seguindo a referência. Não há modo de preservação, e o pipeline nunca testa as duas
políticas ao mesmo tempo.

### Consequência mensurável

Compostos que diferem **apenas** pela marcação isotópica colapsam na mesma
identidade e são contabilizados como duplicatas exatas. Um estudo que use padrões
internos deuterados precisa saber disto antes de deduplicar: a fração do dataset
afetada é o número a reportar.

### Alternativas rejeitadas

- **Preservar isótopos** — rejeitada: exigiria contornar `get_isotope_parent_mol`,
  violando a [D-01](#d-01--motor-químico-central).
- **Coluna auxiliar com a identidade isotópica** — não rejeitada, apenas adiada:
  faz sentido se e quando o corpus contiver compostos marcados em número relevante.

---

## D-12 — Entradas sem carbono

**Status:** Provisória · 2026-09-07 · reavaliação com o corpus real

### Contexto

Entradas puramente inorgânicas atravessam o pipeline e são aprovadas como estrutura
curada:

| entrada | resultado | `exclude_flag` |
| --- | --- | :---: |
| `[Na+].[Cl-]` | `[Cl-].[Na+]` aprovado | não |
| `[Cu+2]` | `[Cu+2]` aprovado | sim |
| `CC(=O)[O-].[Cu+2].CC(=O)[O-]` | inalterado, aprovado | sim |
| `CCO.CC(=O)C` (só solventes) | ambos aprovados | não |

O cloreto de sódio passa porque tanto `Sodium` quanto `Chloride` constam de
`salts.smi`: dispara o guard "todos os componentes são sal, então mantém tudo" da
[D-06](#d-06--entidades-multicomponentes). Não é acidente da implementação — é o
comportamento da referência.

Isso nunca havia sido decidido. Para um pipeline de moléculas pequenas orgânicas,
NaCl aprovado como composto curado é ruído; mas rejeitá-lo no motor divergiria da
referência, contra a [D-01](#d-01--motor-químico-central).

### Decisão

O **motor não muda**: a semântica da referência é preservada e entradas sem carbono
continuam atravessando. A exclusão vira **critério de elegibilidade opcional**,
`require_carbon`, no mesmo lugar e com a mesma natureza dos cortes de peso e de
átomos pesados ([D-10](#d-10--critérios-de-elegibilidade)) — escopo do estudo, não
química.

**Padrão: desligado.** Ligar por omissão mudaria silenciosamente o resultado de
lotes já processados. Quem quiser o dataset estritamente orgânico opta por ele:

```
curation --input entrada.csv --out-dir saida --require-carbon
```

### Impacto no schema

- `EligibilityCriteria.require_carbon: bool = False`
- `rejection_code` admite `ERR_NO_CARBON`

### Consequência mensurável

Número de entradas sem carbono no corpus, com e sem o critério ligado. Se a fração
for material, esta ADR deve ser promovida a Aceita com padrão invertido — e essa
promoção é, ela própria, um resultado da caracterização do dataset.

### Alternativas rejeitadas

- **Rejeitar no motor** — rejeitada: divergiria da referência sem justificativa
  química, apenas de escopo.
- **Ligar por padrão** — rejeitada por ora: muda resultado sem evidência sobre a
  composição do corpus real.

---

## Rastreabilidade

Cada registro processado carrega:

| Campo | Origem |
| --- | --- |
| `rdkit_version` | `rdkit.__version__` |
| `csp_version` | `chembl_structure_pipeline.__version__` |
| `pipeline_version` | versão do pacote próprio |
| `policy_hash` | SHA-256 deste arquivo |

O `policy_hash` é o que torna duas execuções comparáveis. Alterar qualquer ADR muda o
hash e, por construção, marca os lotes anteriores como produzidos sob outra política.

## Matriz de dependências entre decisões

| Decisão | Depende de | Interage com |
| --- | --- | --- |
| D-02 | D-01 (não há como desativar sem violar) | D-05 (mede a perda) |
| D-03 | D-01 | D-07, D-10 (excluídos seguem filtrados) |
| D-04 | — | D-07 (auxiliar, não identidade) |
| D-06 | D-01 | D-10 (soma dos componentes) |
| D-09 | — | D-01 (ordem: instrumentar antes, rejeitar depois) |
| D-10 | D-06 | D-03, D-12 |
| D-11 | D-01 (não há como preservar sem violar) | D-07 (isotopólogos colapsam) |
| D-12 | D-01, D-06 (guard "tudo é sal") | D-10 (mesmo estágio) |

## Registro de revisões

| Data | Alteração |
| --- | --- |
| 2026-09-07 | Criação do documento com D-01 a D-10 |
| 2026-09-07 | D-04: acrescentada a evidência de perda de estereoquímica na canonicalização tautomérica, medida na ablação da Fase 5 |
| 2026-09-07 | D-11: política isotópica registrada (já estava implementada, não documentada) |
| 2026-09-07 | D-12: contrato para entradas sem carbono, com critério opcional `require_carbon` |
