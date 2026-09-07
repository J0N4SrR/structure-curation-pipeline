# Deploy no Streamlit Community Cloud

## Configuração

| campo | valor |
| --- | --- |
| Repository | `J0N4SrR/structure-curation-pipeline` |
| Branch | `main` |
| Main file path | `streamlit_app.py` |

`streamlit_app.py` é um shim: o layout `src/` faz o Streamlit colocar
`src/curation/` no `sys.path` em vez de `src/`, e sem o shim o import
`from curation.app import main` falha.

## Limitação conhecida: sem imagens das estruturas

O app roda no Cloud **sem renderização 2D**. A causa não está neste repositório.

O wheel do RDKit empacota cairo, freetype e libpng em `rdkit.libs`, mas
`rdMolDraw2D` linka contra bibliotecas X11 do sistema:

```
libXrender.so.1 -> libxrender1
libX11.so.6     -> libx11-6
libXext.so.6    -> libxext6
```

O container do Cloud não as tem, e **não é possível instalá-las**: a imagem base
mistura fontes Debian *bullseye* e *trixie*, e o `Release` do bullseye-security
está expirado, então `apt-get update` retorna não-zero e o Streamlit aborta o
passo de dependências inteiro:

```
E: Release file for .../bullseye-security/InRelease is expired
❗️ installer returned a non-zero exit code
```

Um `packages.txt` com os nomes corretos falha do mesmo jeito — o erro é anterior
à instalação. Pior: o passo abortado impede o redeploy, deixando o app preso em
código antigo.

Por isso `packages.txt` foi removido. O app importa `Draw` sob `try/except` e,
sem backend gráfico, exibe um aviso e mantém SMILES em texto. Todas as demais
funcionalidades — funil, exclusões, rastreabilidade, proveniência, exportações —
funcionam normalmente.

Para reavaliar: recriar `packages.txt` com `libxrender1`, `libx11-6`, `libxext6`
e verificar se o `apt-get update` do Cloud passa a retornar zero. Sem comentários
no arquivo — o Cloud o alimenta ao apt via `xargs`, que trataria `#` como nome de
pacote.

## Ambiente

O Cloud usa **Python 3.14** e não permite trocar a versão depois da criação do
app. O `.devcontainer` usa 3.11 e o desenvolvimento local, 3.12. As três variantes
carimbam `versions.python` diferente na proveniência de cada lote. O RDKit está
fixado na mesma versão nas três, então o resultado químico deve coincidir — mas se
dois lotes divergirem, é o primeiro campo a comparar.

## Para rodar localmente com imagens

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[app]"
streamlit run streamlit_app.py
```
