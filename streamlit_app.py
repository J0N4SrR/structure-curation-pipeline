"""Ponto de entrada para o Streamlit Community Cloud.

O pacote vive em ``src/curation`` (layout src). O Streamlit executa
``streamlit run streamlit_app.py`` a partir da raiz do repositório e coloca no
``sys.path`` o diretório **do script**, não ``src/`` — sem este shim, o import
``from curation.app import main`` falha no deploy.

Manter o shim fino é deliberado: nenhuma lógica aqui, apenas resolução de caminho.
Em desenvolvimento local, com o pacote instalado (``pip install -e .``), rodar
``streamlit run src/curation/app.py`` continua funcionando normalmente.
"""

import sys
from pathlib import Path

SRC = Path(__file__).resolve().parent / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from curation.app import main  # noqa: E402  (import após ajuste do sys.path)

main()
