"""Componentes de interface customizados.

Ficam separados de :mod:`curation.app` porque tem ciclo de vida proprio: cada um
serve arquivos estaticos e conversa com o Streamlit por mensagem, nao por chamada
de funcao.
"""
