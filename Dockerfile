FROM python:3.12-slim
WORKDIR /app

# Instalar dependências
COPY requirements.txt .
RUN pip install -r requirements.txt

# Copiar o código da aplicação
COPY . .

# Comando correto para iniciar a sua app Streamlit
CMD ["streamlit", "run", "streamlit_app.py", "--server.port=8000", "--server.address=0.0.0.0"]