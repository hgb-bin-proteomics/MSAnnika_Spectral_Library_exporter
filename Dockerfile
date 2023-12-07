# Dockerfile with MS Annika Spectral Library exporter
# author: Micha Birklbauer
# version: 1.0.0

FROM python:3.12.0

LABEL maintainer="micha.birklbauer@gmail.com"

RUN mkdir app
COPY requirements.txt app
COPY config.py app
COPY create_spectral_library.py app
COPY gui/streamlit_util.py app
COPY gui/streamlit_app.py app
WORKDIR app
RUN pip install -r requirements.txt
RUN pip install streamlit

CMD  ["streamlit", "run", "streamlit_app.py", "--server.maxUploadSize", "5000"]
