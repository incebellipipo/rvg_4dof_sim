FROM ubuntu:22.04

COPY requirements.txt /tmp/requirements.txt

RUN apt update && apt upgrade -yq ;\
    apt install -yq \
        python3-pip \
        python3-venv ;\
    pip install -r /tmp/requirements.txt ;\
    apt clean