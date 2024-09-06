FROM python:3-slim
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y dist-upgrade && \
    apt-get clean && apt-get autoremove && rm -rf /var/lib/apt/lists/*



