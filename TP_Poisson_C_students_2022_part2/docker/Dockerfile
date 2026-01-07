FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get install -y \
    build-essential \
    gcc \
    make \
    vim \
    libc6-dev \
    libblas-dev \
    liblapacke-dev \
    libatlas-base-dev \
    gnuplot

WORKDIR /app

ENV HOSTNAME=ambre

CMD ["tail", "-f", "/dev/null"]