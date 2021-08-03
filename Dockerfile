FROM ubuntu:focal

ENV DEBIAN_FRONTEND=noninteractive
ENV LC_ALL C.UTF-8

RUN apt-get update
RUN apt-get install -y gcc-8 g++-8

ADD https://github.com/vdemichev/DiaNN/releases/download/1.8/diann_1.8.deb .
RUN apt-get install ./diann_1.8.deb
RUN ln -sf /usr/diann/1.8/diann-1.8 /usr/bin/diann
