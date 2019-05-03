FROM python:3-slim

LABEL MAINTAINER "vsande@cimne.upc.edu"
LABEL AUTHOR "vsande@cimne.upc.edu"
LABEL URL "http://fempar.org"

SHELL ["/bin/bash", "-c"]

RUN sed -i '/jessie-updates/d' /etc/apt/sources.list \ 
    && apt-get -qq -y update \
    && apt-get -qq -y install --no-install-recommends git gcc \
    && pip --no-cache-dir install ford \
################################################ \
# Clean packages \
################################################ \
    && apt-get -y clean \
    && apt-get -y autoclean \
    && apt-get -y autoremove \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /var/tmp/*


