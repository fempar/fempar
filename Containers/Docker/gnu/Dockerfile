FROM debian:buster

LABEL MAINTAINER "vsande@cimne.upc.edu"
LABEL AUTHOR "vsande@cimne.upc.edu"
LABEL URL "http://fempar.org"

SHELL ["/bin/bash", "-c"]

################################################ \
# Common dependencies \
################################################ \
RUN apt-get -y update \
    && apt-get -y install --no-install-recommends wget git make cmake vim gcc g++ gfortran python metis libmetis-dev openmpi-bin openmpi-common libopenmpi-dev flex bison gnupg apt-transport-https ca-certificates zlib1g zlib1g-dev valgrind \
    && INSTALL_ROOT=/opt \
################################################ \
# Install MKL \
################################################ \
    && PACKAGE=intel-mkl \
    && VERSION=64bit-2019.3-062 \
    && URL=https://apt.repos.intel.com \
    && KEY=GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB \
    && ROOT_DIR=/tmp \
    && wget $URL/intel-gpg-keys/$KEY -O $ROOT_DIR/$KEY \
    && apt-key add $ROOT_DIR/$KEY \
    && echo "deb $URL/mkl all main" > /etc/apt/sources.list.d/$PACKAGE.list \
    && apt-get -y update \
    && apt-get -y install --no-install-recommends $PACKAGE-$VERSION \
    && rm -rf $ROOT_DIR/$KEY \
    && echo "source /opt/intel/mkl/bin/mklvars.sh intel64" >> /etc/bash.bashrc \
    && source /etc/bash.bashrc \
################################################ \
# Clean packages \
################################################ \
    && apt-get -y clean \
    && apt-get -y autoclean \
    && apt-get -y autoremove \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /var/tmp/*

################################################
# Create fempar user and group, and workdir
################################################
RUN groupadd -r fempar -g 1000 \
    && useradd --no-log-init -r -g fempar -u 1000 fempar \
    && mkdir -p /data \
    && chmod 0777 /data

WORKDIR /data
USER fempar:fempar

################################################
# Export environment variables
################################################
ENV MKL_THREADING_LAYER=GNU



