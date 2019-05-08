FROM gcc:9.1.0

LABEL MAINTAINER "vsande@cimne.upc.edu"
LABEL AUTHOR "vsande@cimne.upc.edu"
LABEL URL "http://fempar.org"

SHELL ["/bin/bash", "-c"]

################################################ \
# Common dependencies \
################################################ \
# This repo fail in gcc:6.2.0 \
RUN apt-get -qq -y update \
    && apt-get -qq -y install --no-install-recommends wget git make cmake vim metis libmetis-dev flex bison gnupg apt-transport-https ca-certificates zlib1g zlib1g-dev valgrind \
    && apt-get -qq -y install --no-install-recommends dkms infiniband-diags libibverbs* ibacm librdmacm* libmlx4* libmlx5* mstflint libibcm.* libibmad.* libibumad* opensm srptools libmlx4-dev librdmacm-dev rdmacm-utils ibverbs-utils perftest vlan ibutils \
    && INSTALL_ROOT=/opt \
################################################ \
# Install MKL \
################################################ \
    && PACKAGE=intel-mkl \
    && VERSION=64bit-2019.3-062 \
    && URL=https://apt.repos.intel.com \
    && KEY=GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB \
    && ROOT_DIR=/tmp \
    && wget -q $URL/intel-gpg-keys/$KEY -O $ROOT_DIR/$KEY \
    && apt-key add $ROOT_DIR/$KEY \
    && echo "deb $URL/mkl all main" > /etc/apt/sources.list.d/$PACKAGE.list \
    && apt-get -qq -y update \
    && apt-get -qq -y install --no-install-recommends $PACKAGE-$VERSION \
    && rm -rf $ROOT_DIR/$KEY \
    && echo "source /opt/intel/mkl/bin/mklvars.sh intel64" >> /etc/bash.bashrc \
    && source /etc/bash.bashrc \
################################################ \
# Install OpenMPI \
################################################ \
    && OMPI_MAJOR_VERSION=3 \
    && OMPI_MINOR_VERSION=1 \
    && OMPI_PATCH_VERSION=3 \
    && OMPI_VERSION=$OMPI_MAJOR_VERSION.$OMPI_MINOR_VERSION.$OMPI_PATCH_VERSION \
    && TAR_FILE=openmpi-$OMPI_VERSION.tar.gz \
    && URL="https://www.open-mpi.org/software/ompi/v$OMPI_MAJOR_VERSION.$OMPI_MINOR_VERSION/downloads" \
    && ROOT_DIR=/tmp \
    && SOURCES_DIR=$ROOT_DIR/openmpi-$OMPI_VERSION \
    && BUILD_DIR=$SOURCES_DIR/build \
    && mkdir -p $BUILD_DIR \
    && cd $ROOT_DIR \
    && wget -q $URL/$TAR_FILE -O $ROOT_DIR/$TAR_FILE \
    && tar xzf $ROOT_DIR/$TAR_FILE -C $SOURCES_DIR --strip-components=1 \
    && cd $BUILD_DIR \
    && $SOURCES_DIR/configure --enable-shared --with-verbs --enable-mpirun-prefix-by-default --with-hwloc --with-libevent --disable-dlopen --with-pmix --prefix=/usr \
    && make --quiet all install \
    && rm -rf $ROOT_DIR/$TAR_FILE $SOURCES_DIR \
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
ENV PATH $PATH:/usr/bin
ENV LD_LIBRARY_PATH /lib:/lib/x86_64-linux-gnu:/usr/lib:/usr/lib/x86_64-linux-gnu:/usr/local/lib:/usr/local/lib64:$LD_LIBRARY_PATH
ENV INCLUDE $INCLUDE:/usr/include/openmpi



