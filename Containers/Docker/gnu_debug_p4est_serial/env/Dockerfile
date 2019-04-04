FROM debian:buster

LABEL MAINTAINER "vsande@cimne.upc.edu"
LABEL AUTHOR "vsande@cimne.upc.edu"
LABEL URL "http://fempar.org"

SHELL ["/bin/bash", "-c"]

################################################ \
# Common dependencies \
################################################ \
RUN apt-get -y update \
    && apt-get -y install --no-install-recommends wget git make cmake gcc g++ gfortran python openmpi-bin openmpi-common libopenmpi-dev flex bison gnupg apt-transport-https ca-certificates \
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
# Install json-fortran \
################################################ \
    && PACKAGE=json-fortran \
    && VERSION=7.0.0 \
    && JSON_FORTRAN_INSTALL=$INSTALL_ROOT/$PACKAGE/$VERSION \
    && URL=https://github.com/jacobwilliams/json-fortran/archive/$VERSION.tar.gz \
    && ROOT_DIR=/tmp \
    && SOURCES_DIR=/tmp/$PACKAGE-$VERSION \
    && TAR_FILE=$SOURCES_DIR.tar.gz \
    && BUILD_DIR=$SOURCES_DIR/build \
    && wget $URL -O $SOURCES_DIR.tar.gz \
    && tar -xzvf $TAR_FILE -C $ROOT_DIR \
    && mkdir -p $BUILD_DIR $JSON_FORTRAN_INSTALL \
    && cd $BUILD_DIR \
    && cmake -DCMAKE_INSTALL_PREFIX=$JSON_FORTRAN_INSTALL $SOURCES_DIR \
    && make install -j4 \
    && rm -rf $TAR_FILE $SOURCES_DIR \
################################################ \
# Install PETSC \
################################################ \
    && PACKAGE=petsc \
    && VERSION=3.9.2 \
    && PETSC_INSTALL=$INSTALL_ROOT/$PACKAGE/$VERSION \
    && TAR_FILE=$PACKAGE-$VERSION.tar.gz \
    && URL="http://ftp.mcs.anl.gov/pub/petsc/release-snapshots" \
    && ROOT_DIR=/tmp \
    && SOURCES_DIR=$ROOT_DIR/$PACKAGE-$VERSION \
    && BUILD_DIR=$SOURCES_DIR/build \
    && wget $URL/$TAR_FILE -O $ROOT_DIR/$TAR_FILE \
    && mkdir -p $SOURCES_DIR \
    && tar xzvf $ROOT_DIR/$TAR_FILE -C $SOURCES_DIR --strip-components=1 \
    && cd $SOURCES_DIR \
    && ./configure --with-mpi=1 --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
                   -with-blas-lapack-dir=/opt/intel/mkl/lib/intel64 \
                   --download-hypre=https://github.com/LLNL/hypre/archive/v2.14.0.tar.gz \
                   --with-x=0 --with-shared=1 --with-64-bit-indices --with-debugging \
                   --prefix=$PETSC_INSTALL \
    && make \
    && make install \
    && rm -rf $ROOT_DIR/$TAR_FILE $SOURCES_DIR \
################################################ \
# Install P4EST \
################################################ \
    && PACKAGE=p4est \
    && VERSION=2.2 \
    && P4EST_INSTALL=$INSTALL_ROOT/$PACKAGE/$VERSION \
    && TAR_FILE=$PACKAGE-$VERSION.tar.gz \
    && URL="https://github.com/p4est/p4est.github.io/raw/master/release" \
    && ROOT_DIR=/tmp \
    && SOURCES_DIR=$ROOT_DIR/$PACKAGE-$VERSION \
    && BUILD_DIR=$SOURCES_DIR/build \
    && wget $URL/$TAR_FILE -O $ROOT_DIR/$TAR_FILE \
    && mkdir -p $SOURCES_DIR \
    && tar xzvf $ROOT_DIR/$TAR_FILE -C $SOURCES_DIR --strip-components=1 \
    && cd $SOURCES_DIR \
    && ./configure --prefix=$P4EST_INSTALL --enable-debug --without-blas --without-lapack \
    && make \
    && make install \
    && rm -rf $ROOT_DIR/$TAR_FILE $SOURCES_DIR \
################################################ \
# Install QHULL \
################################################ \
    && PACKAGE=qhull \
    && VERSION=2015-src-7.2.0 \
    && QHULL_INSTALL=$INSTALL_ROOT/$PACKAGE/$VERSION \
    && TAR_FILE=$PACKAGE-$VERSION.tgz \
    && URL="http://www.qhull.org/download" \
    && ROOT_DIR=/tmp \
    && SOURCES_DIR=$ROOT_DIR/$PACKAGE-$VERSION \
    && BUILD_DIR=$SOURCES_DIR/build \
    && wget $URL/$TAR_FILE -O $ROOT_DIR/$TAR_FILE \
    && mkdir -p $SOURCES_DIR/build \
    && tar xzvf $ROOT_DIR/$TAR_FILE -C $SOURCES_DIR --strip-components=1 \
    && cd $SOURCES_DIR/build \
    && cmake -G "Unix Makefiles" .. \
    && cmake -DCMAKE_INSTALL_PREFIX=/opt/$PACKAGE/$VERSION .. \
    && make \
    && ctest \
    && make install \
    && rm -rf $ROOT_DIR/$TAR_FILE $SOURCES_DIR \
################################################ \
# Install HDF5 \
################################################ \
    && PACKAGE=hdf5 \
    && MAJOR=1 \
    && MINOR=10 \
    && PATCH=5 \
    && VERSION=$MAJOR.$MINOR.$PATCH \
    && HDF5_INSTALL=$INSTALL_ROOT/$PACKAGE/$VERSION \
    && TAR_FILE=$PACKAGE-$VERSION.tar.gz \
    && URL="https://support.hdfgroup.org/ftp/HDF5/releases/$PACKAGE-$MAJOR.$MINOR/$PACKAGE-$VERSION/src" \
    && ROOT_DIR=/tmp \
    && SOURCES_DIR=$ROOT_DIR/$PACKAGE-$VERSION \
    && BUILD_DIR=$SOURCES_DIR/build \
    && wget $URL/$TAR_FILE -O $ROOT_DIR/$TAR_FILE \
    && mkdir -p $SOURCES_DIR \
    && tar xzvf $ROOT_DIR/$TAR_FILE -C $SOURCES_DIR --strip-components=1 \
    && cd $SOURCES_DIR \
    && ./configure CC=mpicc FC=mpif90 CXX=mpicc --enable-parallel --enable-fortran --prefix=$HDF5_INSTALL \
    && make \
    && make install \
    && rm -rf $ROOT_DIR/$TAR_FILE $SOURCES_DIR \
################################################ \
# Clean packages \
################################################ \
    && apt-get -y clean \
    && apt-get -y autoclean \
    && apt-get -y autoremove 

################################################
# Export paths
################################################    
ENV MKL_THREADING_LAYER=GNU
ENV HDF5_ROOT /opt/hdf5/1.10.5
ENV PATH $PATH:/opt/qhull/2015-src-7.2.0
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/opt/qhull/2015-src-7.2.0/lib
ENV PETSC_DIR /opt/petsc/3.9.2
ENV PETSC_ARCH x86_64
ENV PATH $PATH:/opt/p4est/2.2

