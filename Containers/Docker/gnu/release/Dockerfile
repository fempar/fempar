FROM fempar/fempar-env:gnu

USER root:root

RUN source /opt/intel/mkl/bin/mklvars.sh intel64 \
    && ROOT_DIR=/tmp \
    && INSTALL_ROOT=/opt \
################################################ \
# Install json-fortran \
################################################ \
    && PACKAGE=json-fortran \
    && VERSION=7.0.0 \
    && JSON_FORTRAN_INSTALL=$INSTALL_ROOT/$PACKAGE/$VERSION \
    && URL=https://github.com/jacobwilliams/$PACKAGE/archive/$VERSION.tar.gz \
    && SOURCES_DIR=/tmp/$PACKAGE-$VERSION \
    && TAR_FILE=$SOURCES_DIR.tar.gz \
    && BUILD_DIR=$SOURCES_DIR/build \
    && wget -q $URL -O $SOURCES_DIR.tar.gz \
    && tar xzf $TAR_FILE -C $ROOT_DIR \
    && mkdir -p $BUILD_DIR $JSON_FORTRAN_INSTALL \
    && cd $BUILD_DIR \
    && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$JSON_FORTRAN_INSTALL $SOURCES_DIR \
    && cmake --build . --target install \
    && rm -rf $TAR_FILE $SOURCES_DIR \
################################################ \
# Install UMFPACK \
################################################ \
    && PACKAGE=SuiteSparse \
    && PACKAGE_LOWERCASE=suitesparse \
    && VERSION=5.4.0 \
    && TAR_FILENAME=$PACKAGE-$VERSION.tar.gz \
    && TAR_FILE=$ROOT_DIR/$TAR_FILENAME \
    && SOURCES_DIR=/tmp/$PACKAGE-$VERSION \
    && SUITESPARSE_INSTALL=$INSTALL_ROOT/$PACKAGE_LOWERCASE/$VERSION \
    && URL=http://faculty.cse.tamu.edu/davis/$PACKAGE/$TAR_FILENAME \
    && COMPONENTS="SuiteSparse_config AMD CAMD CCOLAMD COLAMD CHOLMOD UMFPACK" \
    && wget -q $URL -O $TAR_FILE \
    && mkdir -p $SOURCES_DIR \
    && tar xzf $TAR_FILE -C $SOURCES_DIR --strip-components=1 \
    && cd $SOURCES_DIR \
    && for COMPONENT in $COMPONENTS; do ( cd $SOURCES_DIR/$COMPONENT && make --quiet ); done \
    && for COMPONENT in $COMPONENTS; do ( cd $SOURCES_DIR/$COMPONENT && make --quiet install INSTALL=$SUITESPARSE_INSTALL ); done \
    && rm -rf $TAR_FILE $SOURCES_DIR \
################################################ \
# Install PETSC \
################################################ \
    && PACKAGE=petsc \
    && VERSION=3.9.2 \
    && PETSC_INSTALL=$INSTALL_ROOT/$PACKAGE/$VERSION \
    && TAR_FILE=$PACKAGE-$VERSION.tar.gz \
    && URL="http://ftp.mcs.anl.gov/pub/petsc/release-snapshots" \
    && SOURCES_DIR=$ROOT_DIR/$PACKAGE-$VERSION \
    && BUILD_DIR=$SOURCES_DIR/build \
    && wget -q $URL/$TAR_FILE -O $ROOT_DIR/$TAR_FILE \
    && mkdir -p $SOURCES_DIR \
    && tar xzf $ROOT_DIR/$TAR_FILE -C $SOURCES_DIR --strip-components=1 \
    && cd $SOURCES_DIR \
    && ./configure --with-mpi=1 --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
                   -with-blas-lapack-dir=/opt/intel/mkl/lib/intel64 \
                   --download-hypre=https://github.com/LLNL/hypre/archive/v2.14.0.tar.gz \
                   --with-x=0 --with-shared=1 --with-64-bit-indices --with-debugging=0 \
                   --prefix=$PETSC_INSTALL \
    && make --quiet \
    && make --quiet install \
    && rm -rf $ROOT_DIR/$TAR_FILE $SOURCES_DIR \
################################################ \
# Install QHULL \
################################################ \
    && PACKAGE=qhull \
    && VERSION=2015-src-7.2.0 \
    && QHULL_INSTALL=$INSTALL_ROOT/$PACKAGE/$VERSION \
    && TAR_FILE=$PACKAGE-$VERSION.tgz \
    && URL="http://www.qhull.org/download" \
    && SOURCES_DIR=$ROOT_DIR/$PACKAGE-$VERSION \
    && BUILD_DIR=$SOURCES_DIR/build \
    && wget -q $URL/$TAR_FILE -O $ROOT_DIR/$TAR_FILE \
    && mkdir -p $SOURCES_DIR/build \
    && tar xzf $ROOT_DIR/$TAR_FILE -C $SOURCES_DIR --strip-components=1 \
    && cd $SOURCES_DIR/build \
    && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$INSTALL_ROOT/$PACKAGE/$VERSION .. \
    && cmake --build . --target install \
    && rm -rf $ROOT_DIR/$TAR_FILE $SOURCES_DIR \
################################################ \
# Install HDF5 \
################################################ \
    && PACKAGE=hdf5 \
    && MAJOR=1 \
    && MINOR=8 \
    && PATCH=21 \
    && VERSION=$MAJOR.$MINOR.$PATCH \
    && HDF5_INSTALL=$INSTALL_ROOT/$PACKAGE/$VERSION \
    && TAR_FILE=$PACKAGE-$VERSION.tar.gz \
    && URL="https://support.hdfgroup.org/ftp/HDF5/releases/$PACKAGE-$MAJOR.$MINOR/$PACKAGE-$VERSION/src" \
    && SOURCES_DIR=$ROOT_DIR/$PACKAGE-$VERSION \
    && BUILD_DIR=$SOURCES_DIR/build \
    && wget -q $URL/$TAR_FILE -O $ROOT_DIR/$TAR_FILE \
    && mkdir -p $SOURCES_DIR \
    && tar xzf $ROOT_DIR/$TAR_FILE -C $SOURCES_DIR --strip-components=1 \
    && cd $SOURCES_DIR \
    && ./configure CC=mpicc FC=mpif90 CXX=mpicc --enable-parallel --enable-fortran --prefix=$HDF5_INSTALL \
    && make --quiet \
    && make --quiet install \
    && rm -rf $ROOT_DIR/$TAR_FILE $SOURCES_DIR 


USER fempar:fempar

################################################
# Export environment variables
################################################
ENV HDF5_ROOT /opt/hdf5/1.8.21
ENV PETSC_DIR /opt/petsc/3.9.2
ENV PATH $PATH:/opt/hdf5/1.8.21/bin:/opt/suitesparse/5.4.0:/opt/qhull/2015-src-7.2.0
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/opt/hdf5/1.8.21/lib:/opt/suitesparse/5.4.0/lib:/opt/qhull/2015-src-7.2.0/lib


