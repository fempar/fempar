FROM fempar/fempar-env:gnu-debug_p4est-parallel

USER root:root

################################################ 
# Fempar compilation 
################################################ 
RUN source /opt/intel/mkl/bin/mklvars.sh intel64 \
    && PACKAGE=fempar \
    && VERSION=experimental \
    && URL="https://gitlab.com/$PACKAGE/$PACKAGE.git" \
    && ROOT_DIR=/tmp \
    && INSTALL_ROOT=/opt \
    && SOURCES_DIR=$ROOT_DIR/$PACKAGE-$VERSION \
    && BUILD_DIR=$INSTALL_ROOT/$PACKAGE \
    && THIRDPARTY_BUILD_DIR=$INSTALL_ROOT/$PACKAGE-thirdparty \
    && FORTRAN_EXTRA_FLAGS="-DFORTRAN_EXTRA_FLAGS=-fimplicit-none" \
    && git clone --single-branch --branch $VERSION --recursive $URL $SOURCES_DIR \
    && ln -sf $SOURCES_DIR/CMake/CTestConfig_CIMNE.cmake $SOURCES_DIR/CTestConfig.cmake \
    && mkdir -p $BUILD_DIR $THIRDPARTY_BUILD_DIR \
################################################ \
# Build fempar thirdparty libraries \
################################################ \
    && cd $THIRDPARTY_BUILD_DIR \
    && cmake -DCMAKE_BUILD_TYPE=DEBUG $FORTRAN_EXTRA_FLAGS $SOURCES_DIR/ThirdParty \
    && cmake --build . \
################################################ \
# Build fempar library \
################################################ \
    && cd $BUILD_DIR \
    && cmake -DCMAKE_BUILD_TYPE=DEBUG -DFEMPAR_ENABLE_TESTS=ON -DFEMPAR_THIRDPARTY_DIR=$THIRDPARTY_BUILD_DIR -DMPIEXEC_PREFLAGS="--allow-run-as-root -oversubscribe --mca btl_vader_single_copy_mechanism none" $SOURCES_DIR \
    && ctest --timeout 99999 -j8 -V -L PARALLEL_FAST -D ExperimentalUpdate -D ExperimentalStart -D ExperimentalConfigure -D ExperimentalBuild -D ExperimentalTest -D ExperimentalSubmit \
################################################ \
# Clean sources \
################################################ \
    && rm -rf $SOURCES_DIR 

USER fempar:fempar

################################################
# Export paths
################################################
ENV FEMPAR_DIR /opt/fempar
