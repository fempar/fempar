FROM fempar/fempar-env:gnu-debug

USER root:root

################################################ 
# Install P4EST 
################################################ 
RUN PACKAGE=p4est \
    && VERSION=2.2 \
    && INSTALL_ROOT=/opt \
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
    && rm -rf $ROOT_DIR/$TAR_FILE $SOURCES_DIR 

USER fempar:fempar

################################################
# Export paths
################################################
ENV PATH $PATH:/opt/p4est/2.2




