FROM golang:stretch

LABEL MAINTAINER "vsande@cimne.upc.edu"
LABEL AUTHOR "vsande@cimne.upc.edu"
LABEL URL "http://fempar.org"

SHELL ["/bin/bash", "-c"]

RUN apt-get update \
    && apt-get install -y build-essential libssl-dev uuid-dev libgpgme11-dev squashfs-tools libseccomp-dev pkg-config \
    && go get -u github.com/golang/dep/cmd/dep \
################################################ \
# Install Singularity \
################################################ \
    && PACKAGE=singularity \
    && VERSION=3.2.1 \
    && mkdir -p $GOPATH/src/github.com/sylabs \
    && cd $GOPATH/src/github.com/sylabs \
    && wget https://github.com/sylabs/$PACKAGE/releases/download/v${VERSION}/$PACKAGE-${VERSION}.tar.gz \
    && tar -xzf $PACKAGE-${VERSION}.tar.gz \
    && cd ./$PACKAGE \
    && ./mconfig \
    && make -C ./builddir \
    && make -C ./builddir install \
    && mkdir -p $HOME/.$PACKAGE \
    && cp /usr/local/etc/$PACKAGE/remote.yaml $HOME/.$PACKAGE/ \
    && rm -rf $GOPATH/src/github.com/sylabs \
################################################ \
# Clean packages \
################################################ \
    && apt-get -y clean \
    && apt-get -y autoclean \
    && apt-get -y autoremove \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /var/tmp/*




