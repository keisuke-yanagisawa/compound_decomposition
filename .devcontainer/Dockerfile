FROM ubuntu:22.04

SHELL ["/bin/bash", "-c"]

# To suppress interactive installation of boost
ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt upgrade -y
RUN apt install -y make g++ libboost-all-dev \
                   wget cmake git


# environment variables for openbabel-2.4.1
# Note that openbabel-3.0.0 CAN'T be used
ENV OBABEL_VERSION=2.4.1
ENV OPEN_BABEL_ROOT=$HOME/.local/openbabel-${OBABEL_VERSION}
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OPEN_BABEL_ROOT/lib

# installation of openbabel version 2.4.1
RUN wget https://sourceforge.net/projects/openbabel/files/openbabel/${OBABEL_VERSION}/openbabel-${OBABEL_VERSION}.tar.gz/download -O openbabel-${OBABEL_VERSION}.tar.gz
RUN tar -zxvf openbabel-${OBABEL_VERSION}.tar.gz
RUN mkdir -p openbabel-${OBABEL_VERSION}/build 
RUN cd openbabel-${OBABEL_VERSION}/build && \
    cmake -DCMAKE_INSTALL_PREFIX=$OPEN_BABEL_ROOT -DCMAKE_CXX_FLAGS="-Wno-narrowing" ..
RUN cd openbabel-${OBABEL_VERSION}/build && make -j8 install

ENV OBABEL_INSTALL_PATH=$OPEN_BABEL_ROOT

# change permission to read ssh settings if exists
RUN chown root:$USER /root/.ssh/config || true
RUN chmod 600        /root/.ssh/config || true