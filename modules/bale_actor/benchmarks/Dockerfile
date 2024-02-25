FROM ubuntu:18.04
LABEL authors="Akihiro Hayashi <ahayashi@gatech.edu>,Sri Raj Paul <sriraj@gatech.edu>"

RUN apt-get update && apt-get install -y \
    apt-utils \
    autoconf \ 
    git \
    libtool \
    pkg-config \
    vim \
    gcc \
    g++ \
    libmpich-dev \
    mpich \
    unzip

ENV LOCAL=/root/local

WORKDIR /root
RUN git clone https://github.com/ofiwg/libfabric.git
WORKDIR /root/libfabric
RUN git checkout tags/v1.10.0rc3 -b v1.10.0rc3
RUN ./autogen.sh && ./configure --prefix=$LOCAL && make install

WORKDIR /root
RUN git clone https://github.com/Sandia-OpenSHMEM/SOS.git
WORKDIR /root/SOS
RUN git checkout tags/v1.4.5 -b v1.4.5
RUN ./autogen.sh && ./configure --prefix=$LOCAL --with-ofi=$LOCAL --enable-pmi-simple && make install

ENV CC=/root/local/bin/oshcc
ENV CXX=/root/local/bin/oshc++
ENV OSHRUN=/root/local/bin/oshrun

WORKDIR /root
RUN git clone https://github.com/jdevinney/bale.git
WORKDIR /root/bale
RUN git checkout tags/bale-2.1 -b bale-2.1
RUN chmod +x ./install.sh
RUN ./install.sh -s -f

ENV BALE_INSTALL=/root/bale/build_unknown

WORKDIR /root
RUN git clone https://github.com/srirajpaul/hclib
WORKDIR /root/hclib
RUN git fetch && git checkout bale_actor
RUN ./install.sh
ENV HCLIB_ROOT=/root/hclib/hclib-install
RUN cd modules/bale_actor && make

WORKDIR /root/hclib/modules/bale_actor/benchmarks
RUN unzip ../inc/boost.zip -d ../inc/

ENV LD_LIBRARY_PATH=$LOCAL/lib:$BALE_INSTALL/lib:$HCLIB_ROOT/lib:$HCLIB_ROOT/../modules/bale_actor/lib
ENV HCLIB_WORKERS=1

