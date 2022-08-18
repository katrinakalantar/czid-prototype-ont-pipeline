FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

LABEL maintainer="IDseq Team idseq-tech@chanzuckerberg.com"
LABEL description = "Image for consensus genome by metagenomic sequencing with spiked primer enrichment or amplicon sequencing"

ARG fastp_version=0.23.1
ARG FLYE_VER="2.9"

RUN apt-get -qq update && apt-get -qq -y install curl locales zip git make build-essential libz-dev wget \
  jq \
  moreutils \
  curl \
  locales \
  zip \
  unzip \
  httpie \
  zlib1g-dev \
  libfreetype6-dev \
  libhts-dev \
  pkg-config \
  apt-utils \
  minimap2 \
  samtools \
  seqtk \
  python3-dev \
  python3-pip \
  python3-setuptools \
  python3-wheel \
  && locale-gen en_US.UTF-8 && apt-get clean && apt-get autoclean

RUN wget -O fastp http://opengene.org/fastp/fastp.${fastp_version} && \
    chmod a+x ./fastp && \
    cp fastp /bin/fastp

#RUN git clone https://github.com/fenderglass/Flye && \
#    cd Flye && \
#    python setup.py install

RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.3.0/seqkit_linux_386.tar.gz && \
  tar -xvf seqkit_linux_386.tar.gz && \
  rm seqkit_linux_386.tar.gz

RUN wget https://github.com/fenderglass/Flye/archive/${FLYE_VER}.tar.gz && \
  tar -xvf ${FLYE_VER}.tar.gz && \
  rm ${FLYE_VER}.tar.gz && \
  cd Flye-${FLYE_VER} && \
  python3 setup.py build && \
  python3 setup.py install

RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/downloads/centrifuge-1.0.3-beta-Linux_x86_64.zip && \
  unzip centrifuge-1.0.3-beta-Linux_x86_64.zip
ENV PATH="/centrifuge-1.0.3-beta:$PATH"

# copied from https://github.com/chanzuckerberg/czid-workflows/blob/main/short-read-mngs/Dockerfile
RUN apt-get -y update && apt-get install -y build-essential libz-dev git python3-pip cmake

WORKDIR /tmp
RUN git clone https://github.com/morsecodist/diamond
WORKDIR /tmp/diamond
RUN git checkout minimal-mods
WORKDIR /tmp/diamond/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make -j6
RUN mv diamond /usr/local/bin

# ARTIC, Medaka and dependencies -- errors due to...
# Could not install packages due to an EnvironmentError: [Errno 28] No space left on device
#RUN apt-get install -y python3-cffi python3-h5py python3-intervaltree python3-edlib git
#RUN pip3 install ont-fast5-api parasail mappy pyspoa protobuf==3.19.* tensorflow==2.7.1 https://github.com/artic-network/fieldbioinformatics/archive/1.2.1.tar.gz
#RUN pip3 install medaka==1.6.0 --no-deps
#RUN pip3 install git+https://github.com/rzlim08/PyVCF.git