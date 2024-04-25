FROM ubuntu:jammy

RUN apt-get update && \
  DEBIAN_FRONTEND='noninteractive' \
  DEBCONF_NONINTERACTIVE_SEEN='true' \
  apt-get install --yes \
    libnetcdf-dev \
    cdftools \
    gfortran \
    build-essential \
    git \
    curl \
    cmake 

    

RUN useradd \
  --create-home \
  hc_user 

USER hc_user

WORKDIR /home/hc_user

RUN git clone 'https://github.com/geodynamics/hc'; curl -L -O https://github.com/GenericMappingTools/gmt/releases/download/4.5.18/gmt-4.5.18-src.tar.bz2;
ENV GMT4HOME="/home/hc_user/gmt-4.5.18"
RUN tar -xvf gmt-4.5.18-src.tar.bz2; cd gmt-4.5.18; ./configure; make install-gmt; cd ..; cd hc; make all;
ENV PATH="/home/hc_user/hc/bin:${PATH}"
