FROM debian:latest AS base
# Install compilers and libraries
ENV DEBIAN_FRONTEND="noninteractive"
RUN apt-get update && apt-get install -y make gpg emacs-nox \
    gfortran gdb lsb-release gpg-agent wget
RUN wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null  
    RUN echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list
RUN apt-get update && apt-get install -y intel-oneapi-mkl-devel \	       
    && rm -rf /var/lib/apt/lists/*
#
FROM base AS compile
COPY . /usr/src/qcrash
WORKDIR /usr/src/qcrash
SHELL ["/bin/bash", "-c"]
RUN source /opt/intel/oneapi/setvars.sh && make clean && make
#
FROM compile AS run
COPY --from=compile /usr/src/qcrash/tdse /usr/local/bin/