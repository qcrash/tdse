FROM ubuntu:latest AS base
# Install compilers and libraries
ENV DEBIAN_FRONTEND="noninteractive"
RUN apt-get update && apt-get install -y cmake gpg parallel gnuplot ffmpeg bc\
    emacs-nox gosu \
    gfortran gdb lsb-release gpg-agent curl procps
ARG TARGETARCH
ENV TARGETARCH=$TARGETARCH
#
FROM base AS toolchains-amd64
RUN curl https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null  
    RUN echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list
RUN apt-get update && apt-get install -y intel-oneapi-mkl-devel \	       
    && rm -rf /var/lib/apt/lists/*
#
FROM base AS toolchains-arm64
RUN curl -O https://developer.arm.com/-/cdn-downloads/permalink/Arm-Performance-Libraries/Version_25.07/arm-performance-libraries_25.07_deb_gcc.tar
RUN tar xf arm-performance-libraries_25.07_deb_gcc.tar
RUN ./arm-performance-libraries_25.07_deb/arm-performance-libraries_25.07_deb.sh --accept
RUN curl "https://developer.arm.com/packages/arm-toolchains%3Aubuntu-24/noble/Release.key" |  tee /etc/apt/trusted.gpg.d/developer-arm-com.asc
RUN echo "deb https://developer.arm.com/packages/arm-toolchains%3Aubuntu-24/noble/ ./"  | tee /etc/apt/sources.list.d/developer-arm-com.list
RUN apt-get update && apt-get install -y arm-performance-libraries environment-modules \
    && rm -rf /var/lib/apt/lists/*
#
FROM toolchains-${TARGETARCH} AS compile
COPY . /usr/src/qcrash
WORKDIR /usr/src/qcrash
RUN useradd -m docker && chown docker:docker .
SHELL ["/bin/bash", "-c"]
RUN source ./premake.sh && cmake -B build --preset $TARGETARCH -DCMAKE_INSTALL_PREFIX=./install && cmake --build build && cmake --install build
#
FROM compile AS run
COPY --from=compile /usr/src/qcrash/install/bin/tdse /usr/local/bin/
COPY --from=compile /usr/src/qcrash/entrypoint.sh /usr/local/bin/
COPY --from=compile /usr/src/qcrash/premake.sh /usr/local/bin/
#
# RUN useradd -m docker
# WORKDIR /home/docker
ENTRYPOINT ["/usr/src/qcrash/entrypoint.sh"]
CMD ["/bin/bash"]