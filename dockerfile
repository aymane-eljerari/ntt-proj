FROM nvidia/cuda:12.6.2-devel-ubuntu24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    vim \
    clangd \
    lldb \
    libomp-dev \
    && rm -rf /var/lib/apt/lists/*

ENV PATH="/usr/local/cuda/bin:${PATH}"

WORKDIR /opt/workspace

RUN git clone https://github.com/openfheorg/openfhe-development.git openfhe && \
    cd openfhe && \
    cmake -B build -S . \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=/usr/local \
      -DBUILD_SHARED=ON \
      -DBUILD_STATIC=ON \
      -DBUILD_UNITTESTS=OFF \
      -DBUILD_EXAMPLES=OFF \
      -DBUILD_BENCHMARKS=ON \
      -DWITH_OPENMP=ON \
      -DWITH_NATIVEOPT=ON && \
    cmake --build build -j$(nproc) && \
    cmake --install build

CMD ["/bin/bash"]
