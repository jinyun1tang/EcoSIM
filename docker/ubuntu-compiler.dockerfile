# Use a specific version for stability
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# We add 'ca-certificates' to ensure curl/git can verify SSL connections
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    git \
    cmake \
    vim \
    gdb \
    automake \
    autoconf \
    libxml2 \
    curl \
    ca-certificates \
    libcurl4-openssl-dev \
    libxml2-dev \
    libtool \
    python3 \
    python-is-python3 \
    && rm -rf /var/lib/apt/lists/*

RUN gcc --version && g++ --version