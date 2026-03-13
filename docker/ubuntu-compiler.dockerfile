# Use the official Ubuntu image as the base
FROM ubuntu:latest

# Set non-interactive mode to prevent installation prompts
ENV DEBIAN_FRONTEND=noninteractive

# Combine update, install, and cleanup in ONE block to prevent caching errors
# and reduce image size.
RUN apt-get update && apt-get install -y \
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
    libcurl4-openssl-dev \
    libxml2-dev \
    libtool \
    python3 \
    python-is-python3 \
    && rm -rf /var/lib/apt/lists/*

# Verify the installation (Optional - usually better to do this in the container, not build)
RUN gcc --version && g++ --version
