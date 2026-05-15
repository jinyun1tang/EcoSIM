#!/bin/bash

# Define variables for easy updates
IMAGE_NAME="my-compiler"
DOCKERFILE_PATH="docker/ubuntu-compiler.dockerfile"
MOUNT_DIR="/EcoSIM"

# In your run_compiler.sh
# Usage: ./run_compiler.sh --fresh
if [[ "$1" == "--fresh" ]]; then
    CACHE_FLAG="--no-cache"
else
    CACHE_FLAG=""
fi

docker build $CACHE_FLAG -f "$DOCKERFILE_PATH" -t "$IMAGE_NAME" .

echo "--- 1. Building the Docker image ($IMAGE_NAME) ---"

# Check if build was successful
if [ $? -eq 0 ]; then
    echo "--- 2. Starting the container ---"
    docker run -it \
        -v "$(pwd):$MOUNT_DIR" \
        -w "$MOUNT_DIR" \
        "$IMAGE_NAME"
else
    echo "Error: Docker build failed. Container will not start."
    exit 1
fi
