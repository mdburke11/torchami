#!/usr/bin/env bash

BUILD_DIR=build
INSTALL_PREFIX="$PWD/install"

echo "[1/3] Configure"
cmake -S . -B "$BUILD_DIR" \
  -DCMAKE_BUILD_TYPE=Release \
  -DMAKE_PYTAMI=OFF

echo "[2/3] Build"
cmake --build "$BUILD_DIR" --parallel

echo "[3/3] Test"
ctest --test-dir "$BUILD_DIR" --output-on-failure || true

echo "[Install] -> $INSTALL_PREFIX"
cmake --install "$BUILD_DIR" --prefix "$INSTALL_PREFIX"
