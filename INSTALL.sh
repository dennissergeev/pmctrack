#!/usr/bin/env bash
mkdir -p build
cd build
cmake -DPLATFORM=linux ..
make
cd ..
