#!/bin/bash

# Install system dependencies
sudo apt-get update
sudo apt-get install -y docker.io build-essential python3-dev

# Python dependencies
pip3 install --user setuptools cython numpy scipy gepard

# Build KM15 Cython module
python3 setup.py build_ext --inplace

# Build dvcsgen
cd dependencies/dvcsgen
make clean
make
chmod +x dvcsgen
cd ../..