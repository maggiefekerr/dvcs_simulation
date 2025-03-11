#!/bin/bash

# Load JLab modules
module load python/3.9.7
module load gcc/9.3.0

# Python dependencies
pip3 install --user cython numpy scipy gepard

# Build KM15 Cython module
python3 setup.py build_ext --inplace

# Build dvcsgen
cd dependencies/dvcsgen
make clean
make
chmod +x dvcsgen
cd ../..