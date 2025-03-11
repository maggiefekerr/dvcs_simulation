#!/bin/bash

# Load JLab environment modules
source /etc/profile.d/modules.sh
module purge
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

echo "Installation complete! Test with:"
echo "python3 main.py --model km15 --trig 10 --fname test"