#!/bin/bash

# Load JLab environment
source /etc/profile.d/modules.sh
module load clas12

# Build KM15 Cython module
python3 setup.py build_ext --inplace

# Install Python dependencies
pip3 install --user cython numpy scipy gepard

# Build dvcsgen
cd dependencies/dvcsgen
make clean
make
chmod +x dvcsgen
cd ../..

# Set PYTHONPATH
echo "export PYTHONPATH=\$PYTHONPATH:$(pwd)" >> ~/.bashrc
# source cshrc

echo "Installation complete! Test with:"
echo "python3 main.py --model km15 --trig 10 --fname test"