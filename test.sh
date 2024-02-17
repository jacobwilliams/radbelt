#!/bin/bash

#
# basic script to build a conda environment, build the python extension, and test it
#

rm -rf ./env
conda env create --prefix ./env -f ./python/environment.yml
conda activate ./env

cd radbeltpy
chmod +x ./build_python.sh
./build_python.sh
cd ..
python ./test/radbeltpy_test.py
