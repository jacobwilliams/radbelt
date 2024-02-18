#!/bin/bash

#
# basic script to build a conda environment, build the python extension, and test it
#

rm -rf ./env
conda env create --prefix ./env -f ./python/environment.yml
conda activate ./env

cd python
chmod +x ./build_python.sh
./build_python.sh
cd ..

rm -rf ./build-python
mkdir ./build-python
mkdir ./build-python/radbeltpy

cp ./python/*.py ./build-python/radbeltpy
cp ./python/*.so ./build-python/radbeltpy
cp ./python/*.dll ./build-python/radbeltpy
cp ./python/*.dylib ./build-python/radbeltpy
cp -r data ./build-python/radbeltpy


python ./test/radbeltpy_test.py
