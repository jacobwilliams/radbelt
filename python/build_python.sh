#!/bin/bash

#
# Build the python extension module using f2py
#

# this is to eliminate the "ld: warning: could not create compact unwind" on MacOS
export LDFLAGS=-Wl,-no_compact_unwind

f2py -c --backend meson ./radbeltpy.pyf \
    ../src/radbelt_kinds_module.F90 \
    ../src/trmfun.f90 \
    ../src/shellig.f90 \
    ../src/radbelt_module.f90 \
    ../src/radbelt_c_module.f90