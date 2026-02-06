#! /bin/bash -l
# Compile the executables in this directory with cmake.
# You may need to add code to set up your environment.

# for release
cmake ./

# # or for debug
# cmake -DCMAKE_BUILD_TYPE=Release ./

# then make
make $1 all
