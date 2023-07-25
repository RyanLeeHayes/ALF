#! /bin/bash -l
# Compile the executables in this directory with cmake.
# You may need to add code to set up your environment.

cmake ./
make $1 wham
