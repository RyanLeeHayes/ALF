#! /bin/bash -l

cmake ./
make $1 all VERBOSE=1
