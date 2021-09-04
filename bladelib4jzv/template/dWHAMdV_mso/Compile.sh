#! /bin/bash -l

source modules

export CC="gcc"
export CXX="g++"

cmake ./

make $1 wham
