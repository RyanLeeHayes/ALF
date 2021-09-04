#! /bin/bash -l

source modules

export CC="gcc"
export CXX="g++"

cmake -DCMAKE_BUILD_TYPE=Debug ./

make $1 PLM
