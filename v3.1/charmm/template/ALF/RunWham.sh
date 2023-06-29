#! /bin/bash

source ../dWHAMdV_mso/modules

# nvprof -o profile.out ../dWHAMdV/wham $1 Energy Lambda
../dWHAMdV_mso/wham $1 Energy Lambda $2 $3
mkdir multisite
mv G*.dat C.dat V.dat multisite/
