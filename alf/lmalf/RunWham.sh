#! /bin/bash

$( dirname -- "${BASH_SOURCE[0]}" )/whamweight $1 Energy Lambda $2 $3
mkdir multisite
mv G*.dat C.dat V.dat multisite/
