#! /bin/bash

source env-charmm
source ../setupenv # python with alf module

python -c "import alf; alf.postprocess($i,$eqS,$S,$N,$skipE,True)"
