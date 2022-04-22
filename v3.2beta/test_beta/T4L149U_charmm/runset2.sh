#!/bin/bash

export CHARMMDIR=/home/rhaye/CHARMM/chv1/charmm_exe
source $CHARMMDIR/modules
export CHARMMEXEC=$CHARMMDIR/gnu/charmm
source ../setupenv # python with alf module

python -c "import alf; alf.initialize()"
python -c "import alf; alf.runflat(1,100,13000,39000)"
python -c "import alf; alf.runflat(101,110,125000,375000)"
