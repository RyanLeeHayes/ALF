#!/bin/bash

export BLADEDIR=/home/rhaye/BLaDE
source $BLADEDIR/modules
export BLADEEXEC=$BLADEDIR/build/blade
source ../setupenv # python with alf module

python -c "import alf; alf.initialize(engine='blade')"
python -c "import alf; alf.runflat(1,200,13000,39000,engine='blade')"
python -c "import alf; alf.runflat(201,210,125000,375000,engine='blade')"
