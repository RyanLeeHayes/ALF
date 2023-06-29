#! /bin/bash

source env-blade
source ../ALF-v3.2beta-20230131/setupenv # python with alf module

python -c "import alf; alf.postvolume($i,$eqS,$S,$N,$skipE,engine='blade')"
