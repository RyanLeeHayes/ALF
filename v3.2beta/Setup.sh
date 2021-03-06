#! /bin/bash

module load cuda
module load anaconda/3.5.3.0

# First time
# From https://charmm-dev.org/wiki/index.php/Main_Page
python -m venv env-alf
source env-alf/bin/activate
pip install -e .

DIR=`pwd`
# echo "export CHARMM_LIB_DIR=$DIR/gnu/install-ALF/lib" > setupenv
echo "source $DIR/env-alf/bin/activate" > setupenv
