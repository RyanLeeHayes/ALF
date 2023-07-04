#! /bin/bash

python -m venv env-alf
source env-alf/bin/activate
pip install -e .

DIR=`pwd`
echo "source $DIR/env-alf/bin/activate" > setupenv
