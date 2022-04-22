#! /bin/bash

module load anaconda/3.5.3.0

./difference.py
./chargecorrect_dsc.py
./couplings.py
./evaluate.py
