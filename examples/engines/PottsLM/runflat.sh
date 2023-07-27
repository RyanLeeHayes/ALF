#!/bin/bash

python -c "import alf; alf.initialize(engine='bladelib')"
python -c "import alf; alf.runflat(1,200,13000,39000,engine='bladelib')"
python -c "import alf; alf.runflat(201,210,125000,375000,engine='bladelib')"
