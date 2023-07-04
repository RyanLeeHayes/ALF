#!/bin/bash

python -c "import alf; alf.initialize()"
python -c "import alf; alf.runflat(1,200,13000,39000)"
python -c "import alf; alf.runflat(201,210,125000,375000)"
