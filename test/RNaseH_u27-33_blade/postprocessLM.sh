#! /bin/bash

if [ $PHASE == 1 ]; then
  python -c "import alf; alf.SetupDCA($i,$NF,$FREQ)"
elif [ $PHASE == 2 ]; then
  iNF=$SLURM_ARRAY_TASK_ID
  python -c "import alf; alf.FilterDCA($i,$iNF,$NF,$FREQ)"
  python -c "import alf; alf.MomentDCA($i,$iNF,$NF,$FREQ)"
elif [ $PHASE == 3 ]; then
  python -c "import alf; alf.BSMomentDCA($i,$NF,$FREQ)"
elif [ $PHASE == 4 ]; then
  iBS=$SLURM_ARRAY_TASK_ID
  python -c "import alf; alf.LMDCA($i,$iBS,$NF,$FREQ)"
elif [ $PHASE == 5 ]; then
  python -c "import alf; alf.FinishDCA($i,$NF,$FREQ)"
fi
