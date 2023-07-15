"""alf runs adaptive landscape flattening to improve sampling in lambda dynamics simulations

Routines:
There are several levels of routines within alf (the module). The first group
of routiness provides basic functionality for running lambda dynamics
simulations with ALF (the method). Examples of how to use them are highlighted
in various examples/engines subdirectories.

Basic Routines:
 * initialize - sets up the analysis0 directory, variables1.inp file, and other
     objects that must be created before beginning ALF
 * runflat - runs several iterative short simulations to get rough                   approximations for the biases; arguments control how many cycles are used
     and how long these simulations are. Common practice is to run many 100 ps
     simulations to get close, and then further refine biases with several more
     1 ns simulations
 * runprod - used for running longer production simulations, usually with
     several indepdendent trials
 * postprocess - analyzes production simulations from runprod to further refine
     biases for additional production simulations and/or produce a Results.txt
     file containing free energy estimates for dG using the histogram-based
     estimator. Relative ddG may be obtained by comparing these results with dG
     for the same perturbations in another physical ensemble

Additional routines are provided to estimate free energies using the Potts model
estimator rather than the histogram-based estimator. This is a more complicated
estimator useful primarily for systems with many perturbation sites.

Potts Model Routines:
 * SetupDCA - copy some relevant files into the dca directory where the analysis
     is performed. Requires postprocess to have been performed on the same
     simulations first
 * FilterDCA - discretize the continuous lambda trajectories to determine which
     (if any) substituent is dominant at each site at each time
 * MomentDCA - calculate moments, or single site probabilities and joint
     probabilities at pairs of sites that are used for likelihood optimization
 * BSMomentDCA - select random bootstrapped samples for uncertainty estimation
     of Potts model results
 * LMDCA - use likelihood maximization to optimize Potts model
 * PLMDCA - use pseudolikelihood maximization to optimize Potts model
 * FinishDCA - calculate the results using LMDCA or PLMDCA fitted Potts model

For more advanced users, it may be desirable to modify or customize alf using
lower level routines. One may begin by making modifications to existing routines
like runflat and runprod, or one may write new routines from scratch. These
lower level routines perform basic functions important in the basic routines.

Lower Level Routines:
 * GetLambdas WORKING HERE
"""

from alf.SetVars import *
from alf.GetLambdas import *
from alf.GetVolumes import *
from alf.GetEnergy import *
from alf.GetFreeEnergy5 import *
from alf.wham.RunWham import *
from alf.GetSteps import *
from alf.GetVariance import *
from alf.dca.GetVariance import *
from alf.dca.GetModel import *
from alf.util.PlotFreeEnergy import *
from alf.util.GetDCD import *

from alf.initialize import *
from alf.runflat import *
from alf.runprod import *
from alf.postprocess import *
from alf.postvolume import *
from alf.postprocessDCA import *
from alf.utilities import *
