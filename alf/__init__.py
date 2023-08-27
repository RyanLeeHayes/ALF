"""adaptive landscape flattening for lambda dynamics simulations

Routines
--------

There are several levels of routines within alf (the module). They can
be roughly divided into four categories: basic routines, which run alf
out of the box with default settings, Potts model routines, which can
give improved free energy estimates for large numbers of mutating sites,
utility routines, which allow users to examine their results more
closely, and low level routines, which are the routines used by the
basic routines to implement ALF (the method).

Basic Routines
--------------

The first group of routiness provides basic functionality for running
lambda dynamics simulations with ALF. Examples of how to use them are
highlighted in various examples/engines subdirectories. Specifically,
these routines are called by various bash and slurm scripts in those
directories. These routines all take a named argument 'engine' that
specifies the molecular dynamics engine in use. Currently supported MD
engines are 'charmm' (CHARMM with DOMDEC GPU ONLY), 'bladelib' (CHARMM
with BLaDE), and 'blade' (a standalone implementation of BLaDE). The
routines also assume the existence of a directory called 'prep'. The
prep directory should contain a file called alf_info.py, which creates a
python dictionary also called alf_info containing information ALF needs
to know about the system, a file called name.inp, where name is the name
of the system given by a string in alf_info.py, which is a script to set
up the system in the selected engine, a file called nbond.str that sets
up nonbonded interaction parameters and cutoffs for the system, and any
supporting files required by name.inp. Further specifications for the
prep directory can be found in the alf README.md file. Example prep
directories are available in examples/systems 

initialize : sets up the analysis0 directory, variables1.inp file, and
    other objects that must be created before beginning ALF
runflat : runs several iterative short simulations to get rough
    approximations for the biases; arguments control how many cycles are
    used and how long these simulations are. Common practice is to run
    many 100 ps simulations to get close, and then further refine biases
    with several more 1 ns simulations
runprod : used for running longer production simulations, usually with
    several indepdendent trials
postprocess : analyzes production simulations from runprod to further
    refine biases for additional production simulations and/or produce a
    Results.txt file containing free energy estimates for dG using the
    histogram-based estimator. Relative ddG may be obtained by comparing
    these results with dG for the same perturbations in another physical
    ensemble

Potts Model Routines
--------------------

Additional routines are provided to estimate free energies using the
Potts model estimator rather than the histogram-based estimator. This is
a more complicated estimator useful primarily for systems with many
perturbation sites.

SetupDCA : copy some relevant files into the dca directory where the
    analysis is performed. Requires postprocess to have been performed
    on the same simulations first
FilterDCA : discretize the continuous lambda trajectories to determine
    which (if any) substituent is dominant at each site at each time
MomentDCA : calculate moments, or single site probabilities and joint
    probabilities at pairs of sites that are used for likelihood
    optimization
BSMomentDCA : select random bootstrapped samples for uncertainty
    estimation of Potts model results
LMDCA : use likelihood maximization to optimize Potts model
PLMDCA : use pseudolikelihood maximization to optimize Potts model
FinishDCA : calculate the results using LMDCA or PLMDCA fitted Potts
    model

Utility Routines
----------------

The utility routine are useful for inspecting the output of the lambda
dynamics simulations to visualize the free energy profiles as a function
of various chemical coordinates to determine the success of flattening,
to collect simulation trajectories to check for relevant spatial
configuration changes, to check transition rates between chemical end
states, and so on. These routines are located in alf/util, and the more
frequently used routines have been polished and included in the alf
module. Those polished routines are

PlotFreeEnergy5 : plot the free energy profiles as a function of lambda
    from a particular analysis directory. If ALF was successful, these
    profiles should be flat to within 1-2 kcal/mol with no large gaps of
    unsampled lambda values
GetDCD : concatenate the many DCD files from a long production run into
    a single DCD file for each independent trial
GetTrans : calculate alchemical transitions at each site, useful for
    assessing convergence of production simulations

Lower Level Routines
--------------------

For more advanced users, it may be desirable to modify or customize alf
using lower level routines. One may begin by making modifications to
existing routines like runflat and runprod, or one may write new
routines from scratch. These lower level routines perform basic
functions important in the basic routines.

initialize_alf_info : read prep/alf_info.py
InitVars : set starting values for analysis0 and variables1.inp
GetLambdas : for run i, convert binary lambda trajectories in run[i]/res
    to human readable lambda trajectories in analysisi/data. Names of
    files are assumed to be those used in alf/default_scripts. Wraps the
    routine GetLambda to create directories and pass it default
    filenames
GetLambda : the function called by GetLambdas to convert and concatenate
    one or more binary lambda trajectories whose filenames are in a list
    to a single human readable lambda trajectory whose filename is also
    passed in
GetEnergy : determine the alchemical trajectories to be used in
    flattening. In short flattening runs by runflat this includes the
    last five simulations (to give a broader view of alchemical space
    than just the last simulation can provide). In longer production
    runs by runprod, this includes all independent trials. The
    alchemical trajectories are copied into analysisi/Lambda, and the
    energy of each of these simulations with respect to the bias
    potential of each of the others is computed in analysisi/Energy.
    These energies are needed by WHAM or MBAR to combine data from all
    simulations into a single model of the alchemical free energy
    landscape
RunWham : runs the cuda code in alf/wham to combine the multiple
    simulations into a single model of the free energy landscape,
    computes free energy profiles as a function of many reaction
    coordinates. Computes the derivative of each of these profiles with
    respect to changes in each of the bias potential parameters, and
    outputs a C.dat matrix and V.dat vector corresponding to the
    quadratic and linear terms of a linearized objective function that
    penalizes deviations from flat free energy landscapes
GetFreeEnergy5 : inverts the C.dat matrix and multiplies by V.dat to
    minimize the linearized objective function. Saves the changes in the
    b (phi), c (psi), x (chi), and s (omega) coefficients into files
    like analysisi/b.dat
SetVars : adds the old bias values like analysisi/b_prev.dat and the
    changes in the bias values like analysisi/b.dat to give the new bias
    values like analysisi/b_sum.dat, and saves those bias values to
    variablesi+1.inp to be read in by the molecular dynamics engine in
    the next simulation
GetVolumes : compute charge change corrections if a key 'q' exists in
    alf_info. The simulation volumes are read from default filenames in
    run[i]/dcd to compute the discrete solvent correction. Wraps the
    routine GetVolume which does not assume a particular filename
    convention
GetVariance : compute free energies of various chemical states using the
    histogram-based estimator, use independent trials to produce
    bootstrap estimates of statistical uncertainty, and place results in
    analysisi/Result.txt
"""

from alf.SetVars import *
from alf.GetLambdas import *
from alf.GetVolumes import *
from alf.GetEnergy import *
from alf.GetFreeEnergy5 import *
from alf.GetFreeEnergyLM import * # undocumented
from alf.wham.RunWham import *
from alf.GetSteps import *
from alf.GetVariance import *
from alf.dca.GetVariance import *
from alf.dca.GetModel import *
from alf.dca.BootstrapMoments import *
from alf.util.PlotFreeEnergy import *
from alf.util.GetDCD import *
from alf.util.GetTrans import *

from alf.initialize import *
from alf.runflat import *
from alf.runprod import *
from alf.postprocess import *
from alf.postprocessDCA import *
from alf.utilities import *
