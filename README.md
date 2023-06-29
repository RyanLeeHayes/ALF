Installation instructions

1. Go into alf/wham and alf/dca and edit the modules files in each to set up the appropriate compilation environment on your machine, including pointing to CUDA if available. Then run Compile.sh in each to compile necessary C and CUDA libraries.

2. Create a virtual environment, and then run pip install to install ALF in that virtual environment. If you follow the script in Setup.sh, it will create a script setupenv which you can source to reload that virtual environment in the future.

Test cases

There is a test case for CHARMM with domdec (RNaseH_u27-33_blade) and a test case for standalone BLaDE (T4L149U_charmm) in test_beta.

T4L149U_charmm - this runs the unfolded side of point mutations to site 149 in T4 lysozyme. Example scripts of how to use ALF are provided. Edit the first few lines of runset2.sh, runset4.sh, and postprocess.sh to set up the environment to run ALF and CHARMM. In particular, make sure CHARMMEXE is set, and that the setupenv script created during installation are sourced. These scripts can then be submitted to slurm with scripts like subset2.sh, subset4.sh, and subsetP.sh. runset2.sh performs flattening, which optimizes biasing potentials before production simulations. This is almost always necessary to get any results at all. After flattening, runset4.sh runs production simulations. Finally, postprocess.sh evaluates free energy using a histogram based estimator.

-------------------------------------OLD----------------------------------------

The contents of this directory WERE as follows

ALF
Contains python scripts used for flattening. There are a few other useful scripts as well, such as PlotFreeEnergy5.m, which can be copied into an analysis directory and run in matlab to visualize the alchemical free energy profiles (to see if they are satisfactorily flat), and other scripts such as GetTrans.py in util which will count transitions in a trajectory.

analysis0
Each step of flattening produces a "run#" directory, where # is a consecutive integer starting from 1, which contains the output of the CHARMM simulation, and an "analysis#" directory which analyzes that simulation and places modified bias variables in a "variables[#+1].inp" file. "analysis0" contains the starting values of the biases in the files analysis0/b_sum.dat, analysis0/c_sum.dat, analysis0/x_sum.dat, and analysis0/s_sum.dat, which are used as input to the first run. Running "InitVars.py" in this directory will initialize these files with 0's, and create the file variables1.inp.

charmm.sh
This file is sourced by runset2.sh, runset3.sh, and runset4.sh, and loads the charmm environment. It should set the environment variable CHARMMEXEC to point to your CHARMM executable, and load any requisite modules.

dWHAMdV_mso
This directory contains the core of the flattening algorithm, which is written in CUDA to make sure it runs efficiently. You will need to recompile it for your machine. To do so, first execute Clean.sh to remove old cmake garbage, edit modules to load an appropriate cmake, c compiler, and cuda module, and then run Compile.sh.

msld_flat.inp
This is the CHARMM input script used during flattening. On line 20, it streams in a file from the prep directory that does all the system setup, everything else the script does should be system independent. It runs two segments of NPT dynamics, one of @esteps which is discarded for equilibration, and one of @nsteps that is used for flattening.

msld_prod.inp
This is the CHARMM input script used during production. It also streams in a system setup script on line 20. This script is set up to be called multiple times, running 1 ns of sampling each time. The 1 ns segments are then stitched together during postprocessing. This makes the production run much more robust against queue time limits and node failures.

nbond.str
This is a stream file that sets up the nonbonded parameters

ntersiteflat
This controls the treatment of intersite coupling during flattening. "0 1" is default, and means that coupling between sites is not adjusted (0), but coupling free energy profiles are considered (1) during flattening.

ntersiteprod
Same thing as ntersiteflat, except for it's used if a production run is used to reoptimize the biases. "1 1" is the default, because production runs are generally long enough that reliable estimates of the coupling parameters can be made.

postprocess.sh
This is a script to do rudimentary postprocessing of a production run. It also reoptimizes the biases based on that run. If the production run was bad, running with reoptimized biases may help. It requires 5 input environment variables that are set in subsetP.sh: "i" which is the step number of the run, "eqS" which is the number of nanoseconds to discard for equilibration, "S" which is the total number of nanoseconds run (including equilibration), "N" which is the number of independent trials (each trial will be in a directory run#L where # is the step number and L is a lowercase letter indicating the trial), and "skipE" which if set to larger that 1, will only consider frames of that modulus in the analysis (used for culling down the amount of data during very large production runs so you don't run out of memory.)

README
You've already figured out what this one is for.

runset2.sh
Run the first phase of flattening with 25 ps equilibration and 75 ps production. There used to be a runset1.sh that made initial guesses for the biases, but it didn't save enough tim to be worth the extra effort. Line 4 loads a python module. You should load an appropriate module for your cluster. "ini", "iri", and "ifi" are set next and control the flattening run. ini is the first step of flattening, set it to 1. Each run chooses a random restart file from a previous run, and uses it to start the current run. "iri" sets the earliest run that can be used for this purpose. Set it to 1. "ifi" is the final step for flattening and controls how many runs will be used in phase 2. 50 is often sufficient, but I generally set it to 100 to be safe. Arginine is the worst, and typically requires "ifi" to be set to 200 to get close to flat. The purpose of this phase, phase 2, is to get the biases close to the correct values. Usually the starting guesses of 0 are off by dozens of kcal/mol, and you want to run very short simulations to get them close to the correct value so you can feed them in to phase 3. Running PlotFreeEnergy5.m on the final analysis directory is a good way to check if the profiles are flat. Generally you want to see sampling all the way accross the profiles with no large unsampled gaps. It doesn't matter if they're not perfectly flat, you just want everything to be thermodynamically accessible.

runset3.sh
These runs use 0.25 ns of equilibration and 0.75 ns of production. They are longer and more expensive, so you only want to use them once you are close to the correct biases after phase 2, otherwise you just end up doing work you could have done in phase 2 at 10 times the computational cost. These runs refine the biases. Many are tempted to skip phase 3, thinking 10 0.1 ns runs are the same as 1 1.0 ns run, but this is not the case. The extra time allows the system to equilibrate and relax further from the initial structure, giving a much better estimate of the true free energy landscape. 10 steps on phase 3 are generally sufficient unless it is a very complex system with multiple sites or many substituents (7 or more - note: we haven't yet successfully exceeded 9 substituents at a site, so consider that an upper bound). "ini" should be set to "ifi" of phase 2 plus one, "iri" is generally set to "ini" minus 5, and "ifi" is set to "ini" plus 9 for 10 cycles.

runset4.sh
Do a production run. The relevant environment variables are set by subset4.sh. Each of these runset scripts is queued to slurm by a corresponding subset script. If you don't use slurm, we'll probably have to rework these files. Each of these files keeps relaunching CHARMM if it fails. Your signal that is happening is if the run# file gets moved to run#_failed during flattening, or if the run#/output file gets moved to run#/failed/output during production. After flattening is complete, you may also want to run short production runs to further improve the biases. For example, if I intend to run a 100 ns production run, a 1 ns flattening run generally doesn't do a good enough job of optimizing the bias, so I'll run a 5 ns production run, then a 20 ns production run, and finally the full 100 ns production run, working my way up by factors of 4 or 5.

subset2.sh
Run this script to submit runset2.sh to the queue.

subset3.sh
Run this script to submit runset3.sh to the queue.

subset4.sh
Run this script to submit runset4.sh to the queue. The for loop over p launches each of the independent trial with a different letter, I generally use 5. "ini" should be set to 1 plus the previous run number, --array=1-40%1 creates a slurm job array of 40 jobs, but only lets one run at a time, each job will run a consecutive time slice, though if one fails, the next one will go back and clean up after it. "nitt" controls the integer number of ns run per job array element, so the total number of ns run in nitt times the job array length.

subsetAll.sh
This queues runset2.sh, runset3.sh, and runset4.sh all in one command. It's probably better to launch them individually until you understand what's going on.

subsetP.sh
Queues postprocess.sh.

systems
This contains example system setups: T4L149F a folded state setup with side chain perturbations, T4L149U an unfolded state setup with side chain perturbations, u58-64 a heptapeptide with whole residue perturbations, u58-64_bprex u58-64 with biasing potential replica exchange, and u58-64_vbrex u58-64 with variable bias replica exchange.



Here are the contents of systems/T4L149F. They should be copied (recursively) into this directory to run that system. Analogous files will be needed for whatever system you decide to run.

analysis0
Same as analysis0 above, except it has been initialized using InitVars.py. Note InitVars.py depends on the file "nblocks" below.

name
The (case sensitive) name of the system. msld_flat.inp and msld_prod.inp will look for prep/name.inp to stream to set up the system.

nblocks
The total number of blocks in the system (minus the environment). This should be the sum of "nsubs". If you are mutating between three residues at one site and 5 at another site, nsubs is 8.

nbond.str
This streams in the nonbonded parameters. It's system independent, except fftx, ffty, and fftz have to be set to highly factorable integers near the box length.

nbshift
This directory contains matrices for replica exchange. Unless you're using replica exchange, they should all be set to 0's. InitVars.py initializes the directory and sets it to 0's by default.

ncentral
Another replica exchange feature. This is the "central" potential, i.e. the replica you want to flatten. Without replica exchange it should be 0.

nnodes
Determines the number of GPUs to use. Unless your system is over 150 A, it generally runs fastest with just one GPU.

nreps
The number of replicas to use. Without replica exchange, set to 1.

nsubs
The number of substituents at each site. If mutating between three residues at the first site, and five residues at the next site, nsubs should be "3 5".

prep
This file has all the system specific files, including name.inp, and minimized.pdb or minimized.psf, the starting structure. In this case, name is T4L149F, and T4L149F.inp sets up the system. Most of the MSLD stuff happens down in the BLOCK command. This script was designed to be fairly flexible. By changing the first few lines, you can change the mutation site, or add additional mutation sites. Currently it says
set box = 71.6297564
set s1seq1 = nat
set s1seq2 = cys
set s1seq3 = ile
set s1seq4 = met
set s1seq5 = ser
set s1seq6 = thr
set segid = PROT
set resid1 = 149
but to also make mutations to tyrosine at F153, you could add the lines
set s2seq1 = nat
set s2seq2 = tyr
set resid2 = 153
(and of course modify nsubs and nblocks appropriately, and rerun InitVars.py). You should also be able to switch to an entirely different protein system. You may need to change the segids and add a few protonation, disulfide, and capping patches as appropriate, but the setup should be largely similar.

Target.txt
This is not part of the input, this is the approximate answer I get out of the final analysis directory when I run this system. Your result should be similar within statistical error. You should see a "Result.txt" file in analysis# (where # is a production run only), that was generated by ALF/GetVariance.py. During flattening runs analysis#/b_sum.dat should move toward Target.txt.

variables1.inp
A stream file that tells the first CHARMM run what to do. The can also be generated with analysis0/InitVars.py.
