Package Description

alf (the module) runs adaptive landscape flattening (ALF - the method)
to optimize bias potentials in a lambda dynamics molecular dynamics
simulation. These optimized bias potentials significantly improve the
efficiency of lambda dynamics, and are required as a starting point for
most modern lambda dynamics simulations. This README contains
information on installation, examples, algorithms, best practices,
python routines, input specification, bug reporting, and citations.



Installation

Detailed installation instructions are available in INSTALL. Briefly,
you will need to compile code in alf/wham and alf/dca, and then use pip
install to install the alf module into a python virtual environment.



Examples

Detailed instructions for how to run examples are available in
examples/INSTRUCTIONS. Briefly, you should copy a directory with scripts
for your particular engine from examples/engines, and then copy a
directory from examples/systems into that directory and rename it prep.



Algorithms

ALF runs cycles of molecular dynamics followed by flattening to
itteratively improve biasing potentials.

Molecular dynamics in the alf software is performed by an external
molecular dynamics engine which must be compiled and installed
independently. The currently supported engines are:
 * charmm : the CHARMM software package utilizing the DOMDEC GPU ONLY
   command for GPU acceleration
 * bladelib : the CHARMM software package utilizing the BLaDE library
   for GPU acceleration. BLaDE is faster that DOMDEC, but has fewer
   features
 * blade : the standalone BLaDE software package
 * pycharmm : the python pyCHARMM package using CHARMM as a library to
   call BLaDE
These engines may be passed to alf routines to specify the engine,
because some routines have slight differences based on the engine used.
Use of other molecular dynamics engines in CHARMM is possible, including
engines that do not use GPUs, see instructions at the end of INSTALL.
Use of other engines will require more extensive modification of alf.

Sampling will be optimal when the free energy landscape is flat, so
flattening seeks to identify the biases that give flat free energy
profiles along several 1D and 2D reaction coordinate projections. Lambda
trajectories are first copied from the format output by the engine, then
energies of several recent trajectories are computed using all
combinations of biases so that a single model of the energy landscape
may be produced by reweighting samples with WHAM/MBAR. This reweighted
model is typically more robust because it includes information from
several trajectories run at different biases. Free energy profiles are
computed from this reweighted ensemble as a function of each lambda,
transitioning from one lambda to another, and as a function of two
lambdas. A linear model of how changes in any bias parameter will change
any profile is computed, and from this linear model a quadratic penalty
function is computed that penalizes deviations of the free energy
profiles from their average values. This penalty function is optimized
to give changes in each of the bias parameters, which are then used in
the next round of molecular dynamics.



Best Practices

An ideal workflow is shown in the file subsetAll.sh within the
examples/engines directory.

ALF goes through many cycles of sampling the system and then estimating
improved bias parameters. Cycles are indexed by integer. In a particular
cycle, for example 3, simulations will be run and lambda trajectories
recorded in run3, then analysis will be performed in analysis3. The
previous values of the biases b (phi), c (psi), x (chi), and s (omega)
are copied from analysis2/b_sum.dat into analysis3/b_prev.dat, the
estimated changes are saved to analysis3/b.dat, and updated values are
saved to analysis3/b_sum.dat and copied into a format readable by the
engine for the next run in variables4.inp.

The first step is to run alf.initialize to set up the analysis0
directory and variables1.inp required by the first run.

The second step is to begin running flattening with alf.runflat with
many short 100 ps runs. The purpose of this step is to get the biases
close to the correct values. Usually the starting guesses of 0 generated
by alf.initialize are off by dozens to hundreds of kcal/mol, and very
short simulations are desirable to get them close to the correct value
with minimal computational effort. For very bad biases the direction of
the change is quite obvious, but one should not change biases by more
than a few kcal/mol because only the bottom few kcal/mol of basins can
be sampled and the shape of the free energy profiles beyond this is
unknown. Consequently, many short cycles of flattening are run for 100
ps with alf.runflat. Discarding the first quarter of alchemical
simulations (in this case 25 ps) for equilibration typically gives the
best results. The number of cycles required varies by system, and while
a halting procedure has been proposed, it has not been implemented in
this package. These short simulations are cheap, so running extra is not
a significant waste. 50 cycles is often sufficient for small hydrophobic
groups, but generally 100 cycles is safer. Arginine perturbations
typically require 200 cycles to get close to flat. Convergence can be
assessed by running alf.PlotFreeEnergy5 on the final analysis directory
to check if the profiles are flat. Generally you want to see sampling
all the way accross the profiles with no large unsampled gaps. It
doesn't matter if they're not perfectly flat, but everything should be
thermodynamically accessible.

The third step is to run a few more cycles of flattening with
alf.runflat, but using somewhat longer simulations of 1 ns to refine the
biases. Again, the first quarter of the simulations is discarded for
equilibration. Typically 10 cycles is sufficient here.

The fourth step is to begin running production using alf.runprod
followed by analysis with alf.postprocess. Typically 5 duplicate
independent trials are run with alf.runprod in production using the same
biases to estimate statistical variability. Often as a simulation runs
longer, new parts of phase space are visited that shift the balance
between alchemical states. Consequently it is advisable to ramp up to
the final production simulation using rounds of shorter production and
postprocessing to refine the biases, typically by factors of 4-5. Thus
if 100 ns is required to sample a system well, it is advised to run 5
copies of 5 ns production and postprocessing, followed by 5 copies of 20
ns of production and postprocessing, followed by 5 copies of 100 ns of
production and postprocessing. For solvation free energies, 10-20 ns is
sufficient. For well behaved protein point mutations 40 ns is
sufficient, for poorly behaved protein point mutations and several
simultaneous protein mutations 100 ns is sufficient, for massive protein
chemical spaces 400 ns has been required, for easy ligand binding
perturbations 10 ns is sufficient, though for challenging ligands 30-60
ns can be advisable, and some massive ligand chemical spaces have
required more sampling.

alf.postprocess with create a file called Results.txt in the
corresponding analysis directory giving dG, the chemical free energy
change in that ensemble. To obtain the relative ddG, take the difference
between this dG and the dG in another ensemble determined with another
run of ALF.

Charge changing perturbations: charge changing perturbations introduce
significant errors into alchemical simulations. Some of these errors are
due to missing physics, e.g. polarizability, but some are due to factors
that can be accurately calculated. For initially neutral boxes using PME
electrostatics, the largest of these factors is the discrete solvent
correction. alf will estimate this correction and automatically include
it in Results.txt IF a key 'q' is placed in the alf_info dictionary with
a corresponding value that is a python list of the net charge of every
substituent or alchemical group.

Substituents per site: we haven't yet successfully exceeded 9
substituents at a site, because with more substituents, more time is
spent in irrelevant alchemical intermediates, so consider that an upper
bound. Try to limit it to six substituents per site. If you care about
more chemical space, try splitting between multiple systems with a
common reference substituent. If you need more substituents and can't
cut down the system, e.g. you care about more than 6 substituents and
their combinations at several sites, reach out to the developers. There
may be options under development that are not polished enough for public
distribution.

Number of sites and Potts model estimator: in principle lambda dynamics
can include as many sites as desired, but in practice, the histogram
based estimator used by alf.postprocess can't robustly make estimates
for more than 100-500 points in chemical space because every physical
point in chemical space must be visited. The Potts model estimator
assumes that free energies can be decomposed into single site terms that
depend only on the substituent identity at each site and two body terms
that account for coupling between all pairs of sites. Three body and
higher order terms are assumed to be zero. This approximation is
reasonable in many systems and significantly reduces sampling
requirements, making it possible to estimate free energies for tens of
thousands of sequences. To use the Potts model estimator, run
alf.postprocess first, then see the examples in examples/engines/PottsLM
and examples/engines/PottsPLM. PottsLM uses likelihood maximization,
best for systems with less than a million chemical end states, and
PottsPLM uses pseudolikelihood maximization, best for systems with more
than a million chemical end states. These use routines defined in
alf/dca.

Coupling between sites: if you have multiple sites, it is possible these
site may interact strongly, and efficiently sampling your system may
require adding biases that modulate the coupling between sites. The
"ntersite" named argument may be passed to alf.runflat and
alf.postprocess to modulate how alf handles coupling. By default, the
[0,0] value of ntersite ignores coupling. The first element of coupling
determines whether coupling biases are adjusted. A value of 0 means no
intersite biases are included, a value of 1 means c, x, and s intersite
biases are adjusted, and a value of 2 only adjusts c biases because this
significantly reduces the overhead of adaptive landscape flattening
calculations for large numbers of sites, and because intersite x and s
values tend to be very small. The second element controls whether 2D
free energy profiles binning lambda at two different sites are used
during flattening. 1 computes these profiles, 0 does not. The short
alf.runflat flattening runs typically do not sample enough alchemical
space to make robust estimates of the coupling biases, so [0,1] is
recommended during these initial flattening runs if significant coupling
is expected to be present. The longer runs during preliminary production
runs analyzed by alf.postprocess can make better estimates of coupling,
so [2,1] is recommended during these cycles if significant coupling is
present. Note this means that the coupling won't begin to be accounted
for until the second production simulation using biases from the first
is run.



Python Routines

Python routines and the alf module are documented using docstrings. To
read this documentation for the entire module, after installation run

python -c "import alf; help(alf)"

or to read it for a routine like alf.runflat run

python -c "import alf; help(alf.runflat)"

You may also find this documentation for the module in alf/__init__.py
and where the routine is defined in alf/runflat.py, respectively.



Input Specification

alf routines assume the existence of several files within a directory
called prep. This directory contains all the system specific infomation
needed by alf to run lambda dynamics on a specific system. A soft link
or symbolic link will be created inside every single run directory so
that the engine can access it from that directory.


alf_info.py: this file sets up a python dictionary called alf_info for
alf and should include the following keys:
     alf_info={} # initialize dictionary
 * nsubs : a list of the number of alchemical groups at each site. Thus
     if you are considering 6 alternative functional groups on a small
     molecule, nsubs is [6]. If you are considering 3 modifications at
     one site and 7 at a second site, nsubs is [3,7]. If you are
     considering the free energies of mutating a protein from its native
     sequence of methionine to lysine or leucine, nsubs is [3] for the
     native and two mutations. If you are considering alternating
     between the native sequence and one candidate mutation
     combinatorially at 15 sites, nsubs is
     [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2].
 * nblocks : the total number of lambda variables, or the sum of nsubs
     alf_info['nblocks']=np.sum(alf_info['nsubs'])
 * name : the (case sensitive) name of the system. msld_flat.inp and
     msld_prod.inp will look for prep/name.inp to stream to set up the
     system.
 * temp : the temperature to run the system at in Kelvin. To use room
     temperature:
     alf_info['temp']=298.15
 * enginepath : a string with the full absolute path to the lambda
     dynamics engine executable
 * nreps : the number of replicas to use. Without replica exchange, set
     to 1.
     alf_info['nreps']=1
 * ncentral : another replica exchange feature. This is the "central"
     potential, i.e. the replica you want to flatten. Without replica
     exchange it should be 0.
     alf_info['ncentral']=0
 * nnodes : determines the number of GPUs to use. Gains for using
     multiple GPUs are generally modest so unless time is of the essence
     and GPUs are abundant, one GPU generally gives the best throughput.
     alf_info['nnodes']=1
 * q : a list of length nblocks containing the net charge of every
     alchemical group. This is used to compute the discrete solvent
     charge changing correction. If missing, no charge changing
     correction will be applied.

name.inp: This file should include everything from loading force field
parameters to setting up the system, initial conditions, periodic
boundary conditions, and alchemical regions. Any supporting files
required should also be placed in and read from the prep directory.
Setting nonbonded options and running dynamics are controlled by the
scripts (msld_flat.inp and msld_prod.inp) that read name.inp. This file
must also include code to add the biases determined by alf to the
potential, otherwise alf will keep adjusting the biases, but sampling
will not improve. The following code will read the biases set by
alf.SetVars in the BLOCK section of charmm and bladelib engine input
files.

```
calc nbiaspot = 5 * ( @nblocks * ( @nblocks - 1 ) ) / 2
ldbi @nbiaspot
set ibias = 1
set iblock = 0
set si = 1
label loop5
if @si .le. @nsites then
   set jblock = @iblock
   set sj = @si
   label loop5b
   if @sj .le. @nsites then
      set ii = 1
      label loop6
      if @ii .le. @nsubs@@{si} then
         calc ip1 = @ii + 1 + @iblock
         set jj = 1
         if @si .eq. @sj then
            calc jj = @ii + 1
         endif
         label loop7
         if @jj .le. @nsubs@@{sj} then
            calc jp1 = @jj + 1 + @jblock
            if @si .eq. @sj then
               adex @ip1 @jp1
            endif

            ldbv @ibias @ip1 @jp1 6 0.0 @cs@@{si}s@@{ii}s@@{sj}s@@{jj} 0
            calc ibias = @ibias + 1
            ldbv @ibias @ip1 @jp1 10 -5.56 @xs@@{si}s@@{ii}s@@{sj}s@@{jj} 0
            calc ibias = @ibias + 1
            ldbv @ibias @ip1 @jp1 8 0.017 @ss@@{si}s@@{ii}s@@{sj}s@@{jj} 0
            calc ibias = @ibias + 1
            ldbv @ibias @jp1 @ip1 10 -5.56 @xs@@{sj}s@@{jj}s@@{si}s@@{ii} 0
            calc ibias = @ibias + 1
            ldbv @ibias @jp1 @ip1 8 0.017 @ss@@{sj}s@@{jj}s@@{si}s@@{ii} 0
            calc ibias = @ibias + 1
            calc jj = @jj + 1
            goto loop7
         endif
         calc ii = @ii + 1
         goto loop6
      endif
      calc jblock = @jblock + @nsubs@@{sj}
      calc sj = @sj + 1
      goto loop5b
   endif
   calc iblock = @iblock + @nsubs@@{si}
   calc si = @si + 1
   goto loop5
endif
```

nbond.str: this is a stream file in prep that sets up the nonbonded
parameters. PME is recommended for electrostatics with an interpolation
order of 6, an easily factorable grid length close to the box size in
angstroms, and an inverse beta of 0.32 A^-1. Van der Waals force
switching ("vfswitch") is recommended for Lennard-Jones interactions
with a switching radius ("ctonnb") of 9 A, a cutoff radius ("ctofnb") of
10 A, and a neighbor list radius ("cutnb") of 12 A.

Other necessary files can be supplied by the user outside the prep
directory if desired, but default values of these files will be created
if they are not supplied.

msld_flat.inp: This is the CHARMM input script used during flattening.
It streams in a file from the prep directory that does all the system
setup, everything else the script does should be system independent. It
runs two segments of NPT dynamics, one of @esteps which is discarded for
equilibration, and one of @nsteps that is used for flattening. If you
want to run NVT dynamics, for example for the vacuum side of a solvation
free energy calculation, make a copy and modify appropriately. This file
is copied from alf/default_scripts for the requested engine.

msld_prod.inp: This is the CHARMM input script used during production.
It also streams in a system setup script. This script is set up to be
called multiple times, running a 1 ns chunk of sampling each time. The 1
ns chunks are then stitched together during postprocessing. This makes
the production run much more robust against queue time limits and node
failures. If you want to run NVT dynamics, for example for the vacuum
side of a solvation free energy calculation, make a copy and modify
appropriately. This file is copied from alf/default_scripts for the
requested engine.

nbshift: This directory contains matrices for replica exchange. Unless
you're using replica exchange, they should all be set to 0's.
InitVars.py initializes the directory and sets it to 0's by default. For
variable biasing potential replica exchange, it also includes a file
vb.inp that is streamed after replica exchange initialization to adjust
the biasing potentials based of the replica index.

G_imp: This is a directory containing the free energy due to the entropy
of the implicit constraints. For a typical implicit constraint fnex
parameter of 5.5, and 10 or fewer substituents, the alf/G_imp directory
containing precomputed values will be used, but if you choose to use a
higher fnex or add a small barrier between states to increase the
sampling of physical states, you should recompute these G_imp values,
and pass alf.runflat and alf.postprocess the path to this directory.



Bug Reporting

Please report any bugs or desired new features to the authors at
https://github.com/RyanLeeHayes/ALF/issues



Citations

Original ALF paper:
Hayes, R. L.; Armacost, K. A.; Vilseck, J. Z. & Brooks III, C. L.
Adaptive Landscape Flattening Accelerates Sampling of Alchemical Space in Multisite λ Dynamics
Journal of Physical Chemistry B, 2017, 121, 3626-3635
DOI: 10.1021/acs.jpcb.6b09656

Current linearized ALF:
Hayes, R. L.; Vilseck, J. Z. & Brooks III, C. L.
Approaching Protein Design with Multisite λ Dynamics: Accurate and Scalable Mutational Folding Free Energies in T4 Lysozyme
Protein Science, 2018, 27, 1910-1922
DOI: 10.1002/pro.3500

Coupling between sites and Potts model estimator:
Hayes, R. L.; Vilseck, J. Z. & Brooks III, C. L.
Addressing Intersite Coupling Unlocks Large Combinatorial Chemical Spaces for Alchemical Free Energy Methods
Journal of Chemical Theory and Computation, 2022, 18, 2114-2123
DOI: 10.1021/acs.jctc.1c00948

Automatic halting procedure:
Raman, E. P.; Paul, T. J.; Hayes, R. L. & Brooks III, C. L.
Automated, Accurate, and Scalable Relative Protein-Ligand Binding Free Energy Calculations using Lambda Dynamics
Journal of Chemical Theory and Computation, 2020, 16, 7895-7914
DOI: 10.1021/acs.jctc.0c00830

Biasing potential replica exchange:
Biasing Potential Replica Exchange Multisite λ-Dynamics for Efficient Free Energy Calculations
Journal of Chemical Theory and Computation, 2015, 11, 1267-1277
DOI: 10.1021/ct500894k
