# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ALF (Adaptive Landscape Flattening) is a Python package (v3.1.3) that optimizes bias potentials for lambda dynamics molecular dynamics simulations. It iteratively runs MD simulations and flattens free energy landscapes using WHAM/MBAR reweighting and bias optimization. The package combines Python orchestration with GPU-accelerated CUDA/C++, MPI C, and Fortran executables.

## Build Commands

### Compile native executables (requires CUDA, MPI, OpenMP, Fortran, CMake 3.20+)
```bash
cd alf/bin
bash Clean.sh    # remove previous build artifacts
cmake ./
make all
```

### Install Python package
```bash
python -m venv env-alf
source env-alf/bin/activate
pip install -e .
```

Or use `bash Setup.sh` from the repo root which does both venv creation and pip install.

### Verify installation
```bash
python -c "import alf"
```

There is no test suite — verification is done by running examples from `examples/`.

## Architecture

### Language Mix
- **Python** (`alf/*.py`): Orchestration, workflow drivers, analysis logic
- **CUDA/C++** (`alf/bin/loss/*.cu`): GPU-accelerated WHAM, loss functions (linear2018, nonlinear2024, whamweight, impcons)
- **C with MPI** (`alf/bin/dca/*.c`): Potts model DCA executables (Filter, Moment, LM)
- **CUDA/C++** (`alf/bin/dca/PLM.cu`): GPU Potts pseudolikelihood maximization
- **Fortran 90** (`alf/bin/io/*.F90`): Binary lambda trajectory readers (GetLambda, GetSteps)
- **CHARMM scripts** (`alf/default_scripts/*.inp`): MD engine input templates

### Python-to-Native Bridge
`liblinear2018.so` is a shared library loaded at runtime via `ctypes.CDLL` (see `alf/loss.py`). All other native executables are invoked via `subprocess.call`.

### Iterative Cycle Model
ALF operates in numbered cycles `[i]`. Each cycle:
1. Runs MD in `run[i]/` using bias parameters from `variables[i].inp`
2. Analyzes results in `analysis[i]/` (lambda trajectories, WHAM reweighting, bias optimization)
3. Writes updated biases to `analysis[i]/b_sum.dat` and `variables[i+1].inp`

Failed cycles get renamed `run[i]_failed/` and retried. Completed cycles (detected by `analysis[i]/b_sum.dat` existing) are skipped.

### Core Workflow (4 basic routines)
- `alf.initialize()` — creates `analysis0/` and `variables1.inp`
- `alf.runflat()` — short iterative flattening cycles (100ps then 1ns)
- `alf.runprod()` — longer production runs with independent trials
- `alf.postprocess()` — bias refinement and free energy estimation → `Result.txt`

### alf_info Dictionary
Central configuration dict loaded from `prep/alf_info.py` by `utilities.initialize_alf_info()`. Required keys: `nsubs`, `nblocks`, `name`, `temp`, `enginepath`, `nreps`, `ncentral`, `nnodes`. Optional: `q` (charges), `loss`, `bias`, `impcons`, `impconsopt`.

### Engine Abstraction
The `engine` parameter (`'charmm'`, `'bladelib'`, `'blade'`, `'pycharmm'`) controls subprocess invocation, `variables[i].inp` format (via `SetVarsCharmm`/`SetVarsBlade`/`SetVarsPycharmm`), and which default script template is copied.

### Two Analysis Paths
- **linear2018**: WHAM → C.dat (Hessian) + V.dat (gradient) → matrix inversion via `GetFreeEnergy5`
- **nonlinear2024**: WHAM weights → nonlinear iterative optimization via standalone executables

### Potts Model Extension (`alf/dca/`, `alf/bin/dca/`)
For large multi-site systems (>100 chemical states): `SetupDCA` → `FilterDCA` → `MomentDCA` → `LMDCA`/`PLMDCA` → `FinishDCA`. Uses compiled C/CUDA executables.

### Bias Potential Types
- `bcxs2018`: b (linear/phi), c (quadratic/psi), x (skew/chi), s (endpoint/omega); alpha=0.017
- `bcxstu2026`: adds t and u endpoint variants; alpha=0.012

## Key Files

| File | Role |
|------|------|
| `alf/runflat.py` | Flattening loop driver — orchestrates cycles of MD + analysis |
| `alf/runprod.py` | Production run driver |
| `alf/postprocess.py` | Production analysis and bias refinement |
| `alf/loss.py` | Python interface to CUDA loss functions (ctypes for linear2018, subprocess for others) |
| `alf/SetVars.py` | Writes `variables[i].inp` files per engine format |
| `alf/GetEnergy.py` | Computes bias energies for WHAM reweighting |
| `alf/GetFreeEnergy5.py` | Matrix inversion for linear bias optimization |
| `alf/utilities.py` | `initialize_alf_info()` — validates and defaults the alf_info dict |
| `alf/bin/CMakeLists.txt` | Single CMake file building all native targets |

## Dependencies
Python: `numpy`, `scipy`, `MDAnalysis`. Build: CMake 3.20+, CUDA 9.0+, C/C++, Fortran, MPI, OpenMP.
