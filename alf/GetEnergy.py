#! /usr/bin/env python

def GetEnergy(alf_info,Fi,Ff,skipE=1):
  """
  Compute energies for wham/mbar reweighting

  This routine selects several simulations for reweighting using WHAM/MBAR
  It computes the relative bias energy from each simulation along the 
  alchemical trajectory. Different replicas in replica exchange are considered 
  separate simulations. 
  
  This routine should be run inside `analysis[Ff]`. . Biases are taken from the 
  appropriate `analysis[i]` directories, and alchemical trajectories are taken 
  from the output of the routine `alf.GetLambdas`, located in `analysis[i]/data`. 
  The alchemical trajectories from `analysis[Fi]/data` to `analysis[Ff]/data` 
  are copied into `analysis[Ff]/Lambda/Lambda[*].dat`, and the computed bias 
  energies are saved in `analysis[Ff]/Energy/ESim[*].dat`. 

  Energies are computed relative to the central replica bias in the `[Ff]` cycle.

  This routine can be called during flattening or production:
    - During **flattening**, `Ff - Fi = 4` is recommended (including 5 cycles in analysis).
    - During **production**, `Ff - Fi = 0` is recommended (including 1 cycle in analysis).

  The function detects flattening vs. production based on the presence of the 
    `run[Ff]` directory. In production, this directory does not exist; instead, 
    `run[Ff]a` and other independent trial directories exist.

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run
        - `nblocks`: Number of blocks.
        - `nsubs`: List specifying the number of sub-states for each site.
        - `nreps`: Number of replicas.
        - `ncentral`: Index of the central replica.
  Fi : int
      The first cycle of alf to include in analysis (inclusive)
  Ff : int
      The final cycle of alf to include in analysis (inclusive)
  skipE : int, optional
      In longer production runs the number of lambda samples may require
      significant amounts of memory to store and analyze. Only alchemical
      frames with index modulus skipE equal to skipE-1 are analyzed.
      - Default: `1` (analyze all frames).
  """

    import sys, os
    import numpy as np

    NF = Ff - Fi + 1

    nblocks = alf_info['nblocks']
    nsubs = alf_info['nsubs']
    nreps = alf_info['nreps']
    ncentral = alf_info['ncentral']

    try:
        b_shift = np.loadtxt('../nbshift/b_shift.dat')
        c_shift = np.loadtxt('../nbshift/c_shift.dat')
        x_shift = np.loadtxt('../nbshift/x_shift.dat')
        s_shift = np.loadtxt('../nbshift/s_shift.dat')
    except FileNotFoundError as e:
        print(f"Error loading shift files: {e}")
        return

    Lambda = []
    b = []
    c = []
    x = []
    s = []

    for i in range(NF):
        analysis_dir = f'../analysis{Fi+i}'
        data_dir = os.path.join(analysis_dir, 'data')
        
        if not os.path.isdir(data_dir):
            print(f"Warning: Directory {data_dir} not found")
            continue

        # print(f"Checking directory: {data_dir}")

        lambda_files = [f for f in os.listdir(data_dir) if f.startswith('Lambda.') and f.endswith('.dat')]
        lambda_files.sort()

        for lambda_file in lambda_files:
            file_path = os.path.join(data_dir, lambda_file)
            print(f"Loading file: {file_path}")
            try:
                Lambda.append(np.loadtxt(file_path)[(skipE-1)::skipE,:])
                
                # Extract j and k from the filename (assuming format Lambda.j.k.dat)
                j, k = map(int, lambda_file.split('.')[1:3])
                
                b.append(np.loadtxt(os.path.join(analysis_dir, 'b_prev.dat')) + b_shift * (k - ncentral))
                c.append(np.loadtxt(os.path.join(analysis_dir, 'c_prev.dat')) + c_shift * (k - ncentral))
                x.append(np.loadtxt(os.path.join(analysis_dir, 'x_prev.dat')) + x_shift * (k - ncentral))
                s.append(np.loadtxt(os.path.join(analysis_dir, 's_prev.dat')) + s_shift * (k - ncentral))
            except Exception as e:
                print(f"Error loading file {file_path}: {e}")

    if not os.path.isdir('Lambda'):
        os.mkdir('Lambda')
    if not os.path.isdir('Energy'):
        os.mkdir('Energy')

    total_simulations = len(Lambda)
    # print(f"Total simulations: {total_simulations}")

    if total_simulations == 0:
        print("Error: No Lambda files found.")
        return

    E = [[] for _ in range(total_simulations)]

    try:
        for i in range(total_simulations):
            for j in range(total_simulations):
                bi, ci, xi, si = b[i], c[i], x[i], s[i]
                Lj = Lambda[j]
                E[i].append(np.reshape(np.dot(Lj,-bi),(-1,1)))
                E[i][-1] += np.sum(np.dot(Lj,-ci)*Lj,axis=1,keepdims=True)
                E[i][-1] += np.sum(np.dot(1-np.exp(-5.56*Lj),-xi)*Lj,axis=1,keepdims=True)
                E[i][-1] += np.sum(np.dot(Lj/(Lj+0.017),-si)*Lj,axis=1,keepdims=True)
    except Exception as e:
        print(f"Error during energy calculation: {e}")
        return

    for i in range(total_simulations):
        Ei = E[total_simulations-1][i]
        for j in range(total_simulations):
            Ei = np.concatenate((Ei,E[j][i]),axis=1)
        np.savetxt(f'Energy/ESim{i+1}.dat', Ei, fmt='%12.5f')

    for i in range(total_simulations):
        Li = Lambda[i]
        np.savetxt(f'Lambda/Lambda{i+1}.dat', Li, fmt='%10.6f')
