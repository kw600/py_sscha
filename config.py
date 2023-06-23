import os

#archer2 details. Please change the account name!!
taskname = 'scha3' #name of the task. Should be different if several auto_run.py are running at the same time 
account = 'e89-ic_m' #account name

#general
#supercell size in each direction
nq1 = 3 
nq2 = 3
nq3 = 3
T0 = 300
nqirr = 6 #number of irreducible q points
N_prim = 13 #number of atoms in the primitive cell


# configurations for DFT calculations
#number of configurations. Better to be a power of 64
N_config = 256

#minimization
#0.1 for the first few ensembles (with large FC gradients) and 0.5 for the later ones
if os.path.isdir('run_dft2'):
    kong_liu_ratio = 0.5
else: 
    kong_liu_ratio = 0.2

