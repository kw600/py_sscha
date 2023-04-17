import os

#archer2 details. Please change the account name!!
taskname = 'scha1' #name of the task. Should be different if several auto_run.py are running at the same time 
account = 'e89-ic_m' #account name
#number of tasks per node. Please leave a space between '=' and number.
n_node_per_job = 1
nrun_per_node = 16 
one_by_one = True

check_one_by_one = True
check_node = 1

cq = 0
for root, dirs, files in os.walk("."):
    for file in files:
        if "harmonic_dyn" in file:
            cq += 1
nqirr = cq - 1 #number of irreducible q points


#general
nq1 = 2 #supercell size in each direction
nq2 = 2
nq3 = 2
T0 = 100

#step1 initial relaxation
ecutwfc = 90 # The plane-wave wave-function cutoff
ecutrho = 360 # The density wave-function cutoff,
conv_thr = 1e-8 # The convergence for the DFT self-consistency
k_spacing = 0.2 #A^-1 The minimum distance in the Brillouin zone sampling

# configurations for DFT calculations
#number of configurations. Better to be a power of 64
N_config = 512
#maximum population/iteration
maxpop = 20 


#step2 DFT calculations for the ensemble of configurations
ecutwfc_2 = 60
ecutrho_2 = 240
conv_thr_2 = 1e-8


#minimization
min_step_dyn = 0.3
kong_liu_ratio = 0.2
