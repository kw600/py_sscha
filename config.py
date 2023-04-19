import os

#archer2 details. Please change the account name!!
taskname = 'scha1' #name of the task. Should be different if several auto_run.py are running at the same time 
account = 'e89-ic_m' #account name
#number of tasks per node. Please leave a space between '=' and number.
n_node_per_job = 1
nrun_per_node = 32
hour = '02'
one_by_one = True


#general
nq1 = 3 #supercell size in each direction
nq2 = 3
nq3 = 3
T0 = 300

cq = 0
for root, dirs, files in os.walk("."):
    for file in files:
        if f"harmonic_{nq1}{nq2}{nq3}_dyn" in file:
            cq += 1
nqirr = cq - 1 #number of irreducible q points





# configurations for DFT calculations
#number of configurations. Better to be a power of 64
N_config = 1024
#maximum population/iteration
maxpop = 20


#step2 DFT calculations for the ensemble of configurations
ecutwfc_2 = 65
conv_thr_2 = 1e-8


#minimization
#0.1 for the first few ensembles (with large FC gradients) and 0.5 for the later ones
kong_liu_ratio = 0.2