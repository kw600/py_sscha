import os

count = 0
for root, dirs, files in os.walk("."):
    for file in files:
        if "harmonic_dyn" in file:
            count += 1
nqirr = count - 1 #number of irreducible q points


#general
nq1 = 3 #supercell size in each direction
nq2 = 3
nq3 = 3
T0 = 100

#step1 initial relaxation
ecutwfc = 90 # The plane-wave wave-function cutoff
ecutrho = 360 # The density wave-function cutoff,
conv_thr = 1e-8 # The convergence for the DFT self-consistency
k_spacing = 0.2 #A^-1 The minimum distance in the Brillouin zone sampling

# configurations for DFT calculations
N_config = 40
population = 1

#step2 DFT calculations for the ensemble of configurations
ecutwfc_2 = 60
ecutrho_2 = 240
conv_thr_2 = 1e-8


#minimization
min_step_dyn = 0.5
kong_liu_ratio = 0.5