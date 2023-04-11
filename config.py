#general
nq1 = 2 #supercell size in each direction
nq2 = 2
nq3 = 2
nqirr = 3 #number of irreducible q points
T0 = 100

#step1 initial relaxation
ecutwfc = 60, # The plane-wave wave-function cutoff
ecutrho = 240, # The density wave-function cutoff,
conv_thr = 1e-6, # The convergence for the DFT self-consistency
k_spacing = 0.2 #A^-1 The minimum distance in the Brillouin zone sampling

# configurations for DFT calculations
N_config = 10
population = 1

#step2 DFT calculations for the ensemble of configurations
ecutwfc_2 = 40
ecutrho_2 = 160
conv_thr_2 = 1e-6


#minimization
min_step_dyn = 0.5
kong_liu_ratio = 0.5
