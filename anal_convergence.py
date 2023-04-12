import random
from mpi4py import MPI
import random
import numpy as np

def generate_random_integers(N):
	"""
	Returns a list of N random non-repetitive integers ranging from 1 to 120.
	"""
	if N > 120:
		raise ValueError("N cannot be greater than 120")
	
	integers = list(range(0, 120))
	random.shuffle(integers)
	
	return integers[:N]


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0:
	# Core 0 receives the calculated results from other cores
	all_results = []
	for i in range(1, size):
		results = comm.recv(source=i)
		all_results.extend(results)
	print(f"All results: {all_results}")
else:
	# Determine the number of configurations to be generated
	num_configs = int(np.ceil(rank/10))*10
	a = generate_random_integers(num_configs)
	split_mask = np.zeros(120, dtype=bool)
	split_mask[a] = True
	dyn = CC.Phonons.Phonons("harmonic_dyn", nqirr = config.nqirr)
	dyn.Symmetrize()
	dyn.ForcePositiveDefinite()
	ensemble = sscha.Ensemble.Ensemble(dyn, T0 = config.T0, supercell= dyn.GetSupercell())
	ensemble.load("data_ensemble_manual", population = config.population, N = config.N_config)
	ensemble.split(split_mask)
	ensemble.update_weights(dyn, config.T0) # Restore the original density matrix at T = 100 K
	minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
	minimizer.minim_struct = False
	minimizer.min_step_dyn = config.min_step_dyn # Values around 1 are good
	minimizer.kong_liu_ratio = config.kong_liu_ratio # Usually 0.5 is a good value
	minimizer.init()
	minimizer.run()
	results = minimizer.results
	
	# Send the calculated results to core 0
	comm.send(results, dest=0)