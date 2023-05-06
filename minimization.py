import sys,os
import numpy as np
import cellconstructor as CC
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
import config

def collect_data(pop,N=-1):
	index=''
	directory = f"run_dft{pop}"
	if N==-1:
		output_filenames = [f for f in os.listdir(directory) if f.endswith(".pwo")] # We select only the output files
		output_files = [os.path.join(directory, f) for f in output_filenames]
		energies = np.zeros(len(output_files))
	else:
		output_filenames = [N]
		output_files = [os.path.join(directory, f) for f in output_filenames]
		energies = np.zeros((10000,))
	 # We add the directory/outpufilename to load them correctly
	# We prepare the array of energies
	 
	for file in output_files:
			# print(file)
		try:	
			# Get the number of the configuration.
			id_number = int(file.split("_")[-1].split(".")[0]) # The same as before, we need the to extract the configuration number from the filename
			
			# Load the file
			ff = open(file, "r")
			lines = [l.strip() for l in ff.readlines()] # Read the whole file removing tailoring spaces
			ff.close()
			
			# Lets look for the energy (in espresso the first line that starts with !)
			# next is used to find only the first occurrence
			#try:
			energy_line = next(l for l in lines if len(l) > 0 if l.split()[0] == "!")
			#except:
			#	print("Error: no energy found in file {}".format(file))
			# Lets collect the energy (the actual number is the 5th item on the line, but python indexes start from 0)
			# note, also the id_number are saved starting from 1
			energies[id_number - 1] = float(energy_line.split()[4])
			
			# Now we can collect the force
			# We need the number of atoms
			nat_line = next( l for l in lines if len(l) > 0 if l.split()[0] == "number" and l.split()[2] == "atoms/cell" )
			nat = int(nat_line.split()[4])
			#except:
			#	print("Error: no number of atoms found in file {}".format(file))
			# Now allocate the forces and read them
			forces = np.zeros((nat, 3))
			forces_lines = [l for l in lines if len(l) > 0 if l.split()[0] == "atom"] # All the lines that starts with atom will contain a force
			for i in range(nat):
				forces[i, :] = [float(x) for x in forces_lines[i].split()[-3:]] # Get the last three number from the line containing the force
			
			# Now we can take the stress tensor
			stress = np.zeros((3,3))
			# We pick the index of the line that starts with the words total stress
			index_before_stress = next(i for i, l in enumerate(lines) if len(l) > 0 if l.split()[0] == "total" and l.split()[1] == "stress")
			# The stress tensor is located just after it
			for i in range(3):
				index = i + index_before_stress + 1
				stress[i, :] = [float(x) for x in lines[index].split()[:3]]

			# We can save the forces_population1_X.dat and pressures_population1_X.dat files
			force_file = os.path.join(f"ens{pop}", "forces_population{}_{}.dat".format(pop,id_number))
			stress_file = os.path.join(f"ens{pop}", "pressures_population{}_{}.dat".format(pop,id_number))
			np.savetxt(force_file, forces)
			np.savetxt(stress_file, stress)
		except:
			print("Error: something went wrong with file {}".format(file))
			# index=index+file.replace('_','.').split('.')[-2]
	# Now we read all the configurations, we can save the energy file
	energy_file = os.path.join(f"ens{pop}", f"energies_supercell_population{pop}.dat")
	np.savetxt(energy_file, energies)
	# if index=='':
	# 	return True
	# return index

def scha(pop):
	IO_freq = sscha.Utilities.IOInfo()
	IO_freq.SetupSaving(f"minim_info{pop}")

	if pop == 1:
			dyn = CC.Phonons.Phonons(f"harmonic_{config.nq1}{config.nq2}{config.nq3}_dyn", nqirr = config.nqirr)
			dyn.Symmetrize()
			dyn.ForcePositiveDefinite()
	else:
			dyn = CC.Phonons.Phonons(f"dyn_pop{int(pop-1)}_", nqirr = config.nqirr)

	ensemble = sscha.Ensemble.Ensemble(dyn, T0 = config.T0, supercell= dyn.GetSupercell())
	ensemble.load(f"ens{pop}", population = pop, N = config.N_config)
	if pop == 1:
			ensemble.update_weights(dyn, config.T0) # Restore the original density matrix at T = 100 K
	minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
	
	# Ignore the structure minimization (is fixed by symmetry)
	minimizer.minim_struct = False


	# Setup the threshold for the ensemble wasting
	minimizer.kong_liu_ratio = config.kong_liu_ratio # Usually 0.5 is a good value

	# Lest start the minimization
	# minimizer.init()
	minimizer.run(custom_function_post = IO_freq.CFP_SaveAll)
	minimizer.finalize()
	minimizer.dyn.save_qe(f"dyn_pop{pop}_")
	print('Converged?',minimizer.is_converged())
	return minimizer.is_converged()

if __name__ == "__main__":
	pop=int(sys.argv[1])
	collect_data(pop)
	min=scha(pop)
	

	
