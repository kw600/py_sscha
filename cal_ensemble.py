import cellconstructor as CC
import cellconstructor.Structure
import cellconstructor.Phonons
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
import sys,os
import config

def generate_ensemble(pop):
	if pop == 1:
		dyn = CC.Phonons.Phonons("harmonic_dyn", nqirr = config.nqirr)
	else:
		dyn = CC.Phonons.Phonons(f"dyn_pop{pop-1}_", nqirr = config.nqirr)
	
	# Apply the sum rule and symmetries
	dyn.Symmetrize()

	# Flip the imaginary frequencies into real ones
	dyn.ForcePositiveDefinite()
        
	ensemble = sscha.Ensemble.Ensemble(dyn, T0 = config.T0, supercell= dyn.GetSupercell())
	# We generate N randomly displaced structures in the supercell
	ensemble.generate(N = config.N_config)

	ensemble.save(f"ens{pop}", population = pop)
	return ensemble

def generate_dft_input(pop):
	ensemble = generate_ensemble(pop)
	typical_espresso_header = f"""
&control
	calculation = "scf"
	tstress = .true.
	tprnfor = .true.
	disk_io = "none"
	pseudo_dir = "/work/e89/e89/kw2318/pseudo_espresso"
&end
&system
	nat = {ensemble.structures[0].N_atoms}
	ntyp = 2
	ibrav = 0
	ecutwfc = {config.ecutwfc_2}
	ecutrho = {config.ecutrho_2}
&end
&electrons
	conv_thr = {config.conv_thr_2}
	!diagonalization = "cg"
&end

ATOMIC_SPECIES
	Pb 207.2 Pb.upf
	Te 127.6 Te.upf

K_POINTS automatic
1 1 1 0 0 0
"""
	all_scf_files = [os.path.join(f"ens{pop}", f) for f in os.listdir(f"ens{pop}") if f.startswith("scf_")]

	# In the previous line  I am reading all the files inside ens{pop} os.listdir(ens{pop}) and iterating over them (the f variable)
	# I iterate only on the filenames that starts with scf_ 
	# Then I join the directory name ens{pop} to f. In unix it will be equal to ens{pop}/scf_....
	# (using os.path.join to concatenate path assure to have the correct behaviour independently on the operating system

	# We will generate the input file in a new directory
	if not os.path.exists(f"run_dft{pop}"):
		os.mkdir(f"run_dft{pop}")

	for file in all_scf_files:
		number = int(file.split("_")[-1].split(".")[0])
		filename = os.path.join(f"run_dft{pop}", "espresso_run_{}.pwi".format(number))
		
		# We start writing the file
		with open(filename, "w") as f:
			f.write(typical_espresso_header)
			
			# Load the scf_population_X.dat file
			ff = open(file, "r")
			structure_lines = ff.readlines()
			ff.close()
			f.writelines(structure_lines)
	return ensemble

if __name__ == "__main__":
	pop = int(sys.argv[1])
	generate_dft_input(pop)



