import cellconstructor as CC
import cellconstructor.Structure
import cellconstructor.Phonons
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
import sys,os
import config

def generate_ensemble():
	if config.population == 1:
		dyn = CC.Phonons.Phonons("harmonic_dyn", nqirr = config.nqirr)
	else:
		dyn = CC.Phonons.Phonons(f"dyn_pop{config.population-1}_", nqirr = config.nqirr)
	
	# Apply the sum rule and symmetries
	dyn.Symmetrize()

	# Flip the imaginary frequencies into real ones
	dyn.ForcePositiveDefinite()
        
	ensemble = sscha.Ensemble.Ensemble(dyn, T0 = config.T0, supercell= dyn.GetSupercell())
	# We generate N randomly displaced structures in the supercell
	ensemble.generate(N = config.N_config)

	ensemble.save(f"ens{config.population}", population = config.population)
	return ensemble

def generate_dft_input():
	ensemble = generate_ensemble()
	typical_espresso_header = f"""
&control
	calculation = "scf"
	tstress = .true.
	tprnfor = .true.
	disk_io = "none"
	pseudo_dir = "pseudo_espresso"
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
	# We extract the number of atoms form the ensemble and the celldm(1) from the dynamical matrix (it is stored in Angstrom, but espresso wants it in Bohr)
	# You can also read it on the fourth value of the first data line on the first dynamical matrix file (dyn_start_popilation1_1); In the latter case, it will be already in Bohr.

	# Now we need to read the scf files
	all_scf_files = [os.path.join(f"ens{config.population}", f) for f in os.listdir(f"ens{config.population}") if f.startswith("scf_")]

	# In the previous line  I am reading all the files inside ens{config.population} os.listdir(ens{config.population}) and iterating over them (the f variable)
	# I iterate only on the filenames that starts with scf_ 
	# Then I join the directory name ens{config.population} to f. In unix it will be equal to ens{config.population}/scf_....
	# (using os.path.join to concatenate path assure to have the correct behaviour independently on the operating system

	# We will generate the input file in a new directory
	if not os.path.exists(f"run_dft{config.population}"):
		os.mkdir(f"run_dft{config.population}")

	for file in all_scf_files:
		# Now we are cycling on the scf_ files we found.
		# We must extract the number of the file
		# The file is the string "ens{config.population}/scf_population1_X.dat"
		# Therefore the X number is after the last "_" and before the "." character
		# We can split before the string file at each "_", isolate the last part "X.dat"
		# and then split it again on "." (obtaining ["X", "dat"]) and select the first element
		# then we convert the "X" string into an integer
		number = int(file.split("_")[-1].split(".")[0])
		
		# We decide the filename for the espresso input
		# We will call it run_dft{config.population}/espresso_run_X.pwi
		filename = os.path.join(f"run_dft{config.population}", "espresso_run_{}.pwi".format(number))
		
		# We start writing the file
		with open(filename, "w") as f:
			# We write the header
			f.write(typical_espresso_header)
			
			# Load the scf_population_X.dat file
			ff = open(file, "r")
			structure_lines = ff.readlines()
			ff.close()
			
			# Write the content on the espresso_run_X.pwi file
			# Note in the files we specify the units for both the cell and the structure [Angstrom]
			f.writelines(structure_lines)
	return ensemble

if __name__ == "__main__":
    generate_dft_input()
    


