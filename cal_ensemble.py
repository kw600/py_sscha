import cellconstructor as CC
import cellconstructor.Structure
import cellconstructor.Phonons
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
import sys,os
import numpy as np
import config

def generate_ensemble(pop):
	if pop == 1:
		dyn = CC.Phonons.Phonons(f"harmonic_{config.nq1}{config.nq2}{config.nq3}_dyn", nqirr = config.nqirr)
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
&end

&electrons
	conv_thr =  1d-6
&end

ATOMIC_SPECIES
Pb 207.2 Pb_ONCV_PBE-1.2.upf
Te 127.6 Te_ONCV_PBE-1.2.upf
K_POINTS automatic
{int(np.ceil(8/config.nq1))} {int(np.ceil(8/config.nq2))} {int(np.ceil(8/config.nq3))}  0 0 0
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

def generate_vasp_input(pop):
	ensemble = generate_ensemble(pop)
	all_scf_files = [os.path.join(f"ens{pop}", f) for f in os.listdir(f"ens{pop}") if f.startswith("scf_")]
	if not os.path.exists(f"run_dft{pop}"):
		os.mkdir(f"run_dft{pop}")

	for i,file in enumerate(all_scf_files):
		if not os.path.exists(f"run_dft{pop}/vasp{i+1}"):
			os.mkdir(f"run_dft{pop}/vasp{i+1}")
		number = int(file.split("_")[-1].split(".")[0])
		filename = os.path.join(f"run_dft{pop}/vasp{i+1}", "POSCAR".format(number))
		
		# We start writing the file
		with open(filename, "w") as f:
			
			# Load the scf_population_X.dat file
			ff = open(file, "r")
			structure_lines = ff.readlines()
			lattice = np.array([i.split() for i in structure_lines[1:4]],dtype=float)
			atoms   = structure_lines[6:]
			sc_lines = []
			v_lines = []
			sn_lines = []
			for line in atoms:
				if line.startswith('Sc'):
					sc_lines.append(line)
				elif line.startswith('V'):
					v_lines.append(line)
				elif line.startswith('Sn'):
					sn_lines.append(line)
			ff.close()
			f.writelines('ScV6Sn6 2 by 2 by 2\n   1.00000000000000\n')
			f.writelines(structure_lines[1:4])
			f.writelines('   Sc   V    Sn\n     8     48     48\nDirect\n')
			for i in sc_lines:
				coord = np.array(i.split()[1:],dtype=float)
				frac = np.dot(np.linalg.inv(lattice),coord)
				f.write(f'{frac[0]}  {frac[1]}  {frac[2]}\n')
			for i in v_lines:
				coord = np.array(i.split()[1:],dtype=float)
				frac = np.dot(np.linalg.inv(lattice),coord)
				f.write(f'{frac[0]}  {frac[1]}  {frac[2]}\n')
			for i in sn_lines:
				coord = np.array(i.split()[1:],dtype=float)
				frac = np.dot(np.linalg.inv(lattice),coord)
				f.write(f'{frac[0]}  {frac[1]}  {frac[2]}\n')
	return ensemble

if __name__ == "__main__":
	
	pop = int(sys.argv[1])
	generate_vasp_input(pop)



