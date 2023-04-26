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
	conv_thr = {config.conv_thr_2}
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

def write_sub():
	l='{'
	r='}'
	dd='\\'
	s1=f"""#!/bin/bash
# 
# Parallel script produced by bolt
#        Resource: ARCHER2 (HPE Cray EX (128-core per node))
#    Batch system: Slurm
#
# bolt is written by EPCC (http://www.epcc.ed.ac.uk)
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --job-name={config.taskname}
#SBATCH --account={config.account}
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=01:0:0
module load quantum_espresso
module list
srun --distribution=block:block --hint=nomultithread pw.x < espresso.pwi > espresso.pwo
srun --distribution=block:block --hint=nomultithread ph.x < harmonic.phi > harmonic.pho
echo 'JOB DONE'
"""

	s2=f"""#!/bin/bash
# Slurm job options (job-name, compute nodes, job time)
#SBATCH --nodes={config.n_node_per_job}
#SBATCH --job-name={config.taskname}
#SBATCH --account={config.account}
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time={config.hour}:00:00
I=$1

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1

module load quantum_espresso

echo "Starting 64 jobs across 4 nodes each using 16 CPUs"

# Loop over subjobs each using 16 CPUs, running in background

for i in $(seq 1 {config.n_node_per_job*config.nrun_per_node})
do
index=$(((I-1)*{config.n_node_per_job*config.nrun_per_node}+i))
echo "Launching job number $i with index $index"
# Launch subjob overriding job settings as required and in the
# background. Make sure to change the `--mem=` flag to the amount
# of memory required. A sensible amount is 1.5 GiB per task as
# this leaves some overhead for the OS etc.

srun --unbuffered --nodes=1 --ntasks={int(128/config.nrun_per_node)} --tasks-per-node={int(128/config.nrun_per_node)} {dd}
		--cpus-per-task=1 --distribution=block:block --hint=nomultithread {dd}
		--mem={int(250/config.nrun_per_node)}G --exact {dd}
		pw.x < ../espresso_run_${l}index{r}.pwi > ../espresso_run_${l}index{r}.pwo &

done
# Wait for all subjobs to finish
echo "Waiting for all jobs to finish ..."

wait

echo "... all jobs finished"
"""
	s3=f"""#!/bin/bash
# Slurm job options (job-name, compute nodes, job time)
#SBATCH --nodes={config.n_node_per_job}
#SBATCH --tasks-per-node=128
#SBATCH --job-name={config.taskname}
#SBATCH --account={config.account}
#SBATCH --partition=standard
#SBATCH --qos=taskfarm
#SBATCH --time={config.hour}:00:00
I=$1

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1

module load quantum_espresso

echo "Starting {config.nrun_per_node*config.n_node_per_job} jobs across 1 node one by one"


for i in $(seq 1 {config.nrun_per_node*config.n_node_per_job})
do
index=$(((I-1)*{config.nrun_per_node*config.n_node_per_job}+i))
echo "Launching job number $i with index $index"
# Launch subjob overriding job settings as required and in the
# background. Make sure to change the `--mem=` flag to the amount
# of memory required. A sensible amount is 1.5 GiB per task as
# this leaves some overhead for the OS etc.

srun --distribution=block:block --hint=nomultithread {dd}
		pw.x < ../espresso_run_${l}index{r}.pwi > ../espresso_run_${l}index{r}.pwo 

done

echo "... all jobs finished"
"""

	with open('sub_archer2', 'w') as f:
		f.write(s1)
	
	with open('sub_dft', 'w') as f:
		if config.one_by_one==False:
			f.write(s2)
		else:
			f.write(s3)

if __name__ == "__main__":
	write_sub()
	pop = int(sys.argv[1])
	generate_dft_input(pop)



