import ase
from ase.calculators.espresso import Espresso
from ase.visualize import view

PbTe_atoms = ase.io.read("PbTe.cif")

from init_structure import *
import config 

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

def sscha():
	PbTe_atoms = structure("PbTe.cif")
	PbTe_atoms.generate_input_for_initial_relax()
	PbTe_atoms.generate_input_for_initial_phonon()
	write_sub()
if __name__ == "__main__":
    sscha()


