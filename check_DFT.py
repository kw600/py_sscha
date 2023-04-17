import os
import config
import sys
import numpy as np
from minimization import *
# Specify the directory containing the output files
pop=sys.argv[1]

output_dir = f"./run_dft{pop}/"
b=''

for i in range(1, config.N_config+1):
	# Construct the filename for the current index
	filename = f"espresso_run_{i}.pwo"
	# Check if the file exists in the output directory
	if not os.path.exists(os.path.join(output_dir, filename)):
		# If the file does not exist, print an error message
		# print(f"File {filename} does not exist.")
		b+=str(i)+" "


# Loop through all the files in the directory
for filename in os.listdir(output_dir):
	# Check if the file is a text file
	if filename.endswith(".pwo"):
		
		# Open the file and read its contents
		with open(os.path.join(output_dir, filename), "r") as f:
			contents = f.read()
		# Check if the keyword "Job done" is in the file contents
		if 'JOB DONE' not in contents:
			a=filename.replace("_",".")
			a=a.split(".")
			# If the keyword is not found, print the filename
			b+=a[-2]+" "

l='{'
r='}'
dd='\\'
n_node=int(np.ceil(len(b.split())/config.nrun_per_node))

n_job=int(np.ceil(n_node/config.n_node_per_job))
l0=b.split()
print(l0)
nn=1
nrun=int(np.ceil(len(l0)/32))
while len(l0)>0:	
	index=''
	if len(l0)>nrun:
		for j in range(nrun):
			index+=l0.pop()+','
	else:
		for j in range(len(l0)):
			index+=l0.pop()+','
	sub1=f"""#!/bin/bash
# Slurm job options (job-name, compute nodes, job time)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --account={config.account}
#SBATCH --job-name={config.taskname}
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=01:0:0

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1

module load quantum_espresso
I='{index}'

for index in $I
do

echo "Launching job number $index"

# Launch subjob overriding job settings as required and in the
# background. Make sure to change the `--mem=` flag to the amount
# of memory required. A sensible amount is 1.5 GiB per task as
# this leaves some overhead for the OS etc.

srun  --distribution=block:block --hint=nomultithread pw.x < espresso_run_${l}index{r}.pwi > espresso_run_${l}index{r}.pwo 

done
# Wait for all subjobs to finish
echo "Waiting for all jobs to finish ..."


echo "JOB DONE"
"""
	with open(f"./run_dft{pop}/dft_continue{nn}", "w") as f:
			f.write(sub1)
	nn+=1

