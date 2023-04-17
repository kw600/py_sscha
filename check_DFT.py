import os
import config
import sys
import numpy as np
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
		if "JOB DONE" not in contents:
			a=filename.replace("_",".")
			a=a.split(".")
			# If the keyword is not found, print the filename
			b+=a[-2]+" "

l='{'
r='}'
dd='\\'
n_node=int(np.ceil(len(b.split())/config.nrun_per_node))
n_job=n_node
l0=b.split()
nn=1
for i in range(n_job):
	index=''
	if len(l0)>=config.nrun_per_node:
		for j in range(config.nrun_per_node):
			index=index+l0.pop()+" "
	else:
		for j in range(len(l0)):
			index=index+l0.pop()+" "
	# print('missing',index)
	sub=f"""#!/bin/bash
# Slurm job options (job-name, compute nodes, job time)
#SBATCH --nodes=1
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

srun --unbuffered --nodes=1 --ntasks={int(128/config.nrun_per_node)} --tasks-per-node={int(128/config.nrun_per_node)} {dd}
		--cpus-per-task=1 --distribution=block:block --hint=nomultithread {dd}
		--mem={int(200/config.nrun_per_node)}G --exact {dd}
		pw.x < espresso_run_${l}index{r}.pwi > espresso_run_${l}index{r}.pwo &

done
# Wait for all subjobs to finish
echo "Waiting for all jobs to finish ..."

wait

echo "JOB DONE"
"""

	with open(f"./run_dft{pop}/dft_continue{nn}", "w") as f:
		f.write(sub)
	nn+=1