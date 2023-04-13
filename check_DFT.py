import os
import config
# Specify the directory containing the output files
output_dir = f"./run_dft{config.population}/"
b=''

for i in range(1, 129):
	# Construct the filename for the current index
	filename = f"espresso_run_{i}.pwo"
	# Check if the file exists in the output directory
	if not os.path.exists(os.path.join(output_dir, filename)):
		# If the file does not exist, print an error message
		print(f"File {filename} does not exist.")
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
print(b)
sub=f"""#!/bin/bash
# 
# Parallel script produced by bolt
#        Resource: ARCHER2 (HPE Cray EX (128-core per node))
#    Batch system: Slurm
#
# bolt is written by EPCC (http://www.epcc.ed.ac.uk)
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --job-name=s1
#SBATCH --account=e89-ic_m
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=01:0:0
module load quantum_espresso
I='{b}'
for index in $I
do
srun --nodes=1 --ntasks=16 --ntasks-per-node=16 --mem=10240M --distribution=block:block --hint=nomultithread pw.x < espresso_run_${l}index{r}.pwi > espresso_run_${l}index{r}.pwo
done"""

with open(f"./run_dft{config.population}/dft_continue", "w") as f:
	f.write(sub)
f.close()