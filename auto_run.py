import subprocess
import sys, os, time
import config
from minimization import *

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
#SBATCH --nodes=4
#SBATCH --job-name={config.taskname}
#SBATCH --account={config.account}
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=01:0:0
I=$1

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1

module load quantum_espresso

echo "Starting 32 jobs across 4 nodes each using 16 CPUs"

# Loop over 32 subjobs each using 16 CPUs, running in background

for i in $(seq 1 32)
do
   index=$(( (I-1)*32 + i ))
   echo "Launching job number $i with index $index"
   # Launch subjob overriding job settings as required and in the
   # background. Make sure to change the `--mem=` flag to the amount
   # of memory required. A sensible amount is 1.5 GiB per task as
   # this leaves some overhead for the OS etc.

   srun --unbuffered --nodes=1 --ntasks=16 --tasks-per-node=16 {dd}
        --cpus-per-task=1 --distribution=block:block --hint=nomultithread {dd}
        --mem=24G --exact {dd}
        pw.x < espresso_run_${l}index{r}.pwi > espresso_run_${l}index{r}.pwo &

done
# Wait for all subjobs to finish
echo "Waiting for all jobs to finish ..."

wait

echo "... all jobs finished"
"""

	with open('sub_archer2', 'w') as f:
		f.write(s1)
	with open('sub_dft', 'w') as f:
		f.write(s2)

def check_dft(output_dir):
	n=0
	for i in range(1, config.N_config+1):
		# Construct the filename for the current index
		filename = f"espresso_run_{i}.pwo"
		# Check if the file exists in the output directory
		if os.path.exists(os.path.join(output_dir, filename)):
			# If the file does not exist, print an error message
			# print(f"File {filename} does not exist.",config.N_config)
			n+=1
	if n==config.N_config:
		return True
	elif n>0.75*config.N_config:
		print("{config.N_config-n} files are missing.")
		return True
	else:
		# print(n)
		return False


def check_complete(output_dir,key='JOB DONE'):
	b=''
	# Loop through all the files in the directory
	for filename in os.listdir(output_dir):
		# Check if the file is a text file
		if filename.endswith(".pwo"):
			
			# Open the file and read its contents
			with open(os.path.join(output_dir, filename), "r") as f:
				contents = f.read()
			# Check if the keyword "Job done" is in the file contents
			if key not in contents:
				a=filename.replace("_",".")
				a=a.split(".")
				# If the keyword is not found, print the filename
				b+=a[-2]+" "
	if b=='':
		return True, ''
	else:
		return False, b

def checkq():
	subprocess.run(["./queue"])
	c=1
	with open(os.path.join('./', 'q.dat'), "r") as f:
		contents = f.readlines()
	for i in contents:
		if config.taskname in i:
			c+=1
	return c

def DFT(pop):
	current_path = os.getcwd()
	if pop==1:
		if not os.path.exists('./harmonic_dyn0'):
			subprocess.run(["./step1"])
			T=True
		else:
			print('Harminc calculations already done. Skip step1.')
			T=False
		#check if the submitted job is finished
		while T:
			if checkq()==1:
				print("Harmonic calculations done.")
				break
			else:
				print('Waiting for harmonic calculations to be finished...')
				time.sleep(30)

	#run the step2
	subprocess.run(["./step2",str(pop)])

	#check if the DFT calculations are finished
	DFT_path=os.path.join(current_path, f"run_dft{pop}")
	while True:
		if check_dft(DFT_path) and checkq()==1:
			print("DFT calculations done. Check whether results are complete.")
			break
		elif not check_dft(DFT_path) and checkq()==1:
			print('Considering to Resubmit the job')
			print('No job is running.')
			exit()
			# subprocess.run(["./step2",str(pop)])
		else:
			print('Waiting for DFT calculations...')
			time.sleep(30)
		
	while True:
		if check_complete(DFT_path)[0]:
			print("DFT calculations complete. Proceed to minimization.")
			break
		else:
			print(f"DFT calculations with index {check_complete(DFT_path)[1]} incomplete. ")
			try:
				subprocess.run(["./step3",str(pop)])
			except:
				print("Error in submitting incomplete job. Please check the error message.")
				exit()
			while True:
				if checkq()==1:
					break
				else:
					print('Waiting for DFT calculations...')
					time.sleep(30)

if __name__ == '__main__':
	write_sub()
	converge = False
	pop = 1
	while not converge:
		DFT(pop)
		collect_data(pop)
		print('nqirr',config.nqirr)
		converge=scha(pop)
		pop+=1
		if pop>config.maxpop:
			print("Maximum population reached. Please check the results.")
			break



