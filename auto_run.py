import subprocess
import sys, os, time
import config
from minimization import *

def check_dft(output_dir):
	for i in range(1, config.N_config+1):
		# Construct the filename for the current index
		filename = f"espresso_run_{i}.pwo"
		# Check if the file exists in the output directory
		if not os.path.exists(os.path.join(output_dir, filename)):
			# If the file does not exist, print an error message
			print(f"File {filename} does not exist.",config.N_config)
			return False
	return True

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
	with open(os.path.join('./', 'q.dat'), "r") as f:
		contents = f.readlines()
	return len(contents)

def DFT(pop):
	current_path = os.getcwd()
	if pop==1:
		subprocess.run(["./step1"])

		#check if the submitted job is finished
		while True:
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
			print('Re-submit job')
			subprocess.run(["./step2",str(pop)])
		else:
			print('Waiting for DFT calculations...')
			time.sleep(30)
		
	while True:
		if check_complete(DFT_path)[0]:
			print("DFT calculations complete. Proceed to minimization.")
			break
		else:
			print(f"DFT calculations with index {check_complete(DFT_path)[1]} incomplete. ")
			subprocess.run(["./step3",str(pop)])
			while True:
				if checkq()==1:
					break
				else:
					print('Waiting for DFT calculations...')
					time.sleep(30)
	
converge = False
pop = 1
while not converge:
	DFT(pop)
	collect_data(pop)
	print('nqirr',config.nqirr)
	converge=scha(pop)
	pop+=1
	if pop>=config.maxpop:
		print("Maximum population reached. Please check the results.")
		break



