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
			return False

	return True

def check_complete(output_dir):
	b=''
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
				return False, b
	return True, b

def DFT(pop):
	current_path = os.getcwd()
	subprocess.run(["./step1"])

	if os.path.isfile(os.path.join(current_path, "dyn3")):
		print("File dyn3 exists. Proceed to next step.")
	else:
		print('Waiting for dyn3 file...')
		time.sleep(60)

	subprocess.run(["./step2"])

	DFT_path=os.path.join(current_path, f"run_dft{pop}")

	if check(DFT_path):
		print("DFT calculations done. Check whether results are complete.")
	else:
		print('Waiting for DFT calculations...')
		time.sleep(30)


	while True:
		if check_complete(DFT_path)[0]:
			print("DFT calculations complete. Proceed to minimization.")
			break
		else:
			print(f"DFT calculations with index {check_complete(DFT_path)[1]} incomplete. ")
			subprocess.run(["./step3",pop])
			print("Waitinf for DFT calculations.")
			time.sleep(30)

converge = False
pop = 1
while converge:
	DFT(pop)
	collect_data(pop)
	converge=scha(pop)
	pop+=1



