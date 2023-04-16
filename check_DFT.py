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

b='298 503 193 205 308 17 91 455 419 261 214 117 210 381 370 233 7 147 119 446 108 219 429 79 476 95 80 94 140 253 65 484 163 354 225 300 75 146 174 72 504 283 195 247 164 319 143 92 393 473 90 199 178 250 344 367 398 435 463 369 145 201 286 408 318 114 495 508 427 461 451 467 144 360 412 325 194 169 255 333 224 454 151 239 34 296 212 342 330 510 404 5 100 275 82 341 500 76 185 200 486 320 248 460 197 376 93 270 189 118 416 444 441 19 155 192 489 139 141 264 162 230 302 365 272 182 420 290 234 443 284 256 196 40 289 131 437 259 115 483 501 134 198 121 36 70 278 167 211 422 138 353 69 445 310 217 215 280 89 85 485 442 267 343 77 97 203 51 59 396 99 130 317 431 240 68 403 161 368 338 171 232 273 490 24 181 266 183 271 274 314 457 312 244 10 364 311 258 227 268 260 206 156 229 439 223 160 487 125 287 305 226 292 166 361 28 129 347 303 315 35 254 42 362 309 497 176 128 493 470 228 288 245 133 450 188 373 251 482 294 276 480 494 281 263 136 477 30 22 6 423 177 109 307 67 74 218 447 66 172'			
l='{'
r='}'
n_node=int(np.ceil(len(b.split())/config.nrun_per_node))
n_job=int(np.ceil(n_node/4))
l0=b.split()
nn=1
for i in range(n_job):
	index=''
	if len(l0)>=64:
		for j in range(64):
			index=index+l0.pop()+" "
	else:
		for j in range(len(l0)):
			index=index+l0.pop()+" "
	print(index)
	sub=f"""#!/bin/bash
# Slurm job options (job-name, compute nodes, job time)
#SBATCH --nodes=4
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

srun --unbuffered --nodes=1 --ntasks={int(128/config.nrun_per_node)} --tasks-per-node={int(128/config.nrun_per_node)} \
		--cpus-per-task=1 --distribution=block:block --hint=nomultithread \
		--mem={int(200/config.nrun_per_node)}G --exact \
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