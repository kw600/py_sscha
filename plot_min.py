import matplotlib.pyplot as plt
import numpy as np

# Path to the file containing the results
f = np.loadtxt("minim_info.dat")
label=['FC gradient','Kong-Liu N_eff']
count=0
for i in [2,6]:
	plt.figure()
	plt.plot(f[:,i], 'o',label=label[count])
	plt.legend()
	count+=1

plt.show()
