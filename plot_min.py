import matplotlib.pyplot as plt
import numpy as np

# Path to the file containing the results
f = np.loadtxt("p1")
label=['Free energy','FC gradient','Kong-Liu N_eff']
count=0
for i in [1,3,7]:
	plt.figure()
	plt.plot(f[:,0], f[:,i], 'o',label=label[count])
	plt.legend()
	count+=1

plt.show()
