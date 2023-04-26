import ase
from ase.calculators.espresso import Espresso
from ase.visualize import view

PbTe_atoms = ase.io.read("PbTe.cif")

import cellconstructor as CC, cellconstructor.Phonons
import cellconstructor.ForceTensor
# Import the numerical libraries and those for plotting
import numpy as np
import matplotlib.pyplot as plt
import sys, os
# Let us define the PATH in the brilluin zone and the total number of points
PATH = "XGKMWLG"
# PATH = "GK"
N_POINTS = 200
# Here we define the position of the special points
SPECIAL_POINTS = {"G": [0,0,0],
"X": [0.5,0, 0.5],
"L": [.5, .5, .5],
"W": [.5, .25, .75],
"K": [3/8., 3/8., 3/4],
"M": [0,0.5,0.5]}
# The two dynamical matrix to be compared

# PATH='GG'
# HARM_DYN = 'harmonic_444_dyn'
SSCHA_DYN2 = 'dyn_pop3_'
# SSCHA_DYN3 = 'dyn333_pop3_'
# The number of irreducible q points
# i.e., the number of files in which the phonons are stored


# Load the harmonic and sscha phonons
# harmonic_dyn = CC.Phonons.Phonons(HARM_DYN, 8)
sscha_dyn2 = CC.Phonons.Phonons(SSCHA_DYN2, 3)
#sscha_dyn3 = CC.Phonons.Phonons(SSCHA_DYN3, 4)
w_s, pols = sscha_dyn2.DyagDinQ(0)
print ("\n".join(["{:.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in  w_s]))

# Get the band path

qpath, data = CC.Methods.get_bandpath(sscha_dyn2.structure.unit_cell, PATH, SPECIAL_POINTS,N_POINTS)
# print(sscha_dyn2.structure.unit_cell)
xaxis, xticks, xlabels = data # Info to plot correclty the x axis
# print(qpath)
# print('qpath')
# print(data)
# Get the phonon dispersion along the path
#harmonic_dispersion = CC.ForceTensor.get_phonons_in_qpath(harmonic_dyn, qpath)
sscha_dispersion2 = CC.ForceTensor.get_phonons_in_qpath(sscha_dyn2, qpath)
print('sscha0')
print(sscha_dispersion2[:,5])
print('sscha1')
print(sscha_dispersion2[:,5])
#sscha_dispersion3 = CC.ForceTensor.get_phonons_in_qpath(sscha_dyn3, qpath)
# nmodes = harmonic_dyn.struct

nmodes = sscha_dyn2.structure.N_atoms * 3
# Plot the two dispersions
plt.figure(dpi = 150)
ax = plt.gca()
for i in range(nmodes):
# for i in [5]:
	lbl=None
	lblsscha = None
	if i == 0:
		#ax.plot(xaxis, harmonic_dispersion[:,i], color = 'k', ls = 'dashed',label = 'Harmonic')
		ax.plot(xaxis, sscha_dispersion2[:,i], color = 'r', label = '222')
		#ax.plot(xaxis, sscha_dispersion3[:,i], color = 'b', label = '333')
	else:
		#ax.plot(xaxis, harmonic_dispersion[:,i], color = 'k', ls = 'dashed')
		ax.plot(xaxis, sscha_dispersion2[:,i], color = 'r')
		#ax.plot(xaxis, sscha_dispersion3[:,i], color = 'b')
# Plot vertical lines for each high symmetry points
for x in xticks:
	ax.axvline(x, 0, 1, color = "k", lw = 0.4)
	ax.axhline(0, 0, 1, color = 'k', ls = ':', lw = 0.4)
# Set the x labels to the high symmetry points

ax.set_xticks(xticks)
ax.set_xticklabels(xlabels)
ax.set_xlabel("Q path")
ax.set_ylabel("Phonons [cm-1]")
plt.tight_layout()
plt.legend()
plt.savefig("dispersion.png")
plt.show()
