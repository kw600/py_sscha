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
PATH = "GXWXKGL"
N_POINTS = 1000
# Here we define the position of the special points
SPECIAL_POINTS = {"G": [0,0,0],
"X": [0, .5, .5],
"L": [.5, .5, .5],
"W": [.25, .75, .5],
"K": [3/8., 3/4., 3/8.]}
# The two dynamical matrix to be compared
HARM_DYN = 'harmonic_dyn'
SSCHA_DYN2 = 'dyn_pop11_'
SSCHA_DYN3 = 'dyn_pop9_'
# The number of irreducible q points
# i.e., the number of files in which the phonons are stored


# Load the harmonic and sscha phonons
harmonic_dyn = CC.Phonons.Phonons(HARM_DYN, 4)
sscha_dyn3 = CC.Phonons.Phonons(SSCHA_DYN3, 4)
sscha_dyn2 = CC.Phonons.Phonons(SSCHA_DYN2, 3)
# Get the band path
qpath, data = CC.Methods.get_bandpath(harmonic_dyn.structure.unit_cell, PATH, SPECIAL_POINTS,N_POINTS)
xaxis, xticks, xlabels = data # Info to plot correclty the x axis

# Get the phonon dispersion along the path
harmonic_dispersion = CC.ForceTensor.get_phonons_in_qpath(harmonic_dyn, qpath)
sscha_dispersion2 = CC.ForceTensor.get_phonons_in_qpath(sscha_dyn2, qpath)
sscha_dispersion3 = CC.ForceTensor.get_phonons_in_qpath(sscha_dyn3, qpath)
# nmodes = harmonic_dyn.struct

nmodes = harmonic_dyn.structure.N_atoms * 3
# Plot the two dispersions
plt.figure(dpi = 150)
ax = plt.gca()
for i in range(nmodes):
	lbl=None
	lblsscha = None
	if i == 0:
		ax.plot(xaxis, harmonic_dispersion[:,i], color = 'k', ls = 'dashed',label = 'Harmonic')
		ax.plot(xaxis, sscha_dispersion2[:,i], color = 'r', label = '222')
		ax.plot(xaxis, sscha_dispersion3[:,i], color = 'b', label = '333')
	else:
		ax.plot(xaxis, harmonic_dispersion[:,i], color = 'k', ls = 'dashed')
		ax.plot(xaxis, sscha_dispersion2[:,i], color = 'r')
		ax.plot(xaxis, sscha_dispersion3[:,i], color = 'b')
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