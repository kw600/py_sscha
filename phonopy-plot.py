import ase
from ase.calculators.espresso import Espresso
from ase.visualize import view
PbTe_atoms = ase.io.read("/home/e89/e89/kw2318/PbTe.cif")

import cellconstructor as CC, cellconstructor.Phonons
import cellconstructor.ForceTensor

import numpy as np
import matplotlib.pyplot as plt
import sys, os
import phonopy

def plot_phonopy(iplot,PATH,ph):
        global ax
        if nplot==1:
                axs=[ax]
        else:
                axs=ax

        cell=ph[0].primitive.cell
        nmodes = ph[0].primitive.get_number_of_atoms() * 3
        qpath, data = CC.Methods.get_bandpath(cell, PATH, SPECIAL_POINTS,N_POINTS)
        xaxis, xticks, xlabels = data # Info to plot correclty the x axis
        ph_dispersion=[]
        qpath_frac = [np.dot(qpath[i],np.transpose(cell)) for i in range(len(qpath))]
        for i in range(len(ph)):

                ph[i].run_band_structure(np.reshape(qpath_frac,(1,-1,3)),path_connections=[True,True], labels=xlabels[:2])

                band_q_points, band_distances, band_frequencies, band_eigvecs = ph[i].get_band_structure()

                f=np.reshape(band_frequencies,(-1,nmodes))
                ph_dispersion.append(f)



        if iplot>0:
                axs[iplot].yaxis.set_visible(False)
        for i in range(nmodes):
                        LS=['dashed','solid','solid','solid']
                        if i == 0:
                                for n in range(len(PH)):

                                        axs[iplot].plot(xaxis, ph_dispersion[n][:,i], color = COLOR[n],ls=LS[n], label = PH[n])
                        else:
                                for n in range(len(PH)):
                                        axs[iplot].plot(xaxis, ph_dispersion[n][:,i],ls=LS[n], color = COLOR[n])

        for x in xticks:
                        # ax.axvline(x, 0, 1, color = "k", lw = 0.4)
                        axs[iplot].axhline(0, 0, 1, color = 'k', ls = ':', lw = 0.4)
        # axs[iplot].set_xlim(xticks[0], xticks[1])
        # axs[iplot].set_ylim((-10,40))
        axs[iplot].set_xticks(xticks)
        axs[iplot].set_xticklabels(xlabels)
        axs[iplot].set_ylabel("Phonons [meV]")



cm_mev=1/8.066
PATH=['GMKGALHAMLHK','gmkg']
nplot=len(PATH)

N_POINTS = 300
# Here we define the position of the special points
SPECIAL_POINTS = {'G': np.array([0., 0., 0.]),
 'M': np.array([ 0.5, 0,  0. ]),
 'K': np.array([ 0.33333333, 0.33333333,  0.        ]),
 'A': np.array([0. , 0. , 0.5]),
 'L': np.array([ 0.5, 0,  0.5]),
 'H': np.array([ 0.33333333, 0.33333333,  0.5       ]),
 'g': np.array([0., 0., 0.3333333333]),
  'm': np.array([0.5, 0., 0.3333333333]),
  'k': np.array([0.3333333333, 0.333333333, 0.3333333333])}


PH=['./dfpt1/dfpt-e0.1/3-vasp-dft/','./dfpt1/dfpt-e0.1-fermi/3-vasp-dfpt/']
ph=[]
for dir in PH:
    try:
        ph.append(phonopy.load(f"{dir}phonopy.yaml",force_constants_filename=f'{dir}FORCE_CONSTANTS'))
    except:
        ph.append(phonopy.load(f"{dir}phonopy.yaml",force_sets_filename=f'{dir}FORCE_SETS'))
COLOR = ['k', 'b','r']
LS=['solid','solid','solid']

if nplot==1:
        plt.figure(dpi = 150)
        ax = plt.gca()
        plot_phonopy(0,PATH[0],ph)
        ax.legend()
else:
        fig, ax = plt.subplots(1, nplot,dpi = 150)
        for i in range(nplot):
                plot_phonopy(i,PATH[i],ph)
        ax[0].legend()

plt.tight_layout()
plt.savefig("dispersion.png")
plt.show()
print('finish')

