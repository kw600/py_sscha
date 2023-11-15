import ase
from ase.calculators.espresso import Espresso
from ase.visualize import view
# PbTe_atoms = ase.io.read("~/PbTe.cif")

import cellconstructor as CC, cellconstructor.Phonons
import cellconstructor.ForceTensor
import numpy as np
import matplotlib.pyplot as plt
import sys, os, time
import phonopy

def plot_gp(filename,color,LS):
        global ax
        d = np.loadtxt(filename)
        nband = int((len(d[0])-1))
        for i in range(nband):
                if i==0:
                        ax.plot(d[:,0],d[:,i+1],color=color,ls=LS,label=f'{filename}')

                else:
                        ax.plot(d[:,0],d[:,i+1],ls=LS,color=color)

        xticks = [d[0,0],d[100,0],d[200,0],d[-1,0]]
        # ax.set_xlim(xticks[0], xticks[1])
        ax.set_xticks(xticks)
        ax.set_xticklabels([ "$\\Gamma$", "M","K","$\\Gamma$"])
        ax.legend()
        ax.set_ylabel("Phonons [cm-1]")
        # ax.set_ylim(-1,18)
        # plt.show()

def plot_bubble(filename):
        global ax
        d = np.loadtxt(filename,skiprows=3)
        nband = int((len(d[0])-1)/2)
        for i in range(nband):
                if i==0:
                        ax.plot(d[:,0],d[:,i+1]*cm_mev,ls=':',color='k',label='v2')
                        ax.plot(d[:,0],d[:,i+1+nband]*cm_mev,color='r',label='v2+static_bubble')
                else:
                        ax.plot(d[:,0],d[:,i+1]*cm_mev,ls=':',color='k')
                        ax.plot(d[:,0],d[:,i+1+nband]*cm_mev,color='r')
        xticks = [d[0,0],d[-1,0]]
        ax.set_xlim(xticks[0], xticks[1])
        ax.set_xticks(xticks)
        ax.set_xticklabels(["$\\Gamma$", "L"])
        ax.legend()
        ax.set_ylabel("Phonons [meV]");ax.set_ylim(-1,18)
        # plt.show()
# Let us define the PATH in the brilluin zone and the total number of points
def plot(iplot,PATH,sscha_dyn):
        global ax
        if nplot==1:
                axs=[ax]
        else:
                axs=ax
        qpath, data = CC.Methods.get_bandpath(sscha_dyn[0].structure.unit_cell, PATH, SPECIAL_POINTS,N_POINTS)
        xaxis, xticks, xlabels = data # Info to plot correclty the x axis
        print('xaxis',xaxis)
        print('xticks',xticks)
        print('xlabels',xlabels)
        sscha_dispersion = [CC.ForceTensor.get_phonons_in_qpath(sscha_dyn[i], qpath) for i in range(len(sscha_dyn))]
        nmodes = sscha_dyn[0].structure.N_atoms * 3

        if iplot>0:
                axs[iplot].yaxis.set_visible(False)
        for i in range(nmodes):
                        LS=['dashed','dashed','solid'];LS=['solid','solid','solid']
                        if i == 0:
                                for n in range(len(sscha_dyn)):

                                        axs[iplot].plot(xaxis, sscha_dispersion[n][:,i], color = COLOR[n],ls=LS[n], label = SSCHA_DYN[n])
                        else:
                                for n in range(len(sscha_dyn)):
                                        axs[iplot].plot(xaxis, sscha_dispersion[n][:,i],ls=LS[n], color = COLOR[n])

        for x in xticks:
                        # ax.axvline(x, 0, 1, color = "k", lw = 0.4)
                        axs[iplot].axhline(0, 0, 1, color = 'k', ls = ':', lw = 0.4)
        # axs[iplot].set_xlim(xticks[0], xticks[1])
        axs[iplot].set_xticks(xticks)
        axs[iplot].set_xticklabels(xlabels)
        axs[iplot].set_ylabel("Phonons [cm$_{-1}$]")
        # axs[iplot].set_ylim(-1,18)

start=time.time()
cm_mev=1/8.066
PATH=['GMKG']
nplot=len(PATH)

N_POINTS = 200
# Here we define the position of the special points
SPECIAL_POINTS = {"G": [0,0,0],
"X": [0.0,0, 0.5],
"L": [.5, .5, .5],
"W": [.5, .25, .75],
"K": [1/3,1/3,0],
"M": [0.5,0,0]}


read_FORCE_SETS = True
read_FORCE_CONSTANTS = False
read_qe = False

if read_FORCE_SETS:
    p1 = phonopy.load(f"phonopy.yaml",force_sets_filename=f'FORCE_SETS')
    phonon = [p1]
    SSCHA_DYN = ['441']
    sscha_dyn = [CC.Methods.sscha_phonons_from_phonopy(i) for i in phonon]
    for i in range(len(SSCHA_DYN)):
            sscha_dyn[i].save_qe(f"qe_{SSCHA_DYN[i]}_dyn")

if read_FORCE_CONSTANTS:
    p1 = phonopy.load(f"phonopy.yaml",force_constants_filename=f'FORCE_CONSTANTS')
    phonon = [p1]
    SSCHA_DYN = ['441']
    sscha_dyn = [CC.Methods.sscha_phonons_from_phonopy(i) for i in phonon]
    for i in range(len(SSCHA_DYN)):
            sscha_dyn[i].save_qe(f"qe_{SSCHA_DYN[i]}_dyn")
            
if read_qe:
    SSCHA_DYN = ['harmonic_441_dyn']
    sscha_dyn = [CC.Phonons.Phonons(SSCHA_DYN[i], nq[i]) for i in range(len(SSCHA_DYN))]

nq = [4,4]
COLOR = ['k', 'r','b']

if nplot==1:
        plt.figure(dpi = 150)
        ax = plt.gca()
        plot(0,PATH[0],sscha_dyn)
        ax.legend()
else:
        fig, ax = plt.subplots(1, nplot,dpi = 150)
        for i in range(nplot):
                plot(i,PATH[i],sscha_dyn)
        ax[0].legend()

plt.tight_layout()
plt.savefig("dispersion.png")
end=time.time()
plt.title('TiSe2')
print(f'Time: {end-start} seconds')
plt.show()
print('finish')
