import ase
from ase.calculators.espresso import Espresso
from ase.visualize import view

PbTe_atoms = ase.io.read("PbTe.cif")

import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.ForceTensor
import cellconstructor.Structure
import cellconstructor.Spectral

# Import the SCHA modules
import sscha, sscha.Ensemble

import numpy as np


# Here the input information
DATA_DIR = 'ens3'# path to the directory ens_pop#lastpop where the last population is stored
N_RANDOM = 512# number elements in the ensamble
DYN_PREFIX =  'dyn_pop2_'  # dyn mat that generated the last population
FINAL_DYN =   'dyn_pop3_'    # SSCHA dyn mat obtained with the last minimization 
SAVE_PREFIX = 'SnTe.Hessian.dyn'  # Free energy Hessian dynamical matrices
NQIRR = 3
Tg = 100
T =  100
POPULATION = 3# number of last population
INCLUDE_V4 = False # True to include the 4th-order SSCHA FC term to calculate the Hessian 


dyn = CC.Phonons.Phonons(DYN_PREFIX, NQIRR)
final_dyn = CC.Phonons.Phonons(FINAL_DYN, NQIRR)
ens = sscha.Ensemble.Ensemble(dyn, Tg, dyn.GetSupercell())
ens.load(DATA_DIR, POPULATION, N_RANDOM)
ens.update_weights(final_dyn, T)
dyn_hessian = ens.get_free_energy_hessian(include_v4 = INCLUDE_V4,
                                          get_full_hessian = True,
                                          verbose = True)
print("Saving the hessian to {}...".format(SAVE_PREFIX))
dyn_hessian.save_qe(SAVE_PREFIX)
print("Done.")

dyn = CC.Phonons.Phonons(FINAL_DYN,NQIRR) 
supercell = dyn.GetSupercell()
tensor3 = CC.ForceTensor.Tensor3(dyn.structure,
                                dyn.structure.generate_supercell(supercell),
                                supercell)

# Assign the tensor3 values
d3 = np.load("d3_realspace_sym.npy")*2.0
tensor3.SetupFromTensor(d3)

# Center and apply ASR
tensor3.Center()
tensor3.Apply_ASR()

# Print it
tensor3.WriteOnFile(fname="FC3",file_format='D3Q')