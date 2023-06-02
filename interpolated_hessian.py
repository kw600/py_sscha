import ase
from ase.calculators.espresso import Espresso
from ase.visualize import view

PbTe_atoms = ase.io.read("PbTe.cif")

import cellconstructor as CC
import cellconstructor.ForceTensor
import cellconstructor.Structure
import cellconstructor.Spectral

import numpy as np

FINAL_DYN =   'dyn_pop3_'
NQIRR = 3

dyn = CC.Phonons.Phonons(FINAL_DYN,NQIRR)
supercell = dyn.GetSupercell()
tensor3 = CC.ForceTensor.Tensor3(dyn.structure,
                                dyn.structure.generate_supercell(supercell),
                                supercell)
tensor3.SetupFromFile(fname="FC3",file_format='D3Q')
k_grid=[4,4,4]

PATH = "GL"
N_POINTS = 10
SPECIAL_POINTS = {"G": [0,0,0],
"X": [0.5,0, 0.5],
"L": [.5, .5, .5],
"W": [.5, .25, .75],
"M": [0,0.5,0.5]}
SSCHA_DYN2 = 'harmonic_222_dyn'
sscha_dyn2 = CC.Phonons.Phonons(SSCHA_DYN2, 3)
qpath, data = CC.Methods.get_bandpath(sscha_dyn2.structure.unit_cell, PATH, SPECIAL_POINTS,N_POINTS)

CC.Spectral.get_static_correction_along_path(dyn=dyn, 
                                             tensor3=tensor3, 
                                             k_grid=k_grid, 
                                             q_path=qpath, 
                                             filename_st="v2_v2+d3static_freq.dat",
                                             T =100.0,
                                             print_dyn = False)

