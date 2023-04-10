import ase
from ase.calculators.espresso import Espresso
from ase.visualize import view

PbTe_atoms = ase.io.read("PbTe.cif")

import initial_density_matrix
from ensemble import *
from dft_ensemble import *
from minization import *
from init_structure import *
import config as cfg



def sscha():
	PbTe_atoms = structure("PbTe.cif")
	PbTe_atoms.generate_input_for_initial_relax()
	PbTe_atoms.generate_input_for_initial_phonon()

if __name__ == "__main__":
    sscha()


