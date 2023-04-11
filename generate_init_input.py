import ase
from ase.calculators.espresso import Espresso
from ase.visualize import view

PbTe_atoms = ase.io.read("PbTe.cif")

from init_structure import *
import config 



def sscha():
	PbTe_atoms = structure("PbTe.cif")
	PbTe_atoms.generate_input_for_initial_relax()
	PbTe_atoms.generate_input_for_initial_phonon()

if __name__ == "__main__":
    sscha()


