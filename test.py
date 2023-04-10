import ase
from ase.calculators.espresso import Espresso
from ase.visualize import view
PbTe_atoms = ase.io.read("PbTe.cif")

# We can view the structure
view(PbTe_atoms)
