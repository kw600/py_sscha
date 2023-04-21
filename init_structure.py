import ase
from ase.calculators.espresso import Espresso
from ase.visualize import view

import cellconstructor as CC
import cellconstructor.Structure
import cellconstructor.Phonons

import config

def init_structure(filename):
	stru = ase.io.read(filename)
	return stru

def get_primitive(unitcell):
	struct = CC.Structure.Structure()
	struct.generate_from_ase_atoms(unitcell)

	# Define the new unit cell
	new_cell = struct.unit_cell.copy()
	new_cell[0,:] = .5 * struct.unit_cell[0,:] + .5*struct.unit_cell[1,:]
	new_cell[1,:] = .5 * struct.unit_cell[0,:] - .5*struct.unit_cell[1,:]
	new_cell[2,:] = .5 * struct.unit_cell[0,:] + .5*struct.unit_cell[2,:]

	# Apply the new unit cell to the structure
	# And remove duplicated atoms
	struct.unit_cell = new_cell
	struct.fix_coords_in_unit_cell()
	primitive = struct.get_ase_atoms()
	return primitive

class structure():
	def __init__(self, filename):
		
		self.unitcell = init_structure(filename)
		self.primitive = get_primitive(self.unitcell)
		
	def generate_input_for_initial_relax(self):
		pw_input = f"""
&CONTROL
   calculation      = 'vc-relax'
   tstress          = .true.
   tprnfor          = .true.
   pseudo_dir       = '/work/e89/e89/kw2318/pseudo_espresso'
/
&SYSTEM
   ecutwfc          = 65
   ntyp             = 2
   nat              = 2
   ibrav            = 0
/
&ELECTRONS
   conv_thr         = 1e-08
/
&IONS
/
&CELL
/

ATOMIC_SPECIES
Pb 207.2 Pb_ONCV_PBE-1.2.upf
Te 127.6 Te_ONCV_PBE-1.2.upf

K_POINTS automatic
8 8 8  0 0 0

CELL_PARAMETERS angstrom
3.275 3.275 0
3.275 -3.275 0
3.275 0 3.275

ATOMIC_POSITIONS angstrom
Pb 0.0000000000 0.0000000000 0.0000000000
Te 3.27500000000 0.0000000000 0.0000000000
"""
		with open("espresso.phi", "w") as f:
				f.write(pw_input)
	def generate_input_for_initial_phonon(self):
		ph_input = f"""
&inputph
	! the final filename
	fildyn = "harmonic_dyn"
	
	! the q mesh
	ldisp = .true.
	nq1 = {config.nq1} 
	nq2 = {config.nq2} 
	nq3 = {config.nq3} 
	
	! compute also the effective charges and the dielectric tensor
	epsil = .true.
&end
					"""
					# We write the input script and execute the phonon calculation program
		with open("harmonic.phi", "w") as f:
						f.write(ph_input)

