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
		pseudos = {"Pb": "Pb.upf","Te": "Te.upf"}

		# Now we define the parameters for the espresso calculations
		input_params = {"calculation" : "vc-relax", # The type of calculation
					"ecutwfc" : config.ecutwfc, # The plane-wave wave-function cutoff
					"ecutrho": config.ecutrho, # The density wave-function cutoff,
				"conv_thr": config.conv_thr, # The convergence for the DFT self-consistency
				"pseudo_dir" : "pseudo_espresso", # The directory of the pseudo potentials
				"tprnfor" : True, # Print the forces
				"tstress" : True # Print the stress tensor
				}

		k_spacing = config.k_spacing #A^-1 The minimum distance in the Brillouin zone sampling

		# espresso_calc = Espresso(input_data = input_params, pseudopotentials = pseudos, kspacing = k_spacing)
		self.primitive.write('espresso.pwi',input_data = input_params, pseudopotentials = pseudos, kspacing = k_spacing)

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

