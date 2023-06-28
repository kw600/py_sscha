"""
Fourier interpolate QE results to a different supercell size.
============
This script read the QE dynamical matrix and use q2r.x to get the force constant matrix. The force constant matrix is transformed into force.dat using qe_to_f90force.py and generate_R.py. The final QE matrix with different supercell size is obtained by farey2qe.py
"""

import os,shutil

input_file = 'dyn_pop2_'
output_file = 'qe2qe_881_dyn'

# generate delta_prim.dat using generate_R.py
shutil.copy(input_file+'1', "harmonic1")
shutil.copy(input_file+'0', "harmonic0")
os.system("python generate_R.py ./")

# use q2r.x to get the force constant matrix
q2r=f'''
&input
    fildyn = '{input_file}'
    zasr = 'simple'
    flfrc = 'format2.txt'
/
'''
with open('q2r.in', 'w') as f:
    f.write(q2r)
os.system("q2r.x < q2r.in > q2r.out")

# transform the force constant matrix into force.dat
os.system("python qe_to_f90force.py")

# generate the final QE matrix with different supercell size
# target supercell size needs to be consistent the harmonic_xxx_dyn files in the directory
os.system("python farey2qe.py ./")




