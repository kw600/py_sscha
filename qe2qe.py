"""
Fourier interpolate QE results to a different supercell size.
============
This script read the QE dynamical matrix and use q2r.x to get the force constant matrix. The force constant matrix is transformed into force.dat using qe_to_f90force.py and generate_R.py. The final QE matrix with different supercell size is obtained by farey2qe.py
"""

import os,shutil

input_file = 'dyn_pop2_'
output_dir = 'lte_881'
output_file = 'qe2qe_881_dyn'

# generate delta_prim.dat using generate_R.py
shutil.copy(input_file+'1', output_dir+"harmonic1")
shutil.copy(input_file+'0', output_dir+"harmonic0")
os.system(f"python generate_R.py {output_dir}")

# use q2r.x to get the force constant matrix
q2r=f'''
&input
    fildyn = '{input_file}'
    zasr = 'simple'
    flfrc = '{output_dir}/format2.txt'
/
'''
with open('q2r.in', 'w') as f:
    f.write(q2r)
os.system("q2r.x < q2r.in > q2r.out")

# transform the force constant matrix into force.dat
os.system(f"python qe_to_f90force.py {output_dir}")

# generate the final QE matrix with different supercell size
# target supercell size needs to be consistent the harmonic_xxx_dyn files in the directory
os.system(f"python farey2qe.py {output_dir}")




