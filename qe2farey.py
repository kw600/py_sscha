"""
READ QE results and TRANSFORMED into FAREY FORMAT.
============
This script read the QE dynamical matrix and use q2r.x to get the force constant matrix. The force constant matrix is transformed into force.dat using qe_to_f90force.py and the corresponding R information is generated by generate_R.py. 
"""

import os,shutil

input_dir = './332_50K'
input_file = 'ScVSn.332.50.1000.dyn'

# generate delta_prim.dat using generate_R.py
shutil.copy(input_dir+'/'+input_file+'1', input_dir+"/harmonic1")
shutil.copy(input_dir+'/'+input_file+'0', input_dir+"/harmonic0")
os.system(f"python generate_R.py {input_dir}")


# use q2r.x to get the force constant matrix
q2r=f'''
&input
    fildyn = './{input_dir}/{input_file}'
    zasr = 'simple'
    flfrc = '{input_dir}/format2.txt'
/
'''
with open('q2r.in', 'w') as f:
    f.write(q2r)
os.system("q2r.x < q2r.in > q2r.out")

# transform the force constant matrix into force.dat
os.system(f"python qe_to_f90force.py {input_dir}")



# Need to the farey grids for the target supercell and the results can be transformed back to QE format using the following script
# os.system(f"python farey2qe.py {output_dir}")




