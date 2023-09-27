import numpy as np
from ase.io import read, write
np.set_printoptions(precision=20)


input_file = 'POSCAR'
output_file = 'new_POSCAR'

new_cell = '''     9.4829268542312395    0.0000000000000000    0.0000000000000000
    -4.7414634271156197    8.2124555579939091    0.0000000000000000
     0.0000000000000000    0.0000000000000000   18.3532012418615906'''

new_cell = np.array([i.split() for i in new_cell.split('\n')],dtype=float)

#use ase to read the structural information
s = read(input_file)
cell = s.get_cell()
pos = s.get_positions()
frac = s.get_scaled_positions()
atom = s.get_chemical_symbols()

#mannually define the new cell and then change basis to new_cell

new_frac = np.dot(pos,np.linalg.inv(new_cell))


def shift_frac(frac):
    for i in range(len(frac)):
            if frac[i] < -0.0001:
                frac[i] += 1
            elif frac[i] >= 0.9999:
                frac[i] -= 1

    for i in range(len(frac)):
            if frac[i] < -0.0001 or frac[i] > 0.9999:
                frac = shift_frac(frac)
    return frac

new_s = []
for i in range(len(new_frac)):
    f = shift_frac(new_frac[i])
    new_s.append(f)
print('Make sure the atoms in the new cell is',len(new_s))

old_poscar = open(input_file,'r').readlines()
with open(output_file,'w') as f:
    f.write(old_poscar[0])
    f.write(old_poscar[1])
    for i in range(3):
        f.write('%.20f %.20f %.20f\n'%(new_cell[i,0],new_cell[i,1],new_cell[i,2]))
    f.write(old_poscar[5])
    f.write(old_poscar[6])
    for i in new_s:
        f.write('%.20f %.20f %.20f\n'%(i[0],i[1],i[2]))



exit()
#shift all the atoms up
with open('new_POSCAR1','w') as f:
    for i in new_s:
        z = i[2] - 0.3397507*2/3
        if z < -0.0001:
            z+=1
        f.write('%.20f %.20f %.20f\n'%(i[0],i[1],z))

def poscar_to_qe():
    d = open('POSCAR','r').readlines()
    num_atom = np.array(d[6].split(),dtype=int)
    atom_list = [] 
    lattice = np.array([i.split() for i in d[2:5]],dtype=float)
    frac =  np.array([i.split() for i in d[8:8+np.sum(num_atom)]],dtype=float)
    atom = d[5].split()

    for i in range(len(atom)):
        atom_list += [atom[i]]*num_atom[i]
    with open('espresso.in','w') as f:
        for i in range(len(frac)):
            pos = np.dot(frac[i],lattice)
            f.write('%s %0.15f %0.15f %0.15f\n' %(atom_list[i],pos[0],pos[1],pos[2]))




