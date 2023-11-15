import numpy as np
import sys

file = sys.argv[1]

#read the poscar
d = open(file,'r').readlines()
lattice = np.zeros((3,3))
for i in range(3):
    lattice[i]=np.array(d[2+i].split(),dtype=float)
atom = np.array(d[6].split(),dtype=int)
id = []
for i in range(len(atom)):
    id += [i+1]*atom[i]
natom = np.sum(atom)
frac = np.zeros((natom,3))
for i in range(natom):
    frac[i]=np.array(d[8+i].split(),dtype=float)

#genreate the qe output
with open('qe.dat','w') as f1:
    f1.write(f'''Dynamical matrix file
File generated with the CellConstructor by Lorenzo Monacelli
{len(atom)} {natom} 0     1.8897259889999998     0.0000000000000000     0.0000000000000000     0.0000000000000000     0.0000000000000000     0.0000000000000000
Basis vectors\n''')
    for i in range(3):
        f1.write(f'     {lattice[i][0]:>14.9f}{lattice[i][1]:>14.9f}{lattice[i][2]:>14.9f}\n')
    f1.write('''        1  'Sc '  40974.8071860994023154
        2  'V '  46430.3369103196018841
        3  'Sn '  108197.5460994290042436\n''')
    for i in range(natom):
        pos = np.dot(frac[i],lattice)
        f1.write(f'   {i+1}    {id[i]}     {pos[0]:>14.9f}{pos[1]:>14.9f}{pos[2]:>14.9f}\n')