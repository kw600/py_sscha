import numpy as np
import sys
def read_poscar(file='POSCAR'):
    d = open(file,'r').readlines()
    num_atom = np.array(d[6].split(),dtype=int)
    atom_list = []
    lattice = np.array([i.split() for i in d[2:5]],dtype=float)
    frac =  np.array([i.split() for i in d[8:8+np.sum(num_atom)]],dtype=float)
    atom = d[5].split()
    num_atom = np.array(d[6].split(),dtype=int)
    for i in range(len(atom)):
        atom_list += [atom[i]]*num_atom[i]
    return atom_list, lattice, frac

def get_bond(atomA,atomB,file='POSCAR'):
    bond=[]
    atom_list, lattice, frac = read_poscar(file)
    for i in range(len(atom_list)):
        if (atom_list[i] == atomA) or (atom_list[i] in atomA):
            for j in range(i,len(atom_list)):
                if (atom_list[j] == atomB) or (atom_list[j] in atomB):
                    l = np.dot(frac[i]-frac[j],lattice)
                    # remember to change the conditions
                    if (np.linalg.norm(frac[i,:2]-frac[j,:2])<0.001) and (np.linalg.norm(l)<4):
                        print(i,j)
                        bond.append(np.linalg.norm(l))

    return bond

path = sys.argv[1]
bond = get_bond('Y','Sn',file=path)
print(bond)
