import os,sys
import numpy as np

qe = open(sys.argv[1],'r').readlines()
for i in range(len(qe)):
    if 'nat' in qe[i]:
        nat = int(qe[i].split()[-1])
    elif 'ntyp' in qe[i]:
        ntyp = int(qe[i].split()[-1])
    elif 'CELL' in qe[i]:
        lattice = np.zeros((3,3),dtype=float)
        for j in range(3):
            lattice[j] = [float(x) for x in qe[i+j+1].split()]
    elif 'ATOMIC_POSITIONS' in qe[i]:
        element=[]
        for j in range(nat):
            atom = qe[i+j+1].split()[0]
            if atom not in element:
                element.append(atom)
        pos = np.zeros((nat,3),dtype=float)
        pos = []
        for j in range(ntyp):
            pos.append([])
        for j in range(nat):
            for k in range(len(element)):
                if qe[i+j+1].startswith(element[k]):
                    pos[k].append([float(x) for x in qe[i+j+1].split()[1:]])

#check the number of supercll and it may needed to be changed.
n_cell = np.gcd.reduce([len(x) for x in pos])

with open("POSCAR",'w') as f:
    f.write(f"supercell {n_cell}\n")
    f.write("1.0\n")
    for i in range(3):
        f.write(" ".join([str(j) for j in lattice[i]])+"\n")
    f.write(" ".join(element)+"\n")
    f.write(" ".join([str(len(x)) for x in pos])+"\n")
    f.write("Cartesian\n")
    for i in range(len(element)):
        for j in range(len(pos[i])):
            p=pos[i][j]
            f.write("%.14f %.14f %.14f \n" % (p[0],p[1],p[2]) )
