import numpy as np
from ase.io import read, write
np.set_printoptions(precision=20)
#change basis to new_cell
s = read('SPOSCAR')
cell = s.get_cell()
pos = s.get_positions()
frac = s.get_scaled_positions()
atom = s.get_chemical_symbols()

a = cell[0,0];c=cell[2,2]
new_cell = np.array([[a/2,a/2/np.sqrt(3),0],[0,a/np.sqrt(3),0],[0,0,c]])
new_frac = np.dot(pos,np.linalg.inv(new_cell))

new_s = []
for i in range(len(new_frac)):
    [x,y,z] = new_frac[i]
    if x >= 0.9999 or y >= 0.9999 or z >= 0.9999:
        # print(new_frac[i],'outside')
        pass
    elif x < -0.0001 or y < -0.0001 or z < -0.0001:
        # print(new_frac[i],'outside')
        pass
    else:
        # print(new_frac[i],'inside')
        new_s.append(new_frac[i])
print(len(new_s))

with open('new_POSCAR','w') as f:
    for i in range(3):
        f.write('%.20f %.20f %.20f\n'%(new_cell[i,0],new_cell[i,1],new_cell[i,2]))
    for i in new_s:
        f.write('%.20f %.20f %.20f\n'%(i[0],i[1],i[2]))


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




def phonopy_to_vesta(eig,output,factor=1):
    eig[:6] = eig[:6]*1.6249883295498455
    eig=eig*factor
    with open(output,'w') as f:
        f.write('VECTR\n')
        index = 1
        for i in eig:
            f.write('   %i    %0.6f    %0.6f    %0.6f 0\n'%(index,i[0],i[1],i[2]))
            f.write(f'    {index}   0    0    0    0\n')
            f.write(' 0 0 0 0 0\n')
            index += 1
        f.write(' 0 0 0 0 0\n')
        f.write('VECTT\n')
        index = 1
        for i in range(len(eig)):
            f.write(f'   {index}  0.500 255   0   0 0\n')
            index += 1
        f.write(' 0 0 0 0 0\n')