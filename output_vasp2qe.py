import sys
import numpy as np


def read_outcar(label):
    bohr = 0.5291772488820865

    POSCAR = open('POSCAR','r').readlines()
    lattice = np.zeros((3,3),dtype=float)
    for i in range(3):
        lattice[i] = [float(x) for x in POSCAR[i+2].split()]
    lattice /= bohr #convert to bohr
    alat = float(lattice[0,0])


    n_cell = int(POSCAR[0].split()[-1])
    natom_per_type = np.array((POSCAR[6].split()),dtype=int)

    pos = np.zeros((natom_per_type.sum(),3),dtype=float)
    for i in range(natom_per_type.sum()):
        pos[i] = [float(x) for x in POSCAR[8+i].split()[0:3]]
    pos /= bohr #convert to bohr

    element = []
    for i in range(len(natom_per_type)):
        element+=natom_per_type[i]*[POSCAR[5].split()[i]]
    n_tot = natom_per_type.sum()
    print('Please confirm the number of cells is ',n_cell)

    #read the OUTCAR
    d = open('OUTCAR','r').readlines()
    if 'Voluntary' not in d[-1]:
        print('Unfinished job')
        exit()
    #initialize the stress tensor
    stress=np.zeros((3,3))

    for i in range(len(d)):
        if 'energy(sigma->0)' in d[i]: #read the total energy
            energy = float(d[i].split()[-1])/(13.605693012183622) #convert eV to Ry
        elif 'in kB' in d[i]: #read the stress tensor
            s = np.array(d[i].split()[2:],dtype=float)
            stress[0,0]=s[0];stress[1,1]=s[1];stress[2,2]=s[2];stress[0,1]=s[3];stress[1,0]=s[3]
            stress[1,2]=s[4];stress[2,1]=s[4];stress[0,2]=s[5];stress[2,0]=s[5]
        elif 'TOTAL-FORCE' in d[i]:#read the force
            force = np.loadtxt(d[i+2:i+2+n_tot])[:,3:]/25.71103168347908 #convert eV/ang to Ry/bohr
            # pos = np.loadtxt(d[i+2:i+2+n_tot])[:,0:3]/bohr #convert ang to bohr
    stress1 = stress/(29421.02648438959/2*10) #convert to Ry/bohr^3   

    index = [i for i in range(n_tot)]
    I =[]
    start = 0
    end = natom_per_type[0]
    # print(natom_per_type)
    # print(start,end,'start end')
    for i in range(len(natom_per_type)):
        I.append(index[start:end])
        if i == len(natom_per_type)-1:
            break
        start = end
        end += natom_per_type[i+1]
    
    # print(I)

    new_index = []
    for cell in range(n_cell):
        for type in range(len(natom_per_type)):
            i = I[type]
            for n in range(int(natom_per_type[type]/n_cell)):
                # print(f'atom type {type}',int(natom_per_type[type]/n_cell))
                # print(cell,type,n)
                a = i.pop(0)
                new_index.append(a)
    # print(new_index)
    with open(f'../{label}.pwo','w') as f:
        f.write(f'number of atoms/cell      =           {n_tot}\n')
        f.write(f'celldm(1)=   {alat:.6f}  celldm(2)=   0.000000  celldm(3)=   0.000000\n')
        f.write('celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000\n\n')
        f.write('     crystal axes: (cart. coord. in units of alat)\n')
        for i in range(3):
            f.write(f'               a({i+1}) = (   {lattice[i,0]/alat:.16f}  {lattice[i,1]/alat:.16f}  {lattice[i,2]/alat:.16f} )\n')

        f.write('   Cartesian axes\n\n')
        f.write('     site n.     atom                  positions (alat units)\n')
        for i in range(n_tot):
            f.write(f'         {i+1}           {element[new_index[i]]}   tau(   {i+1}) = (   {pos[new_index[i],0]/alat:.16f}  {pos[new_index[i],1]/alat:.16f}   {pos[new_index[i],2]/alat:.16f}  )\n')

        f.write(f'!    total energy              =   {energy} Ry\n\n')

        f.write('     Forces acting on atoms (cartesian axes, Ry/au):\n\n')
        for i in range(n_tot):
            f.write(f'     atom {i+1} type {element[new_index[i]]}   force = {force[new_index[i],0]:.16f} {force[new_index[i],1]:.16f} {force[new_index[i],2]:.16f}\n')

        f.write('\n\n')
        f.write('total   stress  (Ry/bohr**3)                   (kbar)\n')
        for j in range(3):
            f.write(f'  {stress1[j,0]:.8f}    {stress1[j,1]:.8f}    {stress1[j,2]:.8f}        {stress[j,0]:.8f}    {stress[j,1]:.8f}    {stress[j,2]:.8f}\n')
        
        f.write('JOB DONE')
    return label, None
read_outcar(sys.argv[1])
