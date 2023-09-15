import cellconstructor as CC
import cellconstructor.Structure
import cellconstructor.Phonons
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
import sys,os
import numpy as np
import phonopy

np.set_printoptions(precision=20)

def generate_vasp_input(pop):
    """
    READ the scf_populations files
    """
    global harmonic_file, nqirr
    all_scf_files = [os.path.join(f"ens{pop}", f) for f in os.listdir(f"ens{pop}") if f.startswith("scf_")]
    if not os.path.exists(f"run_dft{pop}"):
        os.mkdir(f"run_dft{pop}")

    for ii,file in enumerate(all_scf_files):
        number = int(file.split("_")[-1].split(".")[0])
        if not os.path.exists(f"run_dft{pop}/vasp{number}"):
            os.mkdir(f"run_dft{pop}/vasp{number}")
        filename = os.path.join(f"run_dft{pop}/vasp{number}", "POSCAR".format(number))
        
        # We start writing the file
        with open(filename, "w") as f:
            
            # Load the scf_population_X.dat file
            ff = open(file, "r")
            structure_lines = ff.readlines()
            lattice = np.array([i.split() for i in structure_lines[1:4]],dtype=float)
            lattice = np.transpose(lattice)
            atoms   = structure_lines[6:]
            sc_lines = []
            v_lines = []
            sn_lines = []
            for line in atoms:
                if line.startswith('Sc'):
                    sc_lines.append(line)
                elif line.startswith('V'):
                    v_lines.append(line)
                elif line.startswith('Sn'):
                    sn_lines.append(line)
            ff.close()
            f.writelines('ScV6Sn6 3 by 3 by 3\n   1.00000000000000\n')
            f.writelines(structure_lines[1:4])
            dyn = CC.Phonons.Phonons(harmonic_file, nqirr = nqirr)
            [n1,n2,n3]=dyn.GetSupercell()
            factor=int(n1*n2*n3)
            f.writelines(f'   Sc   V    Sn\n     {factor}     {6*factor}     {6*factor}\nDirect\n')
            for i in sc_lines:
                coord = np.array(i.split()[1:],dtype=float)
                frac = np.dot(np.linalg.inv(lattice),coord)
                f.write(f'{frac[0]}  {frac[1]}  {frac[2]}\n')
            for i in v_lines:
                coord = np.array(i.split()[1:],dtype=float)
                frac = np.dot(np.linalg.inv(lattice),coord)
                f.write(f'{frac[0]}  {frac[1]}  {frac[2]}\n')
            for i in sn_lines:
                coord = np.array(i.split()[1:],dtype=float)
                frac = np.dot(np.linalg.inv(lattice),coord)
                f.write(f'{frac[0]}  {frac[1]}  {frac[2]}\n')

def save_scf(harmonic_file,nqirr,Factor,nmode=0):
    global pop
    if not os.path.exists(f'ens{pop}'):
        os.mkdir(f'ens{pop}')

    dyn = CC.Phonons.Phonons(harmonic_file, nqirr = nqirr)
    super_structure = dyn.structure.generate_supercell(dyn.GetSupercell())
    super = dyn.GenerateSupercellDyn(dyn.GetSupercell())
    freqs, pol_vects = super.DyagDinQ(0)
    print(freqs[:10])
    # Compute the displacemets from the polarization vectors
    _m_ = super.structure.get_masses_array()
    _m_ = np.tile(_m_, (3,1)).T.ravel()

    # Compute the atomic displacements
    atomic_disp = np.einsum("ab, a -> ab", pol_vects, 1 / np.sqrt(_m_) )
    # Normalize the displacements
    atomic_disp[:,:] /= np.tile( np.sqrt(np.sum(np.abs(atomic_disp)**2, axis = 0)), (super.structure.N_atoms * 3, 1))
    dis = atomic_disp[:,nmode]
    dis = dis.reshape((-1,3))
    np.save('333.dis',dis)
    for i,factor in enumerate(Factor):
        tmp = super_structure.copy()
        tmp.coords += dis*factor
        # print("%s/scf_population%d_%d.dat" % (f'ens{pop}', pop, i+1))
        tmp.save_scf("%s/scf_population%d_%d.dat" % (f'ens{pop}', pop, i+1))

def quadratic_from_phonopy(path,nmode=0,factor=1):
    #obtain the eigenvector of the supercell
    p=phonopy.load(f"./{path}/phonopy.yaml",force_constants_filename=f'./{path}/FORCE_CONSTANTS')
    p.run_band_structure([[[0,0,0],[0,0,0.5]]],path_connections=[True], labels=['q','G'],with_eigenvectors=True)
    band_q_points, band_distances, band_frequencies, band_eigvecs = p.get_band_structure()
    eig=np.reshape(band_eigvecs[0][0][:,nmode],(-1,3))
    eig=np.real(eig)
    atom = p.get_unitcell().get_chemical_symbols()
    pos=np.array(p.get_unitcell().positions)
    lattice=p.get_unitcell().get_cell()
    mass=p.get_unitcell().get_masses()
    with open(f'POSCAR{factor}','w') as f1:
        f1.write('''python r332
1.0
8.20039500000000032287 4.73450026071126206517 0.00000000000000000000
0.00000000000000000000 9.46900052142252413034 0.00000000000000000000
0.00000000000000000000 0.00000000000000000000 18.31879999999999952820
Sc V Sn
  6  36  36
Direct
''')
        for i in range(len(pos)):
            # print(pos[i])
            # print(eig[i])
            # print(mass[i])
            # # print()
            new_p = pos[i] +factor*eig[i]/np.sqrt(mass[i])
            # print(new_p,'new')
            # print(np.linalg.inv(lattice))
            new_frac = np.dot(new_p,np.linalg.inv(lattice))
            # print(new_p[i],np.linalg.inv(lattice))
            # print(new_frac)
            f1.write("  %10.20f %10.20f %10.20f \n" % (new_frac[0],new_frac[1],new_frac[2]))

def quadratic_from_dis(output_file,factor=1):
    dis = np.load('dis1.5.npy')
    # print(np.max(dis))
    d = open('POSCAR','r').readlines()
    with open(output_file,'w') as f1:
        for i in range(8):
            f1.write(d[i])
        lattice = np.array([d[i].split() for i in range(2,5)],dtype=float)
        # print(lattice)
        frac = np.loadtxt('POSCAR',skiprows=8)
        for i in range(len(frac)):
            new_frac = frac[i] + np.dot(dis[i]*factor,np.linalg.inv(lattice))
            # print(new_frac)
            f1.write("  %0.20f %0.20f %0.20f \n" % (new_frac[0],new_frac[1],new_frac[2]))

def mannually_displace(prim_cell,prim_pos,eig,supercell=[3,3,3]):
    num_cell = supercell[0]*supercell[1]*supercell[2]
    cell_id = np.zeros((supercell[0],supercell[1],supercell[2]))
    sup_atom = np.zeros((num_cell,len(prim_pos)))
    sup_mass = np.zeros((num_cell,len(prim_pos)))
    sup_cell = np.array([prim_cell[i]*supercell[i] for i in range(3)])
    atom=['Sc']+['V']*6+['Sn']*6
    id = 0
    for i in range(supercell[0]):
        for j in range(supercell[1]):
            for k in range(supercell[2]):
                cell_id[i,j,k] = id
                id += 1
    with open('p2','w') as f1:
        for cell in range(num_cell):
            [i,j,k] = np.where(cell_id==cell)
            print(cell,i,j,k)
            sup_pos = prim_pos+i*prim_cell[0]+j*prim_cell[1]+k*prim_cell[2]
            dis = eig*np.exp(1j*2*np.pi*(i/supercell[0]+j/supercell[1]+k/supercell[2]))
            print(dis)
            p = sup_pos + np.real(dis)
            for l in range(len(sup_pos)):
                    f1.write("%s %10.20f %10.20f %10.20f \n" % (atom[l],p[l,0],p[l,1],p[l,2]))

def generate_castep_ph():
    num_cell = supercell[0]*supercell[1]*supercell[2]
    cell_id = np.zeros((supercell[0],supercell[1],supercell[2]))
    sup_atom = np.zeros((num_cell,len(prim_pos)))
    sup_mass = np.zeros((num_cell,len(prim_pos)))
    sup_cell = np.array([prim_cell[i]*supercell[i] for i in range(3)])
    atom=['Sc']+['V']*6+['Sn']*6
    mass=[44.95591]+[50.9414]*6+[118.70999]*6
    id = 0
    for i in range(supercell[0]):
        for j in range(supercell[1]):
            for k in range(supercell[2]):
                cell_id[i,j,k] = id
                id += 1
    header=''' BEGIN header
 Number of ions         351
 Number of branches     1053
 Number of wavevectors  1
 Frequencies in         cm-1
 IR intensities in      (D/A)**2/amu
 Raman activities in    A**4 amu**(-1)
 Unit cell vectors (A)
       16.40079111   0.            0.        
      -8.20039555    14.20350175   0.        
       0.            0.            27.47820186
 Fractional Co-ordinates
 '''
    header2=''' END header
     q-pt=    1   0 0 0      1
       1      -22
       2      67.917047
       3     179.176230
                        Phonon Eigenvectors
Mode Ion                X                                   Y                                   Z
'''
    with open('p2.phonon','w') as f:
        f.write(header)
        count=1
        for cell in range(num_cell):
            [i,j,k] = np.where(cell_id==cell)
            sup_pos = prim_pos+i*prim_cell[0]+j*prim_cell[1]+k*prim_cell[2]
            sup_f = np.dot(sup_pos,np.linalg.inv(sup_cell))
            for l in range(len(sup_pos)):
                f.write(f'{count} {sup_f[l,0]} {sup_f[l,1]} {sup_f[l,2]} {atom[l]} {mass[l]}\n')
                count+=1
    
        f.write(header2)
        count=1
        for cell in range(num_cell):
            [i,j,k] = np.where(cell_id==cell)
            dis = eig*np.exp(1j*2*np.pi*(i/supercell[0]+j/supercell[1]+k/supercell[2]))
            for l in range(len(sup_pos)):
                f.write("1 %i %.4f %.4f %.4f %.4f %.4f %.4f\n" % (count,dis[l,0],0,dis[l,1],0,dis[l,2],0))
                count+=1


if __name__=="__main__":
    for i in range(0,6):
        factor = i
        output_file = f'displace/vasp{factor}/POSCAR'
        quadratic_from_dis(output_file,factor)

