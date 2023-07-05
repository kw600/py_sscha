import cellconstructor as CC
import cellconstructor.Structure
import cellconstructor.Phonons
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
import sys,os
import numpy as np



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
    
    for i,factor in enumerate(Factor):
        tmp = super_structure.copy()
        tmp.coords += dis*factor
        # print("%s/scf_population%d_%d.dat" % (f'ens{pop}', pop, i+1))
        tmp.save_scf("%s/scf_population%d_%d.dat" % (f'ens{pop}', pop, i+1))
    
if __name__=="__main__":
    pop = 1
    harmonic_file = 'f2q_'
    nqirr = 6
    save_scf(harmonic_file,nqirr,[i*0.1 for i in range(11)])
    generate_vasp_input(pop)
