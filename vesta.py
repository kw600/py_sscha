import cellconstructor as CC
import cellconstructor.Structure
import cellconstructor.Phonons
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
import sys,os
import numpy as np
import phonopy

def write_dis_vesta(displacement,f):
    # with open(output_file,'a') as f:
        f.write('VECTR\n')
        index = 1
        for i in displacement:
            f.write('   %i    %0.6f    %0.6f    %0.6f 0\n'%(index,i[0],i[1],i[2]))
            f.write(f'    {index}   0    0    0    0\n')
            f.write(' 0 0 0 0 0\n')
            index += 1
        f.write(' 0 0 0 0 0\n')
        f.write('VECTT\n')
        index = 1
        for i in range(len(displacement)):
            f.write(f'   {index}  0.500 255   0   0 0\n')
            index += 1
        f.write(' 0 0 0 0 0\n')

def read_eig(nmode):
    p=phonopy.load(f"phonopy.yaml",force_constants_filename=f'FORCE_CONSTANTS')
    p.run_band_structure([[[0,0,0],[0,0,0.5]]],path_connections=[True], labels=['q','G'],with_eigenvectors=True)
    band_q_points, band_distances, band_frequencies, band_eigvecs = p.get_band_structure()
    print(band_frequencies[0][0][:10])
    eig=np.reshape(band_eigvecs[0][0][:,nmode],(-1,3))
    eig=np.real(eig)
    atom = p.get_unitcell().get_chemical_symbols()
    pos=np.array(p.get_unitcell().positions)
    lattice=p.get_unitcell().get_cell()
    mass=p.get_unitcell().get_masses()
    return eig,mass

def generate_supercelldyn():
    try:
        p=phonopy.load(f"phonopy.yaml",force_constants_filename=f'FORCE_CONSTANTS')
    except:
        p=phonopy.load(f"phonopy.yaml",force_constants_filename=f'FORCE_SETS')
    s = CC.Methods.sscha_phonons_from_phonopy(p)
    s1 = s.GenerateSupercellDyn(s.GetSupercell())
    s1.save_qe('Gamma_dyn')

def read_eig_new(nmode):
    f = open('Gamma_dyn1','r').readlines()
    natom=int(f[2].split()[1])
    dis = np.zeros((natom,3))
    for i in range(len(f)):
        key = 'freq ('+'%5i' % (nmode+1)+')'
        if key in f[i]:
            for j in range(natom):
                dis[j] = np.array(f[i+j+1].split()[1:6:2])
    return dis

def generate_dis(eig,mass,factor):
    for i in range(len(eig)):
        eig[i] = eig[i]/np.sqrt(mass[i])
    return eig*factor

def output_vesta(input_file,output_file,nmode=0):
    with open(input_file,'r') as f1:
        d = f1.readlines()

    with open(output_file,'w') as f2:
        for i in range(len(d)):
            if 'VECTR' in d[i]:
                id1 = i
            elif 'SPLAN' in d[i]:
                id2 = i
                break
        write=True
        for i in range(len(d)):
            if i < id1 or i >= id2:
                f2.write(d[i])
            else:
                if write == True:

                    #eig0,mass = read_eig(0)
                    #eig1,mass = read_eig(1)
                    #f0=0;f1=1
                    #eig = (f0*eig0+f1*eig1)/np.sqrt(f1**2+f0**2)
                    dis = read_eig_new(nmode)

                    write_dis_vesta(dis*10,f2)
                    write = False

if __name__=="__main__":
    generate_supercelldyn()
    for nmode in range(4):
        input_file  = '333.vesta'
        output_file = f'333.{nmode}.vesta'

        output_vesta(input_file,output_file,nmode)
