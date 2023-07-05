import numpy as np
import sys
N_config = 512
N_prim=13
N_supercell=3*3*3
pop=sys.argv[1]
f=open(f'./energies_supercell_population{pop}.dat','w');f.close()
for i in range(1,N_config+1):
    if i % 50 == 0:
        print(f'Reading configuration {i}')
    #get the force and stress
    stress=np.zeros((3,3))
    with open(f'../run_dft{pop}/vasp{i}/OUTCAR') as f0:
        d=f0.readlines()
    for ii in range(len(d)):
        if 'TOTAL-FORCE' in d[ii]:
            d1=np.loadtxt(d[ii+2:ii+2+N_prim*N_supercell])
            break
        elif 'in kB' in d[ii]:
                s = np.array(d[ii].split()[2:],dtype=float)
                stress[0,0]=s[0];stress[1,1]=s[1];stress[2,2]=s[2];stress[0,1]=s[3];stress[1,0]=s[3]
                stress[1,2]=s[4];stress[2,1]=s[4];stress[0,2]=s[5];stress[2,0]=s[5]        
    force=d1[:,3:]/25.71104309541616
    stress=stress/(29421.02648438959/2*10)

    #get the total energy

    # with open(f'../run_dft{pop}/vasp{i}/OSZICAR') as ff:
    #     dd=ff.readlines()[-1]
    # energy=float(dd.split()[4])/(27.211396/2)
 
            #d is outcar readlines
    for jj in d:
        if 'average (electrostatic) potential at core' in jj:
            break
        elif 'energy(sigma->0)' in jj:
            energy=float(jj.split()[-1])/(27.211396/2)




    # write the position. Only used for test
    # sc=pos[:1*N_supercell];v=pos[1*N_supercell:7*N_supercell];sn=pos[7*N_supercell:13*N_supercell]
    # sc=[i for i in sc];v=[i for i in v];sn=[i for i in sn]
    # with open(f'./vasp_population{pop}_{i}.dat','w') as f1:
    #     f1.writelines(L)
    #     for j in range(N_supercell):
    #         a=sc.pop(0)
    #         f1.write(f'Sc {a[0]} {a[1]} {a[2]}\n')
    #         for k in range(6):
    #             a=v.pop(0)
    #             f1.write(f'V {a[0]} {a[1]} {a[2]}\n')
    #         for k in range(6):
    #             a=sn.pop(0)
    #             f1.write(f'Sn {a[0]} {a[1]} {a[2]}\n')

    #write the force
    sc=force[:1*N_supercell];v=force[1*N_supercell:7*N_supercell];sn=force[7*N_supercell:13*N_supercell]
    sc=[i for i in sc];v=[i for i in v];sn=[i for i in sn]
    with open(f'./forces_population{pop}_{i}.dat','w') as f2:
        for j in range(N_supercell):
            a=sc.pop(0)
            f2.write(f' {a[0]} {a[1]} {a[2]}\n')
            for k in range(6):
                a=v.pop(0)
                f2.write(f' {a[0]} {a[1]} {a[2]}\n')
            for k in range(6):
                a=sn.pop(0)
                f2.write(f' {a[0]} {a[1]} {a[2]}\n')
    #write the stress
    with open(f'./pressures_population{pop}_{i}.dat','w') as f3:
        for k in range(3):
            f3.writelines(f'{stress[k,0]} {stress[k,1]} {stress[k,2]}\n')
    #write the total energy
    with open(f'./energies_supercell_population{pop}.dat','a') as f4:
        f4.write(f'{energy}\n')





