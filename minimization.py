import sys,os
import numpy as np
import cellconstructor as CC
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
import config

def collect_data(pop,Test=False):
    """
    COLLECT DATA
    ============
    This method collects the data from the DFT calculations. It reads the output files and collect the energies and forces.

    Parameters
    ----------
        pop : int
            This is the index of the ensemble.
    """
    index=''
    directory = f"run_dft{pop}"
   
    output_filenames = [f for f in os.listdir(directory) if f.endswith(".pwo")] # We select only the output files
    output_files = [os.path.join(directory, f) for f in output_filenames]
    energies = np.zeros(len(output_files))

    # We add the directory/outpufilename to load them correctly
    # We prepare the array of energies
     
    for file in output_files:
            # print(file)
        try:	
            # Get the number of the configuration.
            id_number = int(file.split("_")[-1].split(".")[0]) # The same as before, we need the to extract the configuration number from the 
            ff = open(file, "r")
            lines = [l.strip() for l in ff.readlines()] # Read the whole file removing tailoring spaces
            ff.close()
            energy_line = next(l for l in lines if len(l) > 0 if l.split()[0] == "!")
            energies[id_number - 1] = float(energy_line.split()[4])
            
            # Now we can collect the force
            nat_line = next( l for l in lines if len(l) > 0 if l.split()[0] == "number" and l.split()[2] == "atoms/cell" )
            nat = int(nat_line.split()[4])
            #except:
            #	print("Error: no number of atoms found in file {}".format(file))
            # Now allocate the forces and read them
            forces = np.zeros((nat, 3))
            forces_lines = [l for l in lines if len(l) > 0 if l.split()[0] == "atom"] # All the lines that starts with atom will contain a force
            for i in range(nat):
                forces[i, :] = [float(x) for x in forces_lines[i].split()[-3:]] # Get the last three number from the line containing the force
            
            # Now we can take the stress tensor
            stress = np.zeros((3,3))
            # We pick the index of the line that starts with the words total stress
            index_before_stress = next(i for i, l in enumerate(lines) if len(l) > 0 if l.split()[0] == "total" and l.split()[1] == "stress")
            # The stress tensor is located just after it
            for i in range(3):
                index = i + index_before_stress + 1
                stress[i, :] = [float(x) for x in lines[index].split()[:3]]

            # We can save the forces_population1_X.dat and pressures_population1_X.dat files
            force_file = os.path.join(f"ens{pop}", "forces_population{}_{}.dat".format(pop,id_number))
            stress_file = os.path.join(f"ens{pop}", "pressures_population{}_{}.dat".format(pop,id_number))
            np.savetxt(force_file, forces)
            np.savetxt(stress_file, stress)
        except:
            print("Error: something went wrong with file {}".format(file))
            # index=index+file.replace('_','.').split('.')[-2]
    # Now we read all the configurations, we can save the energy file
    energy_file = os.path.join(f"ens{pop}", f"energies_supercell_population{pop}.dat")
    np.savetxt(energy_file, energies)

def collect_vaspdata(pop):
    """
    COLLECT VASP DATA
    ============
    This method collects the data from the VASP calculations. It reads the output files and collect the energies and forces and finally transform them into QE format.

    Parameters
    ----------
        pop : int
            This is the index of the ensemble.
    """
    N_config = config.N_config
    N_prim = config.N_prim
    N_supercell=config.nq1*config.nq2*config.nq3
    f=open(f'./energies_supercell_population{pop}.dat','w');f.close()
    for i in range(1,N_config+1):
        L='''CELL_PARAMETERS angstrom
    10.9338607397434711  0.0000000000000000  0.0000000000000000
    -5.4669303698717355  9.4690011644160794  0.0000000000000000
    0.0000000000000000  0.0000000000000000  18.3188012412753238

    ATOMIC_POSITIONS angstrom
    '''
        stress=np.zeros((3,3))
        with open(f'./run_dft{pop}/vasp{i}/OUTCAR') as f0:
            d=f0.readlines()
        with open(f'./run_dft{pop}/vasp{i}/OSZICAR') as ff:
            dd=ff.readlines()[-1]

        energy=float(dd.split()[4])/(27.211396/2)
        for ii in range(len(d)):
            if 'TOTAL-FORCE' in d[ii]:
                d1=np.loadtxt(d[ii+2:ii+2+N_prim*N_supercell])
                break
            elif 'in kB' in d[ii]:
                    s = np.array(d[ii].split()[2:],dtype=float)
                    stress[0,0]=s[0];stress[1,1]=s[1];stress[2,2]=s[2];stress[0,1]=s[3];stress[1,0]=s[3]
                    stress[1,2]=s[4];stress[2,1]=s[4];stress[0,2]=s[5];stress[2,0]=s[5]
        pos=d1[:,:3];force=d1[:,3:]/25.71104309541616
        stress=stress/(29421.02648438959/2*10)

        #change the atom order from vasp back to QE style
        sc=pos[:1*N_supercell];v=pos[1*N_supercell:7*N_supercell];sn=pos[7*N_supercell:13*N_supercell]
        sc=[i for i in sc];v=[i for i in v];sn=[i for i in sn]
        with open(f'./ens{pop}/vasp_population{pop}_{i}.dat','w') as f1:
            f1.writelines(L)
            for j in range(N_supercell):
                a=sc.pop(0)
                f1.write(f'Sc {a[0]} {a[1]} {a[2]}\n')
                for k in range(6):
                    a=v.pop(0)
                    f1.write(f'V {a[0]} {a[1]} {a[2]}\n')
                for k in range(6):
                    a=sn.pop(0)
                    f1.write(f'Sn {a[0]} {a[1]} {a[2]}\n')
        sc=force[:1*N_supercell];v=force[1*N_supercell:7*N_supercell];sn=force[7*N_supercell:13*N_supercell]
        sc=[i for i in sc];v=[i for i in v];sn=[i for i in sn]
        with open(f'./ens{pop}/forces_population{pop}_{i}.dat','w') as f2:
            for j in range(N_supercell):
                a=sc.pop(0)
                f2.write(f' {a[0]} {a[1]} {a[2]}\n')
                for k in range(6):
                    a=v.pop(0)
                    f2.write(f' {a[0]} {a[1]} {a[2]}\n')
                for k in range(6):
                    a=sn.pop(0)
                    f2.write(f' {a[0]} {a[1]} {a[2]}\n')
        with open(f'./ens{pop}/pressures_population{pop}_{i}.dat','w') as f3:
            for k in range(3):
                f3.writelines(f'{stress[k,0]} {stress[k,1]} {stress[k,2]}\n')
        with open(f'./ens{pop}/energies_supercell_population{pop}.dat','a') as f4:
            f4.write(f'{energy}\n')

    #check the energy file is complete
    with open(f'./ens{pop}/energies_supercell_population{pop}.dat') as f5:
        d=f5.readlines()
    if len(d)!=N_config:
        print('energy file is not complete')
        exit()
    E = np.array(d, dtype=float)
    for i in range(N_config):
        if abs(float(d[i])-np.mean(E))>abs(np.mean(E)*0.1):
            print('energy file is not complete')
            exit()


def scha(pop):
    IO_freq = sscha.Utilities.IOInfo()
    IO_freq.SetupSaving(f"minim_info{pop}")

    if pop == 1:
            dyn = CC.Phonons.Phonons(f"harmonic_{config.nq1}{config.nq2}{config.nq3}_dyn", nqirr = config.nqirr)
            dyn.Symmetrize()
            dyn.ForcePositiveDefinite()
    else:
            dyn = CC.Phonons.Phonons(f"dyn_pop{int(pop-1)}_", nqirr = config.nqirr)

    ensemble = sscha.Ensemble.Ensemble(dyn, T0 = config.T0, supercell= dyn.GetSupercell())
    ensemble.load(f"ens{pop}", population = pop, N = config.N_config)
    if pop == 1:
            ensemble.update_weights(dyn, config.T0) # Restore the original density matrix at T = 100 K
    minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
    
    # Ignore the structure minimization (is fixed by symmetry)
    minimizer.minim_struct = False


    # Setup the threshold for the ensemble wasting
    minimizer.kong_liu_ratio = config.kong_liu_ratio # Usually 0.5 is a good value

    # Lest start the minimization
    # minimizer.init()
    minimizer.run(custom_function_post = IO_freq.CFP_SaveAll)
    minimizer.finalize()
    minimizer.dyn.save_qe(f"dyn_pop{pop}_")
    print('Converged?',minimizer.is_converged())
    return minimizer.is_converged()

if __name__ == "__main__":
    pop=int(sys.argv[1])
    if os.path.exists(f"run_dft{pop}/vasp1"):
        collect_vaspdata(pop)
    else:
        collect_data(pop)
    min=scha(pop)
    

    


def SetupFromTensor_parallel(self,id_cell2,id_cell3,tensor=None):
        """
        Setup the third order force constant form 3rd order tensor defined in the supercell
        
        
        Parameters
        ----------
            unitcell_structure : Structure()
                The structure in the unit cell
            supercell_structure : Structure()
                The supercell structure on which the tensor has been computed
            supercell_size : truple
                The number of supercell along each lattice vector
            tensor : ndarray(size =(3*nat_sc, 3*nat_sc, 3*nat_sc, dtype = np.double)
                The third order tensor
        """


        n_sup = np.prod(self.supercell_size)
        nat = self.unitcell_structure.N_atoms
        supercell_size = self.supercell_size


        nat_sc= n_sup * nat
        n_R = self.n_R

        supercell_structure = self.supercell_structure
        unitcell_structure = self.unitcell_structure



        for index_cell2 in [id_cell2]:
            n_cell_x2,n_cell_y2,n_cell_z2=Methods.one_to_three_len(index_cell2,v_min=[0,0,0],
                                                                   v_len=supercell_size)
            for index_cell3 in [id_cell3]:
                n_cell_x3,n_cell_y3,n_cell_z3=Methods.one_to_three_len(index_cell3,v_min=[0,0,0],
                                                                       v_len=supercell_size)
                #
                total_index_cell = index_cell3 + n_sup * index_cell2
                #
                self.x_r_vector2[:, total_index_cell] = (n_cell_x2, n_cell_y2, n_cell_z2)
                self.r_vector2[:, total_index_cell] = unitcell_structure.unit_cell.T.dot(self.x_r_vector2[:,total_index_cell])
                self.x_r_vector3[:, total_index_cell] =  n_cell_x3, n_cell_y3, n_cell_z3
                self.r_vector3[:, total_index_cell] = unitcell_structure.unit_cell.T.dot(self.x_r_vector3[:, total_index_cell])

                for na1 in range(nat):
                    #
                    for na2 in range(nat):
                        # Get the atom in the supercell corresponding to the one in the unit cell
                        na2_vect = unitcell_structure.coords[na2, :] + self.r_vector2[:, total_index_cell]
                        nat2_sc = np.argmin( [np.sum( (supercell_structure.coords[k, :] - na2_vect)**2) for k in range(nat_sc)])
                        #
                        for na3 in range(nat):
                            # Get the atom in the supercell corresponding to the one in the unit cell
                            na3_vect = unitcell_structure.coords[na3, :] + self.r_vector3[:, total_index_cell]
                            nat3_sc = np.argmin( [np.sum( (supercell_structure.coords[k, :] - na3_vect)**2) for k in range(nat_sc)])
                            #
                            self.tensor[total_index_cell,
                                        3*na1 : 3*na1+3,
                                        3*na2 : 3*na2+3,
                                        3*na3 : 3*na3+3] = tensor[3*na1 : 3*na1 +3,
                                                                3*nat2_sc : 3*nat2_sc + 3,
                                                                3*nat3_sc : 3*nat3_sc + 3]
