import cellconstructor as CC
import cellconstructor.Structure
import cellconstructor.Phonons
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
import sys,os
import numpy as np
import config

def generate_ensemble(pop):
    """
    GENERATE THE ENSEMBLE
    ====================

    This method read the harmonic dynamical matrix and generate the ensemble.

    Parameters
    ----------
        pop : int
            This is the index of the ensemble.
    """
    if pop == 1:
        dyn = CC.Phonons.Phonons(f"harmonic_{config.nq1}{config.nq2}{config.nq3}_dyn", nqirr = config.nqirr)
        # make the harmonic dynamical matrix positive definite. The trial H is defined to be definite positive. 
        dyn.Symmetrize()
        dyn.ForcePositiveDefinite()
    else:
        dyn = CC.Phonons.Phonons(f"dyn_pop{pop-1}_", nqirr = config.nqirr)

    # We generate the N randomly displaced structures
    ensemble = sscha.Ensemble.Ensemble(dyn, T0 = config.T0, supercell= dyn.GetSupercell())
    ensemble.generate(N = config.N_config)

    # We save the ensemble to ens{pop}
    ensemble.save(f"ens{pop}", population = pop)
    return ensemble

def generate_dft_input(pop):
    """
    GENERATE THE QE INPUT FILES
    ===========================
    This method generates the input files for the QE calculations. It reads the scf_ files in the ensemble and generate the input files.

    Parameters
    ----------
        pop : int
            This is the index of the ensemble.
    """
    # make sure this is a new ensemble
    if not os.path.exists(f"run_dft{pop}"):
        os.mkdir(f"run_dft{pop}")
    else:
        print("run_dft{} already exists".format(pop))
        exit()
    ensemble = generate_ensemble(pop)

    # read all the scf_ files in the ensemble
    all_scf_files = [os.path.join(f"ens{pop}", f) for f in os.listdir(f"ens{pop}") if f.startswith("scf_")]

    for file in all_scf_files:
        number = int(file.split("_")[-1].split(".")[0])
        typical_espresso_header = f"""
&control
        calculation = "scf"
        tstress = .true.
        tprnfor = .true.
        disk_io = "none"
        pseudo_dir = "/home/kw600/pseudo"
        outdir='/home/kw600/rds/hpc-work/ps_materialscloud/py_sscha/run_dft2/d{number}'
&end

&system
        nat = {ensemble.structures[0].N_atoms}
        ntyp = 2
        ibrav = 0
        ecutwfc = 45
        ecutrho = 360
        occupations='smearing'
        smearing = "fd"
        degauss = 0.0004
&end
&electrons
conv_thr = 1.D-5
&end
ATOMIC_SPECIES
        Ti  47.867 ti_pbe_v1.4.uspp.F.UPF
        Se 78.96   Se_pbe_v1.uspp.F.UPF
K_POINTS automatic
9 9 1  0 0 0
"""
        filename = os.path.join(f"run_dft{pop}", "espresso_run_{}.pwi".format(number))

        with open(filename, "w") as f:
            f.write(typical_espresso_header)
            ff = open(file, "r")
            structure_lines = ff.readlines()
            ff.close()
            f.writelines(structure_lines)
    return ensemble

def generate_vasp_input(pop):
    """
    GENERATE THE VASP INPUT FILES
    ===========================
    This method generates the input files for the VASP calculations. It reads the scf_ files in the ensemble and generate the POSCAR files. 
    Note the INCAR POTCAR KPOINTS files need to be generated manually.

    Parameters
    ----------
        pop : int
            This is the index of the ensemble.
    """
    # make sure this is a new ensemble
    if not os.path.exists(f"run_dft{pop}"):
        os.mkdir(f"run_dft{pop}")
    else:
        print("run_dft{} already exists".format(pop))
        exit()
    ensemble = generate_ensemble(pop)
    all_scf_files = [os.path.join(f"ens{pop}", f) for f in os.listdir(f"ens{pop}") if f.startswith("scf_")]

    for ii,file in enumerate(all_scf_files):
            #read the index of the structure
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
                    #need to sort the atoms
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
                    #factor is the supercell size.
                    #Note the 1 6 6 needs to be mannual changed based on the system.
                    factor=config.nq1*config.nq2*config.nq3
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
    return ensemble


if __name__ == "__main__":
    
    pop = int(sys.argv[1])
    generate_vasp_input(pop)



