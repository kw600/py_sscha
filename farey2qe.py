import numpy as np
import sys

def generate_dyn(k_abs, C, R, M):
    num_cell, num_atoms, _, _, _ = C.shape
    D_q = np.zeros((num_atoms, 3, num_atoms, 3), dtype=complex)
    for i_atom in range(num_atoms):
        for i_cart in range(3):
            for j_atom in range(num_atoms):
                for j_cart in range(3):
                    for i_cell in range(num_cell):
                        exp_k_dot_r = 0.0
                        r = R[i_cell,i_atom,j_atom,:]
                        num_image = np.count_nonzero(r != None)
                        for i_im in range(num_image):
                            exp_k_dot_r += np.exp(-2j*np.pi*k_abs.dot(R[i_cell,i_atom,j_atom][i_im]))
                        exp_k_dot_r = exp_k_dot_r / len(R[i_cell, i_atom, j_atom])
                        D_q[i_atom, i_cart, j_atom, j_cart] += C[i_cell, i_atom, i_cart, j_atom, j_cart] * exp_k_dot_r /np.sqrt(M[i_atom]*M[j_atom])
    
    return D_q

def generate_dyn_qe(k, C, R, alat):
    """
    GENERATE THE DYNAMICAL MATRIX IN QE UNITS
    ====================

    This method read the force constants C to generate dynamical martrix at q points for QE.

    Parameters
    ----------
        k : numpy array
            This is the q point in reduced coordinates in QE output.
        C : numpy array
            This is the force constants in cartesian coordinates.
        R : numpy array
            This is the supercell in cartesian coordinates.
        alat : float
            This is the lattice constant in bohr.
    """
    #Note that k is in cartesian coordinates and needs to be devided by alat
    num_cell, num_atoms, _, _, _ = C.shape
    k_cart = k / alat

    D_q = np.zeros((num_atoms, 3, num_atoms, 3), dtype=complex)
    for i_atom in range(num_atoms):
        for i_cart in range(3):
            for j_atom in range(num_atoms):
                for j_cart in range(3):
                    for i_cell in range(num_cell):
                        exp_k_dot_r = 0.0
                        r = R[i_cell,i_atom,j_atom,:]
                        num_image = np.count_nonzero(r != None)
                        for i_im in range(num_image):
                            # print('kdot',k_cart.dot(r[i_im]))
                            exp_k_dot_r += np.exp(-2j * np.pi * k_cart.dot(r[i_im]))
                        exp_k_dot_r = exp_k_dot_r / num_image
                        D_q[i_atom, i_cart, j_atom, j_cart] += C[i_cell, i_atom, i_cart, j_atom, j_cart] * exp_k_dot_r
    # 2 3 0 1 is simply switch iatom icart and jatom jcart
    # not sure exactly why this is needed but it works. I guess it is because of the way of QE default format. If the code breaks, check this part.
    D_q = np.transpose(D_q, (2, 3, 0, 1)) * 2  # Hatree to Rydberg
    return D_q

def qpoint(input_path,index):
    with open(f'{input_path}{index}','r') as f:
        lines = f.readlines()
    q=[]
    for i in range(len(lines)):
        if ('Dynamical Matrix in cartesian axes' in lines[i]) or ('Dynamical  Matrix in cartesian axes' in lines[i]):
            q.append(lines[i+2].split()[3:6])
    return np.array(q,dtype=float)

def get_charge(i_path,index):
    with open(f'{i_path}{index}','r') as f:
        lines = f.readlines()
    for i in range(len(lines)):
        if 'Dielectric Tensor' in lines[i]:
            return lines[i:i+17]


def output_new(q,ipath,opath,index):
    global reciprocal_lattice, num_cell, num_atoms, C, R, M, header
    with open(f'{opath}{index}','w') as f:
        f.write(header)
        for i in range(len(q)):
            D = generate_dyn_qe(q[i] , C, R, alat)
            D_real = np.real(D); D_imag = np.imag(D)
            f.write('\n     Dynamical  Matrix in cartesian axes\n')
            f.write('\n')
            f.write(f'     q = ( {q[i][0]} {q[i][1]} {q[i][2]} )\n')
            f.write('\n')
            for i_atom in range(num_atoms):
                for j_atom in range(num_atoms):
                    f.write(f'    {i_atom+1}    {j_atom+1}\n')
                    for i_cart in range(3):
                        for j_cart in range(3):
                            f.write(f'  {D_real[i_atom,i_cart,j_atom,j_cart]:>14.9f}{D_imag[i_atom,i_cart,j_atom,j_cart]:>14.9f}')
                        f.write('\n')
        charge=get_charge(ipath,index)
        f.write('\n')
        if charge:
            for ll in charge:
                f.write(f'{ll}')

def get_header(file):
    d = open(file,'r').readlines()
    for i in range(len(d)):
        if 'cartesian axes' in d[i]:
            return d[:i]

def read_header(header):
    data = header.split('\n')
    ntyp = int(data[2].split()[0]); num_atoms = int(data[2].split()[1]); alat = float(data[2].split()[3])
    lattice = np.zeros((3,3))
    for i in range(3):
            lattice[i]=data[4+i].split()
    return ntyp, num_atoms, alat, lattice


if __name__=='__main__':
    script_dir = sys.argv[1]
    supercell = sys.argv[1].replace('/','').split('_')[-1]                

    header = ''
    harmonic=open(f'./{script_dir}/harmonic_{supercell}_dyn1').readlines()
    for i in range(len(harmonic)):
        if 'Matrix in cartesian axes' in harmonic[i]:
            break
        else:
            header += harmonic[i]

    M = np.loadtxt(script_dir + "/equilibrium.dat", dtype=np.float64, comments=['#', '$', '@'], skiprows=1, usecols=1)
    ibz = np.loadtxt(script_dir + "/ibz.dat").reshape((-1,4))
    ntyp, num_atoms, alat, lattice = read_header(header)
    data = np.loadtxt(script_dir + '/force.dat')
    indices = data[:, :5].astype(int) - 1  # subtract 1 to convert to 0-based indexing
    values = data[:, 5]
    # determine shape of 5D array
    shape = np.max(indices, axis=0) + 1
    # create 5D array and fill with values
    C = np.empty(shape)
    C[indices[:, 0], indices[:, 1], indices[:, 2], indices[:, 3], indices[:, 4]] = -values

    data = np.loadtxt(script_dir + '/delta_prim.dat')
    indices = data[:, :4].astype(int) - 1  # subtract 1 to convert to 0-based indexing
    values = data[:, 4:]
    shape = np.max(indices, axis=0) + 1
    R = np.empty(shape, dtype=list)
    R[indices[:, 0], indices[:, 1], indices[:, 2], indices[:, 3]] = values.tolist()

    for i in range(1,len(ibz)+1):
        i_path=f'./{script_dir}/harmonic_{supercell}_dyn'
        q=qpoint(i_path,i)
        print(np.dot(q,np.transpose(lattice)))
        o_path = f'./{script_dir}/' + 'f2q_'
        output_new(q,i_path,o_path,i)
        print(f'file {o_path}{i} done')
