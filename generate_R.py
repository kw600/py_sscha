import numpy as np
import os, sys

def inv_33(mat):
    return np.linalg.inv(mat)

def min_images_brute_force(a, lat_vec):
    """
    Compute the minimum image vectors b of vector a with respect to the lattice specified by the rows of lat_vec.
    """
    tol = 1e-8
    check_shell = 3
    tol_L2 = tol * np.dot(lat_vec[0, :3], lat_vec[0, :3])

    rec_vec = np.linalg.inv(lat_vec).T
    n = np.floor(np.dot(a, rec_vec[:, :3])).astype(int)

    mag_b_sq = -1.0
    nim = -1
    b = np.zeros((8, 3))

    for i in range(n[0] - check_shell, n[0] + check_shell + 1):
        delta1 = a - float(i) * lat_vec[0, :3]
        for j in range(n[1] - check_shell, n[1] + check_shell + 1):
            delta2 = delta1 - float(j) * lat_vec[1, :3]
            for k in range(n[2] - check_shell, n[2] + check_shell + 1):
                delta3 = delta2 - float(k) * lat_vec[2, :3]
                dist2 = np.dot(delta3, delta3)
                # print(i,j,k,dist2,mag_b_sq)
                if abs(dist2 - mag_b_sq) <= tol_L2:
                    # print('image + 1',nim,dist2)
                    # print(b)
                    nim += 1
                    if nim > 8:
                        raise ValueError("Need to increase maxim parameter.")
                    b[nim] = delta3
                elif dist2 < mag_b_sq or nim == -1:

                    mag_b_sq = dist2
                    nim = 0
                    b[0] = delta3
    if nim < 0:
        raise ValueError("Bug.")
    return b[:nim+1],nim

def read_atom_prim_cart(file):
    with open(file, 'r') as f:
        line = f.readlines()
        atom_prim_cart = np.zeros((len(line),3), dtype=np.float64)
        for i in range(len(line)):
            atom_prim_cart[i,:] = np.array(line.split()[1:], dtype=np.float64)
    return atom_prim_cart

def R1_index():
    global q1, q2, q3
    R=np.zeros((q1*q2*q3,3))
    count=0
    for x0 in range(q1):
        for y0 in range(q2):
            for z0 in range(q3):
                R[count]=np.array([x0+1,y0+1,z0+1])
                count+=1
    return R


def atom_super_cart(jatom,icell):
    global prim_latt_vecs,atom_prim_cart
    i=R1_index()[icell]
    R = (i[0]-1)*prim_latt_vecs[0]+(i[1]-1)*prim_latt_vecs[1]+(i[2]-1)*prim_latt_vecs[2]+atom_prim_cart[jatom]
    return R

def super_cart(icell):
    global prim_latt_vecs,atom_prim_cart
    i=R1_index()[icell]
    R = (i[0]-1)*prim_latt_vecs[0]+(i[1]-1)*prim_latt_vecs[1]+(i[2]-1)*prim_latt_vecs[2]
    return R


# Example usage
if __name__ == "__main__":
    path=str(sys.argv[1])
    os.chdir(path)


    with open('harmonic1','r') as f3:
        lines=f3.readlines()
        alat = float(lines[2].split()[3])
        ntype= int(lines[2].split()[0])
        natom= int(lines[2].split()[1])
        for i in range(len(lines)):
            if 'Basis vectors' in lines[i]:
                index = i
                L = np.loadtxt(lines[i+1:i+4])
                prim_latt_vecs = L * alat
                break

        atom_list=[];mass=[]
        for i in range(7,7+ntype):
            l = lines[i].split()
            atom_list.append(l[1])
            mass.append(l[-1])
    atom_prim_cart=[]
    with open('equilibrium.dat', 'w') as f_atoms:
                    f_atoms.write(f'{natom}\n')
                    count=0
                    for j in range(7+ntype,len(lines)):
                        # print(lines[j])
                        if len(lines[j].split())== 5:
                            atom_prim_cart.append(lines[j].split()[2:5])
                            atom = atom_list[int(lines[j].split()[1])-1].replace("'","")
                            m = mass[int(lines[j].split()[1])-1]
                            f_atoms.write(f'{atom} {m} {alat*float(lines[j].split()[2])} {alat*float(lines[j].split()[3])} {alat*float(lines[j].split()[4])}\n')
                            count+=1
                        else:
                            break

    atom_prim_cart = np.array(atom_prim_cart,dtype=float)*alat
    # print(atom_prim_cart)
    with open('harmonic0','r') as f4:
        lines=f4.readlines()
        q=np.array(lines[0].split(),dtype=int)
        nq = int(lines[1].split()[0])
        Q=np.zeros((nq,3))
        for i in range(nq):
            Q[i]=np.array(lines[i+2].split(),dtype=float)
    with open('grid.dat','w') as ff:
        [q1,q2,q3] = q
        ff.write(f'{q1} {q2} {q3}')
    Q = np.dot(Q, np.transpose(L))
    with open('ibz0.dat','w') as f5:
        # f5.write(f'{nq}\n')
        for i in range(nq):
            f5.write(f'{Q[i,0]:.06f} {Q[i,1]:.06f} {Q[i,2]:.06f}\n')

    # assue prim_latt_vecs[0,:] is the a vector
    with open('lattice.dat','w') as f1:
        for i in range(3):
            f1.write(f'{prim_latt_vecs[i,0]:.12f} {prim_latt_vecs[i,1]:.12f} {prim_latt_vecs[i,2]:.12f}\n')

    basis =len(atom_prim_cart); no_prim_cells=q1*q2*q3
    print('no_prim_cells',no_prim_cells)
    delta_prim = np.zeros((8, 3, basis, basis, no_prim_cells))
    super_latt_vecs=np.array([prim_latt_vecs[0]*q1, prim_latt_vecs[1]*q2, prim_latt_vecs[2]*q3])

    with open('delta_prim_1.dat', 'w') as f:
        for i_atom in range(basis):
        # for i_atom in range(1):
            for j_atom in range(basis):
            # for j_atom in range(4,5):
                # print(basis,'basis')
                delta_r_corr = atom_prim_cart[i_atom] - atom_prim_cart[j_atom]
                for i_cell in range(no_prim_cells):
                # for i_cell in range(1,2):
                    temp = atom_super_cart(j_atom,i_cell) - atom_prim_cart[i_atom]
                    delta_r_ims, no_im_cells = min_images_brute_force(temp, super_latt_vecs)
                    for i_im in range(no_im_cells+1):
                        v=delta_prim[i_im,:,i_atom,j_atom,i_cell]
                        delta_prim[i_im,:,i_atom,j_atom,i_cell] = delta_r_ims[i_im] + delta_r_corr
                        output_str = "{:>12d}{:>12d}{:>12d}{:>12d}{:>21.16f}{:>26.16f}{:>26.16f}\n".format(i_cell+1, i_atom+1, j_atom+1, i_im+1, v[0], v[1], v[2])
                        f.write(output_str)

    with open('delta_prim_1.dat', 'r') as f11:
        lines = f11.readlines()

    def get_index(line):
        columns = line.split()
        index = tuple(map(int, columns[:3]))
        return index

    # Sort the lines based on the index using the key function
    sorted_lines = sorted(lines, key=get_index)

    # Write the sorted lines to a new file
    with open('delta_prim.dat', 'w') as f22:
        f22.writelines(sorted_lines)

    exit()
#   generate R for normal fourier transform

    with open('R_1.dat', 'w') as f:
        for i_atom in range(basis):
            for j_atom in range(basis):
                for i_cell in range(no_prim_cells):
                    no_im_cells=0
                    for i_im in range(no_im_cells+1):
                        v=super_cart(i_cell)
                        output_str = "{:>12d}{:>12d}{:>12d}{:>12d}{:>21.16f}{:>26.16f}{:>26.16f}\n".format(i_cell+1, i_atom+1, j_atom+1, i_im+1, v[0], v[1], v[2])
                        f.write(output_str)

    with open('R_1.dat', 'r') as f11:
        lines = f11.readlines()


    # Sort the lines based on the index using the key function
    sorted_lines = sorted(lines, key=get_index)

    # Write the sorted lines to a new file
    with open('R_2.dat', 'w') as f22:
        f22.writelines(sorted_lines)
