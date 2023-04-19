import numpy as np

import numpy as np

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
					# print('image = 0','reset')
					# print(i,j,k,dist2,mag_b_sq)
					# print(b)
	if nim < 0:
		raise ValueError("Bug.")
	return b[:nim+1],nim

def min_images_brute_force1(a, lat_vec):
    check_shell = 3
    tol = 1e-8
    tol_L2 = tol * np.dot(lat_vec[0, :3], lat_vec[0, :3])
    rec_vec = inv_33(lat_vec)
    rec_vec = rec_vec.T
    n = np.array([
        np.floor(np.dot(a, rec_vec[0, :3])),
        np.floor(np.dot(a, rec_vec[1, :3])),
        np.floor(np.dot(a, rec_vec[2, :3]))
    ], dtype=np.int32)
    delta1 = np.zeros(3)
    delta2 = np.zeros(3)
    delta3 = np.zeros(3)
    mag_b_sq = -1.0
    nim = -1
    b = np.zeros((3, 8), dtype=np.float64)
    for i in range(n[0] - check_shell, n[0] + check_shell + 1):
        delta1[:] = a - float(i) * lat_vec[0, :3]
        for j in range(n[1] - check_shell, n[1] + check_shell + 1):
            delta2[:] = delta1 - float(j) * lat_vec[1, :3]
            for k in range(n[2] - check_shell, n[2] + check_shell + 1):
                delta3[:] = delta2 - float(k) * lat_vec[2, :3]
                dist2 = np.dot(delta3, delta3)
                if abs(dist2 - mag_b_sq) <= tol_L2:
                    nim += 1
                    if nim > 8:
                        raise ValueError("Need to increase maxim parameter.")
                    b[:, nim - 1] = delta3[:]
                elif dist2 < mag_b_sq or nim == -1:
                    mag_b_sq = dist2
                    nim = 1
                    b[:, 0] = delta3[:]
    if nim <= 0:
        raise ValueError("Bug.")
    return b, nim

def calculate_minimum_image_vectors(basis, atom_prim_frac, prim_latt_vecs, super_latt_vecs, grid):
    no_grid_points = np.prod(grid)
    cell_pos_cart = np.zeros((3, no_grid_points), dtype=np.float64)
    atom_prim_cart = np.zeros((3, basis), dtype=np.float64)
    atom_super_cart = np.zeros((3, basis, no_grid_points), dtype=np.float64)
    no_prim_cells = np.zeros((basis, basis, no_grid_points), dtype=np.float64)
    delta_prim = np.zeros((3, 8, basis, basis, no_grid_points), dtype=np.float64)
    i_cell = 0
    for m1 in range(grid[0]):
        for m2 in range(grid[1]):
            for m3 in range(grid[2]):
                i_cell += 1
                if i_cell > no_grid_points:
                    raise ValueError("Found more primitive cells than in supercell.")
                cell_pos_cart[:, i_cell - 1] = np.array([
                    float(m1) * prim_latt_vecs[0, :3],
                    float(m2) * prim_latt_vecs[1, :3],
                    float(m3) * prim_latt_vecs[2, :3]
                ]).sum(axis=0)
    for i_atom in range(basis):
        atom_prim_cart[:, i_atom] = np

def read_atom_prim_cart(file):
	with open(file, 'r') as f:
		line = f.readlines()
		atom_prim_cart = np.zeros((len(line),3), dtype=np.float64)
		for i in range(len(line)):
			atom_prim_cart[i,:] = np.array(line.split()[1:], dtype=np.float64)
	return atom_prim_cart

def R_index():
	global q1, q2, q3
	R=np.zeros((q1*q2*q3,3))
	count=0
	m1=np.median(np.array([i for i in range(q1)]))
	m2=np.median(np.array([i for i in range(q2)]))
	m3=np.median(np.array([i for i in range(q3)]))
	for x in range(0,-q1,-1):
		x0=x+1
		if x0<=0:
			x0+=q1
		if x0>m1:
			x0-=q1
		for y in range(0,-q2,-1):
			y0=y+1
			if y0<=0:
				y0+=q2
			if y0>m2:
				y0-=q2
			for z in range(0,-q3,-1):
				z0=z+1
				if z0<=0:
					z0+=q3
				if z0>m3:
					z0-=q3
				R[count]=np.array([x0,y0,z0])
				count+=1
	return R
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
	# print('index i',i)
	# print('supercell',(i[0]-1)*prim_latt_vecs[0]+(i[1]-1)*prim_latt_vecs[1]+(i[2]-1)*prim_latt_vecs[2])
	R = (i[0]-1)*prim_latt_vecs[0]+(i[1]-1)*prim_latt_vecs[1]+(i[2]-1)*prim_latt_vecs[2]+atom_prim_cart[jatom]
	return R



# Example usage
if __name__ == "__main__":
	# l=5
	# a = np.array([l, l, l])*2
	# lat_vec = np.array([[0,l, l], [l, 0, l], [l, l, 0]])*3
	# b = min_images_brute_force(a, lat_vec)
	# print(b)
	with open('grid.dat','r') as f:
		[q1,q2,q3] = np.array(f.readline().split(), dtype=np.int32)
		
	# assue prim_latt_vecs[0,:] is the a vector
	with open('lattice.dat') as f1:
		prim_latt_vecs = np.zeros((3,3))
		for i in range(3):
			prim_latt_vecs[i,:] = np.array(f1.readline().split(), dtype=np.float64)	

	# read position data in the primitive cell
	with open('atoms_in_primitive_cell.1.dat', 'r') as f:
		line = f.readlines()
		atom_prim_cart_pos = np.zeros((len(line),3), dtype=np.float64)
		for i in range(len(line)):
			atom_prim_cart_pos[i,:] = np.array(line[i].split()[2:], dtype=np.float64)
		atom_prim_cart=np.dot(atom_prim_cart_pos, prim_latt_vecs)
		# print('111',atom_prim_cart)
		# print(atom_prim_cart_pos, prim_latt_vecs)
	basis =len(atom_prim_cart); no_prim_cells=q1*q2*q3
	delta_prim = np.zeros((8, 3, basis, basis, no_prim_cells))
	super_latt_vecs=np.array([prim_latt_vecs[0]*q1, prim_latt_vecs[1]*q2, prim_latt_vecs[2]*q3])

	with open('delta_prim_1.dat', 'w') as f:
		for i_atom in range(basis):
		# for i_atom in range(1):
			for j_atom in range(basis):
			# for j_atom in range(4,5):
				delta_r_corr = atom_prim_cart[i_atom] - atom_prim_cart[j_atom]
				for i_cell in range(no_prim_cells):
				# for i_cell in range(1,2):
					temp = atom_super_cart(j_atom,i_cell) - atom_prim_cart[i_atom]
					delta_r_ims, no_im_cells = min_images_brute_force(temp, super_latt_vecs)
					# print(i_atom, j_atom, i_cell)
					# print('sup_atom',atom_super_cart(j_atom,i_cell))
					# print('atom i',atom_prim_cart[i_atom])
					# print('temp',temp)
					# print(delta_r_ims)
					# print(no_im_cells)
					# exit()
					for i_im in range(no_im_cells+1):
						v=delta_prim[i_im,:,i_atom,j_atom,i_cell]
						delta_prim[i_im,:,i_atom,j_atom,i_cell] = delta_r_ims[i_im] + delta_r_corr
						output_str = "{:>12d}{:>12d}{:>12d}{:>12d}{:>21.16f}{:>26.16f}{:>26.16f}\n".format(i_cell+1, i_atom+1, j_atom+1, i_im+1, v[0], v[1], v[2])
						f.write(output_str)
						# print(delta_prim[i_im,:,i_atom,j_atom,i_cell])
						# exit()
						# f.write(f"{i_cell+1} {i_atom+1} {j_atom+1} {i_im+1} ")
						# np.savetxt(f, delta_prim[i_im,:,i_atom,j_atom,i_cell].reshape(1,-1), fmt='%.18e')
		# print(delta_prim)				
		# exit()
	with open('delta_prim_1.dat', 'r') as f1:
		lines = f1.readlines()

	def get_index(line):
		columns = line.split()
		index = tuple(map(int, columns[:3]))
		return index

	# Sort the lines based on the index using the key function
	sorted_lines = sorted(lines, key=get_index)

	# Write the sorted lines to a new file
	with open('delta_prim_2.dat', 'w') as f2:
		f2.writelines(sorted_lines)	   
	
	#f3=np.loadtxt('delta_prim_2.dat')
	#f4=np.loadtxt('delta_prim.dat')
	#for i in range(len(f3)):
	#	if np.linalg.norm(f3[i]-f4[i])>1e-8:
	#		print(i)
	#		print(f3[i])
	#		print(f4[i])
	#		exit()

