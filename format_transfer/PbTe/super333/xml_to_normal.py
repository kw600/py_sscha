import re
import numpy as np

# Define a helper function to extract the indices from the Format 1 block tag
def extract_indices(tag):
    pattern = tag[1:-1].split('.')[1:]
    # print(pattern)
    # match = re.search(pattern, tag)
    [s, s1, m1, m2, m3] = map(int,pattern)
    return s, s1, m1, m2, m3
# i=0
# Open the input and output files

with open('format1.txt', 'r') as f1:
    f = f1.readlines()
    # Read the input file line by line
for i,line in enumerate(f):
        line = line.strip()
        if line.startswith('<NUMBER_OF_ATOMS>'):
             natom=int(line.replace('<',' ').replace(">"," ").split()[1])
        if line.startswith('<MESH_'):
             [q1,q2,q3]=np.array(f[i+1].split(),int)
            #  print(q1,q2,q3,natom)
             break
nq=q1*q2*q3
FC=np.zeros(((nq+1)*9*natom*natom,4))
cell=np.zeros((9*natom*natom,4))
R_index=np.zeros((nq,3))
count=0
for x in range(1,4):
     for y in range(1,4):
          for atom1 in range(1,natom+1):
               for atom2 in range(1,natom+1):
                    cell[count]=np.array([x,y,atom1,atom2])
                    count+=1
# print(cell[121:220])
count=0
for x in range(q1):
     for y in range(q2):
        for z in range(q3):
               R_index[count]=np.array([x+1,y+1,z+1])
               count+=1
# print(cell[0])
count=0
for i in range(9*natom*natom):
     FC[i*(nq+1)]=cell[count]
     count+=1
# print(FC[0:500:nq+1])
for i,line in enumerate(f):
        # Check if the line starts with the Format 1 block tag
        line=line.strip()
        if line.startswith('<s_s'):
            s, s1, m1, m2, m3 = extract_indices(line)
            r1=np.where(np.linalg.norm(R_index-np.array([m1,m2,m3]),axis=1)==0)[0][0]
            r0=np.where(np.linalg.norm(cell-np.array([1,1,s,s1]),axis=1)==0)[0][0]
            D=np.loadtxt(f[i+2:i+5])
            for x in range(1,4):
                 for y in range(1,4):
                    FC[r0*(nq+1)+r1+1]=np.array([m1,m2,m3,D[y-1,x-1]])
                    r0+=natom**2
           

with open('format2.txt', 'w') as f2:
    f2.write(f'{q1} {q2} {q3}\n')
    for i in range(len(FC)):
    # for i in range(500):
        if i%(nq+1)==0:
            a=f'{int(FC[i,3]):.0f}'
            # print(i,FC[i-10:i+10])
            # print(a)
        else:
            a=f'{FC[i,3]:.16E}'
        f2.write(f'{FC[i,0]:.0f} {FC[i,1]:.0f} {FC[i,2]:.0f} {a}\n')
