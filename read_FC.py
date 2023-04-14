import numpy as np


def check_integer(line):
	integer_count = 0
	for a in line.split():
		try:
			int(a)
			integer_count += 1
		except ValueError:
			pass
	return integer_count

def read_FC(file):
	FC=np.zeros((3,3,3,6,6))
	f=open(file,"r")
	lines=f.readlines()
	for l in lines:
		if check_integer(l)==3 and len(l.split())==3:
			nq=int(l.split()[0])*int(l.split()[1])*int(l.split()[2])
			break
	for l1 in range(len(lines)):
		if check_integer(lines[l1])==4:
			# print(lines[l1])
			[l,k,m,n]=map(int, lines[l1].split())
			m1,n1=(m-1)*3+l,(n-1)*3+k
			for i in range(nq):
				# print(l1,i)
				# print(lines[l1+i+1].split())
				
				[a,b,c]=map(int,lines[l1+i+1].split()[:3])
				# print(m,n,l,k)
				# print(a,b,c,m1-1,n1-1)
				FC[a,b,c,m1-1,n1-1]=lines[l1+i+1].split()[3]
	return FC
					
f=read_FC('fc_harmonic_dyn')
print(f[1,1,1])
D=np.zeros((6,6))
for i in range(1,3):
	for j in range(1,3):
		for k in range(1,3):
			D+=f[i,j,k,:,:]
			# print(f[i,j,k])
D = D/8
print(np.linalg.eig(D))