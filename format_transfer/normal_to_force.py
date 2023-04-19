import numpy as np

def check_integer(lines):
	count=0
	#check there are 3 integers in the lines
	for l in lines:
		try:
			int(l)
			count+=1
		except ValueError:
			pass
	return count


with open('format2.txt', 'r') as f:
	lines = f.readlines()
for i in range(len(lines)):
	l = lines[i].split()
	#check there are 3 integers in the lines
	if len(l) == 3 and check_integer(l)==3:
		[q1,q2,q3] = map(int,l)
		start=i


data = lines[start+1:]
fc=np.zeros((int(len(data)*8/9),6))
f1 = open('format3.txt','w')
Cell=np.zeros((q1,q2,q3),int)
cell=0
for x in range(q1):
	for y in range(q2):
		for z in range(q3):
				cell+=1
				Cell[x,y,z]=cell
				
R_index=np.zeros((q1*q2*q3,3))

# count=0
# for x in range(q1):
# 	for y in range(q2):
# 		for z in range(0,-q3,-1):
# 			z0=z+1
# 			if z0<=0:
# 				z0+=6
# 			R_index[count]=np.array([x+1,y+1,z0])
# 			count+=1

count=0
for x in range(0,-q1,-1):
	x0=x+1
	if x0<=0:
		x0+=q1
	for y in range(0,-q2,-1):
		y0=y+1
		if y0<=0:
			y0+=q2
		for z in range(0,-q3,-1):
			z0=z+1
			if z0<=0:
				z0+=q3
			R_index[count]=np.array([x0,y0,z0])
			count+=1
# print(R_index[:12])

for i in data:
	num=check_integer(i.split())
	if num==4:
		[q1,q2,q3,q4] = map(int,i.split())
	elif num==3:
		[r1,r2,r3] = map(int,i.split()[:3])
		# if r1==1 and r2==1 and r3==1 and 
		index=np.where(np.linalg.norm(R_index-np.array([r1,r2,r3]),axis=1)==0)[0][0]
		if float(i.split()[3])>0:
			f1.write(f'           {index+1}           {q3}           {q1}           {q4}           {q2}  {-1*float(i.split()[3])/2:.16E}\n')
		else:
			f1.write(f'           {index+1}           {q3}           {q1}           {q4}           {q2}   {-1*float(i.split()[3])/2:.16E}\n')
f.close()

# Open the file and read its contents
with open('format3.txt', 'r') as f1:
	lines = f1.readlines()

def get_index(line):
	columns = line.split()
	index = tuple(map(int, columns[:5]))
	return index

# Sort the lines based on the index using the key function
sorted_lines = sorted(lines, key=get_index)

# Write the sorted lines to a new file
with open('format4_1.txt', 'w') as f2:
	f2.writelines(sorted_lines)
if __name__=='__main__':
	pass
	f3=np.loadtxt('format4_1.txt')
	f4=np.loadtxt('force.dat')
	for i in range(len(f3)):
		if np.linalg.norm(f3[i]-f4[i])>1e-8:
			print(i)
			print(f3[i])
			print(f4[i])
			exit()
	

