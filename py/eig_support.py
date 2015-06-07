import re	# regex
import sys	# for grabbing input args
import numpy as np	# array handling
import numpy.linalg as la	# eigenvalues

def parse_hgr(infile_name,delim):
	# Read input file
	infile = open(infile_name,'r')
	lines = infile.readlines()
	infile.close()

	# strip off trailing whitespace
	for i in range(len(lines)):
		lines[i] = lines[i].rstrip()

	# get number of nets and cells
	n = lines[0].split(delim)
	num_nets = int(n[0])
	num_cells = int(n[1])

	# [FIX] Use sparse matrix handling
	Q = np.zeros((num_cells,num_cells)) # initialize laplacian matrix

	## Process the HGR file
	for net in range(1,len(lines)):	# iterate over all nets (ignore first, as that just gives us number of cells and nets)
		line_arr = lines[net].split(delim)	# Break up
		
		if(len(line_arr) > 1):
			net_cells = len(line_arr)
			net_weight = 1/(net_cells-1)
			#print(net_weight)

			for i in range(len(line_arr)):
				i_cell = int(line_arr[i])-1 # -1 since hgr file starts with index 1, and python starts with index 0
				Q[i_cell,i_cell] = Q[i_cell,i_cell] + net_weight*(net_cells-1) # Add up connected net weights on diagonals (constructing D)

				for j in range(i+1,len(line_arr)):
					j_cell = int(line_arr[j])-1
					Q[i_cell,j_cell] = Q[i_cell,j_cell] - net_weight
					Q[j_cell,i_cell] = Q[j_cell,i_cell] - net_weight

	return Q

def partition_1d(Q):
	(vals, vecs) = la.eigh(Q) # Q is guaranteed to be hermitian since it is a real symmetric matrix
	sorted_vals = np.argsort(vals)
	
	val2 = vals[sorted_vals[1]]
	vec2 = vecs[:,sorted_vals[1]]

	partition1d = np.argsort(vec2)

	return partition1d
	
def partition_1d_perturbed(Q,Qp,n):
	(vals, vecs) = la.eigh(Q) # Q is guaranteed to be hermitian since it is a real symmetric matrix
	sorted_vals = np.argsort(vals)

	#vals = np.matrix(vals)
	vecs = np.matrix(vecs)
	
	val2 = vals[sorted_vals[1]]
	vec2 = vecs[:,sorted_vals[1]]

	vecs_to_include = [0] + list(range(2,n))

	vec_approx = np.matrix(np.zeros((np.size(vec2),1)))
	Qp = np.matrix(Qp)

	for i in vecs_to_include:
		vv = vecs[:,sorted_vals[i]]

		#top = vv.T * (Qp * vec2)
		top = Qp * vec2
		top = vv.T*top
		top = float(top)
		bot = val2 - vals[sorted_vals[i]]
		vec_approx = vec_approx + (top*vv)/bot

	partition1d = np.argsort(list(vec_approx.T))
	
	return partition1d

def partition_1d_perturbed_alt(Q,Qp,n):
	(vals, vecs) = la.eigh(Q) # Q is guaranteed to be hermitian since it is a real symmetric matrix
	sorted_vals = np.argsort(vals)
	
	val2 = vals[sorted_vals[1]]
	vec2 = vecs[:,sorted_vals[1]]

	vecs_to_include = [0] + list(range(2,n))

	vec_approx = np.matrix(np.zeros((np.size(vec2),1)))
	Qp = np.matrix(Qp)

	for i in vecs_to_include:
		vv = vecs[:,sorted_vals[i]]

		#top = vv.T * (Qp * vec2)
		top = Qp * vec2
		top = vv.T*top
		top = float(top)
		bot = val2 - vals[sorted_vals[i]]
		vec_approx = vec_approx + (top*vv)/bot

	partition1d = np.argsort(list(vec_approx.T))
	
	return partition1d


def compare_lists(A,B,split_ratio):
	# take in two lists of equal length
	# split them into two partitions based on split ratio
	# return number of mismatched elements between first half of A and either half of B

	# determine index to split at
	split_ind = int(round(split_ratio*len(A)))
	A1 = set(A[0:split_ind])
	B1 = set(B[0:split_ind])
	B2 = set(B[split_ind:])

	first_half_diff = A1.difference(B1)
	sec_half_diff = A1.difference(B2)

	print(A)
	print(A1)
	print(split_ind)
	mismatch_ratio_1 = len(first_half_diff)/len(A1)
	mismatch_ratio_2 = len(sec_half_diff)/len(A2)

	return (mismatch_ratio_1,mismatch_ratio_2)


	

