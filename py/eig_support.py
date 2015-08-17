import numpy as np	# array handling
import scipy.linalg as la	# eigenvalues
import scipy.sparse as sps # Sparse matrices
import scipy.sparse.linalg as spsl
import numpy.ma as ma
import math


class map_1to1:
	# 1:1 Mapping
	def __init__(self):
		self.d = {}
	def add(self, k, v):
		self.d[k] = v
		self.d[v] = k
	def remove(self,k):
		self.d.pop(self.d.pop(k))
	def get(self,k):
		return self.d[k]


def parse_hgr(infile_name, delim=" ", index_offset=0):
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
				i_cell = int(line_arr[i])-index_offset # -1 since hgr file starts with index 1, and python starts with index 0
				Q[i_cell,i_cell] = Q[i_cell,i_cell] + net_weight*(net_cells-1) # Add up connected net weights on diagonals (constructing D)

				for j in range(i+1,len(line_arr)):
					j_cell = int(line_arr[j])-index_offset
					Q[i_cell,j_cell] = Q[i_cell,j_cell] - net_weight
					Q[j_cell,i_cell] = Q[j_cell,i_cell] - net_weight

	return Q

def parse_hgr_sparse(infile_name, delim=" ", index_offset=0):
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


	## Process the HGR file
	# Declare empty arrays to hold index and cell data for COO sparse matrix construction
	# Diagonals
	di = []
	ds = []

	# Adjacency
	a_i = []
	a_j = []
	a_s = []

	# Overall
	qi = []
	qj = []
	qs = []

	for net in range(1,len(lines)):	# iterate over all nets (ignore first, as that just gives us number of cells and nets)
		line_arr = lines[net].split(delim)	# Break up
		
		if(len(line_arr) > 1): # Only count nets that have more than one node
			net_cells = len(line_arr)
			net_weight = 1/(net_cells-1)
			#print(net_weight)

			for i in range(len(line_arr)):
				i_cell = int(line_arr[i]) - index_offset # subtract 1 if using matlab indexing. -1 since hgr file starts with index 1, and python starts with index 0
				#Q[i_cell,i_cell] = Q[i_cell,i_cell] + net_weight*(net_cells-1) # Add up connected net weights on diagonals (constructing D)
				#ds.append(net_weight*(net_cells-1)) # Duplicate ij entries will be accumulated during matrix construction, so we don't need to worry about them right now
				#di.append(i_cell)
				# Diagonal, so i and j are the same here
				ds.append(net_weight*(net_cells-1))
				di.append(i_cell)
				#dj.append(i_cell) # Not necessary since symmetric
				qs.append(net_weight*(net_cells-1))
				qi.append(i_cell)
				qj.append(i_cell)

				for j in range(i+1,len(line_arr)):
					j_cell = int(line_arr[j]) - index_offset
					qi.append(i_cell)
					qj.append(j_cell)
					qs.append(-net_weight)
					
					a_i.append(i_cell)
					a_j.append(j_cell)
					a_s.append(net_weight)

					# Matrix is symmetric!
					qi.append(j_cell)
					qj.append(i_cell)
					qs.append(-net_weight)
					
					a_i.append(j_cell)
					a_j.append(i_cell)
					a_s.append(net_weight)

	#Q = sps.coo_matrix( (qs,(qi, qj)), shape=(num_cells,num_cells) )
	#Q = Q.tocsc()

	D = sps.coo_matrix( (ds, (di, di)), shape=(num_cells,num_cells) )
	A = sps.coo_matrix( (a_s, (a_i, a_j)), shape=(num_cells,num_cells) )
	D = D.tocsc()
	A = A.tocsc()
	Q = D - A

	return (Q, D, A)



def partition_1d(Q, eigenval_cutoff=1e-5, num_eigs = 10):
	# Finding more eigs than we need because sparse eigsh sometimes gives spurious eigs
	(vals, vecs) = spsl.eigsh(Q, k=num_eigs, which='LM', sigma=-1) # Q is guaranteed to be hermitian since it is a real symmetric matrix

	# Sparse eigenvalue solver sometimes gives us spurious small eigs. We need to sift them out before we do any processing
	sorted_vals_raw = np.argsort(vals)
	sorted_vals = []
	sorted_vals.append( np.min(np.abs(vals)) ) # First eigenvalue of the laplacian matrix will always be 0. spe
	for el in sorted_vals_raw:
		if (vals[el] > eigenval_cutoff):
			sorted_vals.append(vals[el])

	if ( len(sorted_vals) == 1):
		print("partition_1d: Error! Not enough valid eigenvalues. Rerun with num_eigs > 10")
	
	val2 = vals[sorted_vals[1]]
	vec2 = vecs[:,sorted_vals[1]]

	partition1d = np.argsort(vec2)
	return (partition1d, sorted_vals, vals)


def calc_cutsize_bipart(adjacency_mat, partition_order, area_balance):
	# Reorder the matrix in partition order
	adj_reord = adjacency_mat[:, partition_order]
	adj_reord = adj_reord[partition_order,:]

	# [FIX] Need to properly work out indexing for array slice
	if (area_balance > 0.5):
		area_balance = 1 - area_balance
		
	start_ind = int(round(area_balance * (len(partition_order) - 1) ))
	stop_ind = (len(partition_order)-1) - start_ind
	
	split_inds = range(start_ind, stop_ind+1)
	
	cutsize_vec = np.zeros(len(split_inds))
	for (idx, split_ind) in enumerate(split_inds):
		# top right or bot left quadrant represent interconnections across tiers
		# if split_idn is index of first element in second partition
		#	then top right quadrant is rows 0:split_ind-1 (i.e. not including split_ind)
		#	and cols split_ind:end (including split_ind and end)
		cutsize_vec[idx] = adj_reord[0:split_ind, split_ind:].sum()

	mincut_ind = split_inds[np.argmin(cutsize_vec)]
	mincut_val = np.min(cutsize_vec)
	
	num_elements = np.shape(adjacency_mat)[0]
	split_ind_arr = np.array(split_inds)
	normcut_vec = cutsize_vec/(split_ind_arr*(num_elements-split_ind_arr))
	norm_mincut_ind = split_inds[np.argmin(normcut_vec)]
	norm_mincut_val = np.min(normcut_vec)
	
	p1_size_frac_vec = np.array(split_inds)/len(partition_order)

	return (mincut_val, mincut_ind, cutsize_vec, norm_mincut_val, norm_mincut_ind, normcut_vec, p1_size_frac_vec)

	
def partition_1d_perturbed(Q, Qp, eigenval_cutoff=1e-5, num_eigs_solve=10, num_eigs_corr=5):
	#(vals, vecs) = la.eigh(Q) # Q is guaranteed to be hermitian since it is a real symmetric matrix
	
	(vals, vecs) = spsl.eigsh(Q, k=num_eigs_solve, which="LM", sigma=-1)	
	sorted_vals = np.argsort(vals)
	
	# Sparse eigenvalue solver sometimes gives us spurious small eigs. We need to sift them out before we do any processing
	sorted_vals_raw = np.argsort(vals)
	sorted_vals = []
	sorted_vals.append( np.min(np.abs(vals)) ) # First eigenvalue of the laplacian matrix will always be 0. spe
	for el in sorted_vals_raw:
		if (vals[el] > eigenval_cutoff):
			sorted_vals.append(vals[el])

	if ( len(sorted_vals) == 1):
		print("partition_1d: Error! Not enough valid eigenvalues. Rerun with num_eigs > 10")


	#vals = np.matrix(vals)
	#vecs = np.matrix(vecs)
	
	val2 = vals[sorted_vals[1]]
	vec2 = vecs[:,sorted_vals[1]]

	vecs_to_include = [0] + list(range(2,num_eigs_corr))

	vec_approx = np.zeros((np.size(vec2),1))
	Qp = np.array(Qp)
	
	for i in vecs_to_include:
		vv = np.zeros( (np.size(vec2),1) )
		vv[:,0] = vecs[:,sorted_vals[i]]
		
		# [FIX] Issue with calculation of top. Qp has no shape, could be the problem
		print( np.shape(vv) )
		print( np.shape(Qp) )
		print( np.shape(vec2) )
		print( Qp )

		#print((vv.T).dot( Qp.dot(vec2) ))
		top = float( (vv.T).dot( Qp.dot(vec2) ) )
		bot = val2 - vals[sorted_vals[i]]
		vec_approx = vec_approx + (top/bot)*vv

	partition1d = np.argsort((vec_approx.T))
	partition1d = list(partition1d[0]) # Converting this back to a list. the array access is because the NP ndarray is 10x1 and
									   # we just want a normal "dimensionless" 10 element list
	
	return partition1d


def construct_component_map(dict_file_name):
	component_map = map_1to1()
	with open(dict_file_name,'r') as infile:
		for line in infile:
			arr = line.split()
			component_map.add(arr[0], int(arr[1]))

	return component_map


def write_partitions(block_order, component_1to1, outfile_name, num_partitions=2):
	num_cells = math.floor(len(block_order)/num_partitions)
	partition_list = []
	start_cell = 0
	stop_cell = start_cell + num_cells
	for ind in range(num_partitions):
		if (ind == num_partitions - 1):
			stop_cell = len(block_order)
		partition_list.append( block_order[start_cell:stop_cell])

		start_cell = stop_cell
		stop_cell = start_cell + num_cells


	ostr = "{0:s} {1:d}\n"
	with open(outfile_name,'w') as outfile:
		for pind, partition in enumerate(partition_list):
			for component_id in partition:
				component_name = component_1to1.get(component_id)
				outfile.write(ostr.format(component_name, pind) )

	return




def compare_lists(A,B,split_ratio):
	# take in two lists of equal length
	# split them into two partitions based on split ratio
	# return number of mismatched elements between first half of A and either half of B

	# determine index to split at
	split_ind = int(round(split_ratio*len(A)))
	A1 = set(A[0:split_ind])
	A2 = set(A[split_ind:])
	B1 = set(B[0:split_ind])
	B2 = set(B[split_ind:])

	first_half_diff = A1.difference(B1)
	sec_half_diff = A1.difference(B2)

	mismatch_ratio_1 = len(first_half_diff)/len(A1)
	mismatch_ratio_2 = len(sec_half_diff)/len(A1)

	return (mismatch_ratio_1,mismatch_ratio_2)

def is_symmetric(m):
    """Check if a sparse matrix is symmetric

    Parameters
    ----------
    m : array or sparse matrix
        A square matrix.

    Returns
    -------
    check : bool
        The check result.

    """
    if m.shape[0] != m.shape[1]:
        raise ValueError('m must be a square matrix')

    if not isinstance(m, coo_matrix):
        m = coo_matrix(m)

    r, c, v = m.row, m.col, m.data
    tril_no_diag = r > c
    triu_no_diag = c > r

    if triu_no_diag.sum() != tril_no_diag.sum():
        return False

    rl = r[tril_no_diag]
    cl = c[tril_no_diag]
    vl = v[tril_no_diag]
    ru = r[triu_no_diag]
    cu = c[triu_no_diag]
    vu = v[triu_no_diag]

    sortl = np.lexsort((cl, rl))
    sortu = np.lexsort((ru, cu))
    vl = vl[sortl]
    vu = vu[sortu]

    check = np.allclose(vl, vu)

    return check
	
