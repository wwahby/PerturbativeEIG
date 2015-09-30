import numpy as np	# array handling
#import scipy.linalg as la	# eigenvalues
import scipy.sparse as sps # Sparse matrices
import scipy.sparse.linalg as spsl
import numpy.linalg as npla
#import numpy.ma as ma
import math
import pylab as pl
import time # For function timing
import os


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
	time_start = time.clock()
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
	time_stop = time.clock()
	time_elapsed = time_stop - time_start

	return (Q, time_elapsed)

def parse_hgr_sparse(infile_name, delim=" ", index_offset=0):
	time_start = time.clock()
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

	time_stop = time.clock()

	time_elapsed = time_stop - time_start;

	return (Q, D, A, time_elapsed)


def get_sorted_eigvals(vals, vecs, eigenval_cutoff):
	# Sparse eigenvalue solver sometimes gives us spurious small eigs. We need to sift them out before we do any processing

	sorted_val_inds = np.argsort(vals)
	sorted_vals = []
	sorted_vals.append( np.min(np.abs(vals)) ) # First eigenvalue of the laplacian matrix will always be 0. spe
	for svi_idx, el in enumerate(sorted_val_inds):
		if (vals[el] >= eigenval_cutoff) and ( np.isreal(vals[el]) and (svi_idx > 0) ): # skip first eigenvalue since we're using that as the zero eig
			real_part = np.real(vals[el])
			sorted_vals.append( real_part  ) # stripping off complex part (which should be zero as per conditional above)

	if ( len(sorted_vals) == 1):
		print("place_1d: Error! Not enough valid eigenvalues. Increase num_eigs and rerun!")
	eig2_ind = int(sorted_vals[1]) # index of second eigenvalue/eigenvector
	val2 = vals[ eig2_ind ]
	vec2 = vecs[:, eig2_ind ]

	return (sorted_vals, eig2_ind, val2, vec2)



def place_1d(Q, eigenval_cutoff=1e-5, num_eigs = 10):
	time_start = time.clock()
	# Finding more eigs than we need because sparse eigsh sometimes gives spurious eigs
	(vals, vecs) = spsl.eigsh(Q, k=num_eigs, which='SM') # Q is guaranteed to be hermitian since it is a real symmetric matrix
	#(vals, vecs) = npla.eigh(Q.todense()) # Q is guaranteed to be hermitian since it is a real symmetric matrix

	(sorted_vals, eig2_ind, val2, vec2) = get_sorted_eigvals(vals, vecs, eigenval_cutoff)

	placement_order_1d = np.argsort(vec2)

	time_stop = time.clock()
	time_elapsed = time_stop - time_start

	# Checking for non-orthnormal eigenvectors
#	shape = np.shape(vecs)
#	num_vecs = shape[1]
#	for vec_ind in range(num_vecs):
#		vec = vecs[:, vec_ind]
#		prod = vec.T.dot( vec )
#		print("Ind: {0:d} \t Val: {1:.3g} \t Prod_err: {2:.3g}".format(vec_ind, vals[vec_ind], 1 - abs(prod) ) )


	return (placement_order_1d, sorted_vals, vals, vecs, time_elapsed)


def calc_cutsize_bipart(adjacency_mat, placement_order, area_balance):
	time_start = time.clock()

	# Reorder the matrix in partition order
	adj_reord = adjacency_mat[:, placement_order]
	adj_reord = adj_reord[placement_order,:]

	# [FIX] Need to properly work out indexing for array slice
	if (area_balance > 0.5):
		area_balance = 1 - area_balance

	start_ind = int(round(area_balance * (len(placement_order) - 1) ))
	stop_ind = (len(placement_order)-1) - start_ind

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

	p1_size_frac_vec = np.array(split_inds)/len(placement_order)
	time_stop = time.clock()

	time_elapsed = time_stop - time_start
	(skew_size_vec, skew_cut_vec) = convert_cutsize_to_skew( p1_size_frac_vec, cutsize_vec)
	(skew_size_vec, skew_normcut_vec) = convert_cutsize_to_skew( p1_size_frac_vec, normcut_vec)

	return (mincut_val, mincut_ind, cutsize_vec, norm_mincut_val, norm_mincut_ind, normcut_vec, p1_size_frac_vec, skew_size_vec, skew_cut_vec, skew_normcut_vec, time_elapsed)


def convert_cutsize_to_skew( size_frac_vec, cut_info_vec):
	# Make sure we have the right number of elements in the skew vector
	if (len(size_frac_vec) % 2 == 0):
		skew_length = len(size_frac_vec)/2
		mid_ind = int(skew_length)
		skew_cut_vec = np.minimum( cut_info_vec[0:mid_ind] , cut_info_vec[-1:mid_ind-1:-1] )
		skew_size_vec = size_frac_vec[0:mid_ind]
	else:
		skew_length = (len(size_frac_vec)-1)/2 + 1
		mid_ind = int(skew_length)-1
		skew_cut_vec = np.minimum( cut_info_vec[0:mid_ind] , cut_info_vec[-1:mid_ind:-1] )
		skew_cut_vec = np.append( skew_cut_vec, cut_info_vec[mid_ind] )
		skew_size_vec = size_frac_vec[0:mid_ind+1]

	# Reorient the vectors
	skew_size_vec = 0.5 - skew_size_vec
	skew_size_vec = np.flipud(skew_size_vec)
	skew_cut_vec = np.flipud(skew_cut_vec)

	return ( skew_size_vec, skew_cut_vec)


def place_1d_perturbed(vals, vecs, Qp, eigenval_cutoff=1e-5, num_eigs_solve=10, num_eigs_corr=5):
	#(vals, vecs) = la.eigh(Q) # Q is guaranteed to be hermitian since it is a real symmetric matrix
	time_start = time.clock()

	#(vals, vecs) = spsl.eigsh(Q, k=num_eigs_solve, which="LM", sigma=-1)
	sorted_vals = np.argsort(vals)

	# Sparse eigenvalue solver sometimes gives us spurious small eigs. We need to sift them out before we do any processing
	(sorted_vals, eig2_ind, val2, vec2) = get_sorted_eigvals(vals, vecs, eigenval_cutoff)

	time_eig_standard = time.clock()

	if ( len(sorted_vals) == 1):
		print("place_1d: Error! Not enough valid eigenvalues. Rerun with num_eigs > 10")


	#vals = np.matrix(vals)
	#vecs = np.matrix(vecs)
	mat_size = np.shape(Qp)
	vec_length = mat_size[0] # should be symmetric
	eig2_ind = int(sorted_vals[1]) # index of second eigenvector/eigenvalue
	val2 = vals[eig2_ind]
	vec2 = np.zeros((vec_length,1))
	vec2[:,0] = vecs[:,eig2_ind]

	vecs_to_include = [0] + list(range(2,num_eigs_corr))

	vec_approx = np.zeros((np.size(vec2),1))
	#Qp = np.array(Qp)

	for i in vecs_to_include:
		vv = np.zeros( (np.size(vec2),1) )
		eigi_ind = int(sorted_vals[i]) # index of i-th eigenvector/eigenvalue
		vv[:,0] = vecs[:, eigi_ind ]

		top = float( (vv.T).dot( Qp.dot(vec2) ) )
		bot = val2 - vals[ eigi_ind ]
		vec_approx = vec_approx + (top/bot)*vv

	placement_order_1d = np.argsort((vec_approx.T))
	placement_order_1d = list(placement_order_1d[0]) # Converting this back to a list. the array access is because the NP ndarray is 10x1 and
									   # we just want a normal "dimensionless" 10 element list
	time_stop = time.clock()
	#time_elapsed_tot = time_stop - time_start
	time_elapsed_eig_standard = time_eig_standard - time_start
	time_elapsed_place1d_perturbed = time_stop - time_eig_standard

	return (placement_order_1d, time_elapsed_place1d_perturbed, time_elapsed_eig_standard)


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


def ecdf(data, down_sampling_step=1, type="end"):

	if (type == "mid"):
		adjustment = 0.5
	elif (type == "end"):
		adjustment = 1
	elif (type == "start"):
		adjustment = 0
	else:
		raise ValueError("Valid types are mid, end, and start")

	data_sorted = np.sort( data )
	ecdf_vec = (np.arange( len(data_sorted) ) + adjustment) / float(len(data_sorted))

	return( ecdf_vec, data_sorted )


def compare_partition_schemes(hgr_filename, perturbed_filename, repetitions=1, num_eigs=20):

	cur_dir = os.getcwd()
	base_dir = os.path.dirname( cur_dir )
	netlist_dir = os.path.join(base_dir, "netlists")

	#hgr_filename = "p2.hgr"
	#perturbed_filename = "p2.hgr"

	hgr = os.path.join( netlist_dir, hgr_filename )
	perturbed = os.path.join( netlist_dir, perturbed_filename )

	#num_eigs = 20
	delim =  " "
	num_partitions = 2
	eigenval_cutoff = 0
	index_offset = 1
	unbalance_factor = 5

	infile_name = hgr # Reconstruct filenames with spaces
	infile_p_name = perturbed # Reconstruct filenames with spaces
	run_perturbed = True

	area_balance = 0.5 - unbalance_factor/100

	# Preallocate time vectors
	time_parse_exact_vec = np.zeros(repetitions)
	time_parse_perturbed_exact_vec = np.zeros(repetitions)
	time_place_exact_vec = np.zeros(repetitions)
	time_place_perturbed_approx_vec = np.zeros(repetitions)
	time_place_perturbed_exact_vec = np.zeros(repetitions)
	time_cutsize_orig_vec = np.zeros(repetitions)
	time_cutsize_perturbed_approx_vec = np.zeros(repetitions)
	time_cutsize_perturbed_exact_vec = np.zeros(repetitions)
	time_construct_laplacian_perturbed_vec = np.zeros(repetitions)
	time_total_exact_vec = np.zeros(repetitions)
	time_total_perturbed_approx_vec = np.zeros(repetitions)
	time_total_perturbed_exact_vec = np.zeros(repetitions)

	# Preallocate cutsize vectors
	mincut_val_orig_vec = np.zeros(repetitions)
	mincut_val_perturbed_exact_vec = np.zeros(repetitions)
	mincut_val_perturbed_approx_vec= np.zeros(repetitions)

	mincut_skew_orig_vec = np.zeros(repetitions)
	mincut_skew_perturbed_exact_vec = np.zeros(repetitions)
	mincut_skew_perturbed_approx_vec = np.zeros(repetitions)

	# Preallocate ratiocut vectors
	ratiocut_val_orig_vec = np.zeros(repetitions)
	ratiocut_val_perturbed_exact_vec = np.zeros(repetitions)
	ratiocut_val_perturbed_approx_vec = np.zeros(repetitions)

	ratiocut_skew_orig_vec = np.zeros(repetitions)
	ratiocut_skew_perturbed_exact_vec = np.zeros(repetitions)
	ratiocut_skew_perturbed_approx_vec = np.zeros(repetitions)


	for rep_ind in range(repetitions):
		print("\nIteration {0:d}\n=====================".format(rep_ind) )
		print("Parsing original graph...")
		(Q, D, A, time_parse_exact) = parse_hgr_sparse(infile_name,delim=delim, index_offset=index_offset)

		print("Constructing optimal 1D placement for original graph...")
		(placement_order, eigvals, raw_eigvals, raw_eigvecs, time_place_exact) = place_1d(Q, eigenval_cutoff = eigenval_cutoff, num_eigs = num_eigs)

		print("Finding cutsize of original system...")
		(mincut_val_exact, mincut_ind, cutsize_vec, normcut_ind, normcut_val, normcut_vec, p1_size_frac_vec, skew_size_vec, skew_cut_vec, skew_normcut_vec, time_cutsize_exact) = calc_cutsize_bipart(A, placement_order, area_balance)


		# Perturbed Solutions
		time_start_construct_laplacian_perturbed = time.clock()
		(Qp_exact, Dp_e, Ap_e, time_parse_perturbed_exact) = parse_hgr_sparse(infile_p_name, delim=delim, index_offset=index_offset)
		Qp = Qp_exact - Q
		time_stop_construct_laplacian_perturbed = time.clock()
		time_construct_laplacian_perturbed = time_stop_construct_laplacian_perturbed - time_start_construct_laplacian_perturbed

		print("Constructing approximate 1D placement for perturbed system...")
		(p1dp, time_place_perturbed_approx, time_eig_standard) = place_1d_perturbed(raw_eigvals, raw_eigvecs, Qp, eigenval_cutoff = eigenval_cutoff, num_eigs_solve = num_eigs)

		print("Constructing optimal 1D placement for perturbed system...")
		(p1dpe, vals_sorted_pe, vals_raw_pe, vecs_raw_pe, time_place_perturbed_exact)  = place_1d(Qp_exact, eigenval_cutoff = eigenval_cutoff, num_eigs = num_eigs)

		print("Finding cutsize of approximately perturbed system...")
		(mincut_val_perturbed_approx, mincut_ind_p, cutsize_vec_p, normcut_ind_p, normcut_val_p, normcut_vec_p, p1_size_frac_vec_p, skew_size_vec_p, skew_cut_vec_p, skew_normcut_vec_p, time_cutsize_perturbed_approx) = calc_cutsize_bipart(Ap_e, p1dp, area_balance) # uses exact adjacency matrix, since that part can be known just based on connectivity, without actually solving eig problem

		print("Finding cutsize of exact solution to perturbed system...")
		(mincut_val_perturbed_exact, mincut_ind_pe, cutsize_vec_pe, normcut_ind_pe, normcut_val_pe, normcut_vec_pe, p1_size_frac_vec_pe, skew_size_vec_pe, skew_cut_vec_pe, skew_normcut_vec_pe, time_cutsize_perturbed_exact) = calc_cutsize_bipart(Ap_e, p1dpe, area_balance)

		# Calc times
		time_total_exact = time_parse_exact + time_place_exact + time_cutsize_exact
		time_total_perturbed_exact = time_parse_perturbed_exact + time_place_perturbed_exact + time_cutsize_perturbed_exact
		time_total_perturbed_approx = time_parse_perturbed_exact + time_construct_laplacian_perturbed + time_place_perturbed_approx + time_cutsize_perturbed_approx

		# Get skew values for mincut and ratiocut for each case
		skew_cut_orig_id = np.argmin( skew_cut_vec)
		skew_ratiocut_orig_id = np.argmin( skew_normcut_vec)
		skew_cut_orig = skew_size_vec[skew_cut_orig_id]
		skew_ratiocut_orig = skew_size_vec[skew_ratiocut_orig_id]

		skew_cut_perturbed_exact_id = np.argmin( skew_cut_vec_pe)
		skew_ratiocut_perturbed_exact_id = np.argmin( skew_normcut_vec_pe)
		skew_cut_perturbed_exact = skew_size_vec[skew_cut_perturbed_exact_id]
		skew_ratiocut_perturbed_exact = skew_size_vec[skew_ratiocut_perturbed_exact_id]

		skew_cut_perturbed_approx_id = np.argmin( skew_cut_vec_p)
		skew_ratiocut_perturbed_approx_id = np.argmin( skew_normcut_vec_p)
		skew_cut_perturbed_approx = skew_size_vec[skew_cut_perturbed_approx_id]
		skew_ratiocut_perturbed_approx = skew_size_vec[skew_ratiocut_perturbed_approx_id]

		# Fill in the various vectors
		time_parse_exact_vec[rep_ind] = time_parse_exact
		time_parse_perturbed_exact_vec[rep_ind] = time_parse_perturbed_exact
		time_place_exact_vec[rep_ind] = time_place_exact
		time_place_perturbed_approx_vec[rep_ind] = time_place_perturbed_approx
		time_place_perturbed_exact_vec[rep_ind] = time_place_perturbed_exact
		time_cutsize_orig_vec[rep_ind] = time_parse_exact
		time_cutsize_perturbed_approx_vec[rep_ind] = time_cutsize_perturbed_approx
		time_cutsize_perturbed_exact_vec[rep_ind] = time_cutsize_perturbed_exact
		time_construct_laplacian_perturbed_vec[rep_ind] = time_construct_laplacian_perturbed
		time_total_exact_vec[rep_ind] = time_total_exact
		time_total_perturbed_approx_vec[rep_ind] = time_total_perturbed_approx
		time_total_perturbed_exact_vec[rep_ind] = time_total_perturbed_exact

		# Fill cutsize vectors
		mincut_val_orig_vec[rep_ind] = mincut_val_exact
		mincut_val_perturbed_exact_vec[rep_ind] = mincut_val_perturbed_exact
		mincut_val_perturbed_approx_vec[rep_ind] = mincut_val_perturbed_approx

		mincut_skew_orig_vec[rep_ind] = skew_cut_orig
		mincut_skew_perturbed_exact_vec[rep_ind] = skew_cut_perturbed_exact
		mincut_skew_perturbed_approx_vec[rep_ind] = skew_cut_perturbed_approx

		# Fill Ratiocut vectors
		ratiocut_val_orig_vec[rep_ind] = normcut_val
		ratiocut_val_perturbed_exact_vec[rep_ind] = normcut_val_pe
		ratiocut_val_perturbed_approx_vec[rep_ind] = normcut_val_p

		ratiocut_skew_orig_vec[rep_ind] = skew_ratiocut_orig
		ratiocut_skew_perturbed_exact_vec[rep_ind] = skew_ratiocut_perturbed_exact
		ratiocut_skew_perturbed_approx_vec[rep_ind] = skew_ratiocut_perturbed_approx

	(time_orig_ecdf, time_orig_sorted) = ecdf(time_total_exact_vec)
	(time_pe_ecdf, time_pe_sorted) = ecdf(time_total_perturbed_exact_vec)
	(time_pa_ecdf, time_pa_sorted) = ecdf(time_total_perturbed_approx_vec)

	(mincut_orig_ecdf, mincut_orig_sorted) = ecdf(mincut_val_orig_vec)
	(mincut_pe_ecdf, mincut_pe_sorted) = ecdf(mincut_val_perturbed_exact_vec)
	(mincut_pa_ecdf, mincut_pa_sorted) = ecdf(mincut_val_perturbed_approx_vec)

	(ratiocut_orig_ecdf, ratiocut_orig_sorted) = ecdf(ratiocut_val_orig_vec)
	(ratiocut_pe_ecdf, ratiocut_pe_sorted) = ecdf(ratiocut_val_perturbed_exact_vec)
	(ratiocut_pa_ecdf, ratiocut_pa_sorted) = ecdf(ratiocut_val_perturbed_approx_vec)


	## Figures ( just plot the last run)
	pl.figure(1)
	pl.hold(True)
	pl.plot(p1_size_frac_vec, cutsize_vec, 'k')
	pl.plot(p1_size_frac_vec_p, cutsize_vec_p, 'r')
	pl.plot(p1_size_frac_vec_pe, cutsize_vec_pe, 'b')
	pl.xlabel('P1 Size Fraction')
	pl.ylabel('Cutsize')

	pl.figure(2)
	pl.hold(True)
	pl.plot(p1_size_frac_vec, normcut_vec, 'k')
	pl.plot(p1_size_frac_vec_p, normcut_vec_p, 'r')
	pl.plot(p1_size_frac_vec_pe, normcut_vec_pe, 'b')
	pl.xlabel('Skew')
	pl.ylabel('P1 Size Fraction')

	pl.figure(3)
	pl.hold(True)
	pl.plot(skew_size_vec, skew_cut_vec, 'k')
	pl.plot(skew_size_vec_p, skew_cut_vec_p, 'r')
	pl.plot(skew_size_vec_pe, skew_cut_vec_pe, 'b')
	pl.xlabel('Skew')
	pl.ylabel('Cutsize')

	pl.figure(4)
	pl.hold(True)
	pl.plot(skew_size_vec, skew_normcut_vec, 'k')
	pl.plot(skew_size_vec_p, skew_normcut_vec_p, 'r')
	pl.plot(skew_size_vec_pe, skew_normcut_vec_pe, 'b')
	pl.xlabel('Skew')
	pl.ylabel('Ratio cut')

	pl.figure(5)
	pl.clf()
	pl.hold(True)
	pl.plot(time_orig_sorted, time_orig_ecdf, 'k')
	pl.plot(time_pe_sorted, time_pe_ecdf, 'b')
	pl.plot(time_pa_sorted, time_pe_ecdf, 'r')
	pl.xlabel('Time(s)')
	pl.ylabel('ECDF')

	pl.figure(6)
	pl.clf()
	pl.hold(True)
	pl.plot(mincut_orig_sorted, mincut_orig_ecdf, 'k')
	pl.plot(mincut_pe_sorted, mincut_pe_ecdf, 'b')
	pl.plot(mincut_pa_sorted, mincut_pe_ecdf, 'r')
	pl.xlabel('Mincut')
	pl.ylabel('ECDF')

	pl.figure(7)
	pl.clf()
	pl.hold(True)
	pl.plot(ratiocut_orig_sorted, ratiocut_orig_ecdf, 'k')
	pl.plot(ratiocut_pe_sorted, ratiocut_pe_ecdf, 'b')
	pl.plot(ratiocut_pa_sorted, ratiocut_pe_ecdf, 'r')
	pl.xlabel('Ratio cut')
	pl.ylabel('ECDF')



	print("{0:>16s}\t{1:<16.3g} \n {2:>16s}\t{3:<16.3g} \n {4:>16s}\t{5:<16.3g} \n ".format("t_par_exact", time_place_exact, "t_par_per_ex", time_place_perturbed_exact, "t_par_per", time_place_perturbed_approx) )
	print("{0:>32s}\t{1:<16s}\n{2:>32s}\t{3:<16.3g} \n {4:>32s}\t{5:<16.3g} \n {6:>32s}\t{7:<16.3g} \n {8:>32s}\t{9:<16.3g}".format("Case", "Total Time (s)", "Original", time_total_exact, "Perturbed (exact)", time_total_perturbed_exact, "Perturbed (approx)", time_total_perturbed_approx, "Perturbed (approx, !parse)", time_total_perturbed_approx - time_parse_perturbed_exact) )

