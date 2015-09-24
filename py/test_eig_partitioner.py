import os
import eig_support as es

cur_dir = os.getcwd()
base_dir = os.path.dirname( cur_dir )
netlist_dir = os.path.join(base_dir, "netlists")

hgr_filename = "industry3.hgr"
perturbed_filename = "industry3_add05.hgr"

hgr = os.path.join( netlist_dir, hgr_filename )
perturbed = os.path.join( netlist_dir, perturbed_filename )

num_eigs = 20
delim =  " "
num_partitions = 2
eigenval_cutoff = 1e-3
index_offset = 1
unbalance_factor = 45

infile_name = hgr # Reconstruct filenames with spaces
infile_p_name = perturbed # Reconstruct filenames with spaces
run_perturbed = True

area_balance = 0.5 - unbalance_factor/100

(Q, D, A, time_parse_exact) = es.parse_hgr_sparse(infile_name,delim=delim, index_offset=index_offset)

(partition_order, eigvals, raw_eigvals, raw_eigvecs, time_partition_exact) = es.partition_1d(Q, eigenval_cutoff = eigenval_cutoff, num_eigs = num_eigs)


print("Finding cutsize of system...")
(mincut_val, mincut_ind, cutsize_vec, normcut_ind, normcut_val, normcut_vec, p1_size_frac_vec, time_cutsize_exact) = es.calc_cutsize_bipart(A, partition_order, area_balance)


if (run_perturbed):

	(Qp_exact, Dp_e, Ap_e, time_parse_perturbed) = es.parse_hgr_sparse(infile_p_name, delim=delim, index_offset=index_offset)
	Qp = Qp_exact - Q

	(p1dp, time_partition_perturbed, time_eig_standard) = es.partition_1d_perturbed(raw_eigvals, raw_eigvecs, Qp, eigenval_cutoff = eigenval_cutoff, num_eigs_solve = num_eigs)

	(p1dpe, vals_sorted_pe, vals_raw_pe, vecs_raw_pe, time_partition_perturbed_exact)  = es.partition_1d(Qp_exact, eigenval_cutoff = eigenval_cutoff, num_eigs = num_eigs)

	print("Finding cutsize of perturbed system...")
	(mincut_val_p, mincut_ind_p, cutsize_vec_p, normcut_ind_p, normcut_val_p, normcut_vec_p, p1_size_frac_vec_p, time_cutsize_perturbed) = es.calc_cutsize_bipart(Ap_e, p1dp, area_balance) # uses exact adjacency matrix, since that part can be known just based on connectivity, without actually solving eig problem

	print("Finding cutsize of exact solution to perturbed system...")
	(mincut_val_pe, mincut_ind_pe, cutsize_vec_pe, normcut_ind_pe, normcut_val_pe, normcut_vec_pe, p1_size_frac_vec_pe, time_cutsize_perturbed_exact) = es.calc_cutsize_bipart(Ap_e, p1dpe, area_balance)

	## Figures
	pl.figure(1)
	pl.hold(True)
	pl.plot(p1_size_frac_vec, cutsize_vec, 'k')
	pl.plot(p1_size_frac_vec_p, cutsize_vec_p, 'r')
	pl.plot(p1_size_frac_vec_pe, cutsize_vec_pe, 'b')

	pl.figure(2)
	pl.hold(True)
	pl.plot(p1_size_frac_vec, normcut_vec, 'k')
	pl.plot(p1_size_frac_vec_p, normcut_vec_p, 'r')
	pl.plot(p1_size_frac_vec_pe, normcut_vec_pe, 'b')

	print("{0:>16s}\t{1:<16.3g} \n {2:>16s}\t{3:<16.3g} \n {4:>16s}\t{5:<16.3g} \n ".format("t_par_exact", time_partition_exact, "t_par_per_ex", time_partition_perturbed_exact, "t_par_per", time_partition_perturbed) )
#		print("{0:^16.3g}\t{1:^16.3g}\t{2:^16.3g}".format(time_partition_exact, time_partition_perturbed_exact, time_partition_perturbed) )
