import numpy as np # array handling
import eig_support as es
import argparse
import pylab as pl


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("hgr", help="Name of input HGR file", nargs='+') # '+' so we can reconstruct filenames with spaces
	parser.add_argument("--perturbed", "-p", help="Name of the perturbed HGR file", nargs='+')
	parser.add_argument("--num_eigs", "-n", type=int, help="Number of eigenvectors to use for perturbation approximation", default=10)
	parser.add_argument("--delim", "-d", help="Delimeter used in netlist files. Must be a single character", default = " ")
	parser.add_argument("--num_partitions", help="Number of partitions to create", type=int, default=2)
	parser.add_argument("--eigenval_cutoff", help="Minimum value of a valid eigenvalue. Used to sift out spurious eigs from sparse EIG solver", type=float, default=1e-5)
	parser.add_argument("--index_offset", help="Value of minimum index in netlist graph file (if nonzero).", type=int, default=0)
	parser.add_argument("--unbalance_factor", help="Maximum deviation from equally-split partition. Values range from 0 (50-50 split) to 50 (100-0 split). Defaults to 5 (45/55 split).", type=float, default=5.0 )
	args = parser.parse_args()

	infile_name = ' '.join(args.hgr) # Reconstruct filenames with spaces
	if (args.perturbed):
		infile_p_name = ' '.join(args.perturbed) # Reconstruct filenames with spaces
		run_perturbed = True
	else:
		infile_p_name = "Invalid_name"
		run_perturbed = False

	unbalance_factor = args.unbalance_factor
	area_balance = 0.5 - unbalance_factor/100
	
	(Q, D, A) = es.parse_hgr_sparse(infile_name,delim=args.delim, index_offset=args.index_offset)
	(partition_order, eigvals, raw_eigvals) = es.partition_1d(Q, eigenval_cutoff = args.eigenval_cutoff, num_eigs = args.num_eigs)
	
	print("Finding cutsize of system...")
	(mincut_val, mincut_ind, cutsize_vec, normcut_ind, normcut_val, normcut_vec, p1_size_frac_vec) = es.calc_cutsize_bipart(A, partition_order, area_balance)

	

	if (run_perturbed):
		
		(Qp_exact, Dp_e, Ap_e) = es.parse_hgr_sparse(infile_p_name, delim=args.delim, index_offset=args.index_offset)
		Qp = Qp_exact - Q

		p1dp = es.partition_1d_perturbed(Q, Qp, eigenval_cutoff = args.eigenval_cutoff, num_eigs_solve = args.num_eigs)
		p1dpe = es.partition_1d(Qp_exact)
		
		print("Finding cutsize of perturbed system...")
		(mincut_val_p, mincut_ind_p, cutsize_vec_p, normcut_ind_p, normcut_val_p, normcut_vec_p, p1_size_frac_vec_p) = es.calc_cutsize_bipart(Ap_e, p1dp, area_balance) # uses exact adjacency matrix, since that part can be known just based on connectivity, without actually solving eig problem
		print("Finding cutsize of exact solution to perturbed system...")
		(mincut_val_pe, mincut_ind_pe, cutsize_vec_pe, normcut_ind_pe, normcut_val_pe, normcut_vec_pe, p1_size_frac_vec_pe) = es.calc_cutsize_bipart(Ap_e, p1dpe, area_balance)

		pl.fig(1)
		pl.hold()
		pl.plot(p1_size_frac_vec, cutsize_vec, 'color', 'k')
		pl.plot(p1_size_frac_vec_p, cutsize_vec_p, 'color', 'r')
		pl.plot(p1_size_frac_vec_pe, cutsize_vec_pe, 'color', 'b')
		
		pl.fig(2)
		pl.hold()
		pl.plot(p1_size_frac_vec, normcut_vec, 'color', 'k')
		pl.plot(p1_size_frac_vec_p, normcut_vec_p, 'color', 'r')
		pl.plot(p1_size_frac_vec_pe, normcut_vec_pe, 'color', 'b')
		
#		print("Exact perturbed partition\n==============\n" + str(p1dpe))
#		print()
#		print("Approx perturbed partition\n==============\n" + str(np.array(p1dp)))
#		print()
#
#
#		(r1,r2) = es.compare_lists(list(p1dp),list(p1dpe),0.5)
#		print("{0:>6.4f}\t{1:6.4f}".format(r1,r2))


if ( __name__ == "__main__"):
	main()
