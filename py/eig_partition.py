import numpy as np # array handling
import eig_support as es
import argparse


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("hgr", help="Name of input HGR file", nargs='+') # '+' so we can reconstruct filenames with spaces
	parser.add_argument("--perturbed", "-p", help="Name of the perturbed HGR file", nargs='+')
	parser.add_argument("--num_eigs", "-n", type=int, help="Number of eigenvectors to use for perturbation approximation", default=10)
	parser.add_argument("--delim", "-d", help="Delimeter used in netlist files", default = " ")
	args = parser.parse_args()

	infile_name = ' '.join(args.hgr) # Reconstruct filenames with spaces
	if (args.perturbed):
		infile_p_name = ' '.join(args.perturbed) # Reconstruct filenames with spaces
		run_perturbed = True
	else:
		infile_p_name = "Invalid_name"
		run_perturbed = False

	num_eigs = args.num_eigs
	delimiter = args.delim	# input file delimeter. Must be a single character

	Q = es.parse_hgr_sparse(infile_name,delim=delimiter, index_offset=1)
	p1d = es.partition_1d(Q)

	#print("Exact partition\n==============\n" + str(p1d))
	#print()

	# Require parsed netlist
	#block_map = es.construct_component_map(infile_name + ".dict")
	#es.write_partitions(p1d, block_map, infile_name + ".part")


	if (run_perturbed):
		
		Qp_exact = es.parse_hgr(infile_p_name,delim)
		Qp = Qp_exact - Q
		p1dp = es.partition_1d_perturbed(Q,Qp,num_eigs)
		p1dpe = es.partition_1d(Qp_exact)


		print("Exact perturbed partition\n==============\n" + str(p1dpe))
		print()
		print("Approx perturbed partition\n==============\n" + str(np.array(p1dp)))
		print()


		(r1,r2) = es.compare_lists(list(p1dp),list(p1dpe),0.5)
		print("{0:>6.4f}\t{1:6.4f}".format(r1,r2))


if ( __name__ == "__main__"):
	main()
