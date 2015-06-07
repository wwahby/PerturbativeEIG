import re # regex
import sys # for grabbing input args
import numpy as np # array handling
import eig_support as es
import argparse


#def main():
#	parser = argparse.ArgumentParser()
#	parser.add_argument("infile", type=float, help="Name of input HGR file")
#	parser.add_argument("k", type=float, help="Rent constant for the design")
#	parser.add_argument("p", type=float, help="Rent exponent for the design")
#	parser.add_argument("tiers", type=int, help="Number of tiers in the 3D stack")
#	args = parser.parse_args()

delim = " "	# input file delimeter. Must be a single character

# Grab command line arguments
#infile_name = sys.argv[1]
infile_name = "p2.hgr"
infile_p_name = "p2_p.hgr"

Q = es.parse_hgr(infile_name,delim)

#np.savetxt("testout.txt",Q,fmt='%.4f',delimiter='\t',newline='\n')

p1d = es.partition_1d(Q)

print("Exact partition\n==============\n" + str(p1d))
print()

#if(len(sys.argv) > 2):
	#infile_p_name = sys.argv[2]


n = 10
Qp_exact = es.parse_hgr(infile_p_name,delim)
Qp = Qp_exact - Q
p1dp = es.partition_1d_perturbed(Q,Qp,n)
p1dpe = es.partition_1d(Qp_exact)


print("Exact perturbed partition\n==============\n" + str(p1dpe))
print()
print("Approx perturbed partition\n==============\n" + str(np.array(p1dp)))
print()


(r1,r2) = es.compare_lists(list(p1dp),list(p1dpe),0.5)
print("{0:>6.4f}\t{1:6.4f}".format(r1,r2))
