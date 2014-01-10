# characterize_netlist <infile>
# Quick and dirty script to determine net size distribution
# Takes in a space-delimited netlist
#	First line: "<# nets> <# nodes>"
#	Subsequent lines: "<first node in net> <second node in net> <etc...>"
# Outputs a comma-delimited file listing the number of nets with each net size.
# Filename: <infile>.histo
#	Format: "<# nodes in net>,<# nets>"
# Outputs a comma-delimited file listing the number of nets with each net size.
# Filename: <infile>.norm
#	Format: "<# nodes>,<fraction of nets>"

import re # regex
import sys # for grabbing input args
import numpy as np # for arrays

delim = " "	# input file delimeter. Must be a single character

# Grab command line arguments
infile_name = sys.argv[1]
outfile_name_histo = infile_name + ".histo"
outfile_name_percent = infile_name + ".norm"

# Read input file
infile = open(infile_name,'r')
lines = infile.readlines()
infile.close()

for i in range(len(lines)):
	lines[i] = lines[i].rstrip()
	

# get number of nets and cells
n = lines[0].split(delim)
num_nets = int(n[0])
num_nodes = int(n[1])



histogram = np.zeros(num_nodes)	# Could have one net with all the nodes!

for i in range(1,len(lines)):
	exploded = lines[i].split(delim)
	num_nodes_this_net = len(exploded)
	histogram[num_nodes_this_net-1] = histogram[num_nodes_this_net-1] + 1 #-1 because index 0 corresponds to nets with one node


histogram = np.trim_zeros(histogram,'b') # remove all the trailing zero elements

# Output the raw histogram data
outfile = open(outfile_name_histo,'w')
for i in range(len(histogram)):
	num_nodes_this_net = i + 1
	outfile.write(str(num_nodes_this_net) + "," + str(histogram[i]) )
	outfile.write('\n')

outfile.close()

# Output the normalized histogram data
outfile = open(outfile_name_percent,'w')
cumulative_nets = 0
for i in range(len(histogram)):
	num_nodes_this_net = i + 1
	cumulative_nets = cumulative_nets + histogram[i]
	net_fraction = cumulative_nets/num_nets

	outfile.write( str(i)  + "," + str(net_fraction) )
	outfile.write('\n')

outfile.close()
				


