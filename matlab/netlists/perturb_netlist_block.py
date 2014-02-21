# perturb_hgr <infile> <outfile> <perturbation fraction (between 0 and 1>
# Quick and dirty script to generate random additive perturbations
# to arbitrary netlists.

import re # regex
import sys # for grabbing input args
#import numpy as np # array handling
import random

delim = " "	# input file delimeter. Must be a single character

# initialize PRNG with system time as the seed
random.seed()

# Grab command line arguments
infile_name = sys.argv[1]
outfile_name = sys.argv[2]
#perturbation_chance = float(sys.argv[3])
num_to_perturb = float(sys.argv[3])

# Read input file
infile = open(infile_name,'r')
lines = infile.readlines()
infile.close()

# Remove extra whitespace
for i in range(len(lines)):
	lines[i] = lines[i].rstrip()
	

# get number of nets and cells
n = lines[0].split(delim)
num_nets = int(n[0])
num_cells = int(n[1])

outfile = open(outfile_name,'w')
num_perturbed = 0

start_net = random.randint(1,num_cells-num_to_perturb-1)

for i in range(len(lines)):
	elements = lines[i].split(delim)
	num_nodes_this_net = len(elements)

	# for now we will exclude nets with only two nodes, as those would have only one node after deletion
	# [FIX] Update the PEIG code to work with single-node nets

	if ( (i > start_net+1) and (num_perturbed < num_to_perturb) ):
		add_or_del = random.randint(0,1)

		if(add_or_del == 0): # add some random nodes
			num_to_add = random.randint(1,3)
			new_net = elements

			for j in range(0,num_to_add):
				node_to_add = random.randint(1,num_cells-1)
				lines[i] = lines[i] + " " + str(node_to_add)

			#lines[i] = new_net
			num_perturbed = num_perturbed + 1

		elif ((num_nodes_this_net > 2) and (add_or_del == 1) ): # Don't want to mess up the first line! It's a header
			dice = 0 # always delete if we get here
			perturbation_chance = 1
			add_or_del = random.randint(0,1)

			if (dice <= perturbation_chance) :
				del_node = random.randint(0,num_nodes_this_net-1)
				
				new_net = ""
				for j in range(del_node):
					new_net = new_net + elements[j] + " "

				for j in range(del_node+1,num_nodes_this_net):
					new_net = new_net + elements[j] + " "

				new_net = new_net.rstrip() # we'll have an extra space at the end
				lines[i] = new_net
				num_perturbed = num_perturbed + 1

	outfile.write(lines[i])
	outfile.write('\n')
				


