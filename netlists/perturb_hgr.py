import re # regex
import sys # for grabbing input args
#import numpy as np # array handling
import random

# Quick and dirty script to generate random additive perturbations
# to arbitrary netlists.
# USAGE
# perturb_hgr <infile> <outfile> <perturbation fraction (between 0 and 1>
delim = " "	# input file delimeter. Must be a single character

# initialize PRNG with system time as the seed
random.seed()

# Grab command line arguments
infile_name = sys.argv[1]
outfile_name = sys.argv[2]
perturbation_chance = float(sys.argv[3])

# Read input file
infile = open(infile_name,'r')
lines = infile.readlines()
infile.close()

for i in range(len(lines)):
	lines[i] = lines[i].rstrip()
	

# get number of nets and cells
n = lines[0].split(delim)
num_nets = int(n[0])
num_cells = int(n[1])

outfile = open(outfile_name,'w')

for i in range(len(lines)):
	dice = random.random()
	if (dice < perturbation_chance) and (i > 1):
		new_node = random.randint(1,num_cells+1)
		lines[i] = lines[i] + " " + str(new_node)

	outfile.write(lines[i])
	outfile.write('\n')
				


