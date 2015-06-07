import re # regex
import sys # for grabbing input args
import numpy as np # array handling

delim = "\t"	# input file delimeter. Must be a single character

# Grab command line arguments
infile_name = sys.argv[1]

# Read input file
infile = open(infile_name,'r')
lines = infile.readlines()
infile.close()

# get number of nets and cells
n = lines[0].split(delim)
num_nets = int(n[0])
num_cells = int(n[1])

#print("num_nets " + str(num_nets) + "\tnum_cells " + str(num_cells) )

# [FIX] Use sparse matrix handling
Q = np.zeros((num_cells,num_cells)) # initialize laplacian matrix

## Process the HGR file
for net in range(1,len(lines)):	# iterate over all nets (ignore first, as that just gives us number of cells and nets)
	line_arr = lines[net].split(delim)	# Break up
	
	if(len(line_arr) > 1):
		net_cells = len(line_arr)
		net_weight = 1/(net_cells-1)

		for i in range(len(line_arr)):
			i_cell = int(line_arr[i])-1 # -1 since hgr file starts with index 1, and python starts with index 0
			Q[i_cell,i_cell] = Q[i_cell,i_cell] + net_weight # Add up connected net weights on diagonals (constructing D)

			for j in range(i,len(line_arr)):
				j_cell = int(line_arr[j])-1
				Q[i_cell,j_cell] = Q[i_cell,j_cell] - net_weight
				Q[j_cell,i_cell] = Q[j_cell,i_cell] - net_weight

print(str(Q))
