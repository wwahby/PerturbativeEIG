# Generate 2D mesh

import sys # for grabbing input args
delim = " "	# input file delimeter. Must be a single character


# Grab command line arguments
outfile_name = sys.argv[1]
width = int(sys.argv[2])
	

# get number of nets and cells
num_nets = 3*(width-1)*(width**2)
num_cells = width**3

outfile = open(outfile_name,'w')

header_str = str(num_nets) + delim + str(num_cells) + '\n'
outfile.write(header_str)


for plane in range(width):
	print(str(plane))
	longstr = ""
	for row in range(width):
		for col in range(width):
			node = plane*(width**2) + row*width + col + 1
			node_below = node + width
			node_right = node + 1
			node_under = node + width**2
			
			if(col < width-1): # print horizontal net
				hnet_str = str(node) + delim + str(node_right) + '\n'
				longstr = longstr + hnet_str

			if(row < width-1): # print vertical net
				vnet_str = str(node) + delim + str(node_below) + '\n'
				longstr = longstr + vnet_str

			if(plane < width-1): # print below net
				bnet_str = str(node) + delim + str(node_under) + '\n'
				longstr = longstr + bnet_str

	outfile.write(longstr)
