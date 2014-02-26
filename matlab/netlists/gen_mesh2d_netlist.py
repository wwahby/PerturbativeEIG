# Generate 2D mesh

import sys # for grabbing input args
delim = " "	# input file delimeter. Must be a single character


# Grab command line arguments
outfile_name = sys.argv[1]
width = int(sys.argv[2])
	

# get number of nets and cells
num_nets = 2*(width-1)*width
num_cells = width**2

outfile = open(outfile_name,'w')

header_str = str(num_nets) + delim + str(num_cells) + '\n'
outfile.write(header_str)


for row in range(width):
	print(str(row))
	longstr = ""
	for col in range(width):
		node = row*width + col + 1
		node_below = node + width
		node_right = node + 1
		
		
		if(col < width-1): # print horizontal net
			hnet_str = str(node) + delim + str(node_right) + '\n'
			longstr = longstr + hnet_str
			#outfile.write(hnet_str)

		if(row < width-1): # print vertical net
			vnet_str = str(node) + delim + str(node_below) + '\n'
			longstr = longstr + vnet_str
			#outfile.write(vnet_str)
	outfile.write(longstr)



