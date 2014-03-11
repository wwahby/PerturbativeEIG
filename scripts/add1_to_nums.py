# Read in space-delimited list of numbers and add one to each of them
# Useful for converting numerical netlists to matlab format (can't have 0s)

import sys # for input handling

infile_name = sys.argv[1]
outfile_name = sys.argv[2]

infile = open(infile_name,'r')
lines = infile.readlines()
infile.close()

outfile = open(outfile_name,'w')

for lind in range(1,len(lines)): # skip first line -- just header info
	line = lines[lind]
	strlist = line.split(' ')
	newlist = []
	
	for eind in range(len(strlist)):
		newlist.append( int(strlist[eind]) + 1 )
	
	newline = ' '.join(map(str,newlist))
	lines[lind] = newline

for line in lines:
	outfile.write(line)
	outfile.write('\n')

outfile.close()	
