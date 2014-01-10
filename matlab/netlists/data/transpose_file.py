# transposes a simple array
# assumes comma delimited
# arg1: filename
# arg2: delimeter (optional - if left out script assumes comma delimited)

import re
import sys
import itertools


infile_name = sys.argv[1]
delim = '\t' # default delimiter
outfile_name = infile_name + ".trans"

if (len(sys.argv) > 2):
	delim = sys.argv[2]

infile = open(infile_name,'r')
outfile = open(outfile_name,'w')

data = []
num_rows = 0;
num_cols = 0;

# read the file in
lines = infile.readlines()
infile.close()

# strip off any trailing whitespace and EOLs
for i in range(len(lines)):
	lines[i] = lines[i].rstrip()
	lines[i] = lines[i].split(delim)


# Hacky way to transpose lists in python (v3.0 and later only)
# Shouldn't screw up if lists aren't all the same length
trans = list(map(list,itertools.zip_longest(*lines)))

for line in trans:
	outfile.write( delim.join(line) + "\n")








