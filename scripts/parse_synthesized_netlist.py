# Parses a netlist synthesized with Cadence RTL Compiler 10.10
# parse_synthesized_netlist.py <INFILE> <OUTFILE>
# Outputs a stripped netlist file
#	First line has #nets #cells
# All subsequent lines represent nets, and are space-delimited
#	lists of cells connected to that net
# Before parsing the netlist, the following must be done
#	1. Remove spurious linebreaks
#	2. Remove random backslashes
#	3. Remove wire bus indices
#	See run_these_before_parsing.txt for vim regexes for these

import re
import sys # for grabbing input args

# Get arguments
infile_name = sys.argv[1]
outfile_name = sys.argv[2]

infile = open(infile_name,'r')
lines = infile.readlines()
infile.close()

#wire_reg = re.compile('\s+wire\s+(.+);')
#input_reg = re.compile('\s+input\s+(.+);')
#output_reg = re.compile('\s+output\s+(.+);')

wire_io_reg = re.compile('\s+(wire|output|input)\s+(.+);')
component_reg = re.compile('  \S+ \S+(\s+)?\(') # [FIX] MORE STUFF GOES HERE

wire_dict = dict()
connection_dict = dict()
component_dict = dict()
wire_ind_max = 0

for line_ind in range(len(lines)):

	wire_found = wire_io_reg.match(lines[line_ind])

	if (wire_found):
		wire_string = wire_found.group(2)
		wire_list = re.split(', ',wire_string)

		for wirename in wire_list:
			if wirename not in wire_dict:
				wire_dict.update({wirename : wire_ind_max})
				wire_ind_max = wire_ind_max + 1

	



