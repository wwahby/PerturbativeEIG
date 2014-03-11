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
delim = ' ' # delimiter for outfile

infile = open(infile_name,'r')
lines = infile.readlines()
infile.close()

#wire_reg = re.compile('\s+wire\s+(.+);')
#input_reg = re.compile('\s+input\s+(.+);')
#output_reg = re.compile('\s+output\s+(.+);')

wire_io_reg = re.compile('\s+(wire|output|input)\s+(.+);')
wire_io_bus_reg = re.compile('\s+(wire|output|input)(\s*)\[(\d+):(\d+)\](\s*)(.+);') #  wire [31:0] s1_s9_2;
component_reg = re.compile('  \S+ (\S+)(\s*)\((\..+)\);')
connected_reg = re.compile('.+\((\S+)\)')

# Dictionaries to keep track of the wire names and their corresponding indices
wire_dict = dict()
connection_dict = dict()
component_dict = dict()

# index of the last element
wire_ind_max = -1
component_ind_max = -1

# this list will hold all the nodes that are connected to each wire
wire_connection_list = []

for line_ind in range(len(lines)):

	wire_found = wire_io_reg.match(lines[line_ind])
	bus_found = wire_io_bus_reg.match(lines[line_ind])
	component_found = component_reg.match(lines[line_ind])

	
	# if we have a multi-wire bus, need to do some additional work
	if (bus_found):
		wire_string = bus_found.group(6)
		wire_list = re.split(', ',wire_string)
		bus_max_ind = int(bus_found.group(3))
		bus_min_ind = int(bus_found.group(4))

		for wirename in wire_list:
			for bus_ind in range(bus_min_ind,bus_max_ind+1):
				busname = wirename + '[' + str(bus_ind) + ']'
				
				if busname not in wire_dict:
					wire_ind_max = wire_ind_max + 1
					wire_dict.update({busname : wire_ind_max}) # add new entry to the dictionary
					wire_connection_list.append([]) # extend the wire connection list

	elif (wire_found):
		wire_string = wire_found.group(2)
		wire_list = re.split(', ',wire_string)

		for wirename in wire_list:
			if wirename not in wire_dict:
				wire_ind_max = wire_ind_max + 1
				wire_dict.update({wirename : wire_ind_max}) # add new entry to the dictionary
				wire_connection_list.append([]) # extend the wire connection list
	
	elif (component_found):
		component_name = component_found.group(1)
		connection_string = component_found.group(3)
		connection_list_a = re.split(', ',connection_string)

		# First, check to see if we've already seen this component before
		# if not, add it to the dictionary
		if component_name not in component_dict:
			component_ind_max = component_ind_max + 1
			component_dict.update({component_name : component_ind_max})
			

		for terminal in connection_list_a:
			terminal_found = connected_reg.match(terminal)
			if terminal_found:
				wire_connected = terminal_found.group(1)
				component_ind = component_dict[component_name] # grab the component index from the dictionary
				wire_ind = wire_dict[wire_connected] # we're assuming that all of the wires were listed in the header, so if we've gotten this far we don't need to worry about adding new wires to the dictionary
				wire_connection_list[wire_ind].append(component_ind)



num_nets = len(wire_connection_list)
num_cells = len(component_dict)

# write new file
outfile = open(outfile_name,'w')
outfile.write(str(num_nets) + delim + str(num_cells) + '\n') # print header info
for net_ind in range(num_nets):
	if (len(wire_connection_list[net_ind]) > 1):
		line_out = delim.join(map(str,wire_connection_list[net_ind]))
		outfile.write(line_out + '\n')

