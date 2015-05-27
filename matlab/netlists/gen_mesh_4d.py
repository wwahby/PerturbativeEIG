import numpy as np

def get_pt_ind(x, y, z, w, dim_size):
	pt_ind = x + y*dim_size + z*dim_size**2 + w*dim_size**3 + 1
	return pt_ind


def get_pt_ind_vec(vec, dim_size):
	return get_pt_ind(vec[0], vec[1], vec[1], vec[2], dim_size)


def gen_mesh_4d(dim_size, outfile_name):
	# Generates CLOSED mesh. edges leaving outer nodes are truncated
	outfile = open(outfile_name, 'w')

	num_nodes = dim_size**4
	num_edges = 4*(dim_size-1)*(dim_size)**3
	outfile.write(str(num_edges) + " " + str(num_nodes) + "\n")

	for w in range(dim_size):
		for z in range(dim_size):
			for y in range(dim_size):
				for x in range(dim_size):
					cur_pt_ind = x + y*dim_size + z*dim_size**2 + w*dim_size**3 + 1

					if (w < dim_size - 1):
						new_pt_ind = x + y*dim_size + z*dim_size**2 + (w+1)*dim_size**3 + 1
						outfile.write(str(cur_pt_ind) + " " + str(new_pt_ind) + "\n")

					if (z < dim_size - 1):
						new_pt_ind = x + y*dim_size + (z+1)*dim_size**2 + w*dim_size**3 + 1
						outfile.write(str(cur_pt_ind) + " " + str(new_pt_ind) + "\n")

					if (y < dim_size - 1):
						new_pt_ind = x + (y+1)*dim_size + z*dim_size**2 + w*dim_size**3 + 1
						outfile.write(str(cur_pt_ind) + " " + str(new_pt_ind) + "\n")

					if (x < dim_size - 1):
						new_pt_ind = (x+1) + y*dim_size + z*dim_size**2 + w*dim_size**3 + 1
						outfile.write(str(cur_pt_ind) + " " + str(new_pt_ind) + "\n")


def gen_mesh_4d_vec(dim_size, outfile_name):
	outfile = open(outfile_name, 'w')

	for w in range(dim_size):
		for z in range(dim_size):
			for y in range(dim_size):
				for x in range(dim_size):
					cur_pt = np.array([x, y, z, w])

					if (w < dim_size - 1):
						new_pt = np.array([x, y, z, w+1])
						outfile.write(str(cur_pt) + " " + str(new_pt) + "\n")

					if (z < dim_size - 1):
						new_pt = np.array([x, y, z+1, w])
						outfile.write(str(cur_pt) + " " + str(new_pt) + "\n")

					if (y < dim_size - 1):
						new_pt = np.array([x, y+1, z, w])
						outfile.write(str(cur_pt) + " " + str(new_pt) + "\n")

					if (x < dim_size - 1):
						new_pt = np.array([x+1, y, z, w])
						outfile.write(str(cur_pt) + " " + str(new_pt) + "\n")


def gen_mesh_4d_open(dim_size, outfile_name):
	# Generates OPEN mesh. Edges leaving outer nodes wrap around to other side
	outfile = open(outfile_name, 'w')

	num_nodes = dim_size**4
	num_edges = 4*(dim_size)*(dim_size)**3
	outfile.write(str(num_edges) + " " + str(num_nodes) + "\n")

	for w in range(dim_size):
		for z in range(dim_size):
			for y in range(dim_size):
				for x in range(dim_size):
					cur_pt_ind = x + y*dim_size + z*dim_size**2 + w*dim_size**3 + 1

					if (w < dim_size - 1):
						new_pt_ind = x + y*dim_size + z*dim_size**2 + (w+1)*dim_size**3 + 1
					else:
						new_pt_ind = x + y*dim_size + z*dim_size**2 + 0*dim_size**3 + 1
					outfile.write(str(cur_pt_ind) + " " + str(new_pt_ind) + "\n")

					if (z < dim_size - 1):
						new_pt_ind = x + y*dim_size + (z+1)*dim_size**2 + w*dim_size**3 + 1
					else:
						new_pt_ind = x + y*dim_size + 0*dim_size**2 + w*dim_size**3 + 1
					outfile.write(str(cur_pt_ind) + " " + str(new_pt_ind) + "\n")

					if (y < dim_size - 1):
						new_pt_ind = x + (y+1)*dim_size + z*dim_size**2 + w*dim_size**3 + 1
					else:
						new_pt_ind = x + 0*dim_size + z*dim_size**2 + w*dim_size**3 + 1
					outfile.write(str(cur_pt_ind) + " " + str(new_pt_ind) + "\n")

					if (x < dim_size - 1):
						new_pt_ind = (x+1) + y*dim_size + z*dim_size**2 + w*dim_size**3 + 1
					new_pt_ind = 0 + y*dim_size + z*dim_size**2 + w*dim_size**3 + 1
					outfile.write(str(cur_pt_ind) + " " + str(new_pt_ind) + "\n")
