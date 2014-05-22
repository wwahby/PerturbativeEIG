import glob
import sys # for grabbing input args
import os # hacky way to run the other python scripts


scriptname = sys.argv[1] # the script to run
extension = sys.argv[2] # the files to run it on

file_list = glob.glob(extension)


for fname in file_list:
	exstring = scriptname + " " + fname
	os.system(exstring)

