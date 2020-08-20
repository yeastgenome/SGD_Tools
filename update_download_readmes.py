#!/usr/bin/python
import sys, getopt, os
import re
import subprocess

## Script to copy/update READMEs from pinot to bun (download server) that are listed in a specific file:
##

def main(argv):

	try:
		opts, args = getopt.getopt(argv,"hi:d:",["ifile=", "dir="])
	except getopt.GetoptError:
		print 'update_download_readmes.py -i <inputfile> '
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-i", "--ifile"):
			inputfile = arg
		if opt in ("-d", "--dir"):
			readme_dir = arg
	
	file_obj = open(inputfile, 'r')
#	file_data = open(file).readlines()
	col_count = 0
	
	for directory in file_obj:
		# check to see if directory exists #
		dir = directory.rstrip()

		print "LINE: " + dir
					
		if (os.path.exists(readme_dir + dir)):
			file_to_copy = readme_dir + dir + "/README"
			copy_to_dir = "edith@bun:/share/bun/www-data/html/expression/microarray/" + dir
			
			print "file to copy:" + file_to_copy + "|"
			print "copy to: " + copy_to_dir + "|"
			
			subprocess.call("scp -q "+ file_to_copy + " " + copy_to_dir, shell=True)

if __name__ == "__main__":
   main(sys.argv[1:])