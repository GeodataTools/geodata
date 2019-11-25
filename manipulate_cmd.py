# -*- coding: utf-8 -*-
##	manipulate_cmd.py

##	Manipulate variables in netCDF files -- Command line version

from __future__ import absolute_import

import argparse

# initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("--delete", "-d", type=str, help="variables to delete")
parser.add_argument("--wind", "-w", action="store_true", help="keep wind related variables")
parser.add_argument("--solar", "-s", action="store_true", help="keep solar related variables")
parser.add_argument("--output", "-o", type=str, help="output filename")
parser.add_argument("filename", type=str, help="input filename")

# read arguments from the command line
args = parser.parse_args()
if args.output == None:
	args.output = args.filename

# Printing the arguments:
if args.delete:
    print("Delete variables: %s" % args.delete)
print ("Wind: %s" % args.wind)
print ("Solar: %s" % args.solar)
print ("Input filename: %s" % args.filename)
print ("Output filename: %s" % args.output)

# Call manipulate functions
