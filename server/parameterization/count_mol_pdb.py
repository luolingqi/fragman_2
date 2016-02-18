#!/usr/bin/python
# Scott Mottarella
# 10/06/11

#Counts the number of models in a pdb file

import sys
from optparse import OptionParser

USAGE = "Usage: count_mol_pdb.py pdb"
parser = OptionParser(usage=USAGE)
(options, args) = parser.parse_args()

if len(args)!=1:
	print USAGE
	sys.exit(1)

pdb = open(args[0])
cend = 0
cendmdl = 0
for line in pdb:
	if line.strip()=="END":
		cend += 1
	elif line.strip()=="ENDMDL":
		cendmdl += 1

if cend==cendmdl:
	print cend
elif cend>cendmdl:
	print cend
	#print "Number of END("+str(cend)+") greater than ENDMDL("+str(cendmdl)+")"
else:
	print cendmdl
	#print "Number of ENDMDL("+str(cendmdl)+") greater than END("+str(cend)+")"
