#!/usr/bin/python
# Scott Mottarella
# 10/06/11

#Script takes one pdb file and a number.  Prints single pdb at index of number from given pdb.

import sys
from optparse import OptionParser

USAGE = "Usage: extract_pdb.py pdb number\n"
parser = OptionParser(usage=USAGE)
(options, args) = parser.parse_args()

f = open(args[0])
ex = []
count = 1
for line in f:
	if line.strip()=="ENDMDL" or line.strip()=="END":
		ex.append(line)
		if str(count)==args[1]:
			for l in ex:
				print l.strip()
			break
		else:
			ex = []
			count += 1
			continue
	else:
		ex.append(line)
