#!/usr/bin/python
# Scott Mottarella
# 10/06/11

#Takes in a pdb file and a number to assign the output
#Then preps the pdb to become a probe for mapping with assigned number
# Checks atom names to be unique
# Adds all probes to be chain X
# Makes residue name to be P##
# renames output file to p##.pdb

import sys
from optparse import OptionParser
import re
import os.path
import string

USAGE = "Usage: probe_prep.py [options] pdb probe_num"
parser = OptionParser(usage=USAGE)
parser.add_option("-l", "--letter", action="store", type="string", dest="letter", default="p")
parser.add_option("-m", "--map", action="store_true", dest="mapnames")
(options, args) = parser.parse_args()

if len(args)!=2:
	print USAGE
	sys.exit(1)

if int(args[1])>100:
	print "This is too many probes to be run at once"
	sys.exit(2)

f = open(args[0])
pdb = f.readlines()
f.close()
OUT = open(options.letter+str(args[1]).rjust(2,"0")+".pdb", "w")
if options.mapnames:
	MAPFILE = open(options.letter+str(args[1]).rjust(2,"0")+".mapnames.csv", "w")
names = []
resseqid = -1
for i in range(len(pdb)):
	line = pdb[i].strip()
	if line[:6]=="ATOM  " or line[:6]=="HETATM":
		#Get unique residue seuence id
		if resseqid==-1:
			resseqid = line[22:26]
		name = line[12:16].strip()
		#Check name
		if re.search("^[a-zA-Z]{1,2}$", line[76:78].strip()):
			tmp = 1
			newname = line[76:78].strip()+str(tmp)
			while newname in names:
				tmp += 1
				newname = (line[76:78].strip())+str(tmp)
			names.append(newname)
		elif re.search("^[a-zA-Z]{1,2}$", name):
			tmp = 1
			newname = name+str(tmp)
			while newname in names:
				tmp += 1
				newname = name+str(tmp)
			names.append(newname)
		else:
			tmp = 1
			s = re.search("^[a-zA-Z]{1,2}", name)
			if s:
				newname = s.group()+str(tmp)
				while newname in names:
					tmp += 1
					newname = s.group()+str(tmp)
			else:
				newname = name[0]+str(tmp)
				while newname in names:
					tmp += 1
					newname = name[0]+str(tmp)
			names.append(newname)
		OUT.write("ATOM  "+line[6:11]+" "+newname.rjust(4," ").upper()+line[16:17]+options.letter.upper()+str(args[1]).rjust(2,"0")+" X"+str(resseqid)+line[26:72]+"X"+line[73:]+"\n")
		if options.mapnames:
			MAPFILE.write(newname.rjust(4," ")+"\t"+line[6:27]+"\n")
	elif line[:6]=="TER   ":
		OUT.write("TER   "+line[6:17]+options.letter.upper()+str(args[1]).rjust(2, "0")+line[22:])
	else:
		OUT.write(line+"\n")
OUT.close()
if options.mapnames:
	MAPFILE.close()
