#!/usr/bin/python
# Scott Mottarella
# 08/24/11

#Takes in one sdf file and a number
#Returns only the molecule at the number index
#returns nothing if index is too large

import sys
from rdkit import Chem

sdf = Chem.SDMolSupplier(sys.argv[1])
idx = int(sys.argv[2])

if idx>0 and idx<=len(sdf):
	print Chem.MolToMolBlock(sdf[idx-1])+"$$$$"
