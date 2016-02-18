#!/usr/bin/python

#Scott Mottarella
#Passed a SMILES string, creates 3D coordinates and prints mol block
#capable of taking any number of SMILES strings on the command line (space delimited)

import sys
from rdkit import Chem
from rdkit.Chem import AllChem

smi = sys.argv[1:]

for s in smi:
	mol = Chem.MolFromSmiles(s)
	if mol is None: continue
	#AllChem.Compute2DCoords(mol)
	AllChem.EmbedMolecule(mol)
	try:
		AllChem.UFFOptimizeMolecule(mol)
	except:
		continue
	print Chem.MolToMolBlock(mol)
	print "$$$$"
