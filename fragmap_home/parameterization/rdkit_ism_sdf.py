#!/usr/bin/python
from sys 	import argv
from sys	import exit
from rdkit 	import Chem
from rdkit.Chem	import AllChem

if len(argv) < 3:							# check no of arguments else dies
	exit()

filename = argv[1]
molename = argv[2]

f_handle = open  (filename)
ism_strs = f_handle.read ()
ism_list = ism_strs.split()

for i in ism_list:

	m = Chem.MolFromSmiles(i)					# Constructs a molecule from a SMILES string, OBJECT is NOT indexable
	ism_mdls = Chem.AddHs(m)
	fmcharge = AllChem.GetFormalCharge(ism_mdls)			# Returns the formal charge for the molecule
	print fmcharge

	AllChem.EmbedMultipleConfs(ism_mdls, 1)				# Dist geometry to obtain multiple sets of coordinates for a molecule

	for i in range(ism_mdls.GetNumConformers()):
		AllChem.UFFOptimizeMolecule(ism_mdls, confId=i)		# Force-field optimization of models within ism_mdls

	pdb_mdls = Chem.SDWriter(molename + '.sdf')			# Writes to file: nam_list[j] + '.sdf'

	for i in range(ism_mdls.GetNumConformers()):
		pdb_mdls.write(ism_mdls, confId=i)			# Actual writing
