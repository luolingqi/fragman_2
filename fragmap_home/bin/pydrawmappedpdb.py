from sys import argv
from pymol import cmd, util

if len(argv) < 2:
	print "usage: " + argv[0] + "pdbfile\n"
	quit()

pdbfile = argv[1]


cmd.set('antialias', 1)
cmd.set('depth_cue', 0)
cmd.set('ray_opaque_background', 0)
cmd.set('ray_trace_fog', 0)

cmd.load(pdbfile)

cmd.as('sticks', 'cross*')
cmd.as('surface', 'protein')
cmd.color('wheat', 'protein')

cmd.orient('protein')
cmd.png('mappedpdb', width=600, height=450, dpi=72, ray=1)
