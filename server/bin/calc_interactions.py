import os
import sys
import shutil

if len(sys.argv) < 3:
    sys.exit();
else:
    my_base_dir = sys.argv[1];
    my_complex_file = sys.argv[2];

os.environ['HOME']=my_base_dir

print 'Processing:', my_complex_file

f = open(my_complex_file, 'r') 

pdb = my_complex_file.split('.pdb')[0]

try:
    # Get rid of multiple NMR models and the alternate location in both ATOM and CONNECT records
    atompdb = '%s.atoms' %pdb
    f2 = file(atompdb, 'w')
    header = ['ATOM  ', 'TER   ', 'HETATM']
    bad_numbers = []
    for line in f:
            if line.startswith('MODEL  ') and line[10:14].strip() == '2':
                    break
            elif line[:6] in header:
                    if line[16] not in ['A', ' ']:
                            bad_numbers.append(int(line[6:11].strip()))
                    else:
                            f2.write(line)
            elif line[:6] == 'CONECT' and line[6:11].strip() != '':
                    if int(line[6:11].strip()) not in bad_numbers:
                            f2.write(line)
    f2.flush()
    f2.close()

    # Clean up the pdb file using clean.f
    cleaninp = file('%s.inp' %pdb, 'w')
    print >> cleaninp, atompdb
    cleaninp.close()
    print '    Running Clean'
    os.system(my_base_dir + os.sep + 'hbplus' + os.sep + 'clean < %s.inp > %s.clnlog' %(pdb, pdb))
    os.remove('%s.inp' %pdb)

    #-----------#
    # H B A D D #
    #-----------#
    print '    Running HBAdd...'
    os.system(my_base_dir + os.sep + 'ligplot' + os.sep + 'hbadd ' + atompdb + ' ' + my_base_dir + os.sep + 'het_dictionary_mod.txt > hbadd.log')

    if os.path.isfile('hbplus.rc'):
            rc = '-f hbplus.rc'
    else:
            rc = ''
            error('hbadd', 'blah', pdb)
    if os.path.isfile('%s.new' %pdb):
            new = '%s.new' %pdb
    else:
            new = ''
            error('clean', 'blah',pdb)

    #-------------#
    # H B P L U S #
    #-------------#
    print '    Running HBPlus'
    shutil.copyfile(my_base_dir + os.sep + 'ligplot' + os.sep + 'ligplot.prm', 'ligplot.prm')

    #Calculate Non-Bonded
    os.system(my_base_dir + os.sep + 'hbplus' + os.sep + 'hbplus -L %s -h 2.9 -d 5 -N %s %s > nnb.log' %(rc, new, atompdb))
    #Clean out self interactions 
    interactions = open(pdb + '.nnb', 'r')
    output = open(pdb + '.nnb.out', 'w')

    for line in interactions:
            if line[10] == line[30] and ( line[49] != 'H' and line[50] != 'H' ):
                    continue
            elif line[16:19] == 'HOH' and line[36:39] == 'HOH':
                    continue
            else:
                    output.write(line)

    interactions.close()
    os.remove(pdb + '.nnb')
    output.close()

    #Calculate H Bonds	    
    os.system(my_base_dir + os.sep + 'hbplus' + os.sep + 'hbplus -L %s -h 2.7 -d 3.35 %s %s > hhb.log' %(rc, new, atompdb))
    interactions = open(pdb + '.hhb', 'r')
    output = open(pdb + '.hhb.out', 'w')

    for line in interactions:
            if line[10] == line[30] and ( line[49] != 'H' and line[50] != 'H' ):
                    continue
            elif line[16:19] == 'HOH' and line[36:39] == 'HOH':
                    continue
            else:
                    output.write(line)

    interactions.close()
    os.remove(pdb + '.hhb');
    output.close()

    if not new == '':
            os.remove(new)
    os.remove(atompdb)

except: 
    pass
f.close()
print 'done succesfully!'
