#!/usr/bin/python

#turn atomic contacts to residue contacts

f=open('atomic_contacts.txt')
rec_total=0.0
res_count={}
resns={}
for l in f:
    resi,chain,resn,atom,count=l.strip().split('\t')
    if chain not in res_count:
        res_count[chain]={}
        resns[chain]={}
    if resi not in res_count[chain]:
        res_count[chain][resi]=0
        resns[chain][resi]=resn
    res_count[chain][resi]+=int(count)
    rec_total+=int(count)
f.close()

fw=open('residue_contacts.txt','w')

for chain in sorted(res_count.keys()):
    for resi in sorted(res_count[chain].keys(),key=int):
        #print >>fw, '%s\t%s\t%s\t%d\t%.3f' % (resi, chain, resns[chain][resi], res_count[chain][resi], res_count[chain][resi]/rec_total*100)
        print >>fw, '%s\t%s\t%s\t%d' % (resi, chain, resns[chain][resi], res_count[chain][resi])
fw.close()
