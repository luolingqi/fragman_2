#!/bin/bash
# usage: create_pymol_session.sh fftmap.pse
# needs to be run in ftmap job directory

echo "cd probes/cluster
load ../rec/1rec.pdb, protein
as cartoon, protein
" >tmp.pml
cd probes/cluster

for cc in crosscluster*pdb;do
	echo "load $cc"
done >>../../tmp.pml

#for extended probes
if [ -f ../probes-ext ];then
	moli=1
	for mol in `cut -c-2 ../probes-ext|sort|uniq`;do
		confi=1
		for conf in `ls list_$mol??\_?_uniq|cut -c8-9`;do
			for clus in `cat list_${mol}${conf}_?_uniq`;do
				echo "load $clus.pdb"
				#echo "disable $clus"
				echo "group mol$moli-conf$confi, $clus"
			done
			let confi=$confi+1
		done
		echo "group molecule_$moli, mol$moli-conf*"
		echo "as sticks, molecule_$moli"
		let moli=$moli+1
	done
fi >>../../tmp.pml

cd ../..

echo "as sticks, crosscluster*">>tmp.pml
echo "orient" >>tmp.pml
echo "cd ../.." >>tmp.pml
echo "save $1" >>tmp.pml
echo "bg_color white" >>tmp.pml
echo "png result.png, width=750, height=600, dpi=72, ray=1" >>tmp.pml

pymol -qc tmp.pml
rm tmp.pml
