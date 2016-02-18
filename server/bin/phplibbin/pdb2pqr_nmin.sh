#!/bin/bash

#defaults
prm=~/prms/trunk/parm.prm
rtf=~/prms/trunk/pdbamino.rtf
pdb2pqr=pdb2pqr.py
clean=1

args=$(getopt -o "" -l rtf:,prm:,ph:,smod:,pdb2pqr:,dont-clean,xplor-psf -- "$@")
eval set -- "$args"
while [ ! -z "$1" ]
do
    case "$1" in
        --rtf) rtf=$2;;
        --prm) prm=$2;;
        --ph) ph=$2;;
        --smod) smod=$2;;
        --pdb2pqr) pdb2pqr=$2;;
        --dont-clean) clean=0;;
        --xplor-psf) xplor_psf=1;;
        --) shift; break;;
    esac
    shift
done

if [ $# -lt 2 ]; then
    echo "usage: $0 [options] PDB ch1 ch2 ..." >&2
    exit 1
fi

nmd_opts="--rtf=$rtf --prm=$prm"
if [ -n "$smod" ]; then
    nmd_opts="$nmd_opts --smod=$smod"
fi
if [ "$clean" -eq 0 ]; then
    nmd_opts="$nmd_opts --dont-clean"
fi
pdbnmd.pl $nmd_opts $@
if [ $? -ne 0 ]; then 
    exit $?
fi

base=`basename $1 .pdb`
pdb="$base.pdb"

waters_name="$base.waters.pdb"
nmin_name="${base}_nmin.pdb"
nmin_psf_name="${base}_nmin.psf"
com_name="${base}_nmin.waters.pdb"
postpqr_pqr="$base.charmm.nmin.waters.nodebump.pqr"
postpqr_pdb="$base.charmm.nmin.waters.nodebump.pdb"
postpqr_mol_pdb="$base.charmm.nmin.waters.nodebump.mol.pdb"
postpqr_ngen="$base.charmm.nmin.waters.nodebump_ngen.pdb"
postpqr_nmin="$base.charmm.nmin.waters.nodebump_nmin.pdb"
postpqr_ngen_psf="$base.charmm.nmin.waters.nodebump_ngen.psf"
postpqr_ngen_fixed="$base.charmm.nmin.waters.nodebump_ngen.fixed.pdb"

grep "^HETATM" $pdb | awk '$4 == "HOH"' > $waters_name


cat $nmin_name $waters_name > $com_name

pdb.renumberAtoms.pl $com_name

#nodebump doesn't move the atoms unless hydrogen bond network moves them
#chain keeps the chain id
pdb2pqr_options="--nodebump --mol_charmm_pdb --chain --ff=CHARMM --ffout=LIBMOL"
if [ -n "$ph" ]; then
    pdb2pqr_options="$pdb2pqr_options --with-ph=$ph"
fi
$pdb2pqr $pdb2pqr_options $com_name $postpqr_pqr
if [ $? -ne 0 ]; then 
    exit $?
fi

mv $postpqr_mol_pdb $postpqr_pdb

pdbprep.pl $postpqr_pdb
if [ $? -ne 0 ]; then 
    exit $?
fi

#use psfgen to order the atoms the normal way, as pdb2pqr can reorder the atoms
#also removes all the nonpolar hydrogens added and gets rid of our waters
shift
chains=$@
nmd_opts="$nmd_opts --dont-minimize"
if [ -n "$xplor_psf" ]; then
    nmd_opts="$nmd_opts --xplor-psf"
fi
pdbnmd.pl $nmd_opts $postpqr_pdb $chains
if [ $? -ne 0 ]; then 
    exit $?
fi

tmp1=`mktemp`
tmp2=`mktemp`

#choose atom lines, take characters 22-55
#this is from the chain id to the end of x,y,z coordinates
#thus, we ignore renamed residues and the radius/charge columns added by pdb2pqr
grep '^ATOM' $postpqr_ngen | cut -c 22-55 > $tmp1
grep '^ATOM' $nmin_name    | cut -c 22-55 > $tmp2

#we are comparing our previous minimized to our new generated pdb post pdb2pqr
#to see added atoms, replaced hydrogens, flipped residues, etc.  We choose the lines
#that did not change to be fixed.  We need to add back in the atom numbers and so on for the fixed file
sdiff -l $tmp1 $tmp2 | nl | grep "(" | awk '{ printf "ATOM  %4d\n", $1}' > $postpqr_ngen_fixed

rm $tmp1
rm $tmp2

nmin $postpqr_ngen_psf $prm $rtf $postpqr_ngen $postpqr_ngen_fixed 1000
if [ $? -ne 0 ]; then 
    exit $?
fi

mv $postpqr_nmin $nmin_name
mv $postpqr_ngen_psf $nmin_psf_name
if [ -n "$xplor_psf" ]; then
    postpqr_ngen_xplor_psf="$base.charmm.nmin.waters.nodebump_ngen_xplor.psf"
    nmin_xplor_psf_name="${base}_nmin_xplor.psf"
    mv $postpqr_ngen_xplor_psf $nmin_xplor_psf_name
fi
if [ "$clean" -ne 0 ]; then
    rm $base.charmm.*
fi
