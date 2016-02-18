#!/usr/bin/perl -w

%amino=("G", "GLY", "A", "ALA","P", "PRO", "V", "VAL", "L", "LEU", "I", "ILE", "M", "MET", "F", "PHE", "Y", "TYR", "W", "TRP", "S", "SER", "T", "THR", "C", "CYS", "N", "ASN", "Q", "GLN", "K", "LYS", "H", "HIS", "R", "ARG", "D", "ASP", "E", "GLU", "!", "XXX");


if(@ARGV<3)
{
print "usage:$0 pdbfile #residue_to_keep outfile\n";
exit;
}

$inputfile=shift;
$remainnum=shift;
$outfile=shift;

system "grep ATOM $inputfile >temp.pdb";

$tempdsspfile="$inputfile.dssp";
system "dsspcmbi temp.pdb >$tempdsspfile";

open(OUTPUT, ">$outfile") or die;

$logfile="$tempdsspfile.log";
open(LOG, ">$logfile") or die;

#Determine Chains
open(INPUT, "<temp.pdb") or die;
@inputcontent=<INPUT>;
close(INPUT);

#$chaincount=0;
#$chaintemp="?";

#foreach(@inputcontent)
#{
#$data=$_;

#if($data=~/ATOM\s+\d+/)
#{
#$chainid=substr($data,21,1);
#if($chainid ne $chaintemp)
#{
#$chaintemp=$chainid;
#$chain[$chaincount]=$chainid;
#$chaincount++;

#}
#}
#}



#Find terminal region

open(DSSP, "<$tempdsspfile") or die;

while($s=<DSSP>)
{
if($s=~/#\s+RESIDUE\s+AA/)
{last;}
}

@dsspcontent=<DSSP>;
close(DSSP);
$numdsspcontent=@dsspcontent;

for($i=0;$i<$numdsspcontent;$i++)
{
if($dsspcontent[$i]=~/(\d+)\s+\d+/)
{$dssphead=$1;}

$struc=substr($dsspcontent[$i], 16,1);

if($struc eq "H" ||$struc eq "E" ||$struc eq "G"||$struc eq "I"||$struc eq "T")
{last;}


}

for($i=$numdsspcontent-1;$i>-1;$i--)
{
if($dsspcontent[$i]=~/(\d+)\s+\d+/)
{$dssptail=$1;}

$struc=substr($dsspcontent[$i], 16,1);
if($struc eq "H" ||$struc eq "E" ||$struc eq "G"||$struc eq "I"||$struc eq "T")
{last;}
}

$dsspheadfinal=$dssphead-$remainnum;
$dssptailfinal=$dssptail+$remainnum;

foreach(@inputcontent)
{
$inputline=$_;
$inputresname=substr($inputline,17, 3);
$inputchain=substr($inputline,21,1);
$inputresi=substr($inputline, 22, 5);
$inputresi=~s/\s+//g;
$lock=0;	

	foreach(@dsspcontent)
	{
	$dsspline=$_;
	$dsspid=substr($dsspline, 0, 5);
	$dsspid=~s/\s+//g;
	$dsspresid=substr($dsspline, 5, 5);
	$dsspresid=~s/\s+//g;
	$dsspchain=substr($dsspline, 11, 1);
	$dsspresname=substr($dsspline, 13, 1);
	if($inputresname eq $amino{$dsspresname} && $inputchain eq $dsspchain && $inputresi eq $dsspresid)
	{
	$lock=1;
		if($dsspid>$dsspheadfinal-1 && $dsspid<$dssptailfinal+1)		
		{ print (OUTPUT $inputline); }
		else
		{print (LOG "chain $inputchain res $inputresname resi $inputresi chopped off due to head or tail\n");}
	last;
	}
	}
if($lock==0)
{
	if(!($inputresname eq "GLY"|| $inputresname eq "ALA" || $inputresname eq "PRO"|| $inputresname eq "VAL"|| $inputresname eq "LEU"|| $inputresname eq "ILE"|| $inputresname eq "MET" || $inputresname eq "PHE"|| $inputresname eq "TYR"|| $inputresname eq "TRP"|| $inputresname eq "SER" || $inputresname eq "THR"|| $inputresname eq "CYS"|| $inputresname eq "ASN" || $inputresname eq "GLN" || $inputresname eq "LYS"|| $inputresname eq "HIS"|| $inputresname eq "ARG"|| $inputresname eq "ASP"|| $inputresname eq "GLU"))
	{ print (OUTPUT $inputline);}
	else
	{ print (LOG "chain $inputchain res $inputresname resi $inputresi chopped off due dssp\n"); }

}


}

close(OUTPUT);
close(LOG);

unlink 'temp.pdb';
