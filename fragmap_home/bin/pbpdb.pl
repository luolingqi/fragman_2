#! /usr/bin/perl -w
##########################################################
# pdbchm.pl
# structure charmm minimization on input pdb.
##########################################################
use strict;
use File::Basename;
use Getopt::Long;

my $home = $ENV{'HOME'};
my $debug=0;

# preset variables

if (scalar @ARGV < 5)
{
	print_short_args ();
	exit;
}

my $rtf=shift @ARGV;
my $par=shift @ARGV;
my $wdir=shift @ARGV;
my $pdb=shift @ARGV;
my $mode=shift @ARGV;
my $oscript="$wdir$$.inp";
$pdb=lc($pdb);

# fileparse
my ($pdbname, $pdbdir, $pdbsuffix) = fileparse ($pdb, qr/\.[^.]*/);
my $pdbpre = $pdbdir . $pdbname;
if ($pdbsuffix ne ".pdb")
{
	print STDERR "error: $pdb must have ext .pdb\n";
	exit 1;
}
if (! -e "$wdir$pdbname.psf")
{
	print STDERR "error: $wdir$pdbname.psf file not found\n";
	exit 1;
}

my @chains = @ARGV;
my ($cha, $fil);

open OSCRIPT, "> $oscript" or die "couldn't open file: $oscript\n";
print OSCRIPT "* process protein\n";
print OSCRIPT "* multiple chains with breaks\n";
print OSCRIPT "*\n";
print OSCRIPT "\n";

print OSCRIPT "! RTF AND PARAM FILEiS\n";
print OSCRIPT "open read card unit 2 name $rtf\n";
print OSCRIPT "read rtf card unit 2\n";
print OSCRIPT "close unit 2\n";
print OSCRIPT "open read card unit 2 name $par\n";
print OSCRIPT "read param card unit 2\n";
print OSCRIPT "close unit 2\n";
print OSCRIPT "\n";

print OSCRIPT "! read psf and pdb\n";
print OSCRIPT "open unit 1 read card name $wdir$pdbname.psf\n";
print OSCRIPT "read psf card unit 1\n";
print OSCRIPT "\n";

my $dcel=0.4; # grid spacing (mesh size)
my $epsw=80.0; # dielectric constant of solvent (water)
my $epsp=4.0; # dielectric constant of protein interior
my $conc=0.1; # salt concentration
my $dist=12.0; # distance from the edge of a molecule to the boundary
my $wrad=1.4; #  solute radius to shift dielectric boundary
my $wradsa=1.6; #  solute radius to shift dielectric boundary for sa atoms
my $brad=5.75; # "big" sphere radius to define cavities

print OSCRIPT "open unit 1 read card name $wdir$pdbname.pdb\n";
print OSCRIPT "read coor pdb unit 1\n";
print OSCRIPT "\n";
#print OSCRIPT "coor orient norot\n";
if($mode eq "new")
{
	print OSCRIPT "coor surf acce rpro $brad \n";
	print OSCRIPT "scalar wmain sign\n";
	print OSCRIPT "define ones select prop wmain gt 0.0 show end\n";
	print OSCRIPT "define twos select ones .or. (.BONDED. ones) show end\n";
	print OSCRIPT "define ones select twos .or. (.BONDED.  twos) show end\n";
	print OSCRIPT "scalar wmain set 1.0 select ones end\n";
	print OSCRIPT "scalar wmain mult -$wrad\n";
	print OSCRIPT "scalar wcomp =  radius\n";
	print OSCRIPT "scalar wcomp add $wrad\n";
	print OSCRIPT "scalar wmain sum wcomp\n";
}
elsif($mode eq "newtors")
{
	print OSCRIPT "scalar wmain =  radius\n";
	print OSCRIPT "scalar wmain mult 0.0  select .not. ( TYPE CA .OR. TYPE N .OR. TYPE C .OR. TYPE O ) end\n";
	print OSCRIPT "coor weight surf acce rpro $brad \n";
	print OSCRIPT "scalar wmain sign\n";
	print OSCRIPT "define ones select prop wmain gt 0.0 show end\n";
	print OSCRIPT "define twos select ones .or. (.BONDED. ones) show end\n";
	print OSCRIPT "define threes select twos .or. (.BONDED. twos) show end\n";
	print OSCRIPT "define ones select threes .or. (.BONDED.  threes) show end\n";
	print OSCRIPT "scalar wmain set 1.0 select ones end\n";
	print OSCRIPT "scalar wmain mult -$wradsa\n";
	print OSCRIPT "scalar wcomp =  radius\n";
	print OSCRIPT "scalar wcomp add $wrad\n";
	print OSCRIPT "scalar wmain sum wcomp\n";
	print OSCRIPT "scalar wmain set 0.1 select ones end\n";
}
else
{
   print OSCRIPT "scalar wmain =  radius\n";
}

my $osuffix = "_accs_$mode"; # output file suffix
print OSCRIPT "\n";
print OSCRIPT "coor stat\n";
print OSCRIPT "open unit 1 write card name $wdir$pdbname$osuffix.pdb\n";
print OSCRIPT "write coor pdb unit 1\n";
print OSCRIPT "\n";

print OSCRIPT "! number of grid points in each direction\n";
print OSCRIPT "calc nclx  = int ( ( ?xmax - ?xmin + 2.0 * $dist ) / $dcel ) + 1\n";
print OSCRIPT "calc ncly  = int ( ( ?ymax - ?ymin + 2.0 * $dist ) / $dcel ) + 1\n";
print OSCRIPT "calc nclz  = int ( ( ?zmax - ?zmin + 2.0 * $dist ) / $dcel ) + 1\n";
print OSCRIPT "calc xc = ( ?xmax + ?xmin ) / 2.0\n";
print OSCRIPT "calc yc = ( ?ymax + ?ymin ) / 2.0\n";
print OSCRIPT "calc zc = ( ?zmax + ?zmin ) / 2.0\n";
print OSCRIPT "\n";

print OSCRIPT "set format (f15.5)\n";
print OSCRIPT "\n";

print OSCRIPT "open unit 3 write form name $wdir$pdbname$osuffix.pot\n";
print OSCRIPT "open unit 2 write form name $wdir$pdbname$osuffix.kap\n";
print OSCRIPT "\n";

print OSCRIPT "PBEQ\n";

print OSCRIPT "prnlev 5\n";
print OSCRIPT "scalar wmain statistics select .not. type H* end\n";
print OSCRIPT "define check select (.not type H* ) .and. ( property wmain .eq. 0.0 ) show end\n";
print OSCRIPT "if ?nsel ne 0  stop  !some heavy atom have a zero radius\n";
print OSCRIPT "SOLVE NCLX \@nclx NCLY \@ncly NCLZ \@nclz XCEN \@xc YCEN \@yc ZCEN \@zc DCEL $dcel -\n";
print OSCRIPT "EPSW $epsw EPSP $epsp CONC $conc\n";
print OSCRIPT "\n";

print OSCRIPT "WRITE PHI KCAL/MOL UNIT 3\n";
print OSCRIPT "WRITE FKAPPA2 UNIT 2\n";
print OSCRIPT "\n";

print OSCRIPT "RESET\n";
print OSCRIPT "END\n";
print OSCRIPT "\n";

print OSCRIPT "close unit 3\n";
print OSCRIPT "close unit 2\n";
print OSCRIPT "stop\n";
close OSCRIPT;
my $out = $oscript."\.out";
system "charmm27 < $oscript>$out";
#system "rm $oscript";
#system "rm $out";



sub print_short_args
{
	printf ("usage: %s  rtf par workdir pdbfile mode(orig/new)\n", $0);
}
