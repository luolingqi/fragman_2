#! /usr/bin/perl -w
##########################################################
# apbspdb.pl
# generape pb potential using apbs
##########################################################
use strict;
use File::Basename;
use Getopt::Long;

my $home = $ENV{'HOME'};
my $debug=0;

# preset variables

if (scalar @ARGV < 7)
{
	print_short_args ();
	exit;
}

my $rtf=shift @ARGV;
my $par=shift @ARGV;
my $wdir=shift @ARGV;
$wdir .= '/'; #add training slash to directories
my $pdb=shift @ARGV;
my $mode=shift @ARGV;
my $lnpbe=shift @ARGV;
my $chgm=shift @ARGV;

my $imode=-1;
$mode=lc($mode);
if($mode eq "orig"){$imode=0;}
if($mode eq "new"){$imode=1;}
if($mode eq "newtors"){$imode=2;}
if($imode eq -1)
{
   print "ERROR: unknown mode of operation $mode\n";
   exit 0;
}



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

my $dcel=0.4; # grid spacing (mesh size)
my $epsw=80.0; # dielectric constant of solvent (water)
my $epsp=4.0; # dielectric constant of protein interior
my $conc=0.1; # salt concentration
my $dist=12.0; # distance from the edge of a molecule to the boundary
my $wrad=1.4; #  solute radius to shift dielectric boundary

my $cmd="pbinpgen $wdir$pdbname.psf $par $rtf $wdir$pdbname.pdb $imode out.pqr >dimens";
print "$cmd\n";
system $cmd;
#system "pbinpgen $wdir$pdbname.psf $par $rtf $wdir$pdbname.pdb $imode out.pqr >dimens";

open DIM, "< dimens" or die "couldn't open file dimens: $!\n";
my @dimlines = <DIM>;
chomp @dimlines;
close DIM;
my $dx=0;
my $dy=0;
my $dz=0;
for my $line (@dimlines)
{
    if(substr ($line, 0, 10) eq "DIMENSIONS")
    {

        my @loo=split(/ /, $line);
        $dx=$loo[1];
        $dy=$loo[2];
        $dz=$loo[3];
    }
}

if($dx eq 0 || $dy eq 0 || $dz eq 0)
{
  print "ERROR: cannot read dimensions $dx $dy $dz\n";
  exit 0;
}

$dx=$dx+2*$dist;
$dy=$dy+2*$dist;
$dz=$dz+2*$dist;
my $ndx=int($dx/$dcel)+1;
my $ndy=int($dy/$dcel)+1;
my $ndz=int($dz/$dcel)+1;
my $j=0;
while($j*32+1 < $ndx){$j++;}
$ndx=$j*32+1;
$j=0;
while($j*32+1 < $ndy){$j++;}
$ndy=$j*32+1;
$j=0;
while($j*32+1 < $ndz){$j++;}
$ndz=$j*32+1;

open OSCRIPT, "> $oscript" or die "couldn't open file: $oscript\n";

print OSCRIPT "   read\n";
print OSCRIPT "       mol pqr out.pqr\n";
print OSCRIPT "   end\n";
print OSCRIPT "   elec name out\n";
print OSCRIPT "       mg-manual\n";
print OSCRIPT "       bcfl sdh\n";
print OSCRIPT "       chgm $chgm\n";
print OSCRIPT "       dime $ndx $ndy $ndz\n";
print OSCRIPT "       gcent mol 1\n";
print OSCRIPT "       grid $dcel $dcel $dcel\n";
print OSCRIPT "       ion charge 1 conc $conc radius 2.0\n";
print OSCRIPT "       ion charge -1 conc $conc radius 2.0\n";
print OSCRIPT "       $lnpbe\n";
print OSCRIPT "       mol 1\n";
print OSCRIPT "       pdie $epsp\n";
print OSCRIPT "       sdie $epsw\n";
print OSCRIPT "       srfm mol\n";
print OSCRIPT "       chgm spl0\n";
if($mode eq "newtors")
{
	print OSCRIPT "       srad 0.0\n";
}
else
{
	print OSCRIPT "       srad $wrad\n";
}
print OSCRIPT "       swin 0.3\n";
print OSCRIPT "       sdens 10.0\n";
print OSCRIPT "       temp 298.15\n";
print OSCRIPT "       calcenergy total\n";
print OSCRIPT "       calcforce no\n";
print OSCRIPT "       write pot dx $wdir$pdbname\_apbs_pot_$mode\n";
print OSCRIPT "       write kappa dx $wdir$pdbname\_apbs_kappa_$mode\n";
print OSCRIPT "   end\n";
print OSCRIPT "   quit\n";
close OSCRIPT;
system "apbs $oscript";

#convert dx file from atomic units of kT to kcal/mol
my $dx_file = "$wdir$pdbname\_apbs_pot_$mode.dx";
system "dx_scale 0.5922 $dx_file $dx_file";

#unlink $oscript, 'dimens', 'out.pqr', 'io.mc';



sub print_short_args
{
	printf ("usage: %s  rtf par workdir pdbfile mode(orig/new/newtors) linearity(lpbe/npbe) charge_distribution (spl0/spl2) \n", $0);
}
