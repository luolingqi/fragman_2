#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  pdb_psf_concat.pl
#
#        USAGE:  ./pdb_psf_concat.pl psf1 pdb1 psf2 pdb2 [psf3 pdb3]
#
#  DESCRIPTION:  Joins 2 or more pdbs and psfs
#
#      OPTIONS:  --rtf rtf file
#                --prefix prefix for output pdb and psf
#                --psfgen psfgen executable
#                --wdir working directory
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  David Hall (mn), drhall@bu.edu
#      COMPANY:  Structural Bioinformatics Lab, Boston University
#      VERSION:  1.0
#      CREATED:  07/09/2009 04:35:39 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;
use Fatal qw(open close);
use File::Glob ':glob';
use File::Basename;
use Getopt::Long;

my $psfgen = "psfgen";
my $ertf   = "~/prms/trunk/pdbamino.rtf";
my $wdir   = "";
my $prefix = "";
my @tofix;
my $nfix       = 0;
my $mod;
my $modid;
my $lastmodid  = "nan";

GetOptions(
    "psfgen=s"      => \$psfgen,
    "rtf=s"         => \$ertf,
    "prefix=s"      => \$prefix,
    "wdir=s"        => \$wdir,
);

if (   scalar @ARGV < 4 ||
     ( scalar @ARGV % 2 ) != 0 ) {
    print STDERR "usage: $0 [options] psf1 pdb1 psf2 pdb2 [psf3 pdb3]\n";
    exit 1;
}


#expand ~ in rtf file name
my $rtf     = bsd_glob($ertf, GLOB_TILDE);
my $oscript = "$wdir$$.inp";

open OSCRIPT, '>', $oscript;
print OSCRIPT "# RTF FILE\n";
print OSCRIPT "topology $rtf\n";

my $construct_prefix = 1;
if ( $prefix ne "" ) {
    $construct_prefix = 0;
}

my @pdbs;
print OSCRIPT "# READ IN PSF, PDB\n";
while ( scalar @ARGV ) {
    my $psf = bsd_glob( shift @ARGV, GLOB_TILDE );
    my $pdb = bsd_glob( shift @ARGV, GLOB_TILDE );
    if ( $construct_prefix == 1) {
        my $base = basename( $pdb, '.pdb' );
        $prefix .= $base;
    }
    print OSCRIPT "readpsf $psf\n";
#    print OSCRIPT "coordpdb $pdb\n";

    push @pdbs, $pdb;
}

print OSCRIPT "# WRITE NEW PSF, PDB\n";
print OSCRIPT "writepsf charmm $wdir$prefix.psf\n";
#print OSCRIPT "writepdb $wdir$prefix.pdb\n";
close OSCRIPT;
my $out = $oscript . "\.out";
system "$psfgen < $oscript>$out";

system "cat @pdbs > $wdir$prefix.pdb";
system "pdb.renumberAtoms.pl $wdir$prefix.pdb";
