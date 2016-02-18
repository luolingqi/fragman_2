#! /usr/bin/perl -w

##########################################################
# charmm.min.1.pl
# -------------------
# copyright : (C) 2006 by Ryan Brenke
# email : rbrenke@gmail.com
# -------------------
# Description:
# Parses charmm script template and performs single
# structure charmm minimization on input pdb.
##########################################################

use strict;
use File::Basename;

my $home = $ENV{'HOME'};

if (scalar @ARGV < 5)
{
	print "usage: $0 RTF PRM TMPLT N(10)STEPS PDB\n";
	exit;
}

my $rtf=shift @ARGV;
my $prm=shift @ARGV;
my $tmplt=shift @ARGV; # charmm template script
my $nsteps=shift @ARGV; # number of minimization steps
my $rec=shift @ARGV;

# create names
my ($recname, $recdir, $recsuffix) = fileparse ($rec, qr/\.[^.]*/);
my $recpre = $recdir . $recname;
my $recmin = $recpre . ".min" . $recsuffix;
# get segid
open REC, "< $rec" or die "couldn't open file $rec: $!\n";
my @reclines = <REC>;
my $recline = $reclines[0];
my $segidrec = substr ($recline, 72, 4);
close REC;

my $oscript="$$.in";

open TMPLT, "< $tmplt" or die "couldn't open file: $tmplt\n";
my @lines = <TMPLT>;
close TMPLT;

open OSCRIPT, "> $oscript" or die "couldn't open file: $oscript\n";

for my $line (@lines)
{
	$line =~ s/_NSTEPS_/$nsteps/gs;
	$line =~ s/_PRM_/$prm/gs;
	$line =~ s/_REC_/$rec/gs;
	$line =~ s/_RECMIN_/$recmin/gs;
	$line =~ s/_RTF_/$rtf/gs;
	$line =~ s/_SEGID_/$segidrec/gs;
	print OSCRIPT $line;
}

close OSCRIPT;

system "charmm27 < $oscript";
system "rm $oscript";
