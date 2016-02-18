#! /usr/bin/perl
use strict;
use warnings;
use File::Basename;

if (scalar @ARGV < 2)
{
	print "usage: $0 pdbfile style\n";
	exit;
}

my $pdbfile = shift @ARGV;
my $style = shift @ARGV;

# create pymol script
my $oscript = "$$.pml";
open OSCRIPT, '>>', $oscript or die "couldn't open file $oscript: $!\n";

print OSCRIPT "set antialias, 1\n";
print OSCRIPT "set depth_cue, 0\n";
print OSCRIPT "bg_color white\n";
print OSCRIPT "set ray_trace_fog, 0\n";

print OSCRIPT "load $pdbfile, pdbfile\n";
if ($style eq "rec" || $style eq "prot")
{
	print OSCRIPT 'util.cbss("pdbfile","cyan","magenta","salmon")' . "\n";
}
elsif ($style eq "lig" || $style eq "hep")
{
	print OSCRIPT 'util.cbss("pdbfile","red","yellow","green")' . "\n";
}

print OSCRIPT "show_as cartoon, pdbfile\n";

# show heparin molecule as sticks
if ($style eq "hep")
{
	print OSCRIPT "show_as sticks, resn A00\n";
}

print OSCRIPT 'cmd.orient("visible",animate=-1)' . "\n";
print OSCRIPT "png $style.png, width=160, height=120, dpi=72, ray=1\n";

close OSCRIPT;

my $cmd = "pymol -qc $oscript";
(system $cmd) == 0 or die "system $cmd failed: $!\n";
$cmd = "optipng -quiet -o7 $style.png 2>/dev/null";
(system $cmd) == 0 or die "system $cmd failed: $!\n";

unlink($oscript);
