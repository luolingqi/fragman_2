#! /usr/bin/perl
use strict;
use warnings;
use File::Basename;

if (scalar @ARGV < 3)
{
	print "usage: $0 REC LIG MODELI\n";
	exit;
}

my $rec = shift @ARGV;
my $lig = shift @ARGV;
my $modeli = shift @ARGV;

# create pymol script
my $oscript = "$$.pml";
open OSCRIPT, '>>', $oscript or die "couldn't open file $oscript: $!\n";

print OSCRIPT "set antialias, 1\n";
print OSCRIPT "set depth_cue, 0\n";
print OSCRIPT "bg_color white\n";
print OSCRIPT "set ray_trace_fog, 0\n";

print OSCRIPT "load $rec, rec\n";
print OSCRIPT 'util.cbss("rec","cyan","magenta","salmon")' . "\n";

print OSCRIPT "load $lig, lig\n";
print OSCRIPT 'util.cbss("lig*","red","yellow","green")' . "\n";
print OSCRIPT "show_as cartoon\n";

print OSCRIPT 'cmd.orient("visible",animate=-1)' . "\n";
print OSCRIPT "png model.$modeli.png, width=320, height=240, dpi=72, ray=1\n";

close OSCRIPT;


my $cmd = "pymol -qc $oscript";
(system $cmd) == 0 or die "system $cmd failed: $!\n";
$cmd = "optipng -quiet -o7 model.$modeli.png 2>/dev/null";
(system $cmd) == 0 or die "system $cmd failed: $!\n";

unlink($oscript);
