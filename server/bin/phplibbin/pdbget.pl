#! /usr/bin/perl
use strict;
use warnings;
use File::Basename;

my $prefix = "http://www.rcsb.org/pdb/files/"; # pdb prefix

if (scalar @ARGV < 1)
{
	print "usage: $0 pdbid\n";
	exit;
}

my $pdbid = shift @ARGV;
my $ofile = $pdbid;
if (scalar @ARGV > 0)
{
	$ofile = shift @ARGV;
}
my ($ofilename, $ofiledir, $ofilesuffix) = fileparse ($ofile, qr/\.[^.]*/);
my $suffix = ".pdb";
if ($ofilesuffix ne $suffix)
{
	$ofile .= $suffix;
}

my $wget_file = $prefix . $pdbid . '.pdb';
my $cmd = "lwp-download '$ofile' '$wget_file'";
#my $cmd = "curl -s --compressed -f -o '$ofile' '$wget_file'";
#print "$cmd\n";
system $cmd;
#for some reason, have to bit shift $? to get legit
#status.  Perl is weird
exit($? >> 8);
