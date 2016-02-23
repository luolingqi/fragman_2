#! /usr/bin/perl -w
use strict;

if (scalar @ARGV < 1)
{
	print "usage: $0 FILE:CHARMM.RTF STRING:RESI.name\n";
	exit;
}

my $ifile = shift @ARGV;
my $resi = shift @ARGV;

open (IFILE, "< $ifile") or die "couldn't open $ifile\n";
my @lines = <IFILE>;
close IFILE;
chomp @lines;


# readstate -
# 0: not in RESI $resi
# 1: in RESI $resi
# 2: in RESI $resi, and found ATOM line
my $readstate = 0;

for my $line (@lines)
{
	next if ($line =~ m/^!/); # comment

	if ($readstate == 0 && $line =~ m/^RESI\s+($resi)\s+/)
	{
		$readstate = 1;
	}
	elsif (($readstate == 1 || $readstate == 2) && $line =~ m/^ATOM\s+/)
	{
		print "$line\n";
		$readstate = 2;
	}
	elsif ($readstate == 2 && ($line !~ m/^ATOM\s+/ && $line !~ m/^GROU/))
	{
		last;
	}
}
