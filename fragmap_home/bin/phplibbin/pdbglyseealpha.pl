#!/usr/bin/perl
use strict;
use warnings;

if (scalar @ARGV < 2)
{
	print "usage: $0 maskfile outfile\n";
	exit;
}

my $maskfile = shift @ARGV;
my $outfile  = shift @ARGV;

open MASK, "<", $maskfile or die "Failed opening maskfile $maskfile\n";
my @lines = <MASK>;
close MASK;

open OUT, ">", $outfile or die "Failed opening outfile $outfile\n";
foreach (@lines)
{
	next unless /^ATOM/;
	#blank out letter before (deleting alternate location)
	#for libmol compatibility
	substr($_, 16, 5, ' GLY ');
	substr($_, 12, 4, ' CA ');

	print OUT $_;
}

close OUT;

