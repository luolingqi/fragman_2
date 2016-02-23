#!/usr/bin/perl

if (scalar @ARGV < 2)
{
        print "usage: $0 pdbfile msfile\n";
        exit;
}


my @atoms, @residues, @chains, @lines;

open PDB, "<",  $ARGV[0] or die "Failed Opening PDB\n";
while (<PDB>)
{
	chomp;
	if (/^ATOM/)
	{
		push @atoms, substr $_, 12, 4;
		push @chains, substr $_, 73, 1;
		push @residues, substr $_, 22, 4;
	}
	s/\s//g for @atoms;
}
close PDB;

#read ms file into an array so the file can be overwritten
open MS,  "<",  $ARGV[1] or die "Failed Opening MS\n";
while (<MS>)
{
	if (/^ATOM/) {
		push @lines, ($_);
	}
}
close MS;

open MS,  ">",  $ARGV[1] or die "Failed Writing to MS\n";
foreach (@lines)
{
	substr($_, 21, 1) = shift @chains;
	substr($_, 22, 4) = shift @residues;
	print MS $_;
}

close MS;
