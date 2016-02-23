#!/usr/bin/perl

if (scalar @ARGV < 2)
{
	print "usage: $0 pdbfile maskfile [chain-res list]\n";
	exit;
}

my $pdbfile = shift @ARGV;
my $maskfile = shift @ARGV;

my @chainres = @ARGV;

open MASK, ">", $maskfile or die "Failed opening maskfile $maskfile\n";
open PDB, "<", $pdbfile or die "Failed opening PDB $pdbfile\n";
while (<PDB>)
{
	next unless /^ATOM/;
	$chain    = substr $_, 73, 1;
	$residues = substr $_, 22, 4;
        $residues =~ s/\s//g;
        $linecr   = $chain . '-' . $residues;

        $mask = 0;
	foreach $globalcr (@chainres)
        {
		$mask = 1 if lc($linecr) eq lc($globalcr);
	}
	
	print MASK $_ if $mask == 1;
}

close PDB;
close MASK;

