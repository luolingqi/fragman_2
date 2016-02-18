#!/usr/bin/perl -w
use strict;

my ($pdb, $new, @pdblines, $pdbline, @words, $part1, $part2, $i, $tp, $no, $cat, $het);

$pdb = shift(@ARGV);
$new = "$pdb.new";
$het = "HETATM";

open (PDB, "<$pdb") || die "cannot open $pdb: !\n";
chomp(@pdblines = <PDB>); close PDB;

open (OUT, ">$new") || die "cannot open $new: !\n";

$i = 1;
foreach $pdbline(@pdblines)
{
	@words = split(/\s+/, $pdbline);
	$tp = $words[0];
	$no = $words[2];
########if (defined $no)
########{
########	$no = substr($no, 0, 1);			# takes only first letter of atom name representative of atom element
########}

	if (($tp eq "ATOM") or ($tp eq "HETATM"))
	{
		$part1 = substr($pdbline, 7,  4);
		$part2 = substr($pdbline,17, 61);

		$cat = "$no$i";
		printf (OUT "%6s %4s %4s %61s\n", $het,$part1,$cat,$part2);
		
		$i++;
	}
	else
	{
		printf (OUT "%s\n", $pdbline);
	}
#	else
#	{
#		die "check pdb file after rdkit(.sdf) and openbabel(.pdb): !\n";
#	}
}

close OUT;
system "mv $new $pdb";
