#!/usr/bin/perl
use strict;
use warnings;
use Fatal qw(open close);


if ( scalar @ARGV < 2 ) {
	print STDERR "Usage: pdb.renumberResidues.pl pdb_file renumbered\n";
}

my $pdb = shift @ARGV;
my $ofile = shift @ARGV;

my $line = '';

my $lastResNum = 0;
my $counter = 0;
my $currentResNum;
my $insert = ' ';

open my $PDB, '<',  $pdb;
open my $OUT, '>',  $ofile;
while (<$PDB>) {
	if ( /^(HETATM)|^(ATOM)/ ) {

		$currentResNum = substr($_, 22, 5);
		if ($currentResNum ne $lastResNum) {
			$counter++;
			if ($counter == 1000) {
				$counter = 1;
				$insert = 'A';
			}
		}
		$lastResNum = $currentResNum;

		$counter = sprintf("%4i", $counter);
		substr($_, 22, 5) = $counter.$insert;
	}
	print $OUT $_;
}

close($PDB);
close($OUT);
