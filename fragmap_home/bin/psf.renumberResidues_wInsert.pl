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
my $inside = 0;
my $resi;

open my $PDB, '<',  $pdb;
open my $OUT, '>',  $ofile;
while (<$PDB>) {
	if (/NBOND/) {
		$inside = 0;
	}
	if ( $inside == 1 && length($_) > 17 ) {

		$currentResNum = substr($_, 14, 5);
		if ($currentResNum ne $lastResNum) {
			$counter++;
			if ($counter == 1000) {
				$counter = 1;
				$insert = 'A';
			}
		}
		$lastResNum = $currentResNum;

		$resi = sprintf("%i", $counter);
		$resi = sprintf("%-5s", $resi.$insert);
		substr($_, 14, 5) = $resi;
	}
	print $OUT $_;
	if (/NATOM/) {
		$inside = 1;
	}
}

close($PDB);
close($OUT);
