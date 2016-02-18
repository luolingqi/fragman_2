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
my $lastChain  = ' ';
my $currentChain;
my $counter = 0;
my $currentResNum;
my $insert = ' ';

open my $PDB, '<',  $pdb;
open my $OUT, '>',  $ofile;
while (<$PDB>) {
	if ( /^(HETATM)|^(ATOM)/ ) {

		$currentResNum = substr($_, 22, 5);
        $currentChain  = substr($_, 21, 1);
		if ($currentResNum ne $lastResNum) {
			$counter++;
			if ($counter == 1000) {
				$insert = 'A';
			}
			if ($counter == 2000) {
				$insert = 'B';
			}

			my $curnum = substr($currentResNum,0,4);
			my $lastnum = substr($lastResNum,0,4);
			my $delta = $curnum -$lastnum;
			if ($delta > 1 || $delta < 0) {
				$counter++;
				if ($counter == 1000) {
					$insert = 'A';
				}
				if ($counter == 2000) {
					$insert = 'B';
				}
			}
		}
		if ($currentChain ne $lastChain) {
			$counter = 1;
			$insert = ' ';
		}
		$lastResNum = $currentResNum;
		$lastChain  = $currentChain;

		$counter = sprintf("%4i", $counter);
		substr($_, 22, 5) = $counter.$insert;
	}
	print $OUT $_;
}

close($PDB);
close($OUT);
