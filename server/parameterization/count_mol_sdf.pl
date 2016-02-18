#!/usr/bin/perl -w
use strict;

#Scott Mottarella
#Passed an sdf file, prints number of molecules in that file
#based on number of molecule termination lines ('$$$$')

my $count=0;
open IN, "$ARGV[0]" or die $!;
while (<IN>) {
	if ($_ =~ /^\s*\$\$\$\$\s*$/) {
		$count++;
	}
}
close IN;
print "$count\n";
