#!/usr/bin/perl
use strict;
use warnings;
use Fatal qw(open close);

#model list
open my $MODEL, '>>', 'modellist';

for my $file (@ARGV)
{

	open my $PDB, '<', $file;
	
	my $count = 0;
	my $OUT;  #output file

	while (<$PDB>)
	{

		#open a file for each model
		if ( /^MODEL/ )
		{
			$count++;
			open $OUT, '>', $file.$count;
			print $MODEL $file.$count, "\n";
		}

		print $OUT $_ if /^ATOM/;

		#close the file at the end of the model
		close $OUT if /^ENDMDL/;

	}
	close $PDB;
}
close $MODEL;


