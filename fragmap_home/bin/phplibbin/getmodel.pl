#!/usr/bin/perl
use strict;
use warnings;

use Fatal qw(open close);

if (scalar @ARGV < 1)
{
	print STDERR "usage: $0 inpdb [model] [outpdb]\n";
	exit 1;
}

my $pdb = $ARGV[0];

my $model = defined $ARGV[1] ? $ARGV[1] : 1;

my $ofile = defined $ARGV[2] ? $ARGV[2] : '/dev/stdout';

die 'Model must be 1 or greater' if $model < 1;

open my $FILE, '<', $pdb;
#enter slurp mode (see perlvar)
local $/;

my @pieces = split(/\nMODEL/, <$FILE>);
if ($#pieces == 0) {
        @pieces = split(/\nENDMDL/, $pieces[0]);
	unshift @pieces, ' ';
}
close $FILE;

die "PDB has $#pieces models, model $model requested" if ($model >= @pieces);

substr($pieces[$model], 0, 0, "\n") unless (" " eq substr($pieces[$model], 0, 1) );

open my $OFILE, '>', $ofile;
print $OFILE 'MODEL', $pieces[$model], "\n";
close $OFILE;
