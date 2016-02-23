#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  fake_msur.pl
#
#        USAGE:  ./fake_msur.pl  IN_FILE OUT_FILE
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  David Hall (mn), drhall@bu.edu
#      COMPANY:  Structural Bioinformatics Lab, Boston University
#      VERSION:  1.0
#      CREATED:  03/18/2010 05:21:32 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;
use Fatal qw(open close);
use 5.010;

if (scalar @ARGV < 2) {
    say STDERR "usage $0 IN_PDB OUT_MS";
    exit 1;
}


my $in_file = shift;
my $out_file = shift;

open my $IN, '<', $in_file;
open my $OUT, '>', $out_file;

while (<$IN>) {
    next unless /^ATOM/;

    substr($_, 56) = " 1.0\n";


    print $OUT $_;

}

close $IN;
close $OUT;
