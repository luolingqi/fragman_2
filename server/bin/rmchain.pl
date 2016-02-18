#!/usr/bin/perl 
use strict;
use warnings;

if ( @ARGV < 1 ) {
    print "usage:$0 file\n";
    exit;
}

my $file = shift;
open FILE, '<', $file  or die;
my @lines = <FILE>;
close FILE;

open OUT, '>', $file;

foreach (@lines) {

    if ( $_ =~ /^ATOM/ || $_ =~ /^HETATM/ ) {
        substr($_, 21, 1) = ' ';

    }
    print OUT $_;

}

close OUT;
