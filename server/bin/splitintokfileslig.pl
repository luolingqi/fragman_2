#!/usr/bin/perl 
use strict;
use warnings;

if ( @ARGV < 2 ) {
    print "usage:$0 pdblist k\n";
    exit;
}
my $pdblist = shift;
my $k       = shift;

open LIST, '<', $pdblist  or die;
my @pdbs = <LIST>;
close LIST;

chomp @pdbs;

foreach my $file (@pdbs) {
    open FILE, '<', $file or die;
    my @content = <FILE>;
    close(FILE);

    my $filelength = @content;

    my $bin = int( $filelength / $k );

    for ( my $l = 0; $l < $k; $l++ ) {

        my $fileout = $file . $l;
        open OUT, '>', $fileout;
        my $start  = $l * $bin;
        my $ceil   = ( $l + 1 ) * $bin;

        for ( my $j = $start; $j < $ceil; $j++ ) {
            print( OUT $content[ $j ] );
        }
        close(OUT);

    }
}
