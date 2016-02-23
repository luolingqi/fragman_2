#!/usr/bin/perl
use strict;
use warnings;
use Fatal qw(open close);

my $file = $ARGV[0];

open FILE, '<', $file;
my @line = <FILE>;
close FILE;

open FILE, '>', $file;
foreach my $line (@line) {
    chomp $line;
    if ( $line =~ /^ATOM/ ) {
        my $chain =
          ( ' ' eq substr( $line, 73, 1 ) )
          ? substr( $line, 72, 1 )
          : substr( $line, 73, 1 );
        substr( $line, 21, 1 ) = $chain;
    }
    print FILE "$line\n";
}
close FILE;
