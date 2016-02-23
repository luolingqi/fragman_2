#!/usr/bin/perl 

use strict;
use warnings;

if ( @ARGV < 1 ) {
    print "usage:$0 nnb.out\n";
    exit;
}

my $file = shift;

my $outfile = "$file.rawextract";

open FILE, '<', $file  or die;

my %count;
my %name;

while (<FILE>) {
	
	 my $line = $_;

    if (!(  (      ( substr( $line, 10, 1 ) eq "-" )
                && ( substr( $line, 30, 1 ) eq "-" )
            )
            || (   ( substr( $line, 10, 1 ) ne "-" )
                && ( substr( $line, 30, 1 ) ne "-" ) )
        )
        )
    {
        if ( substr( $line, 10, 1 ) ne "-" ) {
            my $resi = substr( $line, 11, 4 );
            $resi =~ s/^0+//;
            my $chain = substr( $line, 10, 1 );

            my $resn = substr( $line, 16, 3 );

            #$count[$resi]=$count[$resi]+1;
            #$name[$resi]=$resn;

            my $tag = "$resi\t$chain";

            $count{$tag}++;
            $name{$tag}  = $resn;

            next;
        }

        if ( substr( $line, 30, 1 ) ne "-" ) {
            my $resi = substr( $line, 31, 4 );
            $resi =~ s/^0+//;
            my $chain = substr( $line, 30, 1 );
            my $resn  = substr( $line, 36, 3 );

            #$count[$resi]=$count[$resi]+1;
            #$name[$resi]=$resn;

            my $tag = "$resi\t$chain";

            $count{$tag}++;
            $name{$tag}  = $resn;

            next;
        }

    }
}

my @key1 = keys(%count);

open OUTFILE, '>', $outfile  or die;

foreach (@key1) {
    print OUTFILE "$_\t$name{$_}\t$count{$_}\n";
}
