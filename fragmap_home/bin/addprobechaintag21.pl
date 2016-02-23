#!/usr/bin/perl

if ( @ARGV < 2 ) {
    print "usage:$0 file chain\n";
    exit;
}

my $file      = shift;
my $recchains = shift;
$recchains =~ s/'//g;
my @recchain = split( /\s+/, $recchains );

my $chainstring = "ZYXWVUTSRQPONMLKJIHGFEDCBA";
foreach (@recchain) {
    if ( length($_) == 1 ) {
        $chainstring =~ s/$_//i;
    }
}


open FILE, '<', $file or die;
my @lines = <FILE>;
close FILE;

open OUT, '>', $file;

my $lock = 0;
my $chain;
my $cluster;

foreach (@lines) {
    if ( $_ =~ /HEADER crosscluster.(\d+).\d+.pdb/ ) {
        $cluster = $1;
        $cluster =~ s/^0//g;
        $cluster =~ s/^0//g;
        $chain = substr( $chainstring, $cluster, 1 );
        $lock = 1;
        print OUT $_;
        next;
    }

    if ( $_ =~ /END/ ) {
        $lock = 0;
        print OUT $_;
        next;
    }

    if ( $lock == 1 ) {
        if ( $_ =~ /^ATOM/ ) {
				substr($_, 21, 1, $chain);	
            print OUT $_;
            next;
        }
    }

    print OUT $_;

}

$output = substr $chainstring, 0, $cluster+1;
print OUT "#PROTEINTAG $recchains\n";
print OUT "#PROBETAG $output\n";
close OUT;
