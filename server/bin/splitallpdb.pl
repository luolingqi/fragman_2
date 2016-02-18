#!/usr/bin/perl 
use strict;
use warnings;

if ( @ARGV < 2 ) {
    print "usage:$0 listallpdb complexs3\n";
    exit;
}

my $list      = shift;
my $complexs3 = shift;

open LIST, '<', $list or die;
my @pdbs = <LIST>;
close LIST;

chomp @pdbs;

open COMPLEX, '<', $complexs3 or die;
my @probes = <COMPLEX>;
close COMPLEX;

chomp @probes;

#in pdb file, probes are characterized
#by the last 3 letters in their ids
$_ = uc( substr( $_, -3 ) ) for @probes;

open SLIST, '>', 'splitlist' or die;

foreach my $pdb (@pdbs) {

    foreach my $probe (@probes) {
        my $outfile = $pdb . $probe;
        $outfile =~ s/\.//g;
        print SLIST "$outfile\n";
        system "grep $probe $pdb >$outfile";

    }
}

close SLIST;
