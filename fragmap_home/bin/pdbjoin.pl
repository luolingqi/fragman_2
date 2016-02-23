#! /usr/bin/perl
##########################################################
# pdbnmd.pl
# generate input pdb/psf.
##########################################################
use strict;
use warnings;
#use Fatal qw(open close);
use File::Basename;
use Getopt::Long;

my $home = $ENV{'HOME'};

my $wdir       = "";
my $osuffix = '_join';
GetOptions(
    "wdir=s"        => \$wdir,
    "osuffix:s"     => \$osuffix,
);

if ( scalar @ARGV < 2 ) {
    exit 1;
}
my $pdb    = shift @ARGV;
my @chains = @ARGV;

# check name
my ( $pdbname, $pdbdir, $pdbsuffix ) = fileparse( $pdb, qr/\.[^.]*/ );
if ( $pdbsuffix eq ".pdb" ) {
    $pdb = substr( $pdb, 0, length($pdb) - 4 );
}

# deal with ? chain id
my $ncho = -1;
my @nk;
my $mi;
foreach my $cho (@chains) {
    $ncho = $ncho + 1;
    if ( $cho eq "?" ) {
        @nk = ( @nk, $ncho );
    }
}

if (@nk) {
    for ( $mi = @nk - 1 ; $mi >= 0 ; $mi-- ) {
        splice( @chains, $nk[$mi], 1 );
    }

    my @filss = glob("$wdir$pdb-?.*.pdb");
    chomp @filss;
    @filss = sort { substr( $a, -9 ) cmp substr( $b, -9 ) } @filss;
    foreach my $lin (@filss) {
        my $nu = length($lin) - 10;
        my $no = length( $wdir . $pdb . "-" );
        $lin = substr( $lin, $no, $nu - $no + 1 );
    }
    my %chainh = map { $_, 1 } @chains;    # get unique elements
    for (@filss) {
        next if $chainh{$_}++;
        push @chains, $_;
    }
}



my $cha;

open my $OUT, '>', "$wdir$pdb$osuffix.pdb";

foreach my $cho (@chains) {
    $cho = lc($cho);
    my @files = glob("$wdir$pdb-$cho.*.pdb");
    chomp @files;
    foreach my $fil (@files) {
		open my $IN, '<', $fil;
		while (<$IN>) {
			print $OUT $_ unless /END/;
		}
		close($IN);
    }
}
close($OUT);
