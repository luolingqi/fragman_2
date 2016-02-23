#!/usr/bin/perl
use strict;
use warnings;

my $pdbfiles = shift @ARGV;
open PDBFILES, '<', $pdbfiles;
my @files = <PDBFILES>;
chomp @files;
close PDBFILES;

foreach my $file (@files) {
    open INFILE, '<', $file;
    my @data = <INFILE>;
    chomp @data;
    close INFILE;
    open OUTFILE, '>', $file;

    foreach my $line (@data) {
        if (   ( $line =~ /CCN/ )
            || ( $line =~ /ACN/ )
            || ( $line =~ /EOH/ )
            || ( $line =~ /_NC/ )
            || ( $line =~ /HBE/ )
            || ( $line =~ /BNZ/ )
            || ( $line =~ /TBU/ )
            || ( $line =~ /CHX/ )
            || ( $line =~ /PES/ )
            || ( $line =~ /EHN/ )
            || ( $line =~ /IPH/ )
            || ( $line =~ /MOH/ )
            || ( $line =~ /IOH/ )
            || ( $line =~ /URE/ )
            || ( $line =~ /EOL/ )
            || ( $line =~ /DMF/ )
            || ( $line =~ /ACM/ ) )
        {
            $line =~ s/ATOM  /HETATM/;
            print OUTFILE "$line\n";
        }
        else {
            print OUTFILE "$line\n";
        }
    }
    close OUTFILE;
}

