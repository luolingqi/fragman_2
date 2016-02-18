#!/usr/bin/perl


use strict;
#use Data::Dumper;$Data::Dumper::Indent = 1;
#use Storable;
#use Getopt::Long;
#require ("./PDB.pl");

my $lastResNum;
my $counter = 0;
my $lastResNum;
my $lastAtomNum;
my $currentResNum;
my $currentResName;
my $currentAtomNum;
my $counter = 0;


if (scalar(@ARGV) < 1 || scalar(@ARGV) != 1) {
  print("Usage: pdb.renumberResidues.pl pdb_file\nThis will renumber the file 'pdb_file'\n");
}

#my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);

my @files = @ARGV;

my $line = '';

for (my $i = 0; $i < scalar(@files); $i++) {
  my $file = $files[$i].".tmp";
  $counter = 0;
  $lastResNum = 0;
  open(PDB1, $files[$i]) or die("$!: $files[$i]");
  open(TMP, "> $file") or die($!);
  while (<PDB1>) {
    if ($_ =~ /^ATOM/) {
      $line = $_;
      $currentResName = substr($_, 17, 3);
      $currentResNum = substr($_, 22, 5);
     # print "$currentResName\t$currentResNum\n";
      $currentAtomNum = substr($_, 6, 5);
      $currentResNum =~ s/\s+//g;
      $currentAtomNum =~ s/\s+//g;
      if ($currentResNum ne $lastResNum || $currentAtomNum != ($lastAtomNum + 1)) {
	 $counter++;
      }
      $lastResNum = $currentResNum;
      $lastAtomNum = $currentAtomNum;

      $counter = sprintf("%4i", $counter);
      substr($_, 22, 5) = $counter." ";
      #print $_;
      print TMP $_;
    }
  }
  #print the last one
#   $counter++;
#   substr($line, 22, 4) = $counter;
#   #print $_;
#   print TMP $line;

  close(PDB1);
  close(TMP);
  system("mv $file $files[$i]");
}
