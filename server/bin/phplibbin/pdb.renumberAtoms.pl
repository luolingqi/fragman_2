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
  print("Usage: pdb.renumberAtoms.pl pdb_file\nThis will renumber the file 'pdb_file'\n");
}


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
      #       $currentResName = substr($_, 17, 3);
      $currentResNum = substr($_, 22, 4);
      $currentAtomNum = substr($_, 6, 5);
      $currentResNum =~ s/\s+//g;
      $currentAtomNum =~ s/\s+//g;
      #      if ($currentResNum != $lastResNum || $currentAtomNum != ($lastAtomNum + 1)) {
      $counter++;
      #      }
      #      $lastAtomNum = $currentAtomNum;


      if ($counter >= 99900 && $lastResNum != $currentResNum) {
	$counter = 1;
      }

      $counter = sprintf("%5i", $counter);
      substr($_, 6, 5) = $counter;
      #print $_;
      print TMP $_;

      $lastResNum = $currentResNum;

    } else {
      # print "";

	$counter = 0 if ($_=~/HEADER/);

      if(!($_=~/END/))
      {
      print TMP $_;
      }
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
