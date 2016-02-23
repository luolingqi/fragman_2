#!/usr/bin/perl -w


#===============================================================================#
# Here, Spencer (Mon Jan 09 19:58:40 2006)				        #
# 									        #
# You need to make the script look at the name of the current residue,	        #
# then open up the file that contains the connect records.  I might want        #
# to add another level to the has so that separate hashes are created	        #
# for every residue found.						        #
#===============================================================================#


use Data::Dumper;
use strict;

#require("$ENV{'CSMAPPERL'}/PDB.pl");
require("PDB.pl");
if (!@ARGV || scalar(@ARGV) != 2) {
  print "Usage: addConnect.pl pdb /path/to/rtf_directory\n"; exit(1);
}

#my $rtfFile = $ARGV[0];
my $pdbFile = $ARGV[0];
my $rtfDir = $ARGV[1];

my $tmp;

#my %rtfHash;
my @tmpA;
#my %pdbHash;
my $r;


my $resName;

# Renumber the pdb
# my $command = "pdb.renumberAtoms.pl $pdbFile";
# system $command;

# Get the unique residues in the pdb
$tmp = PDB::getSequenceThreeLetter($pdbFile, 0);
my @sequence = split(/\s+/, $tmp);
my %seen = ();
my @uniqueSequence = grep { ! $seen{$_} ++ } @sequence;


use vars qw($lastResNum $currentResNum $currentResName $lastResName $rtfFile $rtfFileNoOne $firstResNum);
$lastResNum = -999;
$lastResName = "";
my @pdbArray;
my $counter = 0;

open (PDB, "$pdbFile") or die($!);
while (<PDB>) {
  if ($_ =~ /^ATOM/ || $_ =~ /^HET/) {
    $pdbArray[$counter] = $_;
    $counter++;
  }
}
close(PDB);

chomp(@pdbArray);

# Make a hash for the pdb where the keys are the resname.resnum.
my %masterPdbHash;

# Once you have the hash, take a unique list of all the residues
# inside the pdb and then create a hash for each of these residues of
# the rtf.

for (my $j = 0; $j < scalar(@pdbArray); $j++) {
  $currentResNum = PDB::getResNum($pdbArray[$j]);
  $currentResName = PDB::getResName($pdbArray[$j]);

  my %pdbHash;
  $pdbHash{PDB::getAtomName($pdbArray[$j])} = PDB::getAtomNum($pdbArray[$j]);

  if (!exists($masterPdbHash{$currentResName.$currentResNum})) {
    $masterPdbHash{$currentResName.$currentResNum} = {};
    $masterPdbHash{$currentResName.$currentResNum}->{'atomHash'} = [];
    $masterPdbHash{$currentResName.$currentResNum}->{'atomHash2'} = {};
    $masterPdbHash{$currentResName.$currentResNum}->{'residue'} = $currentResName;
    $masterPdbHash{$currentResName.$currentResNum}->{'resNum'} = $currentResNum;
  }

  $masterPdbHash{$currentResName.$currentResNum}->{'atomHash2'}->{PDB::getAtomName($pdbArray[$j])} = PDB::getAtomNum($pdbArray[$j]);

  push(@{$masterPdbHash{$currentResName.$currentResNum}->{'atomHash'}}, %pdbHash);
}


# Now that we have the pdbhash, create the rtfhash

my %masterRtfHash;


for (my $i = 0; $i < scalar(@uniqueSequence); $i++) {
  my %rtfHash;
#  $rtfFile = "$ENV{'CSMAP'}/probes/base/1".lc($uniqueSequence[$i]).".rtf";
  $rtfFile = "$rtfDir/1".lc($uniqueSequence[$i]).".rtf";
  $rtfFileNoOne = "$rtfDir/".lc($uniqueSequence[$i]).".rtf";
  open (RTF, "$rtfFile") or open(RTF, "$rtfFileNoOne") or die("can't open rtf file $rtfFile or $rtfFileNoOne".$!);
  while (<RTF>) {
    if ($_ =~ /^BOND/) {
      @tmpA = split (/\s+/, $_);
      for (my $j = 1; $j < scalar(@tmpA); $j+=2) {
	push (@{$rtfHash{$tmpA[$j]}}, $tmpA[$j+1]);
	push (@{$rtfHash{$tmpA[$j+1]}}, $tmpA[$j]);
      }
    }
  }
  $masterRtfHash{$uniqueSequence[$i]} = \%rtfHash;
  close (RTF);
}


# HERE SPENCER(Tue Jan 10 23:32:07 2006)

# You need do dereference the rtfhash until you get to the array for a
# particular atom.  Once you find that, then print the CONECT list.

#print Dumper($masterPdbHash{PS132}->{'residue'}); exit;


my @keysArray = sort { $masterPdbHash{$a}->{'resNum'} <=> $masterPdbHash{$b}->{'resNum'} } (keys %masterPdbHash);

#print Dumper(@keysArray); exit;

open (PDB, ">>$pdbFile") or die($!);
#foreach my $keys (keys %masterPdbHash) {
foreach my $keys (@keysArray) {
  if (exists($masterPdbHash{$keys}->{'residue'})) {
    print PDB "REMARK ".$keys."\n";
    foreach my $atom (keys %{$masterPdbHash{$keys}->{'atomHash2'}}) {
      for (my $i= 0; $i < scalar(@{$masterRtfHash{$masterPdbHash{$keys}->{'residue'}}->{$atom}}); $i++) {
	my $partnerAtom = $masterRtfHash{$masterPdbHash{$keys}->{'residue'}}->{$atom}->[$i];
	if ($masterPdbHash{$keys}->{'atomHash2'}->{$atom} < $masterPdbHash{$keys}->{'atomHash2'}->{$partnerAtom}) {
# 	  printf PDB ("%5i", $masterPdbHash{$keys}->{'atomHash2'}->{$atom});
# 	  printf PDB ("%5i", $masterPdbHash{$keys}->{'atomHash2'}->{$partnerAtom});

	  print PDB "CONECT";
	  printf PDB ("%5i", $masterPdbHash{$keys}->{'atomHash2'}->{$atom});
	  printf PDB ("%5i\n", $masterPdbHash{$keys}->{'atomHash2'}->{$partnerAtom});
	} 
      }
    }
  }
}
close(PDB);




#print Dumper(%masterPdbHash);


