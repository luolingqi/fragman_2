package PDB;
use Data::Dumper;
# vars

%aa1to3 = (
	   'G' => 'GLY', 'A' => 'ALA', 'V' => 'VAL', 'L' => 'LEU', 'I' => 'ILE',
	   'S' => 'SER', 'C' => 'CYS', 'T' => 'THR', 'M' => 'MET', 'P' => 'PRO',
	   'F' => 'PHE', 'Y' => 'TYR', 'H' => 'HIS', 'K' => 'LYS', 'R' => 'ARG',
	   'D' => 'ASP', 'E' => 'GLU', 'N' => 'ASN', 'Q' => 'GLN', 'W' => 'TRP'
	  );
%aa3to1 = reverse %aa1to3;

# routines for doing pdb stuff

sub new { 
  my $class = shift; return bless {}, $class; 
}


sub extractResidueOneFile {
  # vars
  my $pdb = $_[0];
  my $res = $_[1];
  my $chain = $_[2]; # optional ??
  my $outname = $_[3];

  $pdb =~ m#([\w\.]*)\.pdb#i;

  my $pdbID = $1;

  $res = uc($res);
  $res =~ /([A-Z0-9]{3})([0-9]+)/;
  my $resName = $1;
  my $resNum = $2;

  if (length($resName) == 1) {
    $resName = $aa1to3{$resName};
  }

  my $found = 0;
  my $pwd = `pwd`;
  chomp ($pwd);
  # end vars

  # work
  open(PDB, $pdb) or die($!.": $pdb");
#  my $outname = "$pwd/$pdbID.".$aa3to1{$resName}."$resNum".".pdb";
  open(OUT, ">>$outname") or die($!.": >$outname");
  while (<PDB>) {
    if ($_ =~ /^ATOM/) {
      print getResNum($_)."\n";
      if ((getResNum($_) == $resNum) && (getResName($_) eq $resName)) {
	$found = 1;
	if ($chain) {
	  if (getChain($_) eq $chain) {
	    print OUT $_;
	    next;
	  } # else continue below
	} else {
	  print OUT $_;
	}
      } elsif ($found == 1 && getResNum($_) != $resNum) {
	print "closing\n";
	last;
      }
    }
  }
  close(PDB);
  close(OUT);
  if ($found == 1) {
    # return the path to the new file
    return $outname;
  } else {
    system("rm $outname");
    return 0;
  }
}

sub extractResiduesOneFile {
  # vars
  my $pdb = $_[0];
  my $residues = $_[1];
  my $chain = $_[2]; # optional ??
  my $outname = $_[3];

  $pdb =~ m#([\w\.]*)\.pdb#i;

  my $pdbID = $1;


#   if (length($resName) == 1) {
#     $resName = $aa1to3{$resName};
#   }

  my $found = 0;
  my $pwd = `pwd`;
  chomp ($pwd);
  my ($res, $resName, $resNum);
  my $counter = 0;

  my %tracker;

  for (my $i = 0; $i < scalar(@residues); $i++) {
    $res = uc(@{$residues}[$i]);
    $res =~ /([A-Z0-9]{3})([0-9]+)/;
    $resName = $1;
    $resNum = $2;
    $tracker{$resName.$resNum} = 0;
  }

  # end vars


  my $currentRes = 0;
  my $lastRes = 0;

  open(PDB, $pdb) or die($!.": $pdb");
#  my $outname = "$pwd/$pdbID.".$aa3to1{$resName}."$resNum".".pdb";
  open(OUT, ">>$outname") or die($!.": >$outname");
  while (<PDB>) {
    if ($_ =~ /^ATOM/) {
      $currentRes = getResNum($_);
      if ($lastRes) {
	if ($lastRes != $currentRes && $counter == scalar(@{$residues})) {
	  last;
	}
      }
      for (my $i = 0; $i < scalar(@{$residues}); $i++) {
	$res = uc(@{$residues}[$i]);
	$res =~ /([A-Z0-9]{3})([0-9]+)/;
	$resName = $1;
	$resNum = $2;
	if ((getResNum($_) == $resNum) && (getResName($_) eq $resName)) {
	  if ($currentRes != $lastRes) {
	    $tracker{$resName.$resNum} = 1;
	    $counter++;
	    #print $counter."\n";
	  }
	  $found = 1;
	  if ($chain) {
	    if (getChain($_) eq $chain) {
	      print OUT $_;
	      next;
	    }			# else continue below
	  } else {
	    print OUT $_;
	  }
	}
      }
      $lastRes = $currentRes;
    }
  }
  close(PDB);
  close(OUT);
  if ($found == 1) {
    # return the path to the new file
    return $outname;
  } else {
    system("rm $outname");
    return 0;
  }
}


sub extractResidue {
  # vars
  my $pdb = $_[0];
  my $res = $_[1];
  my $chain = $_[2]; # optional ??

  $pdb =~ m#([\w\.]*)\.pdb#i;

  my $pdbID = $1;

  $res = uc($res);
  $res =~ /([A-Z0-9]{3})([0-9]+)/;
  my $resName = $1;
  my $resNum = $2;

  if (length($resName) == 1) {
    $resName = $aa1to3{$resName};
  }

  my $found = 0;
  my $pwd = `pwd`;
  chomp ($pwd);
  # end vars

  # work
  open(PDB, $pdb) or die($!.": $pdb");
  my $outname = "$pwd/$pdbID.".$aa3to1{$resName}."$resNum".".pdb";
  open(OUT, ">$outname") or die($!.": >$outname");
  while (<PDB>) {
    if ($_ =~ /^ATOM/) {
      if ((getResNum($_) == $resNum) && (getResName($_) eq $resName)) {
	$found = 1;
	if ($chain) {
	  if (getChain($_) eq $chain) {
	    print OUT $_;
	    next;
	  } # else continue below
	} else {
	  print OUT $_;
	}
      }
    }
  }
  close(PDB);
  close(OUT);
  if ($found == 1) {
    # return the path to the new file
    return $outname;
  } else {
    system("rm $outname");
    return 0;
  }
}

sub listChains  {
  # Sub takes a pdb as an argument and returns an array of the chains.
  my @chains;
  my $pdb = $_[0];
  open(FH, $pdb) or die($!);
  my $lastChain = '';
  while (<FH>) {
    if ($_ =~ /^ATOM/) {
      my $currentChain = substr($_, 21, 1);
      if ($lastChain ne $currentChain && $currentChain =~ /[A-Z]/) {
	push(@chains, $currentChain);
      }
      $lastChain = $currentChain;
    }
  }
  close(FH);
  return(@chains);
}


sub extractChain {
  my $pdb = $_[0];
  my $chain = $_[1];
  $pdb =~ m#([\w\.]*)\.pdb#i;
  my $output = "$1.$chain.pdb";
  open (OUT, ">$output") or die($!);
  open (PDB, $pdb) or die("Couldn't open the pdb $pdb: $!\n");
  while (<PDB>) {
    if (substr($_, 0, 7) =~ /^ATOM/) {
      my $currentChain = substr($_, 21, 1);
      if ($currentChain eq $chain) {
	print OUT $_;
      }
    }
  }
  print OUT "TER";
  close (PDB);
  close(OUT);
}





sub mutateResidue {
  # vars
  my $pdb = $_[0];
  my $res = $_[1];
  my $chain = $_[2]; # optional ??

  $pdb =~ m#([\w\.]*)\.pdb#i;

  my $pdbID = $1;

  $res = uc($res);
  $res =~ /([A-Z0-9]{3})([0-9]+)/;
  my $resName = $1;
  my $resNum = $2;

  if (length($resName) == 1) {
    $resName = $aa1to3{$resName};
  }

  my $found = 0;
  my $pwd = `pwd`;
  chomp ($pwd);
  # end vars

  # work
  open(PDB, $pdb) or die($!.": $pdb");
  my $outname = "$pwd/$pdbID.".$aa3to1{$resName}."$resNum"."A.pdb";
  open(OUT, ">$outname") or die($!.": >$outname");
  while (<PDB>) {
    if ($_ =~ /^ATOM/) {
      if ((getResNum($_) == $resNum) && (getResName($_) eq $resName)) {
	if ($chain) {
	  if (getChain($_) ne $chain) {
	    print OUT $_;
	    next;
	  } # else continue below
	}
	$found = 1;
	my $atom = substr($_, 12, 4);
	$atom =~ s/\s//g;
	if ($atom eq 'N' || $atom eq 'H' || $atom eq 'CA' || $atom eq 'CB' || $atom eq 'C' || $atom eq 'O') {
	  substr($_, 17, 3) = "ALA";
	  print OUT $_;
	}
      } else {
	print OUT $_;
      }
    } else {
      print OUT $_;
    }
  }
  close(PDB);
  close(OUT);
  if ($found == 1) {
    # return the path to the new file
    return $outname;
  } else {
    system("rm $outname");
    return 0;
  }
}

sub getSequenceThreeLetter {
  #=====================================================================#
  # Takes a pdb and gives you the sequence of amino acids as one big    #
  # chain
  #=====================================================================#
  my $pdb = $_[0];
  my $convert = $_[1];
  open(PDB, $pdb) or die($!);
  my $currResNum;
  my $lastResNum = 0;
  my $sequence = "";			# return value
  while (<PDB>) {
    if ($_ =~ /^ATOM/) {
      $currResNum = substr($_, 22, 4);
      $currResNum =~ s/\s+//g;
      if ($currResNum != $lastResNum) {
	my $currResName = substr($_, 17, 3);
	$currResName =~ s/\s+//g;
	#print $currResName."\n";
	if ($convert == 1) {
	  $sequence .= $aa3to1{$currResName};
	} else {
	  $sequence .= $currResName." ";
	}
      }
      $lastResNum = $currResNum;
    }
  }
  close(PDB);
  return $sequence;
}

sub getSequence {
  #=====================================================================#
  # Takes a pdb and gives you the sequence of amino acids as one big    #
  # chain
  #=====================================================================#
  my $pdb = $_[0];
  my $convert = $_[1];
  open(PDB, $pdb) or die($!);
  my $currResNum;
  my $lastResNum = 0;
  my $sequence = "";			# return value
  while (<PDB>) {
    if ($_ =~ /^ATOM/) {
      $currResNum = substr($_, 22, 4);
      $currResNum =~ s/\s+//g;
      if ($currResNum != $lastResNum) {
	my $currResName = substr($_, 17, 3);
	$currResName =~ s/\s+//g;
	#print $currResName."\n";
	if ($convert == 1) {
	  $sequence .= $aa3to1{$currResName};
	} else {
	  $sequence .= $currResName;
	}
      }
      $lastResNum = $currResNum;
    }
  }
  close(PDB);
  return $sequence;
}


sub getDistance {
  # Returns all the pairwise distances between two residues in two
  # pdbs, presumably ligand and receptor, but could be used to compare
  # two identical proteins.
  my $pdb1 = $_[0];
  my $res1 = $_[1];	   # should be glu32 or e32 (case insensitive)
  my $pdb2 = $_[2];
  my $res2 = $_[3];	   # should be glu32 or e32 (case insensitive)
  my ($coords1, $coords2);
  if ($_[4] && $_[5]) {
    $coords1 = $_[4];
    $coords2 = $_[5];
  } else {
    $coords1 = getCoordinates($pdb1, $res1);
    $coords2 = getCoordinates($pdb2, $res2);
  }


  # Hashes are the return value and store the pairwise distances as
  # {OE} = { 'C' => 1.2, 'HH' => 1.0 ...}
  my (%res1d, %res2d);
  if (!$coords1 || !$coords2) {
    return 0;
  }

  foreach my $atom1 (keys %{$coords1}) {
    foreach my $atom2 (keys %{$coords2}) {
      my $distance = _calcDistance($$coords1{$atom1}->{'x'}, $$coords2{$atom2}->{'x'}, $$coords1{$atom1}->{'y'}, $$coords2{$atom2}->{'y'},$$coords1{$atom1}->{'z'},$$coords2{$atom2}->{'z'});
      $res1d{$atom1}{$atom2} = $distance;
      $res2d{$atom2}{$atom1} = $distance;
    }
  }
  return (\%res1d, \%res2d);
}


sub getAllCoordinates {
  # Returns all the coordinates for all residues in a pdb as a hash.

  my $pdb1 = $_[0];
  # Hashes are the return value and store the pairwise distances as
  # $coordinates{OE} = { 'C' => 1.2, 'HH' => 1.0 ...}
  my %coordinates;
  open(PDB1, $pdb1) or die("$!: $pdb1");
  while (<PDB1>) {
    if ($_ =~ /^ATOM/) {
      my $chain = getChain($_);
      my $currentResName = substr($_, 17, 3);
      my $currentResNum = substr($_, 22, 4);
      $currentResNum =~ s/\s+//g;
      my $x = substr($_, 30, 8);
      $x =~ s/\s+//g;
      my $y = substr($_, 38, 8);
      $y =~ s/\s+//g;
      my $z = substr($_, 46, 8);
      $z =~ s/\s+//g;
      my $atom = substr($_, 12, 4);
      $atom =~ s/\s+//g;
      $coordinates{$chain}{$currentResName.$currentResNum}{$atom} = { 'x' => $x, 'y' => $y, 'z' => $z };
    }
  }
  close(PDB1);

  return (\%coordinates);
}

sub getCoordinates {
  # Returns all the coordinates for a single residue in a pdb as a hash.

  my $pdb1 = $_[0];
  my $res1 = $_[1];	   # should be glu32 or e32 (case insensitive)
  my ($resName1, $resNum1, $found);
  # Hashes are the return value and store the pairwise distances as
  # $coordinates{OE} = { 'C' => 1.2, 'HH' => 1.0 ...}
  my %coordinates;

  $res1 = uc($res1);
  $res1 =~ /([A-Z0-9]{3})([0-9]+)/;
  $resName1 = $1;
  $resNum1 = $2;
  if (length($resName1) == 1) {
    $resName1 = $aa1to3{$resName1};
  }

  open(PDB1, $pdb1) or die("$!: $pdb1");
  while (<PDB1>) {
    if ($_ =~ /^ATOM/) {
      my $currentResName = substr($_, 17, 3);
      my $currentResNum = substr($_, 22, 4);
      $currentResNum =~ s/\s+//g;
      if ($currentResNum == $resNum1 && $currentResName eq $resName1) {
	$found = 1;
	my $x = substr($_, 30, 8);
	$x =~ s/\s+//g;
	my $y = substr($_, 38, 8);
	$y =~ s/\s+//g;
	my $z = substr($_, 46, 8);
	$z =~ s/\s+//g;
	my $atom = substr($_, 12, 4);
	$atom =~ s/\s+//g;
	$coordinates{$atom} = { 'x' => $x, 'y' => $y, 'z' => $z };
      }
    }
  }
  close(PDB1);

  if (!$found) {
    die ("Cannot get the coordinates for residue $res1 in $pdb1 because it is inavlid.\nPlease check your arguments and try again.\n");
  } else {
    return (\%coordinates);
  }
}


sub runSAS {
  # Runs the sas program and returns the sas for a PDB as a hash.
  # calculates SAS with the newsur program
  my $pdb = $_[0];
  my $mapHashRef = _mapResidueNumbers($pdb);
  system("cp $pdb tmp.pdb");
  system("/usr/local/bin/newsur_oct.x tmp.pdb .tmp.sas 2>&1 > .tmp");
  system("rm fort*");
  system("rm .tmp");
  open(FH, ".tmp.sas") or die($!." .tmp.sas");
  my @sasFile = <FH>;
  for (my $i = 0; $i < scalar(@sasFile); $i++) {
    chomp($sasFile[$i]);
    $sasFile[$i] =~ s/^\s//;
    if ($sasFile[$i] =~ /^[A-Z]/) {
      my ($resName, $sc, $aa, $foo, $index) = split(/\s+/, $sasFile[$i]);
      $$mapHashRef{$index}->{'scSAS'} = $sc;
      $$mapHashRef{$index}->{'aaSAS'} = $aa;
    }
  }
  close(FH);
  return (\$mapHashRef);
}

sub _mapResidueNumbers {
  # Will create a file called protein_name.map which will give you a
  # list of all the residues as they appear in the pdb, so that
  # programs that don't use the residue index numbers listed in the
  # pdb can be mapped to these values.
  my $pdb = $_[0];
  my (%mapHash, $currentResidue, $currentResidueNo, $lastResidueNo, $newResidueNo);
  open(FH, "$pdb") or die($!.": $pdb");
  $pdb =~ /(.*)\.pdb/;
  $newResidueNo = 1;
  $lastResidueNo = 0;
  while (<FH>) {
    if (substr($_, 0, 6) =~ /ATOM/) {
      $chain = getChain($_);
      $currentResidue = substr($_, 17, 3);
      $currentResidueNo = substr($_, 22, 5);
      $currentResidueNo =~ s/\s*//g;
      if ($lastResidueNo ne $currentResidueNo) {
	$mapHash{$newResidueNo} = { 'resNum' => $currentResidueNo, 'resName' => $currentResidue, 'chain' => $chain};
	$newResidueNo++;
      }
      $lastResidueNo = $currentResidueNo;
    }
  }
  close(FH);
  return \%mapHash;
}


sub _calcDistance {
  # returns the 3d distance between 6 sets of numbers: x1, x2, y1, y2,
  # z1, z2

  my ($x1, $x2, $y1, $y2, $z1, $z2) = @_;
  my $distance = sqrt(($x2 - $x1)**2 + ($y2 - $y1)**2 + ($z2 - $z1)**2);
  return sprintf("%.4f", $distance);
}


sub bindingFreeEnergy {
  # Runs the fast.x program to calculate the binding free energy for
  # two pdbs.  Will return the location of the output files for the
  # first pdb and the second pdb in that order.
  my ($pdb1, $pdb2) = @_;
  my $pwd = `pwd`;  chomp ($pwd);
  my $tmp1 = ".lig.tmp";
  my $tmp2 = ".rec.tmp";
  $pdb1 =~ m#([\w\.]*)\.pdb#i;
  my $out1 = "$pwd/$1.energy";
  $pdb2 =~ m#([\w\.]*)\.pdb#i;
  my $out2 = "$pwd/$1.energy";

  # remove the OCT1, HT1, etc...
  open(TMPL, ">$tmp1") or die($!);
  open(LPDB, "<$pdb1") or die($!.": $pdb1");
  while (<LPDB>) {
    # substitute OCT1 for O
    # sub HT1 for H
    # delete HT2, HT3, OCT2
    if ($_ =~ /^ATOM/ && substr($_, 12, 4) =~ /^\sHT|^OCT/) {
      if (substr($_, 12, 4) =~ /^\sHT1/) {
	$_ =~ s/HT1/H  /;
	print TMPL $_;
      } elsif (substr($_, 12, 4) =~ /^OCT1/) {
	$_ =~ s/OCT1/ O  /;
	print TMPL $_;
      }
    } else {
      print TMPL $_;
    }
  }
  close(LPDB);
  close(TMPL);

  open(RPDB, "<$pdb2") or die($!);
  open(TMPR, ">$tmp2") or die($!);
  while (<RPDB>) {
    # substitute OCT1 for O
    # sub HT1 for H
    # delete HT2, HT3, OCT2
    if ($_ =~ /^ATOM/ && substr($_, 12, 4) =~ /^\sHT|^OCT/) {
      if (substr($_, 12, 4) =~ /^\sHT1/) {
	$_ =~ s/HT1/H  /;
	print TMPR $_;
      } elsif (substr($_, 12, 4) =~ /^OCT1/) {
	$_ =~ s/OCT1/ O  /;
	print TMPR $_;
      }
    } else {
      print TMPR $_;
    }
  }
  close(RPDB);
  close(TMPR);

  

  #   my $command = "/usr/local/bin/fast.x $pdb2 $pdb1 1 -5 -5 -5 0.0 > $out1 2>&1";
  #   my $command2 = "/usr/local/bin/fast.x $pdb1 $pdb2 1 -5 -5 -5 0.0 > $out2 2>&1";

  my $command = "/usr/local/bin/fast.x $tmp2 $tmp1 1 -5 -5 -5 0.0 > $out1 2>&1";
  my $command2 = "/usr/local/bin/fast.x $tmp1 $tmp2 1 -5 -5 -5 0.0 > $out2 2>&1";

  my $pid = fork;
  if (not defined $pid) {
    die "couldn't fork: $!";
  }
  if ($pid) {
    #    sleep 30;
    waitpid($pid, 0);	 # this seems to be working -- doublecheck....
    kill 'INT', $pid;
  } else {
    exec("$command; $command2");
  }
  return ($out1, $out2);
}

sub bindingFreeEnergy2 {
  # This version needs to keep the c and n termini intact and instead
  # give the length of each successive chain.  These will be entered
  # as the arguments for fast.x

  # Runs the fast.x program to calculate the binding free energy for
  # two pdbs.  Will return the location of the output files for the
  # first pdb and the second pdb in that order.
  my ($pdb1, $pdb2) = @_;
  my $pwd = `pwd`;  chomp ($pwd);
  $pdb1 =~ m#([\w\.]*)\.pdb#i;
  my $out1 = "$pwd/$1.energy";
  $pdb2 =~ m#([\w\.]*)\.pdb#i;

  # create an index of the residue numbers, and where each chain ends
  my $lastResNo = 0;
  my $numResInCurrentChain = 0;
  my $totalRes;
  open(LPDB, "<$pdb1") or die($!.": $pdb1");
  my (@numResInLig, @numResInRec);
  while (<LPDB>) {
    my $currentRes = substr($_, 17, 3);
    my $currentResNo = substr($_, 22, 5);

    if ($_ =~ /^TER/) {
      push (@numResInLig, $totalRes);
      $numResInCurrentChain = 0;
      $lastResNo = 0;
      next;
    } else {
      $currentResNo =~ s/\s//g;
      if ($lastResNo != $currentResNo) {
	$numResInCurrentChain++;
	$totalRes++;
      }
      $lastResNo = $currentResNo;
    }
  }
  close(LPDB);

  $lastResNo = 0;
  $numResInCurrentChain = 0;
  $totalRes = 0;

  open(RPDB, "<$pdb2") or die($!);
  while (<RPDB>) {
    my $currentRes = substr($_, 17, 3);
    my $currentResNo = substr($_, 22, 5);
    if ($_ =~ /^TER/) {
      push (@numResInRec, $totalRes);
      $numResInCurrentChain = 0;
      $lastResNo = 0;
      next;
    } else {
      $currentResNo =~ s/\s//g;
      if ($lastResNo != $currentResNo) {
	$numResInCurrentChain++;
	$totalRes++;
      }
      $lastResNo = $currentResNo;
    }
  }
  close(RPDB);

  if (scalar(@numResInLig) > 3 || scalar(@numResInRec) > 3) {
    die ("Won't work because number of chains in ligand or receptor too high");
  }

  my $command;

  if (scalar(@numResInLig) > 1) {
    if (!$numResInLig[1]) {
      $numResInLig[1] = -5;
    }
    $command = "/usr/bin/nohup nice /usr/local/bin/fast.x $pdb1 $pdb2 1 -5 $numResInLig[0] $numResInLig[1] 0 > $out1 2>&1";
  } elsif (scalar(@numResInRec) > 1) {
    if (!$numResInRec[1]) {
      $numResInRec[1] = -5;
    }
    $command = "/usr/bin/nohup nice /usr/local/bin/fast.x $pdb2 $pdb1 1 -5 $numResInRec[0] $numResInRec[1] 0 > $out1 2>&1";
  } else {
    # default for one chain in each ligand and receptor.
    $command = "/usr/bin/nohup nice /usr/local/bin/fast.x $pdb1 $pdb2 1 -5 -5 -5 0 > $out1 2>&1";
  }
  my $pid = fork;
  if (not defined $pid) {
    die "couldn't fork: $!";
  }
  if ($pid) {
    #    sleep 30;
    waitpid($pid, 0);	 # this seems to be working -- doublecheck....
    kill 'INT', $pid;
  } else {
    exec("$command");
  }

  # do a little cleanup

  my $tmp = 'tmp.energy';
  system("rm fort.*");
  return ($out1);
}

sub energyDifference {
  # Takes two *.energy files and returns the difference (delta G)
  my $WT = $_[0];
  my $MUT = $_[1];
  open(WT, $WT) or die($!.": $WT");
  my @wt = <WT>;
  close(WT);

  open(MUT, $MUT) or die($!."$MUT");
  my @mut = <MUT>;
  close(MUT);

  $wt[scalar(@wt) - 1] =~ s/^\s+//g;
  $wt[scalar(@wt) - 1] =~ s/\s+/,/g;
  my @wtG = split(/,/, $wt[scalar(@wt) - 1]);
  $mut[scalar(@mut) - 1] =~ s/^\s+//g;
  $mut[scalar(@mut) - 1] =~ s/\s+/,/g;
  my @mutG = split(/,/, $mut[scalar(@mut) - 1]);

  my $wtSum = $wtG[0] + $wtG[1];
  my $mutSum = $mutG[0] + $mutG[1];
  my $delta = ($mutSum - $wtSum);
  $delta = sprintf("%.1f", $delta);
  return $delta;
}

sub runCharmm {
  my $pdb = $_[0];
  my $pid = fork;
  if (not defined $pid) {
    die "couldn't fork: $!";
  }
  if ($pid) {
    waitpid($pid, 0);	 # this seems to be working -- doublecheck....
    kill 'INT', $pid;
  } else {
    exec("/usr/bin/nohup nice -10 /usr/local/bin/pdb.addHtoPdb.pl $pdb > .out");
  }

  open (FH, ".out");
  my @return = <FH>;
  close(FH);
  return $return[0];
}


# below parse the current line and return information like residue
# name, chain, x, y, z, etc...

sub getResNum {
  # returns the residue number of current atom
  my $tmp = substr($_[0], 22, 4);
  $tmp =~ s/\s+//g;
  return $tmp;
}

sub getResName {
  return substr($_[0], 17, 3);
}

sub getX {
  return substr($_[0], 30, 8);
}

sub getY {
  return substr($_[0], 38, 8);
}

sub getZ {
  return substr($_[0], 46, 8);
}

sub getAtomName {
  my $tmp = substr($_[0], 12, 4);
  $tmp =~ s/\s+//g;
  return $tmp;
}

sub getAtomNum {
  my $tmp = substr($_[0], 6, 5);
  $tmp =~ s/\s+//g;
  return $tmp;
}

sub getChain {
  return substr($_[0], 21, 1);
}

sub getOccupancy {
  my $tmp = substr($_[0], 54, 6);
  $tmp =~ s/\s+//g;
  return $tmp;
}



return 1;
