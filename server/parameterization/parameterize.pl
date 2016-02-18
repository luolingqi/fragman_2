#!/usr/bin/perl -w
use strict;

#Scott Mottarella
#12/07/11

#Perl script for running parameterization server
#Takes one file from user with one valid SMILES per line, current limit 10 total lines

#REQUIRED SCRIPTS/PROGRAMS MUST BE IN YOUR PATH!
#1. rdkit_3d_from_smi.py
#2. babel (or obabel, but there are conflicts with confab)
#3. confab (http://code.google.com/p/confab/)
#4. count_mol_pdb.py
#5. extract_pdb.py
#6. probe_prep.py
#7. run_multiple.pl *Requires the rest of the parameterization server
#8. addchain.pl
#9. msur

#Export variables
#$ENV{RDBASE} = "/home/tanggis/src/rdkit";
#$ENV{LD_LIBRARY_PATH} = "/home/tanggis/src/rdkit/lib";
#$ENV{PYTHONPATH} = "/home/tanggis/src/rdkit";
$ENV{PATH} = "$ENV{PATH}:/home/ftmap/bin";
my $param_root="/home/ftmap/parameterization";
$ENV{PATH} = "$param_root:$ENV{PATH}";
$ENV{PATH} = "$param_root/antechamber-1.27/exe:$ENV{PATH}";
#$ENV{PATH} = "$param_root/openbabel-2.2.3/bin:$ENV{PATH}";
$ENV{PATH} = "$param_root/sarnoff:$ENV{PATH}";
$ENV{ACHOME} = "$param_root/antechamber-1.27";

#Allow for command line options (may be none, but allows for future options to be added)
use Getopt::Long;
my $filename = "parameters";

GetOptions("name=s" => \$filename);

#Check for one input and that it is a file that exists
#It's on the user to check that the file is of the correct format
if (@ARGV!=1){
	print "Usage: parameterize.pl file\n";
	exit 0;
}else{
	unless (-e $ARGV[0]){
		print "Usage: parameterize.pl file\n";
		exit 0;
	}
}

#Read in input SMILES and check against limit for maximum SMILES
open(FILE, $ARGV[0]) or die "Could not open $ARGV[0]: $!";
my @tmp = <FILE>;
close FILE;
unless (@tmp <= 10){
	print "Current limit: 10 SMILES per run\n";
	exit 0;
}
my (@smi, @charge);
foreach (@tmp){
	if ($_ =~ /^([+|-]?\d+)[\ \t]*(.*)$/){
		push @charge,$1;
		push @smi,$2;
	}else{
		print "\"$_\" does not match format \"charge SMILES\". Continuing...\n";
		next;
	}
}
chomp @smi;

#Create a new working directory and change to it
if (-e "$filename"){
	print "Could not create directory $filename. File already exists!\n";
	exit 0;
}else{
	system "mkdir $filename";
}
chdir "$filename" or die "Could not change to directory $filename: $1";
mkdir "probes";
mkdir "probes/rec";

#my $create = 0;
#my $filename = 0;
#while (!$create){
#	if (-e "$name$filename"){
#		$filename += 1;
#	}else{
#		system "mkdir $name$filename";
#		$create = 1;
#	}
#}
#chdir "$name$filename" or die "Could not change to directory $name$filename: $!";


#Iterate over the SMILES, generate conformers and parameterize
for (my $i = 1; $i <= @smi; $i++){
	my $name = sprintf("%02d", $i);
	chomp $smi[$i-1];
	`echo \"BEGIN SMILES $name $smi[$i-1]\" >> ../$filename.log`;
	mkdir "$name";
	chdir $name;
	system "rdkit_3d_from_smi.py \"$smi[$i-1]\" > $name.sdf ";
	#system "echo \"$smi[$i-1]\" > $name.smi";
	#system "babel --gen3d -h $name.smi $name.sdf";
	if (`wc -l $name.sdf | awk '{print \$1}'`==0){
		`echo \"SMILES $name $smi[$i-1] is not a valid SMILES. Proceeding with next SMILES...\" >> ../probes/$filename.out`;
		`echo \"SMILES $name $smi[$i-1] is not a valid SMILES. Proceeding with next SMILES...\" >> ../../$filename.log`;
		chdir "..";
		next;
	}
	system "confab -a $name.sdf $name.confab.pdb >>../../$filename.log";
	my $confs = `count_mol_pdb.py $name.confab.pdb`;
	my $start = 1;
	if ($confs > 99){
		chomp $confs;
		`echo \"Current limit on number of conformers generated is 99.  SMILES $name $smi[$i-1] generated $confs and will be skipped.\" >> ../probes/$filename.out`;
		`echo \"Current limit on number of conformers generated is 99.  SMILES $name $smi[$i-1] generated $confs and will be skipped.\" >> ../../$filename.log`;
		chdir "..";
		next;
	}elsif ($confs > 1){
		$start = 2;
	}
	chomp $confs;
	my $success = 0;
	for (my $j = $start; $j <= $confs; $j++){
		#Make the count start at number 1 instead of 2 for consistent file naming
		my $number = $j;
		if ($start == 2){
			$number = $j-1;
		}
		my $conf = sprintf("%02d", $number);
		mkdir "$name$conf";
		chdir "$name$conf";
		system "extract_pdb.py ../$name.confab.pdb $number > $name.confab.$number.pdb";
		system "probe_prep.py -l ".chr($i+96)." $name.confab.$number.pdb $number";
		system "echo 1 > count1";
		system "run_multiple.pl pdb bcc pre $charge[$i-1] ".chr($i+96).sprintf("%02d",$number).".pdb 2>&1 1>>../../../$filename.log";
		system "cat 0001/".chr($i+96).sprintf("%02d",$number)."-bcc_nmin.pdb | sed \"s/ X0 / X  /\" > 1lig.pdb";
		system "cat 0001/".chr($i+96).sprintf("%02d",$number)."-bcc_nmin.psf | sed \"s/ X0 / X  /\" > 1lig.psf";
		system "cat 0001/".chr($i+96).sprintf("%02d",$number)."-bcc_nmin_xplor.psf | sed \"s/ X0 / X  /\" > 1lig_xplor.psf";
		system "addchain.pl 1lig.pdb";
		system "head -n 1025 /home/ftmap/prms/atoms.0.0.4.prm.ms.3cap+0.5ace_hphobe3.Hr0rec >> atoms_msur.ats";
		system "grep -v \"^\$\" 0001/files_components/atoms_subatom | awk '{print \$1,\"\t\",\$2,\"\t\",\$3,\"\t\",\"-1\",\"\t\",\$5,\"\t\",\$6}' >> atoms_msur.ats";
		system "tail -n 8 /home/ftmap/prms/atoms.0.0.4.prm.ms.3cap+0.5ace_hphobe3.Hr0rec >> atoms_msur.ats";
		system "msur -p atoms_msur.ats 1lig.pdb 1lig.ms 2>&1 1>>../../../$filename.log";
		if (`wc -l 1lig.pdb | awk '{print \$1}'`==0){
			system "echo \"Parameterization of SMILES $name conformer $conf ($smi[$i-1]) did not succeed.\" >> ../../../$filename.log";
			chdir "..";
			next;
		}
		$success = 1;
		mkdir "../../probes/$name$conf";
		system "cp 1lig.pdb 1lig.psf 1lig_xplor.psf 1lig.ms 0001/atoms_new.ats 0001/pdbamino_new.rtf 0001/parm_new.prm ../../probes/$name$conf";
		system "ln -s ../rec/1rec.1sidehp3.ms ../../probes/$name$conf";
		system "ln -s ../rec/1rec_for_min.pdb ../../probes/$name$conf";
		system "ln -s ../rec/1rec_for_min_xplor.psf ../../probes/$name$conf";
		system "ln -s ../rec/1rec.charmm.pdb ../../probes/$name$conf";
		system "ln -s ../rec/1rec.pdb ../../probes/$name$conf";
		system "ln -s ../rec/1rec.psf ../../probes/$name$conf";
		system "ln -s ../rec/1rec_xplor.psf ../../probes/$name$conf";
		system "ln -s ../rec/mask.pdb ../../probes/$name$conf";
		system "ln -s ../rec/mask.ms ../../probes/$name$conf";
		system "ln -s ../rec/pb.dx ../../probes/$name$conf";
		system "ln -s ../rec/pb.pot ../../probes/$name$conf";
		system "echo $name$conf >> ../../probes/probes";
		chdir "..";
	}
	if ($success){
		if ($start == 2){
			$confs--;
		}
		system "echo \"Parameterization of SMILES $name $smi[$i-1] COMPLETE with $confs conformations.\" >> ../probes/$filename.out";
	}else{
		system "echo \"Parameterization of SMILES $name $smi[$i-1] FAILED.\" >> ../probes/$filename.out";
	}
	chdir "..";
}
