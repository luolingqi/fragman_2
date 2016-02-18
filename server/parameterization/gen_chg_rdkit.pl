#!/usr/bin/perl -w
use strict;use Getopt::Long;
use v5.18;

my ($filename, $chartype, $net_char, @pdblines, @inp_lines, @gau_lines, $k, $hydrogen, $name, $method, $m, @charge);
my ($ism_line, $count_pos, $count_neg, $count_net);

if (scalar @ARGV < 7){short_args();exit;}							# SERVER MOD

$chartype = '';
$hydrogen = '';
$count_net= '';											# SERVER MOD
GetOptions ('chartype=s'=>\$chartype, 'hydrogen=s'=>\$hydrogen, 'count_net=s'=>\$count_net);
$filename = shift @ARGV;

# system "babel -imol2 $filename.mol2 -osmi $filename.smi";

$net_char = $count_net;

system "antechamber -i $filename.mol2 -fi mol2 -o $filename.gau -fo gcrt";			# produce gaussian input file

open (GAU, "$filename.gau") || die "couldn't open $filename.gau: $!\n";
chomp(@gau_lines = <GAU>); close GAU;

system "mv $filename.gau $filename.gau_orig";

open (OUT, ">$filename.gau") || die "cannot write $filename.gau: $!\n";

for ($k=0; $k<@gau_lines; $k++)
{
	if ($k == 6) {print OUT "$net_char 1\n";}
	else {print OUT "$gau_lines[$k]\n";}
}
close OUT;

given($chartype)
{
	when ("rsp") {$m = "r"; $method = "resp";}

	when ("bcc") {$m = "a"; $method = "bcc";}
}

system "gjftogms.sh <$filename.gau> $filename.inp -$m -d";

open (INP, "$filename.inp") || die "couldn't open file $filename.inp: $!\n";
chomp(@inp_lines = <INP>); close INP;

system "cp $filename.inp $filename.inp_orig";

if ($method eq "resp" or $method eq "bcc")							# SERVER MOD
{
	open (OUT, ">$filename.inp") || die "cannot write file $filename.inp: $!\n";
	for ($k=0; $k<@inp_lines; $k++)
	{
        	if ($k == 1)
        	{
                	print OUT " \$SYSTEM MEMORY=70000000 \$END\n";
			print OUT "$inp_lines[$k]\n";
        	}
        	else
        	{
                	print OUT "$inp_lines[$k]\n";
        	}
	}
	close OUT;
}

system "rungms $filename.inp >& $filename.log";							# GAMESS computes electrostatic potentials needed for AM1-
												# -BCC or RESP charges
system "mv $filename.inp $filename-gms.inp";

system "gmsto$method.sh <$filename.log> $filename.ac";						# converts .log .ac file

system "antechamber -i $filename.ac -fi ac -o $filename-$chartype.mol2 -fo mol2";		# produces final .mol2 file with charges

system "antechamber -i $filename.ac -fi ac -o $filename-$chartype.pdb -fo pdb";			# produces final .pdb file

system "charmmgen -i $filename.ac -f ac -o $filename";						# produces rtf and prm files using charmmgen

#-----------------------------------------------------------------------------------		# replaces mol in .pdb file with ligand name

open (PDB, "$filename-$chartype.pdb") || die "couldn't open file $filename-$chartype.pdb: $!\n";
chomp(@pdblines = <PDB>); close PDB;

open (OUT, ">$filename-$chartype.pdb") || die "couldn't open file $filename-$chartype.pdb: $!\n";

$name = uc($filename);

for ($k=0; $k<@pdblines; $k++)
{
	if (substr($pdblines[$k],17,3) =~ "mol")
	{
		substr($pdblines[$k],17,3) = $name;
                print OUT "$pdblines[$k] $net_char\n";
	}
	else
	{
		print OUT "$pdblines[$k]\n";
	}
}
close (OUT) || die "can't close $filename-$chartype.pdb: $!";

#----------------------------------------------------------------------------------------------------------------------

sub short_args
{
        printf ("usage  :       $0 options filename\n");
        printf ("argvs  :       --chartype bcc, rsp\n");
        printf ("               --hydrogen pre, abs\n");
	printf ("               --count_net 	   \n");					#SERVER MOD
        printf ("\n");
}

sub count_cha
{
	$ism_line = shift(@_);
	$count_pos = ($ism_line =~ tr/+//);
	$count_neg = ($ism_line =~ tr/-//);
	$count_net = $count_pos - $count_neg;
}
