#!/usr/bin/perl -w
use strict;

my ($rtf, @rtf, $prm, @prm, $pdb, @pdb, $par, @par, @tmp, $tmp);
my ($i, $j, $k, $l, $m, $s, $rtf_0, $rtf_1, $rem, $flp, $dum);
my ($ats, @ats, $sub, @sub, $mass_switch);

#usage: $0 pdbamino.rtf parm.prm atoms lig.rtf lig.prm atoms_subatom

if (scalar @ARGV < 6)
{
	short_args ();
	exit;
}

$pdb = shift (@ARGV); read_files($pdb); @pdb = @tmp;
$par = shift (@ARGV); read_files($par); @par = @tmp;
$ats = shift (@ARGV); read_files($ats); @ats = @tmp;
$rtf = shift (@ARGV); read_files($rtf); @rtf = @tmp;
$prm = shift (@ARGV); read_files($prm); @prm = @tmp;
$sub = shift (@ARGV); read_files($sub); @sub = @tmp;

open OUT, "> pdbamino_new.rtf" || die "cannot open file pdbamino_new.rtf\n";

$i = 0;
$j = 0;
$s = 1;
$pdb[$i] =~ /(\w+|\W)/;
	
while ($1 ne "END")
{
	if ($1 eq "\n" and $s == 1)
	{
		print_segm ($j, "MASS");	# print_segm is only for printing concatenating MASS lines, seg2 is for all remaining lines
		$s = 2;
	}

	print OUT $pdb[$i];
	
	$i++;		
	$pdb[$i] =~ /(\w+|\W)/;
}

for($k = $j + 1; $k < @rtf; $k++)
{
	print OUT $rtf[$k];
}                

print OUT "END\n";

close OUT;

#----------------------------------------------------------------------------------------------------------------------

open OUT, "> parm_new.prm" || die "cannot open file parm_new.prm\n";

$i = 0;
$j = 0;
$s = 1;
$par[$i] =~ /(\w+|\W)/;

while ($1 ne "END")
{
        if ($1 eq "THETAS")
        {
                print_seg2 ("BOND"); #this section from lig_final.rtf is now being appended
        }

	if ($1 eq "PHI")
        {
                print_seg2 ("ANGLE");
        }

 	if ($1 eq "IMPHI")
        {
                print_seg2 ("DIHEDRAL");
        }

	if ($1 eq "NBONDED")
        {
                print_seg2 ("IMPHI");
        }

	if ($1 eq "!")
	{	
		$par[$i+1] =~ /(\w+|\W)/;
		
		if ($1 eq "NBFIX")
        	{
			print_seg2 ("NONBONDED");
        	}

	}

        print OUT $par[$i];

        $i++;
        $par[$i] =~ /(\w+|\W)/;
}

print OUT "END\n";

close OUT;

#----------------------------------------------------------------------------------------------------------------------

open OUT, "> atoms_new.ats" || die "cannot open file atoms_new.ats\n";

$i = 0;
$j = 0;
$ats[$i] =~ /(\w+|\W)/;

while ($1 ne "END")
{
	if ($1 eq "\n")
	{
		$ats[$i+1] =~ /(\#\w+|\#\W)/;

#		if ($1 eq "HYDROGEN")
		if ($1 eq "#PAIRWISE")
        	{
			while ($1 ne "\n")
        		{
                        #	print OUT $sub[$j];

				my $tmp = $sub[$j];
				$tmp =~ s/^\s+//;
				printf(OUT "atom        $tmp");
                	
                		$j++;
                		$sub[$j] =~ /(\w+|\W)/;
			}
        	}
	}
	
	print OUT $ats[$i];
	
	$i++;
	$ats[$i] =~ /(\w+|\W)/;
}

print OUT "END";
close OUT;

#----------------------------------------------------------------------------------------------------------------------

sub short_args
{
	printf ("usage: $0 pdbamino.rtf parm.prm atoms lig.rtf lig.prm atoms_subatom\n", $0);
}

sub read_files
{
	$tmp = shift(@_);
	@tmp = ();

	open(TMP, "$tmp") || die "cannot open file $tmp: $!\n";

	if (-s "$tmp") {}
	else
	{
		print "$tmp not properly created...exiting\n"; exit;
	}

        @tmp = <TMP>; 
	close TMP;
}

sub print_segm
{
	$j = shift(@_);
	$l = shift(@_);
	$rtf[$j] =~ /(\w+|\W)/;

	$mass_switch = 0;
	while ($1 ne "\n")
	{
		if ($1 eq "$l")
		{
			print OUT $rtf[$j];
			$mass_switch = 1;
		}

		$j++;
		$rtf[$j] =~ /(\w+|\W)/;
	}
	if ($mass_switch != 1)
	{
		print "$rtf not properly created...exiting\n"; exit;
	}
}

sub print_seg2
{
	$l = shift(@_);
	$flp = 0; $m = 0;

	foreach $j(@prm)
	{
		$j =~ /(\w+|\W)/;

		if ($1 eq $l)
		{
			$flp = "1";
			next;
		}
		if ($flp eq "1")
		{
			if ($l eq "NONBONDED" and $m < 3)
			{
				$m++; next;
			}

			print OUT $j;
		}
		if ($1 eq "\n")
		{
			$flp = "0";
		}
	}
}
