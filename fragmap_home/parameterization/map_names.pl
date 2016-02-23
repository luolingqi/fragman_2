#!/usr/bin/perl -w
use strict;
use v5.18;

my ($pdb, @pdblines, $pdbline, $min, @minlines, $minline, @words_min, @words_pdb, $seg1, $seg2);
my (%name_table, $rtf, $sat, @rtflines, $rtfline, @satlines, $satline, @words_rtf, @words_sat);
my ($x_old, $y_old, $z_old, $hcounter, $x_new, $y_new, $z_new);

if (scalar @ARGV < 4) {short_args(); exit;}
$pdb = shift (@ARGV);					# original pdb file / May2809 uses lig.mol2 instead because of openbabel adding hydrogen
$min = shift (@ARGV);
$rtf = shift (@ARGV);
$sat = shift (@ARGV);

open (PDB, "$pdb") || die "couldn't open $pdb: $!\n";
chomp(@pdblines=<PDB>); close PDB;

open (MIN, "$min") || die "couldn't open $min: $!\n";
chomp(@minlines=<MIN>); close MIN;

open (RTF, "$rtf") || die "couldn't open $rtf: $!\n";
chomp(@rtflines=<RTF>); close RTF;

open (SAT, "$sat") || die "couldn't open $sat: $!\n";
chomp(@satlines=<SAT>); close SAT;

open OUT, ">$min.rename" || die "couldn't write $min.rename: $!\n";

$hcounter = 99;
foreach $minline(@minlines)
{
	@words_min = split(/\s+/, $minline);

	if ($words_min[0] eq "ATOM")
	{
		foreach $pdbline(@pdblines)
		{
			$pdbline =~ s/^\s+//;
			@words_pdb = split(/\s+/, $pdbline);

			if (defined $words_pdb[6]) 	# modified 102009, originally checks for z-coordinate, but modified to check for y-coordinate
			{
				$words_pdb[2] =~ /([-+]?[0-9]*\.?[0-9]{3})/;
                        	$x_old = $1;
                        	$words_pdb[3] =~ /([-+]?[0-9]*\.?[0-9]{3})/;
                        	$y_old = $1;
                        	$words_pdb[4] =~ /([-+]?[0-9]*\.?[0-9]{3})/;
                        	$z_old = $1;

				$words_min[5] =~ /([-+]?[0-9]*\.?[0-9]{3})/;
                                $x_new = $1;
                                $words_min[6] =~ /([-+]?[0-9]*\.?[0-9]{3})/;
                                $y_new = $1;
                                $words_min[7] =~ /([-+]?[0-9]*\.?[0-9]{3})/;
                                $z_new = $1;

				my $length = length($words_min[2]);     # modified 102009
				if (length($words_min[2]) > 4)
				{
					$words_min[4] =~ /([-+]?[0-9]*\.?[0-9]{3})/;
                                	$x_new = $1;
                                	$words_min[5] =~ /([-+]?[0-9]*\.?[0-9]{3})/;
                                	$y_new = $1;
                                	$words_min[6] =~ /([-+]?[0-9]*\.?[0-9]{3})/;
                                	$z_new = $1;
				}

				if ((abs($x_old - $x_new) lt 0.006 ) and (abs($y_old - $y_new) lt 0.006) and (abs($z_old - $z_new) lt 0.006))
				{
					$seg1 = substr($minline, 0,11);
					$seg2 = substr($minline,17,39);

				########if ($words_pdb[1] eq "H")
				########{
				########	$words_pdb[1] = "H$hcounter";
				########	print "Hydrogen added by Openbabel: H$hcounter\n";
				########	$hcounter = $hcounter - 1;
				########}
			
					if (length($words_min[2]) > 4)
                                	{
                                        	$words_min[2] = substr($words_min[2], 0, 4);
                                	}
	
					$name_table{$words_min[2]} = $words_pdb[1];

					printf (OUT "%11s %4s %39s\n", "$seg1", "$words_pdb[1]", "$seg2");
				}
			}
		}
	}
	else
	{
		print OUT "$minline\n";
	}
}

close OUT;

open OUT, ">$rtf.rename" || die "couldn't write $rtf.rename: $!\n";

foreach $rtfline(@rtflines)
{
	@words_rtf = split(/\s+/, $rtfline);

	if (!defined $words_rtf[0]) {print OUT "\n";}
        else
        {
        given ($words_rtf[0])
        {
		when ("ATOM")
		{
			printf (OUT "%4s %-4s %4s    %9s\n", "$words_rtf[0]", "$name_table{$words_rtf[1]}", "$words_rtf[2]", "$words_rtf[3]");
		}
		when ("BOND")
		{		
			$seg2 = substr($rtfline, 25, 50);
			printf (OUT "%4s %-4s  %-4s         %25s\n", "$words_rtf[0]", "$name_table{$words_rtf[1]}", "$name_table{$words_rtf[2]}", "$seg2");
		}
		when ("ANGL")
		{
			$seg2 = substr($rtfline, 22, 28);
			printf (OUT "%4s %-4s  %-4s  %-4s %28s\n", "$words_rtf[0]", "$name_table{$words_rtf[1]}", "$name_table{$words_rtf[2]}", "$name_table{$words_rtf[3]}" ,"$seg2");

		}
		when ("DIHE")
		{
			$seg2 = substr($rtfline, 29, 21);
			printf (OUT "%4s %-4s  %-4s  %-4s  %-4s %21s\n", "$words_rtf[0]", "$name_table{$words_rtf[1]}", "$name_table{$words_rtf[2]}", "$name_table{$words_rtf[3]}", "$name_table{$words_rtf[4]}", "$seg2");
		}
		when ("IMPH")
		{
			printf (OUT "%4s %-4s  %-4s  %-4s  %-4s\n", "$words_rtf[0]", "$name_table{$words_rtf[1]}", "$name_table{$words_rtf[2]}", "$name_table{$words_rtf[3]}", "$name_table{$words_rtf[4]}");
		}
		default	{print OUT "$rtfline\n";}
        }
	}
}

close OUT;

open OUT, ">$sat.rename" || die "couldn't write $sat.rename: $!\n";	# subatoms

foreach $satline(@satlines)
{
	if ($satline ne "")
	{
		$seg1 = substr($satline,  0, 23);
		$seg2 = substr($satline, 29, 24);

		$satline =~ s/^\s+//;
		@words_sat = split(/\s+/, $satline);
	
		if (defined $words_sat[2])
		{
			printf (OUT "%-23s %-4s %-24s\n", "$seg1", "$name_table{$words_sat[2]}", "$seg2");
		}
	}
}

print OUT "\n";
close OUT;

#-----------------------------------------------------------------------------------------------------------------------

sub short_args
{
	printf ("usage: $0 lig.mol2 lig-bcc-x.0000.pdb lig_final.rtf atoms_subatom\n", $0);
}
