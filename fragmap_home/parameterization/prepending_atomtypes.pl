#!/usr/bin/perl -w
use strict;
use v5.18;

my ($cha, $uni, $prp, $rtf, $prm, @rtflines, $rtfline, @elements, $number, @prmlines, $prmline);
my ($bond, $angl, $dihe, $imph, $nbnd, $tag, $atm_1u, $atm_2u, $atm_3u, $atm_4u, $n, $segment);
my ($y1, $y2, $y3, $y4);
my ($imphi);

if (scalar @ARGV < 4) {short_args (); exit;}

$cha = $ARGV[0];
$uni = $ARGV[1];
$prp = $ARGV[2];
$rtf = $ARGV[3];
$prm = $ARGV[3];

open (RTF, "$rtf\_unify.rtf") || die "couldn't open file $rtf-unify.rtf: $!\n";
chomp(@rtflines = <RTF>); close RTF;

open (OUT, ">$rtf\_final.rtf") || die "couldn't open file $rtf\_final.rtf: $!\n";

$number = 201;
foreach $rtfline (@rtflines)
{
	$rtfline =~ s/^\s+//;
	#$_ = $rtfline;
        @elements = ();
        @elements = split(/\s+/, $rtfline);
	
	if (defined $elements[0])
        {
		given ($elements[0])
		{
			when ( "MASS")
			{
				printf (OUT "%-7s %2s %-3s    %9.6f\n", "$elements[0]", "$number", "Y$elements[2]", "$elements[3]");
				$number++;
			}

			when ( "ATOM")
                	{
                                printf (OUT "%-4s %-3s   %-4s   %9.6f\n", "$elements[0]", "$elements[1]", "Y$elements[2]", "$elements[3]");
    			}
			default
			{
				print OUT "$rtfline\n";
			}
		}
	}
	else
	{
		print OUT "$rtfline\n";
	}
}

close (OUT) || die "can't close $rtf\_final.rtf: $!";

open (PRM, "$prm.prm") || die "couldn't open file $prm.prm: $!\n";
chomp (@prmlines = <PRM>); close PRM;

open (OUTU, ">$prm\_unify.prm") || die "couldn't open file $prm\_unify.prm: $!\n";
open (OUTF, ">$prm\_final.prm") || die "couldn't open file $prm\_final.prm: $!\n";

$bond = 0; $angl = 0; $dihe = 0; $imph = 0; $nbnd = 0; $n = 0;
foreach $prmline (@prmlines)
{
        $prmline =~ s/^\s+//;													# remove leading space
        @elements = ();
        @elements = split(/\s+/, $prmline);

        if (defined $elements[0])
        {
		$tag = $elements[0];

		if 	($bond eq "1") {$tag = "bond_lines";}
		elsif	($angl eq "1") {$tag = "angl_lines";}
                elsif   ($dihe eq "1") {$tag = "dihe_lines";}
                elsif   ($imph eq "1") {$tag = "imph_lines";}
                elsif   ($nbnd eq "1") {$tag = "nbnd_lines";}

		given ($tag)
		{
			 when ("BOND") 
			{
				$bond = "1";
				print OUTU "$prmline\n";
                                print OUTF "$prmline\n";
			}
			 when ("bond_lines")
			{	
				$atm_1u = uc($elements[0]);
				$atm_2u = uc($elements[1]);
				printf (OUTU "%-4s %-4s   %5s       %5s\n", "$atm_1u", "$atm_2u", "$elements[2]", "$elements[3]");
				printf (OUTF "%-4s %-4s   %5s       %5s\n", "Y$atm_1u", "Y$atm_2u", "$elements[2]", "$elements[3]");	
			}	

			 when ("ANGLE")
			{
				$angl = "1";
				print OUTU "$prmline\n";
                                print OUTF "$prmline\n";
			}
			
			 when ("angl_lines")
			{
				$atm_1u = uc($elements[0]);
                                $atm_2u = uc($elements[1]);
				$atm_3u = uc($elements[2]);
				printf (OUTU "%-4s %-4s %-4s    %5s     %5s\n", "$atm_1u", "$atm_2u", "$atm_3u", "$elements[3]", "$elements[4]");
				printf (OUTF "%-4s %-4s %-4s    %5s     %5s\n", "Y$atm_1u", "Y$atm_2u", "Y$atm_3u", "$elements[3]", "$elements[4]");
			}
		
			 when ("DIHEDRAL")
			{
				$dihe = "1";
				print OUTU "$prmline\n";
                                print OUTF "$prmline\n";
			}

			 when ("dihe_lines")
			{
				$atm_1u = uc($elements[0]);
                                $atm_2u = uc($elements[1]);
                                $atm_3u = uc($elements[2]);
				$atm_4u = uc($elements[3]);

				$segment = substr($prmline,16,70);
				
				if ($atm_1u eq "X") {$y1 = "";} else {$y1 = "Y";}
                                if ($atm_2u eq "X") {$y2 = "";} else {$y2 = "Y";}
                                if ($atm_3u eq "X") {$y3 = "";} else {$y3 = "Y";}
                                if ($atm_4u eq "X") {$y4 = "";} else {$y4 = "Y";}
				printf (OUTU "%-4s %-4s %-4s %-4s     %-4s %-2s %5s\n", "$atm_1u", "$atm_2u", "$atm_3u", "$atm_4u", 								    "$elements[4]", "$elements[5]", "$elements[6]");
                                printf (OUTF "%-4s %-4s %-4s %-4s     %-4s %-2s %5s\n", "$y1$atm_1u", "$y2$atm_2u", "$y3$atm_3u", "$y4$atm_4u", 						    "$elements[4]", "$elements[5]", "$elements[6]");
			}

			 when ("IMPHI")
			{
				$imph = "1";
				print OUTU "$prmline\n";
                                print OUTF "$prmline\n";
			}

			 when ("imph_lines")
			{
				$atm_1u = uc($elements[0]);
                                $atm_2u = uc($elements[1]);
                                $atm_3u = uc($elements[2]);
                                $atm_4u = uc($elements[3]);

				if ($atm_1u eq "X") {$y1 = "";} else {$y1 = "Y";}
				if ($atm_2u eq "X") {$y2 = "";} else {$y2 = "Y";}
				if ($atm_3u eq "X") {$y3 = "";} else {$y3 = "Y";}
				if ($atm_4u eq "X") {$y4 = "";} else {$y4 = "Y";}

				if ($elements[4] < 10) {$imphi = 50;}	
				
				printf (OUTU "%-4s %-4s %-4s %-4s     %5s %-2s %5s\n", "$atm_1u", "$atm_2u", "$atm_3u", "$atm_4u", 								    "$imphi$elements[4]", "$elements[5]", "$elements[6]");
                                printf (OUTF "%-4s %-4s %-4s %-4s     %5s %-2s %5s\n", "$y1$atm_1u", "$y2$atm_2u", "$y3$atm_3u", "$y4$atm_4u", 						    	    "$imphi$elements[4]", "$elements[5]", "$elements[6]");
			}

			 when ("NONBONDED")
			{
				$nbnd = "1";
				print OUTU "$prmline\n";
                                print OUTF "$prmline\n";
			}

			 when ("nbnd_lines")
			{
				if ($n < 3)
				{
					print OUTU "$prmline\n";
					print OUTF "$prmline\n";
					$n += 1;
				}
				else
				{					
					$atm_1u = uc($elements[0]);

					printf (OUTU "%-4s     %-4s    %-7s       %-5s   %-4s  %-4s  %-4s\n", "$atm_1u", "$elements[1]", "$elements[2]", 						    "$elements[3]", "$elements[4]", "$elements[5]", "$elements[6]");
					printf (OUTF "%-4s     %-4s    %-7s       %-5s   %-4s  %-4s  %-4s\n", "Y$atm_1u", "$elements[1]", "$elements[2]", 						    "$elements[3]", "$elements[4]", "$elements[5]", "$elements[6]");
				}
			}
			default
        		{
				$bond = 0; $angl = 0; $dihe = 0; $imph = 0; $nbnd = 0;
                		print OUTU "$prmline\n";
                		print OUTF "$prmline\n";
        		}	
		}
	}
	else
	{
		$bond = 0; $angl = 0; $dihe = 0; $imph = 0; $nbnd = 0;
		print OUTU "$prmline\n";
		print OUTF "$prmline\n";
	}
}

close (OUTU) || die "can't close $prm\_unify.prm: $!";
close (OUTF) || die "can't close $prm\_final.prm: $!";
#----------------------------------------------------------------------------------------------------------------------

sub short_args
{
        printf ("usage	:	$0 options filename\n");
        printf ("argvs	:	--chartype bcc, rsp\n");
	printf ("	:	--rtfunify filename\n");
	printf ("	:	--prepfile filename\n");
        printf ("\n");
}
