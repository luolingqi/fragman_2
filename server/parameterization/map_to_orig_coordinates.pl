#!/usr/bin/perl -w
use strict;

my($nac, $new, $mol, @nac_lines, @new_lines, $new_line, @words_new, @words_nac, $i, $j, $coor, $seg1, $seg2);

if (scalar @ARGV < 3)
{
        short_args (); 
	exit;
}

$nac = shift (@ARGV);
$new = shift (@ARGV);
$mol = shift (@ARGV);

open (NAC, "$nac") || die "couldn't open file $nac: $!\n";
chomp(@nac_lines = <NAC>);
close NAC;

open (NEW, "$new") || die "couldn't open file $new: $!\n";
chomp(@new_lines = <NEW>);
close NEW;

open OUT, "> $mol\_new.pdb" || die "cannot open file $mol\_new.pdb\n";

$j = 0;
foreach $new_line(@new_lines)
{
	@words_new = split(/\s+/, $new_line);

	if ($words_new[0] eq "ATOM")
	{
		for($i = $j ; $j < @nac_lines; $j++)				# $j++ is not working because of 'last'
		{
			@words_nac = split(/\s+/, $nac_lines[$j]);

			if ($words_nac[0] eq "ATOM")
			{
				$coor = substr($nac_lines[$j], 30, 24);
				$seg1 = substr($new_line, 0, 30);
				$seg2 = substr($new_line, 54, 2);
				print OUT "$seg1$coor$seg2\n";
		
				$j++;
				last;
			}
		}
	}
	else
	{
		print OUT "$new_line\n";
	}
}

close OUT;

#----------------------------------------------------------------------------------------------------------------------
sub short_args
{
       printf ("usage: $0 tmp_gftogms/nac  molename-chartype.pdb molename\n", $0);
}


