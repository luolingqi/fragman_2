#! /usr/bin/perl -w
use strict;
use v5.18;

my ($rtf, $prm, $gaf, $num, @rtflines, $rtfline, @elements, $elements, %atm_type, @sublines, $subline);
my (%gaf_type, @gaflines, $gafline, @gaf_line, @rtf_words, @prmlines, @prm_line, $prmline);

if (scalar @ARGV < 4) {short_args(); exit;}

$rtf = lc($ARGV[0]);
$prm = lc($ARGV[1]);
$gaf = lc($ARGV[2]);
$num = lc($ARGV[3]);

open (RTF, "$rtf") || die "couldn't open file $rtf: $!\n";
chomp (@rtflines = <RTF>); close RTF;

foreach $rtfline(@rtflines)
{
        $rtfline =~ s/^\s+//;                                                                                                   # remove leading space
        @elements = ();
        @elements = split(/\s+/, $rtfline);

        if (defined $elements[0])
        {
                given ($elements[0])
                {
                        when ("ATOM")
                        {
				$atm_type{$elements[1]} = $elements[2];
			}
		}
	}
}

open (GAF, "$gaf") || die "couldn't open file $gaf: $!\n";
chomp (@gaflines = <GAF>); close GAF;

foreach $gafline(@gaflines)
{
	@gaf_line = split /\s+/, $gafline;
        $gaf_type{$gaf_line[0]} = $gaf_line[3];
}


system "touch temp";
system "charmm.rtf,.prm,map.subatom2prm.atom,.subatom.pl $rtf $prm temp 1 $num > tmp_atoms1";
system "rm temp";

open (SUB, "tmp_atoms1") || die "couldn't open file tmp_atoms1: $!\n";
chomp (@sublines = <SUB>); close SUB;

open (OUT, ">subatom.inp") || die "cannot create subatom.inp: $!";

foreach $subline(@sublines)
{
	$subline =~ s/^\s+//;
	@rtf_words = split /\s+/, $subline;

	if (exists $gaf_type{"$atm_type{$rtf_words[2]}"})
	{
		print OUT "$rtf_words[1] $rtf_words[2] $gaf_type{\"$atm_type{$rtf_words[2]}\"}\n";
	}

}	
close (OUT) || die "can't close subatom.inp: $!";

system "charmm.rtf,.prm,map.subatom2prm.atom,.subatom.pl $rtf $prm subatom.inp 1 $num > tmp_atoms2";

open (PRM, "tmp_atoms2") || die "couldn't open file tmp_atoms: $!\n";
chomp (@prmlines = <PRM>); close PRM;

open (OUT, ">atoms_subatom") || die "cannot create atoms_subatom: $!";

foreach $prmline(@prmlines)
{
        @prm_line = split /\s+/, $prmline;
	printf (OUT "        %-5s   %-4s    %-2s      %-2s      %-5s  %6s\n", "$prm_line[1]", "$prm_line[2]", "$prm_line[3]", "$prm_line[4]", "$prm_line[5]"		    , "$prm_line[6]");
}

print OUT "\n";

close (OUT) || die "can't close atoms_subatom: $!";

#----------------------------------------------------------------------------------------------------------------------
sub short_args
{
        printf ("usage  :       $0 lig\_unify.rtf lig.prm gaff2mol2 100\n");
}
