#!/usr/bin/perl -w
use strict;
use v5.18;

my ($rtf, @rtflines, $rtfline, @elements, $identifier, $hx, $cx, $l, $carbon, $hydrogen);
my (@h_removed, $h_recorded, %charge_table, %bond_table, $swit, $h_safety, $final_charge);

if (scalar @ARGV < 1) {short_args (); exit;}
$rtf = shift(@ARGV);

open (RTF, "$rtf.rtf") || die "couldn't open file $rtf.rtf: $!\n";
chomp(@rtflines = <RTF>); close RTF;

open OUT, ">$rtf\_unify.rtf" || die "couldn't write $rtf\_unify.rtf: $!";

$l = 0;
$final_charge = 0;
foreach $rtfline(@rtflines)
{
	$rtfline =~ s/^\s+//;                                   			                                        # remove leading space
	@elements = ();
        @elements = split(/\s+/, $rtfline);
	$h_safety = 1;

        if (defined $elements[0])
        {
        	given ($elements[0])
		{
			when ( "ATOM")
			{
				$charge_table{$elements[1]} = sprintf "%.6f", $elements[3];

				$cx = $elements[1];             # atom name
				$hx = $elements[2];             # atom type

				if ($hx =~ m/h\d|hc|ha|hx/)
				{
        				$h_removed[$l] = $cx;
        				$l++;
				}
			}

			when ( "BOND")
                	{
				$carbon=$elements[1]; $hydrogen=$elements[2];

				if ($hydrogen =~ m/(H\d'?)/)
				{
					if ($carbon =~ m/(C\d'?)/)
					{
						$charge_table{$carbon} = $charge_table{$carbon} + $charge_table{$hydrogen};

						foreach $h_recorded(@h_removed)
						{
							if ($hydrogen eq $h_recorded)
							{
								$h_safety = "0"; last;
							}
						}
						if ($h_safety eq "1")
						{
							$h_removed[$l] = $hydrogen;
							$l++;
							print "C-H Hydrogen missed: $hydrogen";
						}
					}
					else
                                	{
                                        	print "Hydrogen not removed: $hydrogen\n";
                                	}
				}		
			}
		}
	}
}

foreach $carbon(keys %charge_table)
{
	$h_safety = 1;
	foreach $h_recorded(@h_removed)
	{
		if ($carbon eq $h_recorded)
		{
			$h_safety = "0"; last;
		}
	}
	if ($h_safety eq "1")
	{
		$final_charge = $final_charge + $charge_table{$carbon};
		printf ("%-3s %6.3f\n", "$carbon", "$charge_table{$carbon}\n");
	}
}

printf ("%s %6.3f\n", "Overall charge after unifying:", "$final_charge");

foreach $rtfline(@rtflines)
{
	$rtfline =~ s/^\s+//;
	@elements = ();
	@elements = split(/\s+/, $rtfline);

	if (!defined $elements[0]) {print OUT "\n";}
	else
	{
	given ($elements[0])
	{
			$swit = 1;
			when ( "MASS")
                        {
				if ($elements[2] !~ m/h\d|hc|ha|hx/)
				{
					$hx = uc($elements[2]);
					printf (OUT "%-7s %2s %-2s     %9s\n", "$elements[0]", "$elements[1]", "$hx", "$elements[3]");
				}
                        }

			when ("ATOM")
			{
				$cx = $elements[1];
				search_chs($elements[1], \@h_removed);
                                if($swit eq "1")
                                {
                                        $hx = uc($elements[2]);
					printf (OUT "%-4s %-3s   %-3s    %9.6f\n", "$elements[0]", "$elements[1]", "$hx", "$charge_table{\"$cx\"}");
                                }
			}

			when ("BOND")
			{	
				search_chs($elements[1], \@h_removed);
                                search_chs($elements[2], \@h_removed);

                                if($swit eq "1"){print (OUT "$rtfline\n");}
        		}

			when ( "ANGL") {searchXXXX ($rtfline);}
                        when ( "DIHE") {searchXXXX ($rtfline);}
                        when ( "IMPH") {searchXXXX ($rtfline);}

                        default
                        {
                                if ($elements[0] eq "RESI")
                                {
                                        $rtf = uc($rtf);
					printf (OUT "%4s %3s  %6.3f\n", "RESI", "$rtf", "$final_charge");
                                }
                                elsif ($elements[0] eq "99")
                                {
                                        printf (OUT "   %6s\n", "$rtfline");
                                }
                                else
                                {
                                        print OUT "$rtfline\n";
                                }
                        }
		}
		}
	}


close (OUT) || die "can't close rtf_unify.rtf: $!";

print "Hydrogens removed   : @h_removed\n";

#----------------------------------------------------------------------------------------------------------------------

sub short_args
{
        printf ("usage: $0 name [small molecule w/o extension; require name.rtf and name.prm]\n", $0);
}

sub search_chs
{
	my $c = $_[0]; my @h = @{$_[1]};
	my $h;

	foreach $h(@h)
	{
		if($c eq $h){$swit = 0; last;}
	}
}

sub searchXXXX
{
	my $line = $_[0];
		
	search_chs($elements[1], \@h_removed);
	search_chs($elements[2], \@h_removed);
	search_chs($elements[3], \@h_removed);
	search_chs($elements[4], \@h_removed);

	if($swit eq "1"){print (OUT "$line\n");}
}

# $rtfline =~ /(\w+'?)\s+(\w+'?)\s+(\w+'?)\s+([-+]?[0-9]*\.?[0-9]+)/;                            # for future reference	
