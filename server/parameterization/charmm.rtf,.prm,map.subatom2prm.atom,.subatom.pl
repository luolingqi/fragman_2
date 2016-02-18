#! /usr/bin/perl -w
use strict;

if (scalar @ARGV < 5)
{
	print "usage: $0 FILE:CHARMM.RTF FILE:CHARMM.PRM FILE:MAP.SUBATOM BOOL:IFILE_ZEROBASED? INT:START_ATOMN\n";
	exit;
}

my $file_rtf = shift @ARGV;
my $file_prm = shift @ARGV;
my $file_map = shift @ARGV;
my $zerobased = shift @ARGV;
my $start_atomn = pop @ARGV;

my @maps = `cat $file_map`;
my %maph;
for my $map (@maps)
{
	my @map_words = split /\s+/, $map;
	my $key = $map_words[0] . ' ' . $map_words[1];
	$maph{$key} = $map_words[2];
}

for my $key (sort keys %maph)
{
	my $val = $maph{$key};
#print "$key => $val\n";
}

my $cmd = "grep ^RESI $file_rtf | " . "awk '{print " . '$2' . "}'";
my @resis = `$cmd`;
chomp @resis;

my @atom_lines = ();
my $atomn = $start_atomn;
for my $resi (@resis)
{
	my $major_type = $resi;

	my @rtf_lines = `charmm.rtf.resi.atom.pl $file_rtf $major_type`;
	for my $rtf_line (@rtf_lines)
	{
		my @rtf_words = split /\s+/, $rtf_line;
		my $minor_type = $rtf_words[1];
		my $minor_type_alt = $rtf_words[2];
		my $charge = $rtf_words[3];

		my $prm_line = `charmm.prm.nbonded.type.pl $file_prm $minor_type_alt`;
		my @prm_words = split /\s+/, $prm_line;
		my $radius = $prm_words[3];

		my $key = "$major_type $minor_type";
		my $subatomn = -1;
		if (exists $maph{$key})
		{
			$subatomn = $maph{$key};
			if (! $zerobased)
			{
				$subatomn--; # start from zero
			}
		}
		else
		{
			if ($minor_type !~ m/^H/)
			{
#				print "warning: subatomn of $key not defined\n";
				$subatomn--;
			}
		}

		my $atom_line = sprintf ("%6d%6s%-10s%-8s%4d%4s%8.3f%12.3f", $atomn, '', $major_type, $minor_type, $subatomn, '', $radius, $charge);
	print "$atom_line\n";

		$atomn++;
	}
}
