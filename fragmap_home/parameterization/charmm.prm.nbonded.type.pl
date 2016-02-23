#! /usr/bin/perl -w
use strict;

if (scalar @ARGV < 1)
{
	print "usage: $0 FILE:CHARMM.PRM STRING:atomtype\n";
	exit;
}

my $ifile = shift @ARGV;
my $target = shift @ARGV; # target atom type

open (IFILE, "< $ifile") or die "couldn't open $ifile\n";
my @lines = <IFILE>;
close IFILE;
chomp @lines;

# readstate -
# 0: not in NONBONDED
# 1: in NONBONDED
my $readstate = 0;

my $bestat = ''; # best match atom type
my $bestline = ''; # best match line, paired with $bestat
my $i;
for ($i = 0; $i < scalar @lines; $i++)
{
	my $line = $lines[$i];
	next if ($line =~ m/^!/); # comment

	if ($readstate == 0 && ($line =~ m/^NBONDED\s+/ || $line =~ m/^NONBONDED\s+/))
	{
		while ($line =~ m/-$/) # line continuation
		{
			$i++;
			$line = $lines[$i];
		}
		$readstate = 1;
	}
	elsif ($readstate == 1 && ($line =~ m/^NBFIX/ || $line =~ m/^VOLUME/))
	{
		last;
	}
	elsif ($readstate == 1 && $line !~ m/^\s*$/)
	{	
		$line=~s/^\s+//;
		my @words = split /\s+/, $line;
		my $atomtype = $words[0];
#		print "$line\t";
#		print "$atomtype\n";
		# get perl regex from charmm regex:
		my $perl_atomtype = '^' . $atomtype;
		$perl_atomtype =~ s/\*/\\w*/g;
		$perl_atomtype =~ s/%/\\w/g;
		$perl_atomtype .= '$';
	#	print $perl_atomtype;
		if ($target =~ m/$perl_atomtype/)
		{
			$atomtype =~ s/%/+/g; # allow for easy lt string comparison
			if ($bestat eq '')
			{
				$bestat = $atomtype;
				$bestline = $line;
			}
			elsif ($atomtype gt $bestat)
			{
				$bestat = $atomtype;
				$bestline = $line;
			}
		}
	}
}

print "$bestline\n";
