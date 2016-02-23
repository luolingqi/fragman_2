#!/usr/bin/perl

if (scalar @ARGV < 1)
{
        print "usage: $0 msfile attraction [chain-res list]\n";
        exit;
}

my $msfile = shift @ARGV;

#set attraction and make sure it's fixed point with one decimal
my $attraction  = shift @ARGV;
$attraction = sprintf('%.1f', $attraction);

my @chainres = @ARGV;

my @lines;

#read ms file into an array so the file can be overwritten
open MS,  "<",  $msfile or die "Failed Opening MS\n";
while (<MS>)
{
	if (/^ATOM/) {
		push @lines, ($_);
	}
}
close MS;

open MS,  ">",  $msfile or die "Failed Writing to MS\n";
foreach $line (@lines)
{
	#skip the line if the atom is buried
	if ('0.0' ne substr $line, 57, 3)
	{
		$linecr = substr $line, 21, 5;
		$linecr =~ s/\s//g;
		$linecr = substr($linecr, 0, 1) . '-' . substr($linecr, 1);
		
		foreach $globalcr (@chainres)
		{
			substr($line, 57, 3) = $attraction if lc($linecr) eq lc($globalcr);
		}
		
	}
	
	
	print MS $line;
}

close MS;
