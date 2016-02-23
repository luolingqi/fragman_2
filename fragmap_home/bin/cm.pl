#! /usr/bin/perl -w

$file = $ARGV[0];

$sumx = 0;
$sumy = 0;
$sumz = 0;
$count = 0;
$i=0;

open(FILE, "< $file") || die "Can't open $file\n";
while($line = <FILE>)
{
	if ($line =~ /HETATM/ || $line =~ /ATOM/)
	{
		$x1[$i] = substr($line, 30, 8);
                $y1[$i] = substr($line, 38, 8);
                $z1[$i] = substr($line, 46, 8);
                $x1[$i] =~ s/\s+//g;
                $y1[$i] =~ s/\s+//g;
                $z1[$i] =~ s/\s+//g;
	        
		$sumx += $x1[$i];
                $sumy += $y1[$i];
                $sumz += $z1[$i];
                $count++;
		$i++;
	}
}
close FILE;

$x0 = $sumx/$count;
$y0 = $sumy/$count;
$z0 = $sumz/$count;

printf ("%3.3f\t%3.3f\t%3.3f\n", $x0,$y0,$z0);
