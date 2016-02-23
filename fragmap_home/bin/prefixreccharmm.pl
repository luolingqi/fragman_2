#! /usr/bin/perl -w

if(@ARGV<2)
{

print "Usage:$0 1rec.ms outfile\n";

exit;
}

$file=shift;
$outfile=shift;

open(FILE, "<$file") or die;
@liglines=<FILE>;
close(FILE);


open(OFILE, ">$outfile") or die;



for my $line (@liglines)
{
        next if ($line !~ /^ATOM/);

        my $prefix = substr ($line, 0, 30);
        my $suffix = substr ($line, 54,4 );

        my $newx = substr ($line, 30, 8);
        my $newy = substr ($line, 38, 8);
        my $newz = substr ($line, 46, 8);
#print "$x, $y, $z\n";
	$text=substr($line, 0, 70);

$text=~s/\n//g;
	$length=length($text);
for($i=0;$i<70-$length;$i++)
{
$text=$text." ";

}
	

        # rotate and translate

  #      printf OFILE ("%30s%8.3f%8.3f%8.3f  %4.1f", $prefix, $newx, $newy, $newz, $suffix);
 #       print (OFILE "   0.0      AAAA\n");

	print OFILE "$text  AAAA\n";
}
print (OFILE "END\n");
close(OFILE);

