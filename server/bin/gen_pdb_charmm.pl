#! /usr/bin/perl -w


if (scalar @ARGV < 4)
{
	print "usage: $0 ROTPRM FTDATFILE PDB_IFILE OUTFILE\n";
	exit;
}

# command line args
my $rotprm = shift @ARGV;
my $datfile = shift @ARGV;
#my $datfile_i = shift @ARGV; # datfile index (0 based)
my $ligand = shift @ARGV;
my $ofile= shift @ARGV;


if ($ligand eq $ofile)
{
	print "PDB_IFILE is identical to PDB_OFILE, exiting\n";
	exit;
}


# read in ftdatfile
open (DATFILE, "< $datfile") or die "could not open $datfile: $!\n";
my @datlines = <DATFILE>;
close DATFILE;
$datanum=@datlines;


if (scalar @datlines < 1)
{
        print "error: $datfile contains no data\n";
        exit;
}
# 

open (ROTFILE, "< $rotprm") or die "could not open $rotprm: $!\n";

my @rotlines = <ROTFILE>;
close ROTFILE;
if (scalar @rotlines < 1)
{
        print "error: $rotprm contains no data\n";
        exit;
}

open (LIG, "< $ligand") or die "couldn't open $ligand: $!\n";
#open (OFILE, "> $ofile") or die "could not open $ofile: $!\n";
my @liglines = <LIG>;
close LIG;
chomp @liglines; # remove newlines
my $modelnum=0;






open(OFILE,">$ofile");

for($i=0;$i<$datanum;$i++)
{
@words=split(/\s+/,$datlines[$modelnum]);
$modelnum=$modelnum+1;
my $roti = $words[0]; # rotation index into rotprm
my $tx = $words[1]; # x translation
my $ty = $words[2]; # y translation
my $tz = $words[3]; # z translation
my $rotline = $rotlines[$roti]; # get rotation roti

$rotline =~ s/^\s+//; # remove leading whitespace
# get values
my $rotseti;
my $a11; my $a12; my $a13;
my $a21; my $a22; my $a23;
my $a31; my $a32; my $a33;
($rotseti, $a11,$a12,$a13,$a21,$a22,$a23,$a31,$a32,$a33) = split(/\s+/, $rotline);

print (OFILE "MODEL     $modelnum\n");
for my $line (@liglines)
{
        next if ($line !~ /^ATOM/);

        my $prefix = substr ($line, 0, 30);
        my $suffix = substr ($line, 54);

        my $x = substr ($line, 30, 8);
        my $y = substr ($line, 38, 8);
        my $z = substr ($line, 46, 8);
#print "$x, $y, $z\n";

        # rotate and translate
        my $newx = $a11*$x + $a12*$y + $a13*$z + $tx;
        my $newy = $a21*$x + $a22*$y + $a23*$z + $ty;
        my $newz = $a31*$x + $a32*$y + $a33*$z + $tz;
	
        printf OFILE ("%30s%8.3f%8.3f%8.3f  %4.1f", $prefix, $newx, $newy, $newz, $suffix);
	print (OFILE "   0.0      BBBB\n");
}
print (OFILE "ENDMDL\n");


}
close(OFILE);





#print "$tx, $ty, $tz\n";


# read in rotprm
#shift @rotlines; # remove rotprm header

#print "$rotline\n";
#print "$a11 $a12 $a13 $a21 $a22 $a23 $a31 $a32 $a33\n";



# read in ligand pdb file and print ofile


# print the ofile
#print "$x, $y, $z\n";

	# rotate and translate


#close OFILE;
