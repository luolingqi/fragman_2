#!/usr/bin/perl
if (scalar @ARGV <3) {
print "usage : $0 workdir complexs clusterdirectory\n";
exit;
}
$workdir=shift @ARGV;
#$referencedir=shift @ARGV;
$complexs=shift @ARGV;
$dir=shift @ARGV;

chdir($workdir);
open (FILE,"<$complexs");
@a=<FILE>;
$l=@a;
close (FILE);
rmdir("$dir");
mkdir("$dir");
for ($i=0;$i<$l;$i++){
chomp $a[$i];
printf("Processing %s\n",$a[$i]);

#my $ll = `ls -t $a[$i]\/ll.*.out | head -n 1`;
my $ll = `ls -t $a[$i]\/sge.jcf.o* | head -n 1`;
chomp $ll;

print "grep Delta $ll |sort -nk2 >$a[$i]\_min.sort", "\n";
system "grep Delta $ll |sort -nk2 >$a[$i]\_min.sort";


open(ENSORT, "<$a[$i]\_min.sort") or die;
@ensortlines=<ENSORT>;
close(ENSORT);


open(ENER, ">$workdir/$dir/$a[$i].ener") or die;
$nmodel=0;
foreach(@ensortlines)
{
	@data=split(/\s+/,$_);
	$nmodel=$nmodel+1;
	$data[5]=~s/;//;
	$enerline="$data[5]\t$data[8]\t$data[11]\n";
	print(ENER $enerline);


}

close(ENER);

print "$nmodel\n";

open(MINPDB, ">$workdir/$dir/$a[$i].min.pdb") or die;


for($j=0;$j<$nmodel;$j++)
{
$x=$j+1;

print (MINPDB "MODEL\t$x\n");

$y=sprintf("%04\d",$j);

$minout="min_out$y.pdb";


#print "$minout\n";

open(MINOUT, "<$a[$i]/$minout") or die "no $minout";
@minlines=<MINOUT>;
close(MINOUT);
#print (MINPDB "TER\n");
foreach(@minlines)
{
	if($_=~/^ATOM/)
	{
		print (MINPDB "$_");
	}
}
print (MINPDB "ENDMDL\n");
}



close(MINPDB);
#system($cmdline);
}
#$cmdline="conv $dir";
#system($cmdline);
