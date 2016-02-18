#!/usr/bin/perl
if (scalar @ARGV <4) {
print "usage : $0 workdir referencedir complexs directory\n";
exit;
}
$workdir=shift @ARGV;
$referencedir=shift @ARGV;
$complexs=shift @ARGV;
$dir=shift @ARGV;

chdir($workdir);
open (FILE,"<$complexs");
@a=<FILE>;
close (FILE);
$l=$#a;
$str="cat rec/1rec*.ener | grep 'EXTERN>' |awk '{print \$4}'";
$rec_elec=`$str`;
$str="cat rec/1rec*.ener | grep 'EXTERN>' |awk '{print \$3}'";
$rec_vdw=`$str`;
chomp $rec_elec;
chomp $rec_vdw;

printf ("Rec self  elec energy:%s %s \n",$rec_vdw,$rec_elec);
rmdir("$dir");
mkdir("$dir");
for ($i=0;$i<=$l;$i++){
chomp $a[$i];
printf("Processing %s\n",$a[$i]);
$str="cat $referencedir/charmm.$a[$i].ener | grep 'EXTERN>' |awk '{print \$4}'";
$lig_elec=`$str`;
chomp $lig_elec;
$str="cat $referencedir/charmm.$a[$i].ener | grep 'EXTERN>' |awk '{print \$3}'";
$lig_vdw=`$str`;
chomp $lig_vdw;
printf ("Lig %s self   energy:%s %s \n",$a[$i],$lig_vdw,$lig_elec);

#$cmdline="seq 0 199 |awk '{print \$1*10+1}' |xargs -i bash -c \"cat $a[$i]/charmmoutput/$a[$i].{}.min.pdb >> $dir/tmp\"";
#$cmdline="find $a[$i]/charmmoutput/ -name '$a[$i].*.min.pdb' | perl -e 'print sort{\$c=\$a;\$d=\$b;\$c=~s/\$a[\$i]\\/charmmoutput\\/\$a[\$i]\\.//;\$d=~s/\$a[\$i]\\/charmmoutput\\/\$a[\$i]\\.//;\$c <=> \$d;} <>' |xargs cat > $dir/tmp";
$cmdline="find $a[$i]/charmmoutput/ -name '$a[$i].*.min.pdb' | perl -e 'print sort{\$c=\$a;\$d=\$b;\$c=~s/$a[$i]\\/charmmoutput\\/$a[$i].\/\/;\$d=~s/$a[$i]\\/charmmoutput\\/$a[$i].\/\/;\$c <=> \$d;} <>' |xargs cat > $dir/tmp";
print $cmdline;
system($cmdline);
system("mv $dir/tmp $dir/$a[$i].min.pdb");
#$cmdline="seq 0 199 |awk '{print \$1*10+1}' |xargs -i bash -c \"cat $a[$i]/charmmoutput/$a[$i].{}.ener | grep 'EXTERN>' >> $dir/$a[$i].en\"";
$cmdline="find $a[$i]/charmmoutput/ -name  '$a[$i].*.ener' | perl -e 'print sort{\$c=\$a;\$d=\$b;\$c=~s/$a[$i]\\/charmmoutput\\/$a[$i].\/\/;\$d=~s/$a[$i]\\/charmmoutput\\/$a[$i].\/\/;\$c <=> \$d;} <>' |xargs grep 'EXTERN>' > $dir/$a[$i].en";
system($cmdline);
print $cmdline;
$cmdline="cat $dir/$a[$i].en |awk '{printf(\"%.3f %s %.3f\\n\", \$3+\$4-($lig_elec)-($rec_elec)-($lig_vdw)-($rec_vdw),\$3-($lig_vdw)-($rec_vdw),\$4-($lig_elec)-($rec_elec))}' >$dir/$a[$i].ener";
print $cmdline."\n";
system($cmdline);
$cmdline="rm $dir/$a[$i].en;";
#system($cmdline);
}
#$cmdline="conv $dir";
#system($cmdline);
