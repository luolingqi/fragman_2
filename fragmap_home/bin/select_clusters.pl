#! /usr/bin/perl -w
use strict;

my ($list_complexs, @complexs, $line_complexs, $user_dist, $cluster);

$user_dist = $ARGV[0];
$list_complexs = $ARGV[1];

open(LIST, "<$list_complexs") || die "Cannot open $list_complexs\n ";
chomp(@complexs=<LIST>); close LIST;

system "ls crosscluster*.pdb | sed \"s/\.pdb//g\" | sort | head -n 6 > list_xclusters";

$cluster="cluster";
foreach $line_complexs(@complexs)
{
	system "ls $line_complexs$cluster*.pdb | sed \"s/\.pdb//g\" | sort > list_$line_complexs";
	system "find_clusters.pl $user_dist list_xclusters list_$line_complexs";
	system "uniq list_$line_complexs\_$user_dist | sort | uniq | head -n 6 > list_$line_complexs\_$user_dist\_uniq";
}
