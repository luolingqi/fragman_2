#! /usr/bin/perl -w
use strict;

# Script arguments: clustering_distance, list_of_consensus_clusters, list_of_user_clusters_for one_ligand

my ($list_xclusters, $list_pclusters, @xclusters, @pclusters, $line_xclusters, $line_pclusters, $user_dist);
my ($COM_xcluster, @COM_xnumbers, $COM_ycluster, @COM_ynumbers, $dist_x, $dist_y, $dist_z, $dist);

$user_dist = $ARGV[0];
$list_xclusters = $ARGV[1];
$list_pclusters = $ARGV[2];

open(LIST_X, "<$list_xclusters") || die "Cannot open $list_xclusters\n";
chomp(@xclusters=<LIST_X>); close LIST_X;

open(LIST_P, "<$list_pclusters") || die "Cannot open $list_pclusters\n";
chomp(@pclusters=<LIST_P>); close LIST_P;

foreach $line_xclusters(@xclusters)
{
	$COM_xcluster = `cm.pl $line_xclusters.pdb`;
	@COM_xnumbers = split(/\s+/, $COM_xcluster);
	
	foreach $line_pclusters(@pclusters)
	{
		$COM_ycluster = `cm.pl $line_pclusters.pdb`;
		@COM_ynumbers = split(/\s+/, $COM_ycluster);
		
		$dist_x = $COM_xnumbers[0] - $COM_ynumbers[0];
		$dist_y = $COM_xnumbers[1] - $COM_ynumbers[1];
		$dist_z = $COM_xnumbers[2] - $COM_ynumbers[2];

		$dist = sqrt(($dist_x) ** 2 + ($dist_y) ** 2 + ($dist_z) ** 2);

		if ($dist <= $user_dist)
		{
			system "echo $line_pclusters >> $list_pclusters\_$user_dist";
		}
	}
}
