#!/usr/bin/php -q
<?php
require '/home/ftmap/public_html/env/appvars.php';
require '/home/ftmap/public_html/env/dbvars.php';


$query = "SELECT id FROM jobs WHERE userid NOT IN (SELECT userid FROM users WHERE email LIKE '%@bu.edu') AND status <> 'l.d' AND time < NOW() - interval '3 months' order by time";

$sth = $dbh->query($query);

while ( $line = $sth->fetch(PDO::FETCH_ASSOC) ) {
	$id = $line['id'];
	$job = new App_Job($line['id']);
	$local_dir = $job->id_dir();
	$cmd = 'rm -r '.escapeshellarg($local_dir);
	echo $cmd;
	system($cmd);
	$job->status = 'l.d';
}


