<?php
$dbname='ftmap';
$dbuser='ftmap';

$dbconn = pg_connect("dbname=$dbname user=$dbuser") or die('Could not connect: ' . pg_last_error());
$dbh = new PDO("pgsql:dbname=$dbname", $dbuser);
$dbh->setAttribute(PDO::ATTR_ERRMODE, PDO::ERRMODE_WARNING);
$dbh->setAttribute(PDO::ATTR_EMULATE_PREPARES, true);
