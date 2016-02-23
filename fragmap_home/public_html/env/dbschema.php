<?php
require("login/isset.l2.php");

/*
   createdb cluspro
*/

/*
	create sequence serial_jobid start 100;
*/

/*
	CREATE TABLE USERS (
		userid serial PRIMARY KEY,
		username varchar(30) UNIQUE,
		password varchar(65),
		email varchar(120) UNIQUE,
		firstname varchar(30),
		lastname varchar(30),
		affiliation varchar(60),
		isloggedin bool DEFAULT '0',
		hasloggedin bool DEFAULT '0',
		timecreated timestamp(0)
	);
*/

/*
	CREATE TABLE JOBS (
		id serial PRIMARY KEY,
		jobname varchar(40),
		userid integer REFERENCES users,
		recname varchar(60) DEFAULT '',
		recchains varchar(20) DEFAULT '',
		ligname varchar(60) DEFAULT '',
		ligchains varchar(20) DEFAULT '',
		errstring varchar(120) DEFAULT '',  -- reported errors for the job
		status varchar(3),
		queueid varchar(40),   -- id in the cluster queue
		nrots int NOT NULL,   -- number of rotations
		time timestamp(0),
		ip cidr
	);
*/

?>
