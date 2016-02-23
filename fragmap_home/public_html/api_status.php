<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'env/path.php';

if ( !isset($_POST) || !isset($_POST['username']) ||
     !isset($_POST['jobid']) || !isset($_POST['sig']) ) {
	header('HTTP/1.0 400 Bad Request');
	echo json_encode( array( 'status' => 'failure', 'errors' => array('Missing Parameters')) );
	exit;
}

//check shared secret
$sub_sig  = $_POST['sig'];
unset($_POST['sig']);
$username = pg_escape_string($_POST['username']);
$jobid = (int) $_POST['jobid'];
$query = "SELECT secret FROM users WHERE username = '$username'";
list($secret) = pg_fetch_row( pg_query($query) );

$keys = array_keys($_POST);

sort($keys, SORT_STRING);
$sig_string = '';
foreach ($keys as $key) {
	$sig_string .= $key . $_POST[$key];
}

$calc_sig = hash_hmac('md5', $sig_string, $secret);

if ($calc_sig !== $sub_sig) {
	header('HTTP/1.0 400 Bad Request');
	echo json_encode( array( 'status' => 'failure', 'errors' => array('Incorrect signiature')) );
	exit;
}

//send thorugh actiondock
$query = "SELECT userid FROM users WHERE username = '$username'";
list($userid) = pg_fetch_row( pg_query($query) );

$query = "
        select userid,status
        from jobs
        where id = $jobid
        ";
$result = pg_query($query) or die('Query failed: ' . pg_last_error());

$line = pg_fetch_array($result, null, PGSQL_ASSOC);
$jobowner=$line['userid'];
$jobstatus=$line['status'];

$query= " select count(*) from jobs where status NOT IN ('l.f', 'l.d', 'l.e', 'r.e') and id<$jobid and id>17342;   ";
$result = pg_query($query) or die('Query failed: ' . pg_last_error());
$line = pg_fetch_array($result, null, PGSQL_ASSOC);
$totalcount=$line['count'];

$liuser = new App_User($userid);

if (! $liuser->ownsJob($jobid) )
{
   header('HTTP/1.0 403 Forbidden');
	echo json_encode( array( 'status' => 'failure', 'errors' => array('Access to Job Denied')) );
   exit;
}

$job = new App_Job($jobid);

echo json_encode( array( 'status' => $job->status, 'count' =>$totalcount ) );
