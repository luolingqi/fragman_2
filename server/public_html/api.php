<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'env/path.php';
//check shared secret
$sub_sig  = $_POST['sig'];
unset($_POST['sig']);
$username = pg_escape_string($_POST['username']);
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
list($_SESSION['userid']) = pg_fetch_row( pg_query($query) );

$env = 'api';

$defaultform = array(
   'server'     => 'scc2', //server option only available to privileged
	'protpdb'     => '',
	'protchains'     => '',
	'jobname'    => '',
   'useprotpdbid'   => '1',
	'pbmode'         => 'newtors', //option only available to privileged
	'probeset'         => 'orig',
);

$_POST = array_merge($defaultform, $_POST);

require 'actiondock.php';

//check for errors
if (isset($errors) && sizeof($errors) > 0) {
	header('HTTP/1.0 400 Bad Request');
	echo json_encode( array( 'status' => 'failure', 'errors' => $errors ) );
} else {
	echo json_encode( array( 'status' => 'success', 'id' => $id ) );
}
