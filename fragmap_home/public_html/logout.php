<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'sessionstart.php';

$liuser->isloggedin = 'f';

unset($_SESSION);
session_destroy();

header('Location: login.php');

?>
