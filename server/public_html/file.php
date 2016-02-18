<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'sessionstart.php';

if (! isset($_GET['jobid']) ||
    ! isset($_GET['filetype']) )
{
   header('HTTP/1.0 404 Not Found');
   echo 'Parameters Missing';
   exit;
}

$jobid = (int) $_GET['jobid'];

if (! $liuser->ownsJob($jobid) )
{
   header('HTTP/1.0 403 Forbidden');
   echo 'Access to File Denied';
   exit;
}

$job = new App_Job($jobid);

$model  = isset($_GET['model'])  ? $_GET['model']  : null;
$coeffi = isset($_GET['coeffi']) ? $_GET['coeffi'] : null;

$file = new App_File($job, $_GET['filetype'], $model, $coeffi);

if (! $file->readable() ) {
   header('HTTP/1.0 403 Forbidden');
   echo 'File Not Readable';
   exit;
}

$out = fopen($file->fullpath(), 'rb');

header('Content-Type: ' . $file->content_type());
header('Content-Length: ' . filesize( $file->fullpath() ) );
header('Content-Disposition: attachment; filename=' . basename( $file->fullpath() ) );

fpassthru($out);

