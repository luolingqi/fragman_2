<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'sessionstart.php';

if (! isset($_GET['job']) ||
    ! isset($_GET['coeffi'])|| 
    ! isset($_GET['nmodels']) )
{
   header('HTTP/1.0 404 Not Found');
   echo 'Parameters Missing';
   exit;
}

$jobid   = (int) $_GET['job'];
$coeffi  = (int) $_GET['coeffi'];
$coeffi  = sprintf ("%03d", $coeffi);
$nmodels = (int) $_GET['nmodels'];

if (! $liuser->ownsJob($jobid) )
{
   header('HTTP/1.0 403 Forbidden');
   echo 'Access to File Denied';
   exit;
}

$job = new App_Job($jobid);

$zipname = tempnam('/tmp', 'CPZip_');

$zip = new ZipArchive;
$zip->open($zipname, ZIPARCHIVE::OVERWRITE);

$modeldir = STORAGE_DIR . "/$job->jobid/";

for ($modeli = 0; $modeli < $nmodels; $modeli++)
{
   $model = sprintf ("%02d", $modeli);
   $model = "model.$coeffi.$model.pdb";
   $zip->addFile($modeldir.$model, "cluspro.$job->jobid/$model");
}

$zip->close();

$nmodels = sprintf ("%02d", $nmodels);
$outfilename = "$job->jobid.$coeffi.$nmodels.zip";

header('Content-Type: application/zip');
header('Content-Length: ' . filesize( $zipname ) );
header('Content-Disposition: attachment; filename=' . $outfilename );

$out = fopen($zipname, 'rb');
fpassthru($out);

fclose($out);

unlink($zipname);
