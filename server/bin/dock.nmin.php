#!/usr/bin/php -q
<?php
require 'env/appvars.php';
require 'env/path.php';
require 'env/dbvars.php';

if ($argc < 3)
{
  echo "usage: $argv[0] id poststring\n";
  exit;
}
echo $argv[2];
echo "\n";

$id = (int) $argv[1];
$job = new App_Job($id);
$_POST = unserialize ($argv[2]);

// set up env
$iddir   = $job->id_dir();
chdir($iddir) or die ("couldn't chdir $iddir\n");

# get rec and lig if given
$recofile = "$iddir/$recorig";
if ($_POST['useprotpdbid'])
{
	PDB::get($_POST[$recpdb], $recofile) or $job->le("$_POST[$recpdb] not found in PDB");
   $job->protname = $_POST[$recpdb] . '.pdb';

	// copy rec, dir "user" should exist
	$src = $recofile;
	$dst = "$iddir/user/" . $_POST[$recpdb] . ".pdb";
	if (! copy($src, $dst) )
	{
		$job->le('Copy of protein failed');
	}
   $job->pdbid = $_POST[$recpdb];
}

// create images of user rec, lig for jobdetail.php
$protname = $job->protname;

chdir("$iddir/user") or die("couldn't chdir $iddir/user\n");
$cmd = "pydrawpdb.pl $protname prot";
$ret0 = system($cmd);

//check for number of models in rec, lig, get first if needed
$cmd1 = "grep -c '^MODEL' $recofile";
echo $cmd1, "\n";
$output1 = system($cmd1);
$cmd2 = "grep -c '^ENDMDL' $recofile";
echo $cmd2, "\n";
$output2 = system($cmd2);
if ( ($output1 > 1) || ($output2 > 1) )
{
	$cmd = "getmodel.pl $recofile 1 $recofile";
	echo $cmd, "\n";
	system($cmd);
//	$job->add_warning('Protein has multiple models, using first.');
}

// pdbprep
$job->status = 'l.p';
PDB::prep($iddir, $recorig);

// pdbnmd
$job->status = 'l.m';
App_PDB::nmin($iddir, $recorig, $_POST['protchains']);

chdir($iddir) or die("couldn't chdir $iddir\n");
if ((! file_exists("$iddir/${recorigpre}_nmin.pdb")) ||
       filesize("$iddir/${recorigpre}_nmin.pdb") === 0)
{
  $cmd = "grep -h 'unknown residue type' $iddir/*.inp.out";
  $output = array();//reset so output only contains stuff from this command
  exec($cmd, $output, $status);
  if ($status == 0) {
    $cmd = "grep -h 'unknown residue type' $iddir/*.inp.out | awk '{print \$4}'";
    echo $cmd, "\n";
   	$output = system($cmd);
    $job->le("Unknown residue $output in receptor.  Please remove.");
  }
  $job->le('Processing failed on receptor');
}
		$ret0 != 0 ||
		$ret2 != 0 
	)
{
	$job->le('jobdetail image create failed');
}

mkdir($job->prms_dir()) or die("couldn't mkdir " . $job->prms_dir() );
copy(PRMS_DIR . '/' . DEFAULT_ATOMS , $job->prms_dir() . '/' . DEFAULT_ATOMS );
copy(PRMS_DIR . '/' . DEFAULT_COEFFS, $job->prms_dir() . '/' . DEFAULT_COEFFS );
copy(PRMS_DIR . '/' . DEFAULT_ROTS  , $job->prms_dir() . '/' . DEFAULT_ROTS );
copy(PRMS_DIR . '/pdbamino.rtf'     , $job->prms_dir() . '/pdbamino.rtf' );
copy(PRMS_DIR . '/parm.prm'         , $job->prms_dir() . '/parm.prm' );
$job->atomsfile  = 'atoms_new.ats';
$job->coeffsfile = DEFAULT_COEFFS;
$job->rotsfile   = DEFAULT_ROTS;

// Advanced Features below
if ( $job->owner()->isPrivileged() && isset($_POST['server']) )
{
   $job->server = $_POST['server'];
}

if ( $job->owner()->isPrivileged() && isset($_POST['pbmode']) )
{
   $job->pbmode = $_POST['pbmode'];
}

if ($job->protmask === 'uploaded')
{
   PDB::glyseealpha("$iddir/mask.pdb");
}

$job->probeset = $_POST['probeset'];

if ( isset( $_POST['skipcharmm']) ) {
	$job->skipcharmm = 't';
}

if ( isset( $_POST['keep_metals']) ) {
	$job->keep_metals = 't';
}

if ( isset( $_POST['ppimode']) ) {
	copy(PRMS_DIR . '/' . PPI_COEFFS, $job->prms_dir() . '/' . PPI_COEFFS );
	$job->coeffsfile = PPI_COEFFS;
}

//uploaded coeffs file
if ($job->uploadcoeff  == 't')
{
	rename($job->id_dir() . '/uploaded_coeff.0.0.4', $job->prms_dir() . '/uploaded_coeff.0.0.4' );
	unlink($job->prms_dir() . '/' . $job->coeffsfile );
	$job->coeffsfile = 'uploaded_coeff.0.0.4';
}

// job is ready to be sent to cluster server
$job->status = 'l.b';
$job->touch();

chdir(ROOT_DIR.'/public_html');
$jobcmd = "jobrunner.php $job->jobid";
exec("$jobcmd 1>$iddir/jobrunner.log &");
//system("$jobcmd 2>&1");
