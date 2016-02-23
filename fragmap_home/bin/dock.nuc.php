#!/usr/bin/php -q
<?php
require 'env/appvars.php';
require 'env/path.php';
require 'env/dbvars.php';

if ($argc < 2)
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

//nucleic acid mode
if ( isset( $_POST['nucleic_acid']) ) {
	$job->nucleic_acid = 't';
}

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

//check if chains are legit
pdbchains($job, $recorigpre, 'protchains');

// check if protein has residues numbered over 10000
$pwd = getcwd();
chdir($iddir);
$cmd = "grep -E '^ATOM' $recorig | cut -c 23 | awk 'NF>0' | wc -l";
echo $cmd, "\n";
$output = array();
exec($cmd, $output);
print_r($output);
chdir($pwd);
if ($output[0] > 0) {
	PDB::join($iddir, $recorig, $_POST['protchains']);

	$pwd = getcwd();
	chdir($iddir);
	$cmd = "pdb.renumberResidues_wInsert_and_break.pl $recorig $recorig.tmp";
	echo $cmd, "\n";
	system($cmd);
	rename("$recorig.tmp", $recorig);
	chdir($pwd);

	PDB::prep($iddir, $recorig);
}

// pdbchm
$job->status = 'l.m';

if ($job->nucleic_acid == 't') {
		  PDB::nmin($iddir, $recorig, $_POST['protchains'], null, PRMS_DIR.'/beta/pdbamino.rtf', PRMS_DIR.'/beta/parm.prm');
} else {
		  App_PDB::nmin($iddir, $recorig, $_POST['protchains']);
}

chdir($iddir) or die("couldn't chdir $iddir\n");
if ((! file_exists("$iddir/${recorigpre}_nmin.pdb")) ||
       filesize("$iddir/${recorigpre}_nmin.pdb") === 0)
{
	$cmd = "grep -h 'was not found' $iddir/*.inp.out";
	$output = array();//reset so output only contains stuff from this command
	exec($cmd, $output, $status);
	if ($status == 0) {
		$cmd = "grep -h 'was not found' $iddir/*.inp.out | awk '{print \$7}'";
		echo $cmd, "\n";
   	$output = system($cmd);
		$job->le("Unknown residue $output' in receptor.  Please remove.");
	}
	$job->le('Processing failed on receptor');
}
if (! copy("${recorigpre}_nmin.pdb", "$rec") )
{
	$job->le('copy (protorig, prot) failed');
}
//copy($recorig, $rec);
// create images of processed rec for jobdetail.php
$cmd = "pydrawpdb.pl $rec prot";
$ret2 = system($cmd);

if (
		$ret0 != 0 ||
		$ret2 != 0 
	)
{
	$job->le('jobdetail image create failed');
}

mkdir($job->prms_dir()) or die("couldn't mkdir " . $job->prms_dir() );
copy(PRMS_DIR . '/beta/' . DEFAULT_ATOMS , $job->prms_dir() . '/' . DEFAULT_ATOMS );
copy(PRMS_DIR . '/' . DEFAULT_COEFFS, $job->prms_dir() . '/' . DEFAULT_COEFFS );
copy(PRMS_DIR . '/' . DEFAULT_ROTS  , $job->prms_dir() . '/' . DEFAULT_ROTS );
copy(PRMS_DIR . '/beta/pdbamino.rtf'     , $job->prms_dir() . '/pdbamino.rtf' );
copy(PRMS_DIR . '/beta/parm.prm'         , $job->prms_dir() . '/parm.prm' );
$job->atomsfile  = 'atoms_new.ats';
$job->coeffsfile = DEFAULT_COEFFS;
$job->rotsfile   = DEFAULT_ROTS;

if ($job->nucleic_acid == 't') {
    $job->coeffsfile = 'coeffs.dna';
    copy(PRMS_DIR . '/beta/coeffs.dna', $job->prms_dir() . '/' . $job->coeffsfile );
}

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
$jobcmd = "jobrunner.nuc.php $job->jobid";
exec("$jobcmd 1>$iddir/jobrunner.log &");
//system("$jobcmd 2>&1");

function pdbchains($job, $name, $chains)
{
  $dir = $job->id_dir();
  $files = glob("$dir/$name*");

  foreach ($files AS $file)
  {
    $file = basename($file);
    echo "$file\n";
    $parts = split('-', $file);
    if (isset($parts[1]))
    {
      $parts = split('\.', $parts[1]);
      $goodchains[] = $parts[0];
    }
  }

  $postedchains = preg_split('/\s+/', trim($_POST[$chains]), -1, PREG_SPLIT_NO_EMPTY);

  foreach($postedchains AS $chain)
  {
    $chain = strtolower($chain);
    if (!in_array($chain, $goodchains) )
    {

      $job->le("Chain $chain not found");
    }
  }
}
