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
echo $iddir;
echo "\n";
chdir($iddir) or die ("couldn't chdir $iddir\n");

# get rec and lig if given
$recofile = glob("$iddir/$recorigpre*")[0];
//$recofile = "$iddir/$recorigpre"; //without extension
/*
if ($_POST['useprotpdbid'])
{
	PDB::get2($_POST[$recpdb], $recofile) or $job->le("$_POST[$recpdb] not found in PDB");
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
*/
// create images of user rec, lig for jobdetail.php
//$protname = $job->protname;
//$ftmapname = $job->ftmapname;

chdir("$iddir/user") or die("couldn't chdir $iddir/user\n");
//$cmd = "pydrawpdb.pl '$protname' prot";
//$ret0 = system($cmd);
// job is undergoing fragment parameterization
$job->status = 'l.p';
$job->touch();
//Paremeterization of fragments
$cmd="atlas_parameterize --gen3D --pH 7.3 $job->fragmentfilename $job->fragmentfilename.json";
$cmd_frag = "python " . BIN_DIR . "/prepare_fragment.py $job->fragmentfilename";
echo $cmd_frag;
echo "\n";
system($cmd_frag . " >> $iddir/dock.log 2>&1",$ret_frag);

//Save and Process the fftmap result file
//Save the protein and consensus cluster objects into individual files

//$cmd_save="sh /vagrant/fragmap_home/bin/save_ftmap.sh";
$cmd_ftmap = "python " . BIN_DIR . "/prepare_ftmap.py fftmap*";
echo "$cmd_ftmap\n";
system($cmd_ftmap . " >> $iddir/dock.log 2>&1", $ret_ftmap);

if (
	$ret_frag != 0 ||
	$ret_ftmap != 0
)
{
	$job->le('fragment parameterization or ftmap preprocessing failed!');
}


// Advanced Features below
if ( $job->owner()->isPrivileged() && isset($_POST['server']) )
{
   $job->server = $_POST['server'];
}

// job is ready to be sent to cluster server
$job->status = 'l.b';
$job->touch();

chdir(ROOT_DIR.'/public_html');
//$jobcmd = "jobrunner.php $job->jobid";
//exec("$jobcmd 1>$iddir/jobrunner.log &");
//system("$jobcmd 2>&1");

/*
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
*/
