<?php
//initialize error messages
$errors = array();

//save arrays $_POST['charge'][] and $_POST['smiles'][] before they are killed by trim()
if (isset($_POST['charge']) ) {
	$charge=$_POST['charge'];
	$smiles=$_POST['smiles'];
}
// save user input in session
$_POST = array_map('trim', $_POST);

// get ip
$ip = $_SERVER['REMOTE_ADDR'];

// check the chains
// rec
$lenmax = 20;
if (strlen ($_POST['protchains']) > $lenmax)
{
	$errors[] = "Protein chains must be fewer than $lenmax characters.";
}
if (! preg_match('/^(\w|\?|h\w{3,4})?(\s(\w|\?|h\w{3,4}))*$/i', $_POST['protchains']) )
{
	$errors[] = 'Protein chains must be white space separated alphanumeric characters.';
}

// check pdb id
// rec
if ( $_POST['useprotpdbid'] )
{
 	if (! preg_match('/^\d\w{3}$/', $_POST['protpdb']) )
  	{
  		$errors[] = 'Protein PDB id must be 4 alphanumeric characters.';
  	}
}
if (isset ($_FILES[$recupfile]) &&
	 ($_FILES[$recupfile]['error'] !== UPLOAD_ERR_OK) &&
    (! $_POST['useprotpdbid']) )
{
   switch($_FILES[$recupfile]['error']) {
      case UPLOAD_ERR_INI_SIZE:
      case UPLOAD_ERR_FORM_SIZE:
         $errors[] = 'PDB file too large.';
         break;
      case UPLOAD_ERR_PARTIAL:
         $errors[] = 'PDB file only partially uploaded.  Please try again';
         break;
      case UPLOAD_ERR_NO_FILE:
         $errors[] = 'Please provide a protein file';
         break;
      case UPLOAD_ERR_NO_TMP_DIR:
      case UPLOAD_ERR_CANT_WRITE:
         //todo we should let ourselves know if this happens
         break;
   }
}

if (isset($_POST['pbmode']))
{
	$modes = array('orig', 'new', 'newtors');
	if (! in_array($_POST['pbmode'], $modes) )
	{
		$errors[] = 'Invalid PB mode';
	}
}

$sets = array('orig', 'aa', 'expanded', 'expanded-orig', 'customize', 'upload', 'customize2');
if (! in_array($_POST['probeset'], $sets) ) {
	$errors[] = 'selected:'.$_POST['probeset'].'Invalid Probe Set';
}
if ($_POST['probeset'] === 'upload'){
	if ($_FILES['upload_pdb']['error'] !== UPLOAD_ERR_OK){
		$errors[] = 'Error uploading PDB';
	}
	if (! preg_match('/^[+|-]?\d+$/', $_POST['upload_chg'])){
		$errors[] = 'No charge provided for uploaded structure';
	}
	/*
	if ($_FILES['upload_psf']['error'] !== UPLOAD_ERR_OK){
		$errors[] = 'Error uploading PSF';
	}
	if ($_FILES['upload_xplor']['error'] !== UPLOAD_ERR_OK){
		$errors[] = 'Error uploading XPLOR';
	}
	if ($_FILES['upload_prm']['error'] !== UPLOAD_ERR_OK){
		$errors[] = 'Error uploading PRM';
	}
	if ($_FILES['upload_rtf']['error'] !== UPLOAD_ERR_OK){
		$errors[] = 'Error uploading RTF';
	}
	if ($_FILES['upload_ats']['error'] !== UPLOAD_ERR_OK){
		$errors[] = 'Error uploading ATS';
	}
	 */
}

if (isset($_POST['callbackurl']) ) {
	$callbackurl = filter_var($_POST['callbackurl'], FILTER_VALIDATE_URL);
	if ($callbackurl === false) {
		$errors[] = 'Callback url not valid';
	}
} else {
	$callbackurl = '';
}

// if no error, process page
if (sizeof($errors) === 0)
{

   // insert job into db
   $result = pg_query("select nextval('jobs_id_seq') as key");
   $row = pg_fetch_array($result, 0);
   $id = $row['key'];
   $lenmax = 20;
   $jobname = trim($_POST['jobname']);
   if ( strlen($jobname) > $lenmax)
   	$jobname = substr($jobname, 0, $lenmax);
   elseif ( $jobname === '')
   	$jobname = $id;
   
   $protchains = pg_escape_string($_POST['protchains']);
   
   $userid = pg_escape_string($_SESSION['userid']);
   $callbackurl = pg_escape_string($callbackurl);
   $query = <<<QUERY
   		insert
   		into jobs
   		(id, jobname, userid, status, nrots, time, ip, protchains, touched, callbackurl)
   		values
   		('$id', '$jobname', '$userid', 'l.q', '$nrots', 'now', '$ip', '$protchains', 'now', '$callbackurl')
QUERY;
   pg_query($query) or die('query failed: ' . pg_last_error());
   

   // redirect to job status/result
	if (!isset($env) || $env != 'api') {
		header("Location: queue.php");
	}
   
   $job = new App_Job($id);
   
   // set up env
   $iddir = STORAGE_DIR . "/$id";
   mkdir($iddir) or die("couldn't mkdir $iddir\n");
   mkdir("$iddir/user") or die("couldn't mkdir $iddir/user\n");
   
   // upload rec and lig if given
   $job->status = 'l.g';
   
   $recofile = "$iddir/$recorig";
   if (isset ($_FILES[$recupfile]) && ($_FILES[$recupfile]['error'] === UPLOAD_ERR_OK) &&
       (! $_POST['useprotpdbid']) )
   {
   	// save a copy of users files
   	if ( is_uploaded_file($_FILES[$recupfile]['tmp_name']) )
   	{
   		$job->protname = $_FILES[$recupfile]['name'];
   
   		$src = $_FILES[$recupfile]['tmp_name'];
   		$dst = "$iddir/user/" . $_FILES[$recupfile]['name'];
   		if (! copy ($src, $dst) )
   		{
            $job->le('rec copy failed');
   		}
   	}
   	upload_file($recupfile, $recofile);
   }

   if ( isset( $_FILES['protmask'] ) && $_FILES['protmask']['error'] === UPLOAD_ERR_OK )
   {
      if ( is_uploaded_file($_FILES['protmask']['tmp_name']) )
      {
         $job->protmask = 'uploaded';
      }
      upload_file('protmask', "$iddir/mask.pdb");
   }
   
   if ( isset( $_FILES['coeff_file'] ) && $_FILES['coeff_file']['error'] === UPLOAD_ERR_OK )
   {
      if ( is_uploaded_file($_FILES['coeff_file']['tmp_name']) )
      {
         $job->uploadcoeff = 't';
      }
      upload_file('coeff_file', "$iddir/uploaded_coeff.0.0.4");
   }
   
   if ( isset( $_FILES['prms_file'] ) && $_FILES['prms_file']['error'] === UPLOAD_ERR_OK )
   {
      if ( is_uploaded_file($_FILES['prms_file']['tmp_name']) )
      {
         $job->uploadprms = 't';
      }
      upload_file('prms_file', "$iddir/uploaded_prms.0.0.6");
   }
   
   /*
   if ( $_POST['probeset'] === 'customize' && $_POST['select_input']=='show_upload_tar' && isset( $_FILES['probes_file'] ) && $_FILES['probes_file']['error'] === UPLOAD_ERR_OK )
   {
	   upload_file('probes_file', "$iddir/uploaded_probes.tar");
	   chdir($iddir);
	   $cmd="tar xf uploaded_probes.tar";
	   system($cmd);
	   system("ls -d probes/???? |cut -c 8-11 >probes/probes-ext");
	   chdir(ROOT_DIR."/public_html");

   } elseif ( $_POST['probeset'] === 'customize' && $_POST['select_input']=='show_input_smiles' && $_POST['smiles_box'] !='' && $_POST['smiles_box'] !='#   SMILES' )*/
   if ( $_POST['probeset'] === 'customize' && $_POST['smiles_box'] !='' && $_POST['smiles_box'] !='#   SMILES' )
   {
	$text = trim($_POST["smiles_box"]);
	$fh = fopen("$iddir/SMILES.txt", "w") or die("can't create SMILES.txt");
	fwrite($fh, $text);
	fclose($fh);

	chdir($iddir);

	system("parameterize.pl SMILES.txt >/dev/null 2>&1");
	rename("parameters/probes", "probes");
	system("tar cf parameters.tar probes");
	system("ls -d probes/???? |cut -c 8-11 >probes/probes-ext");
	chdir(ROOT_DIR."/public_html");
   } /*elseif ( $_POST['probeset'] === 'customize' && $_POST['select_input']=='show_input_by_line' && $charge[0] !='' )
   {
	   $text='';
	   for ($i=0;$i<sizeof($charge);$i++){
		   if ($charge[$i]=='' || $smiles[$i]=='') continue;
		   $text.=trim($charge[$i]).' '.trim($smiles[$i])."\n";
		   }
	   $text = trim($text);
	   $fh = fopen("$iddir/SMILES.txt", "w") or die("can't create SMILES.txt");
	   fwrite($fh, $text);
	   fclose($fh);

	   chdir($iddir);

	   system("parameterize.pl SMILES.txt >/dev/null 2>&1");
	   rename("parameters/probes", "probes");
	   system("tar cf parameters.tar probes");
	   system("ls -d probes/???? |cut -c 8-11 >probes/probes-ext");
	   chdir(ROOT_DIR."/public_html");
   }*/
   if ( $_POST['probeset'] === 'customize2' )
   {

	chdir($iddir);
	mkdir("probes", 0755);
	chdir("probes");
	mkdir("rec", 0755);

	$text = trim($_POST["smiles_str"]);
	$fh = fopen("$iddir/probes/SMILES.ism", "w") or die("can't create SMILES.txt");
	fwrite($fh, $text);
	fclose($fh);

	$chg = trim($_POST['smiles_chg']);
	$cmd = "alt_small_molecule_prepare.pl -msur -number 1 -formlchg $chg -hydrogens -genconfs SMILES.ism";
	system("echo $cmd >> param.log");
	system("$cmd >> param.log 2>&1");

	chdir("..");
	system("tar cf parameters.tar probes");
	system("ls -d probes/???? |cut -c 8-11 >probes/probes-ext");
	chdir(ROOT_DIR."/public_html");
   }
   if ( $_POST['probeset'] === 'upload' )
   {
       chdir($iddir);
       mkdir("probes", 0755);
       chdir("probes");
       mkdir("rec", 0755);

       upload_file('upload_pdb', "1lig.pdb");

       $chg = trim($_POST['upload_chg']);
       $cmd = "alt_small_molecule_prepare.pl -msur -number 1 -hydrogens -formlchg $chg ";
       /*
       if(isset($_POST['upload_hs'])){
	       $cmd .= "-hydrogens ";
       }
	*/
       if(isset($_POST['upload_confs'])){
	       $cmd .= "-genconfs ";
       }
       $cmd .= "1lig.pdb";
       system("echo $cmd >> param.log");
       system("$cmd >> param.log 2>&1");

       chdir("..");
       system("ls -d probes/???? |cut -c 8-11 >probes/probes-ext");
       system("tar cf parameters.tar probes");
       chdir(ROOT_DIR."/public_html");
   }

   // exec dock script
   // serialize post array and escape args for command line passing
   $poststring  = escapeshellarg(serialize ($_POST));
   $dockcmd     = "dock.php $id $poststring";
   //echo "$dockcmd<br>";
   //exec ("$dockcmd > /dev/null &");
   exec("$dockcmd 1>$iddir/dock.log &");
   //system ("$dockcmd 2>&1");
   
	if (!isset($env) || $env != 'api') {
		exit;
	}
}
   
function upload_file ($file, $ofile)
{
   if (! move_uploaded_file($_FILES[$file]['tmp_name'], $ofile) )
 	{
 		die ("file upload move failed\n");
 	}
}
