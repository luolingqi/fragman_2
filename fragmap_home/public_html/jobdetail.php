<?php
require 'env/appvars.php';
require 'sessionstart.php';
require 'env/dbvars.php';

$id = (int) $_GET['job'];

if (! $liuser->ownsJob($id) )
{
   header('HTTP/1.0 403 Forbidden');
   echo 'Access to Job Denied';
   exit();
}  

$page = new App_Page( array('require_login' => true, 'name' => 'Results') );
  
$job = new App_Job($id); 

$page->header();

echo "<div class='results'>";

$job->quickPropertyLoad('jobname', 'ftmapfilename', 'fragmentfilename', 'protchains',
                        'errstring', 'status', 'time', 'nrots');

echo "<h3>$job->jobid: $job->jobname</h3>";

echo "<table class='nice' id='jobproperties'>";
echo '<tr><th>Status</th><td>';
switch ($job->status) {
   case 'l.e':
	   echo "error on local system";
      break;
   case 'r.e':
	   echo "error on remote supercomputer (rec+lig may be too large)";
      break;
   case 'l.f':
   	echo "finished (<a href='models.php?job=$id'>View Map</a>)";
      break;
   case 'r.q':
   	echo "in queue on remote supercomputer";
      break;
   case 'r.r':
   	echo "running on remote supercomputer";
      break;
   default:
	   echo $job->status;
}
echo "</td></tr>";

echo "<tr><th>Submitted</th><td>$job->time</td></tr>";
echo "<tr><th>Errors</th><td>";
echo ($job->errstring === '') ? '(none reported)' : $job->errstring;
echo "</td></tr>";
echo "<tr><th>PDB id</th><td>";
echo ($job->pdbid === '') ? '(unspecified)' : $job->pdbid;
echo "</td></tr>";
echo "<tr><th>Chains</th><td>";
echo ($job->protchains === '') ? '(all)' : $job->protchains;
echo "</td></tr>";

echo "</table>";
echo "<br />";
echo "<br />";

echo "<table class='nice'>";

// user input
echo "<tr>";
echo "<th>";
echo "user input";
echo "</th>";

echo "<th>";
echo "processed input";
echo "</th>";
echo "</tr>";

// pdb files
echo "<tr>";
echo "<td>";
echo $job->user_file_html('prot'), '<br />';
echo "</td>";

echo "<td>";  
echo $job->proc_file_html('prot'), '<br />';
echo "</td>";
echo "</tr>";

// images
echo "<tr>";
echo "<td>";
echo $job->pdb_img('user', 'prot');
echo "</td>";

echo "<td>";
echo $job->pdb_img('proc', 'prot');
echo "</td>";
echo "</tr>";

// chains
/*echo "<tr>";
echo "<td>";
echo "chains:";
echo "<br />";
echo ($job->protchains === '') ? '(all)' : $job->protchains;
echo "</td>";
echo "</tr>";*/


echo "</table>";

echo $job->probes_table_html();

echo "</div>";

echo "<br />";
echo "<br />";
echo "<br />";

$page->footer();
