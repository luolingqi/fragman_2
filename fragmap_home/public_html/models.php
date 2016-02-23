<?php
require 'env/appvars.php';
require 'sessionstart.php';
require 'env/path.php';
require 'env/dbvars.php';

$id = (int) $_GET['job'];
$coeffi = isset($_GET['coeffi']) ? $_GET['coeffi'] : 0;

if ( ! isset($liuser) ) {
  header('Location: login.php?redir='.htmlentities($_SERVER['REQUEST_URI']));
  exit();
}
if (! $liuser->ownsJob($id) )
{
  header('HTTP/1.0 403 Forbidden');
  echo 'Access to Job Denied';
  exit();
}

$page = new App_Page( array( 'require_login' => true, 'name' => 'Results' ) );

$page->header();

$job = new App_Job($id);
echo "<h3><a href='jobdetail.php?job=$id'>Job Details</a>: $job->jobname</h3>";

$file = new App_File($job, 'model_file');
$file1 = new App_File($job, 'model_pse');
echo "Download Map: <a id='modelslink' href='file.php?". $file->http_query() . "'>PDB</a>";
echo " or <a id='modelslink' href='file.php?". $file1->http_query() . "'>PyMol session</a>";
echo "<br/>";
?>
<a href='#jsmol_app' id='load_jsmol'>Click to view results using JSmol</a>
<div id="jsmol_app"></div>
<div id="jsmol_control">
  <script src="../js/jsmol/JSmol.min.nojq.js"></script>
  <script src="../js/jsmol/js/JSmolThree.js"></script>
  <script src="../js/jsmol/js/JSmolGLmol.js"></script>

  <!-- Load pse javascript-->
  <script type='text/javascript'>
   $(document).ready(function() {
     Jmol.setDocument(0);

     $("#jsmol_control").hide();
     var info = {
       background: '#000000',
       height: 600,      // pixels (but it may be in percent)
       width: '100%',
       use: 'WebGL HTML5',         // 'HTML5' or 'Java' (case-insensitive)
       j2sPath: '../js/jsmol/j2s',           // only used in the HTML5 modality
       script: 'load file.php?<?php echo "jobid=".$job->jobid."&filetype=modelb64_pse"; ?>; background white;',
       addSelectionOptions: false,
       debug: false,
       coverImage: 'file.php?<?php echo "jobid=".$job->jobid."&filetype=result_img"; ?>',
       coverTitle: 'Click to view model using JSmol',
       deferApplet: true,  // wait to load applet until click
       deferUncover: true, // wait to uncover applet until script completed
       readyFunction: function() {
         $("#load_jsmol").hide();
         $("#jsmol_control").show();
       }
     };

     $('#load_jsmol').click(function(e) {
       info.deferApplet=false;
       $('#jsmol_app').html(Jmol.getAppletHtml('jsmol', info));
     });
     $('#jsmol_app').html(Jmol.getAppletHtml('jsmol', info));

   });
  </script>

  <script language="JavaScript" type="text/javascript">
   document.write('<b>Show/hide clusters:</b> ');
   Jmol.jmolCheckbox('jsmol','', '', 'All', true, 'cc_all');
   <?php
   // read crossclusters for JSMol
   chdir(STORAGE_DIR . "/$job->id/probes/cluster");
   $output = array();
   $cmd = "ls crosscluster.*.pdb";
   exec($cmd, $output);
   $clusters=array();
   foreach ($output AS $crosscluster)
   {
     $parts = explode('.', $crosscluster);
     $clusters[]=$parts[2];
   }

   $jsarg="'cc_all'";
   foreach ($clusters as $i => $size) {
     $cc=sprintf('__crosscluster_%03d_%03d', $i, $size);
     $size=sprintf('%d',$size);
     echo "Jmol.jmolCheckbox('jsmol','select {$cc};wireframe 0.3', 'select {$cc};wireframe off', '$i ($size)', true, 'cc$i');\n";
     $jsarg.=sprintf(",'cc%s'",$i);
   }
   echo "Jmol.jmolSetCheckboxGroup($jsarg);\n";
   ?>
   document.write('<br/><b>Color protein by:</b> ');
   Jmol.jmolRadioGroup('jsmol',[
     ['select protein;color group', 'Direction'],
     ['select protein;color structure', 'Secondary structures'],
     ['select protein;color temperature', 'b-factor'],
     //['select all;color chain', 'Chain']
   ] );
   document.write('<br/>');
   Jmol.jmolCommandInput('jsmol', 'Execute Jmol script', 30);
   document.write('<br/>');
  </script>
  <br/>
  <b>Rotate: </b> Left drag |
  <b>Zoom: </b> Scroll wheel
</div>

<hr/>

<?php
//prepare data array for google chart
//nonbonded
$lines=file("/data/ftmap/$job->id/nonbonded.$job->id.rawextract", FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);
$js_data_nb='["Sequence","non bonded %", {role: "tooltip"}]';
$cmd = "awk -F'\t' '{sum+=$4} END {print sum}' /data/ftmap/$job->id/nonbonded.$job->id.rawextract";
$output = array();
exec($cmd, $output);
$total_hits=$output[0];
foreach ($lines AS $l) {
  $p = explode("\t", $l);
  $percent=sprintf('%.2f',$p[3]/$total_hits*100);
  $data=array( "$p[1]$p[0]", (float) $percent, "$p[1]/$p[0]/$p[2]: $percent%" );
  $js_data_nb.=",\n".json_encode($data);
}
//echo $js_data;

//hbond
if ($job->nucleic_acid != 't') { // hydrogen bond for nucleic acid/probe has not been written yet
  $lines=file("/data/ftmap/$job->id/hbonded.$job->id.rawextract", FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);
  $js_data_hb='["Sequence","Hbond %", {role: "tooltip"}]';
  $cmd = "awk -F'\t' '{sum+=$4} END {print sum}' /data/ftmap/$job->id/hbonded.$job->id.rawextract";
  $output = array();
  exec($cmd, $output);
  $total_hits=$output[0];
  foreach ($lines AS $l) {
    $p = explode("\t", $l);
    $percent=sprintf('%.2f',$p[3]/$total_hits*100);
    $data=array( "$p[1]$p[0]", (float) $percent, "$p[1]/$p[0]/$p[2]: $percent%" );
    $js_data_hb.=",\n".json_encode($data);
  }
}
?>
<!-- Google Charts API for data plotting -->
<script type="text/javascript" src="https://www.google.com/jsapi"></script>
<script type="text/javascript">
 google.load("visualization", "1", {packages:["corechart"]});
 google.setOnLoadCallback(drawChart);
 function drawChart() {
   var data_nb = google.visualization.arrayToDataTable([ <?php echo "$js_data_nb\n"; ?> ]);
   var options_nb = {
     hAxis: {title: 'Sequence', textStyle: {fontSize: 8}}
   };
   var chart_nb = new google.visualization.ColumnChart(document.getElementById('chart_div_nb'));
   chart_nb.draw(data_nb, options_nb);
   <?php if (isset($js_data_hb)) { ?>
   var data_hb = google.visualization.arrayToDataTable([ <?php echo "$js_data_hb\n"; ?> ]);
   var options_hb = {
     hAxis: {title: 'Sequence', textStyle: {fontSize: 8}}
   };
   var chart_hb = new google.visualization.ColumnChart(document.getElementById('chart_div_hb'));
   chart_hb.draw(data_hb, options_hb);
   <?php } ?>
 }
</script>

<?php
$file = new App_File($job, 'nb_file');
echo "<div id='chart_div_nb'></div>";
echo "<a id='nblink' href='file.php?". $file->http_query() . "'>Nonbonded Interactions</a><br/>";

if ($job->nucleic_acid != 't') { // hydrogen bond for nucleic acid/probe has not been written yet
  $file = new App_File($job, 'hb_file');
  echo "<hr/><div id='chart_div_hb'></div>";
  echo "<a id='hblink' href='file.php?". $file->http_query() . "'>Hbond Interactions</a><br/>";
}

$file = new App_File($job, 'summary_file');
echo "<hr/><a id='summarylink' href='file.php?". $file->http_query() . "'>Probes summary</a><br/>";

?>

<?php
$page->footer();
