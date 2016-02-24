<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'sessionstart.php';
require 'env/path.php';
// init dock session vars

if (isset($_POST['action']))
{
	require 'actiondock.php';
}

$js = '
   var advanced = {
      show: function(){
         $("#advancedoptions").show();
         $("#advanced").attr("value",1);
         $("#advancedtoggle > img").attr("src", "image/triangle.down.gif");
      },
      hide: function(){
         $("#advancedoptions").hide();
         $("#advanced").attr("value",0);
         $("#advancedtoggle > img").attr("src", "image/triangle.right.gif");
      },
      toggle: function(){
         if ( "0" === $("#advanced").attr("value") ) {
            advanced.show();
         } else {
            advanced.hide();
         }
      },
      init: function(){
         if ( "0" === $("#advanced").attr("value") ) {
            advanced.hide();
         } else {
            advanced.show();
         }
      } 
   };
   var prot = {
      showpdb: function(){
         $("#protpdbcode,#showprotfile").show();
         $("#useprotpdbid").attr("value",1);
         $("#showprotpdb,#prot").hide();
      },
      showfile: function(){
         $("#protpdbcode,#showprotfile").hide();
         $("#useprotpdbid").attr("value",0);
         $("#showprotpdb,#prot").show();
      },
      init: function(){
         if ( "0" === $("#useprotpdbid").attr("value") ) {
            prot.showfile();
         } else {
            prot.showpdb();
         }
      }
   };

$(document).ready(function(){
   $("#protpdbcode,#showprotfile").hide();

   advanced.init();
   prot.init();

   $("#showprotpdb").click(prot.showpdb);
   $("#showprotfile").click(prot.showfile);
   $("#advancedtoggle").click(advanced.toggle);
   
})';

$page = new App_Page( array( 'name' => 'Map', 'require_login' => true, 'js' => $js) );

$page->header();

$defaultform = array(
   'server'     => 'scc2', //server option only available to privileged
   'protpdb'     => '',
   'protchains'     => '',
   'jobname'    => '',
   'advanced'      => '0',
   'useprotpdbid'   => '1',
   'pbmode'         => 'newtors', //option only available to privileged
   'probeset'         => 'orig',
   'skipcharmm'  => '',
   'ppimode'  => '',
   'nucleic_acid'  => '',
   'keep_metals'  => '',
);

$form = array_merge($defaultform, $_POST);

?>

<?php
if ( isset($errors) && sizeof($errors) > 0)
{
	echo '<div id="errorbox">';
	echo '<h2>Sorry, the job you submitted contains errors:</h2>';
   echo '<ul>';
   foreach ($errors AS $error)
   {
      echo "<li>$error</li>";
   }
   echo '</ul>';
	echo '</div>';
}

?>

<h2>Fragment Binding Prediction</h2>
<!--<p style="margin: 1em auto;padding:0; width;70%"><strong><font color='red'>Warning: </font></strong>Currently there is a <a href="http://www.bu.edu/tech/about/research/computation/scc/updates/">problem with our cluster</a> which is causing jobs to fail. Please submit jobs after this is resolved.</p>-->
<!--<p style="margin: 1em auto;padding:0; width;70%"><strong><font color='red'>Warning: </font></strong>[Wed Jun 11 11:15:51 EDT] Currently our computing cluster is unstable. Please submit jobs after this issue is resolved. </p>-->

<?php
if ( 'piper' === $liuser->username) {
    echo '<h4 id="embargo">Note: all jobs by non logged in users will be publicly accessible.  Please
              create an account if data is embargoed and needs to remain confidential</h4>';
}

?>
<form id="jobsform" enctype='multipart/form-data' action='home.php' method='post'>
<!--<input type="hidden" name="advanced" id="advanced" value="--><?php //echo $form['advanced']; ?><!--" />-->
<!--<input type="hidden" name="useprotpdbid" id="useprotpdbid" value="--><?php //echo $form['useprotpdbid']; ?><!--" />-->

<!-- job name -->
<div>
   <label for="jobname">Job Name:</label>
   <input type="text" name="jobname" value="<?php echo $form['jobname']; ?>" />
</div>


   <h3 style="padding-top:0.5em;">Upload FTMap Result</h3>
	<div>
		Please upload the FTMap output file in the format of either pdb or pse!
	</div>

   <div>
      <input type="file" name="ftmap" id="ftmap" size="10" />

   </div>

   <!--<div>
      <label for="protchains">Chains:</label> 
      <input type="text" name="protchains" id="protchains" value="<?php echo $form['protchains']; ?>" />


	<div>
		Whitespace separate desired chains. Leave chains blank to use all chains.
	</div>
	</div>-->
	<h3 style="padding-top:0.5em;"> Upload Fragments</h3>
	<div>
		<input type="file" name="fraglist" id="fraglist" size="10" />
	</div>
	<div>
		Please upload a file containing the fragments to dock. The compatible chemical formats include: <br />
		sdf, smi, mol2, pdb
	</div>


<!--	<h3 id="advancedtoggle"><img src="image/triangle.down.gif" alt="" /> Advanced Options</h3>
<div id="advancedoptions" class="yui-g">
	<div style="width:105%;">
      <label for="protmask" style="width:100px;">Protein Mask:</label>
      <input type="file" name="protmask" id="protmask" size="10" />
   </div>
	-->
<!--
<?php /* if ($liuser->isPrivileged() )
		{
			$modes = array('orig', 'new', 'newtors');
			echo '<div style="width:60%;">';
			echo '	<label for="pbmode">P-B Mode:</label>';
			echo '	<select id="pbmode" name="pbmode">';
			foreach ($modes AS $mode)
			{
				echo "<option value='$mode'";
				echo ($form['pbmode'] == $mode) ? ' selected="selected"' : '';
				echo ">$mode</option>";
			}
			echo '	</select>';
			echo '</div>';
		}
 */
?>
-->


<!-- dock -->
<div id="submit">
    <input name="action" type="submit" value="Map" />
</div>
</form>

<?php
   $page->footer();
