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

<h2>Map</h2>
<!--<p style="margin: 1em auto;padding:0; width;70%"><strong><font color='red'>Warning: </font></strong>Currently there is a <a href="http://www.bu.edu/tech/about/research/computation/scc/updates/">problem with our cluster</a> which is causing jobs to fail. Please submit jobs after this is resolved.</p>-->
<!--<p style="margin: 1em auto;padding:0; width;70%"><strong><font color='red'>Warning: </font></strong>[Wed Jun 11 11:15:51 EDT] Currently our computing cluster is unstable. Please submit jobs after this issue is resolved. </p>-->

<?php
if ( 'piper' === $liuser->username) {
    echo '<h4 id="embargo">Note: all jobs by non logged in users will be publicly accessible.  Please
              create an account if data is embargoed and needs to remain confidential</h4>';
}

?>
<form id="jobsform" enctype='multipart/form-data' action='home.php' method='post'>
<input type="hidden" name="advanced" id="advanced" value="<?php echo $form['advanced']; ?>" />
<input type="hidden" name="useprotpdbid" id="useprotpdbid" value="<?php echo $form['useprotpdbid']; ?>" />

<!-- job name -->
<div>
   <label for="jobname">Job Name:</label>
   <input type="text" name="jobname" value="<?php echo $form['jobname']; ?>" />
</div>

<div style="padding-top:1.5em;">
Accepted PDB Input:<br />
20 standard amino acids
</div>

   <h3>Protein</h3>
   <div>
      <input type="file" name="prot" id="prot" size="10" />
      <span id="protpdbcode"><label for="protpdb">PDB ID:</label>
      <input type="text" name="protpdb" id="protpdb" size="4" value="<?php echo $form['protpdb']; ?>" />
      </span>
   </div>
   <div>
      <span id="showprotpdb" class="link">Use PDB ID</span><span id="showprotfile" class="link">Upload PDB</span>
   </div>
   <div>
      <label for="protchains">Chains:</label> 
      <input type="text" name="protchains" id="protchains" value="<?php echo $form['protchains']; ?>" />
   </div>
<div>
     Whitespace separate desired chains. Leave chains blank to use all chains.
</div>
<h3 id="advancedtoggle"><img src="image/triangle.down.gif" alt="" /> Advanced Options</h3>
<div id="advancedoptions" class="yui-g">
   <div style="width:105%;">
      <label for="protmask" style="width:100px;">Protein Mask:</label>
      <input type="file" name="protmask" id="protmask" size="10" />
   </div>
<?php if ($liuser->isPrivileged() )
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
?>
<?php if ($liuser->isPrivileged() )
			{ ?>
	<div style="width:70%">
		<label for="probeset">Probe Set:</label>
		<select id="probeset" name="probeset">
		<?php
			$sets = array('orig'=>'Original', 'aa'=>'Peptides', '1c2a' => 'thrombin', '1172' => 'urokinase', '1173' => '1173', '1174' => '1174');
			foreach ($sets AS $set=>$name)
			{
				echo "<option value='$set'";
				echo ($form['probeset'] == $set) ? ' selected="selected"' : '';
				echo ">$name</option>";
			}
		?>
		</select>
	</div>
	<div class="doublespan">
		<input type="checkbox" name="skipcharmm" id="skipcharmm" <?php echo ($form['skipcharmm']) ? 'checked="checked"' : ''; ?> value="checked" />
		<label for="skipcharmm">Skip Charmm</label>
	</div>
   	<div style="width:105%">
      		<label for="coeff_file" style="width:110px;">Coefficients File:</label>
         	<input type="file" name="coeff_file" id="coeff_file" size="10" />
   	</div>
	<div class="doublespan">
		<input type="checkbox" name="keep_metals" id="keep_metals" <?php echo ($form['keep_metals']) ? 'checked="checked"' : ''; ?> value="checked" />
		<label for="keep_metals">Keep metals and heme</label>
	</div>
<?php } else { ?>
<input type="hidden" name="probeset" value="orig"  />
<input type="hidden" name="pbmode" value="newtors"  />
<?php }
?>
	<div class="doublespan">
		<input type="checkbox" name="ppimode" id="ppimode" <?php echo ($form['ppimode']) ? 'checked="checked"' : ''; ?> value="checked" />
		<label for="ppimode">PPI Mode (for binding hot spots on protein protein interfaces)</label>
	</div>

    <div class="doublespan">
        <input type="checkbox" name="nucleic_acid" id="nucleic_acid" <?php echo ($form['nucleic_acid']) ? 'checked="checked"' : ''; ?> value="checked" />
        <label for="nucleic_acid">Has Nucleic Acid</label>
    </div>
<br/>
</div>

<!-- dock -->
<div id="submit">
    <input name="action" type="submit" value="Map" />
</div>
</form>

<?php
   $page->footer();
