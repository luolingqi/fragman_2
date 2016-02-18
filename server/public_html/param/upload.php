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

function cleartext(element){
	if (element.value == "SMILES" || element.value == "0   "){
		element.value = "";
	}
};

function validate(){
	if (document.getElementById("orig").checked) {
		$("#loading").show();
		return true;
	} else if (document.getElementById("upload").checked) {
		var x=$("#upload_chg").val().trim();
		if (x.length == 0){
			alert("You must provide a formal charge for this molecule.");
			return false;
		}
		$("#loading").show();
	} else {
		var x=$("#smiles_str").val().trim();
		if (x.length == 0){
			alert("You must provide a SMILES string.");
			return false;
		}
		var x=$("#smiles_chg").val().trim();
		if (x.length == 0){
			alert("You must provide a formal charge for this molecule.");
			return false;
		}
		$("#loading").show();
	}
}

$(document).ready(function(){
   $("#loading").hide();

   $("#protpdbcode,#showprotfile").hide();

   advanced.init();
   prot.init();

   $("input:radio[name=probeset]")[0].checked=true;
   $("#customize_probes").hide();
   $("#upload_probes").hide();

   $("#showprotpdb").click(prot.showpdb);
   $("#showprotfile").click(prot.showfile);
   $("#advancedtoggle").click(advanced.toggle);

   $("div#upload_tar").hide();
   $("div#input_by_line").hide();

   $("input:radio[name=select_input]")[0].checked=true;
   $("div#input_smiles").hide();
   $("div#upload_probes").hide();
   $("div#input_by_line").show();

   $("#show_input_smiles").click(function() {
	   $("div#upload_probes").hide();
	   $("div#input_by_line").hide();
	   $("div#input_smiles").show();
	   });
   $("#show_upload_probes").click(function() {
	   $("div#input_smiles").hide();
	   $("div#input_by_line").hide();
	   $("div#upload_probes").show();
	   });
   $("#show_input_by_line").click(function() {
	   $("div#input_smiles").hide();
	   $("div#upload_probes").hide();
	   $("div#input_by_line").show();
	   });

   $(function() {
		   $( "#radio" ).buttonset();
		   });

   $("#btnDel").attr("disabled","disabled");

   $("#btnAdd").click(function() {
                var num     = $(".clonedInput").length; // how many "duplicatable" input fields we currently have
                var newNum  = new Number(num + 1);      // the numeric ID of the new input field being added
 
                // create the new element via clone(), and manipulate it"s ID using newNum value
                var newElem = $("#input" + num).clone().attr("id", "input" + newNum);
 
                // manipulate the name/id values of the input inside the new element
                //newElem.children(":first").attr("id", "name" + newNum).attr("name", "name" + newNum);
                newElem.find("#charge").val("");
                newElem.find("#smiles").val("");
 
                // insert the new element after the last "duplicatable" input field
                $("#input" + num).after(newElem);
 
                // enable the "remove" button
                $("#btnDel").attr("disabled","");
                $("#btnDel").removeAttr("disabled");
 
                // business rule: you can only add 5 names
                if (newNum == 10)
                    $("#btnAdd").attr("disabled","disabled");
            });
 
   $("#btnDel").click(function() {
                var num = $(".clonedInput").length; // how many "duplicatable" input fields we currently have
                $("#input" + num).remove();     // remove the last element
 
                // enable the "add" button
                $("#btnAdd").attr("disabled","");
                $("#btnAdd").removeAttr("disabled");
 
                // if only one element remains, disable the "remove" button
                if (num-1 == 1)
                    $("#btnDel").attr("disabled","disabled");
            });

   $( "#dialog1" ).dialog({
			autoOpen: false,
			show: "blind",
			hide: "explode"
		});
   $( "#dialog2" ).dialog({
			autoOpen: false,
			show: "blind",
			hide: "explode"
		});
   $( "#dialog3" ).dialog({
			autoOpen: false,
			show: "blind",
			hide: "explode"
		});
   $( "#dialog4" ).dialog({
			autoOpen: false,
			show: "blind",
			hide: "explode"
		});

   $( "#opener1" ).click(function() {
		   $( "#dialog1" ).dialog( "open" );
		   return false;
		});
   $( "#opener2" ).click(function() {
		   $( "#dialog2" ).dialog( "open" );
		   return false;
		});
   $( "#opener3" ).click(function() {
		   $( "#dialog3" ).dialog( "open" );
		   return false;
		});
   $( "#opener4" ).click(function() {
		   $( "#dialog4" ).dialog( "open" );
		   return false;
		});
 
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
   'uploadprobes'  => '',
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

<h2>UPLOAD TEST</h2>

<form id="jobsform" enctype='multipart/form-data' action='' onsubmit="return validate()" method='post'>
<input type="hidden" name="advanced" id="advanced" value="<?php echo $form['advanced']; ?>" />
<input type="hidden" name="useprotpdbid" id="useprotpdbid" value="<?php echo $form['useprotpdbid']; ?>" />



<div style='width:auto'>
<!-- JOB NAME -->
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
      <span id="showprotpdb" class="link">Use PDB ID</span>
      <span id="showprotfile" class="link">Upload PDB</span>
   </div>
   <div>
      <label for="protchains">Chains:</label> 
      <input type="text" name="protchains" id="protchains" style="" value="<?php echo $form['protchains']; ?>" />
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
</div>

<label style="float:none" for="probeset">Probes:</label>
<label style="float:none"><input type="radio" id="orig" name="probeset" value="orig" onClick='$("#customize_probes").hide();$("#upload_probes").hide()'>Original</label>
<label style="float:none"><input type="radio" id="customize" name="probeset" value="customize2" onClick='$("#customize_probes").show();$("#upload_probes").hide()'>Provide SMILES</label>
<label style="float:none"><input type="radio" id="upload" name="probeset" value="upload" onClick='$("#upload_probes").show();$("#customize_probes").hide()'>Upload Files</label>
</div>

<div id="customize_probes">

<!--
	<div id="radio" style="width: auto;" >
		<input type="radio" id="show_input_by_line" name="select_input" value="show_input_by_line" /><label style="float:none; text-align:center" for="show_input_by_line" >SINGLE</label>
		<input type="radio" id="show_input_smiles"  name="select_input" value="show_input_smiles"  /><label style="float:none; text-align:center" for="show_input_smiles"  >MULTIPLE</label>
		<input type="radio" id="show_upload_tar"    name="select_input" value="show_upload_tar"    /><label style="float:none; text-align:center" for="show_upload_tar"    >TAR FILE</label>
	</div>
-->


   <!-- <div id="input_by_line" style="width:auto">
	   <table border="0" style="margin:auto">
	   <tr>
	   <th>Formal Charge</th>
	   <th>SMILES
	   <input type="button" id="btnAdd" style="width:22px"value="+" />
	   <input type="button" id="btnDel" style="width:22px"value="-" />
	   </th>
	<a href="#" id='opener1'>?</a>
	<div id="dialog1" style='text-align:justify'title="Input SMILES">
		<p>User can insert SMILES of additional small molecules to be used in mapping. Up to 10
small molecules can be added by clicking the "+" sign</p>
	<p>Format: <br/>Formal_charge &nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp SMILES_string</p>
	<p>Example:<br/>
	0  NCc1cccc(Cl)c1
	</p>
	</div>

	   </tr>
	   <tr id="input1" style="margin-bottom:4px;width:auto" class="clonedInput">
	   <td><input type="text" name="charge[]" id="charge" size="2" style="width:auto"/></td>
	   <td><input type="text" name="smiles[]" id="smiles" /></td>
	   </tr>
	   </table>
   </div> -->

   <div id="input_smiles" style='width:auto'>
<!--	<a href="#" id='opener2'>?</a><br/> 
	<div id="dialog2" style='width:auto' title="Input SMILES">
		<p>Provide SMILES string</p>
	</div>
		<label for="smiles_str">SMILES:</label>
		<input type="text" name="smiles_str" id="smiles_str" value="SMILES" style="" onfocus=cleartext(this) />
		<br /><br />
		<label for="smiles_chg">Charge:</label>
		<input type="text" name="smiles_chg" id="smiles_chg" value="0   " style="" onfocus=cleartext(this) />
-->
	<div id="dialog2" style='width:auto' title="Input Multiple SMILES">
		<p>Copy and paste up to 10 SMILES strings in the following format</p>
	
	<p>Format:</p>
	<p>Formal_charge SMILES_string
	<br>Formal_charge SMILES_string
	<br>......</p>
	<p>Example:<br/>
	0  NCc1cccc(Cl)c1
	</p>
	</div>
	   <textarea id='smiles_box' name="smiles_box" style="width:100%" rows="5" cols="100" onfocus=cleartext(this)></textarea>
   </div>
</div>
<div id="upload_probes">
   <div id="upload_probes" style='width:auto'>
        <div id="dialog2" style='width:auto' title="Upload Probe Files">
                <p>Upload probe PDB file</p>
	</div>
	   <div style="width:105%;">
		<label for="upload_pdb">File:</label>
		<input type="file" name="upload_pdb" id="upload_pdb" size="10" />
		<br /><br />
		<label for="upload_chg">Charge:</label>
		<input type="text" name="upload_chg" id="upload_chg" value="0   " style="" onfocus=cleartext(this) />
		<br/><br />
	<!--<input type="checkbox" name="upload_hs" id="upload_hs" checked>Add hydrogens to model<br/>-->
	<input type="checkbox" name="upload_confs" id="upload_confs">Generate Conformers
	</div>
   </div>

   <!-- <div id="upload_tar" style='width:auto;overflow:hidden'>
	<a href="#" id='opener3'>?</a>
	<br/>
	<input type="file" name="probes_file" id="probes_file" size="10" />
	<div id="dialog3" style='text-align:justify'title="Upload a TAR file" >
	<p>First time users: </p>
	<p>1. Generate small molecule parameters by going to the <a href='parameterization_home.php' target='_blanck'><b>Parameterization Tab</b></a></p> 
	<p>2. Download the files to inspect or edit, TAR the files, and upload here</p>
        </div>

    </div>-->


</div>
<br>
<!-- dock -->
<div id="submit">
    <input name="action" type="submit" value="Map" />
</div>
<img id="loading" src="image/loading9.gif"/><br>
</form>

<?php
   $page->footer();
