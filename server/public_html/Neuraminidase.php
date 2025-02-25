<?php
	session_start();
	$_SESSION["preloginpage"]="examples.php";
//	require("login/isset.php");

require_once('env/dbvars.php');
require_once('env/appvars.php');

	require("$phpdir/header.php");
?>
<link rel='stylesheet' type='text/css' href='http://ajax.googleapis.com/ajax/libs/jqueryui/1.8.16/themes/overcast/jquery-ui.css' />
<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.6.4/jquery.min.js"></script>
<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jqueryui/1.8.16/jquery-ui.min.js"></script>
<script type="text/javascript" src="examples/example.min.js"></script>
<script type='text/javascript'>
<!--
tabfocus('tabExamples');
$(document).ready(function() {
  $('input').val(['cc0','cc2','lig1']);
  $('label[for="cc0"], label[for="cc2"], label[for="lig1"]').addClass('ui-state-active').attr('aria-pressed', true);
});
//-->
</script>
<style>
.ui-button {margin-bottom:0.2em;}
#applet_controls {margin:auto}
</style>


<div class="main">
<h2 style="margin-top:0;padding-top:16px;">Neuraminidase N2</h2>
<p>Mapping of unbound <a href="http://www.pdb.org/pdb/explore/explore.do?structureId=1ivg">PDB 1ivg</a> with <a href="http://en.wikipedia.org/wiki/Sialic_acid">sialic acid</a> from <a href="http://www.pdb.org/pdb/explore/explore.do?structureId=2bat">PDB 2bat</a> superimposed.</p>
<div id="applet_container" style="width:632px;float:left;padding-left:16px;">
	<applet width="600" height="600" id="av" name="av" code="MoleculeViewerApplet.class" archive="examples/MoleculeViewer.jar" style="padding-bottom:16px;">
		<param name="scriptFile" value="examples/neuraminidase.script">
		<img src="examples/Neuraminidase.png" />
		<param name="java_arguments" value="-Djnlp.packEnabled=true">
	</applet>
	<!-- HTML5 some day
            <object type="application/x-java-applet" id="av" name="av">
		<param name="code" value="MoleculeViewerApplet">
		<param name="archive" value="/MoleculeViewer.jar">
		<param name="scriptFile" value="/molecule_viewer.script">
		<param name="java_arguments" value="-Djnlp.packEnabled=true">
		<img src="site_1.png" />
	</object>-->
</div>
<div style="min-height:616px;">
<?php
for ($i=0; $i <=7; $i++) {
	$j = $i+1;
	echo "<input type='checkbox' id='cc$i' value='cc$i'><label for='cc$i'>Consensus Site $j</label><br>";
}
?>
<input type='checkbox' id='lig1' value='lig1'><label for='lig1'>Sialic Acid</label><br>
<table id="applet_controls">
<tr>
<th scope="row">Left Mouse</th>
<td>Rotate</td>
</tr>

<tr>
<th scope="row"> Shift + Left </th>
<td>Scale</td>
</tr>

<tr>
<th scope="row"> Ctrl + Left </th>
<td>Translate</td>
</tr>

<tr>
<th scope="row">+,&minus;</th>
<td>Clipping</td>
</tr>

</table>
</div>
</div>



<?php
	require("$phpdir/footer.php");
