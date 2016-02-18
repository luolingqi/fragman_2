<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'sessionstart.php';

$page = new App_Page( array('name' => 'Examples' ) );
$page->header();
?>
<h2> Examples</h2>

<div class='examples'>
<p align="Justify">

<table class='results' border=0 cellpadding=5 align=center>

<tr>
<th scope="row" style="text-align:left">Angiotensin Converting Enzyme (PDB 1o8a)</th>
<td><a href="examples/ace_example.pdb">PDB File</a></td>
<td><a href="examples/ace_example.pse">Pymol Session</a></td>
<td><a href="show_example.php?example=ace">View Online</a></td>
</tr>
<tr>
<th scope="row" style="text-align:left">Carbonic Anhydrase I (PDB 2cab)</th>
<td><a href="examples/ca_example.pdb">PDB File</a></td>
<td><a href="examples/ca_example.pse">Pymol Session</a></td>
<td><a href="show_example.php?example=ca">View Online</a></td>
</tr>
<tr>
<th scope="row" style="text-align:left">Neuraminidase N2 (PDB 1ivg)</th>
<td><a href="examples/neuraminidase_example.pdb">PDB File</a></td>
<td><a href="examples/neuraminidase_example.pse">Pymol Session</a></td>
<td><a href="show_example.php?example=neuraminidase">View Online</a></td>
</tr>
<tr>
<th scope="row" style="text-align:left">Phospholipase A2 (PDB 1bbc)</th>
<td><a href="examples/phospholipase_example.pdb">PDB File</a></td>
<td><a href="examples/phospholipase_example.pse">Pymol Session</a></td>
<td><a href="show_example.php?example=phospholipase">View Online</a></td>
</tr>
<tr>
<th scope="row" style="text-align:left">Urokinase Type Plasminogen Activator (PDB 2o8t)</th>
<td><a href="examples/urokinase_example.pdb">PDB File</a></td>
<td><a href="examples/urokinase_example.pse">Pymol Session</a></td>
<td><a href="show_example.php?example=urokinase">View Online</a></td>
</tr>
</table>

</div>
<?php
$page->footer();
