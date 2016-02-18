<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'sessionstart.php';

$page = new App_Page( array('name' => 'Papers' ) );
$page->header();
?>
<h3>Publications</h3>

<div class='publications'>

<p align="Justify">
Kozakov D, Grove LE, Hall DR, Bohnuud T, Mottarella SE, Luo L, Xia B, Beglov D, Vajda S.
The FTMap family of web servers for determining and characterizing ligand-binding hot spots of proteins.
 <em>Nature Protocols</em>. 2015 10(5):733-755;&nbsp;&nbsp;<a href='http://structure.bu.edu/sites/default/files/nprot.2015.043.pdf'>pdf</a>
</p>

<p align="Justify">
Ngan CH, Bohnuud T, Mottarella SE, Beglov D, Villar EA, Hall DR, Kozakov D, Vajda S. FTMAP: extended protein mapping with user-selected probe molecules. 
<em>Nucleic Acids Research</em>. 2012 Jul;40(Web Server issue):W271-5.
 <a href='http://structure.bu.edu/sites/default/files/Nucl.%20Acids%20Res.-2012-Ngan-W271-5.pdf'>pdf</a>
</p>

<p align="Justify">
Brenke R, Kozakov D, Chuang GY, Beglov D, Hall D, Landon MR, Mattos C, and Vajda S. Fragment-based identification of druggable 'hot spots' of proteins using Fourier domain correlation techniques.
 <em>Bioinformatics</em>. 2009 Mar 1;&nbsp;&nbsp;<a href='http://structure.bu.edu/sites/default/files/621.pdf'>pdf</a>
</p>

<p align="Justify">
Hall DR, Kozakov D, Vajda S.  2012.  Analysis of Protein Binding Sites by Computational Solvent Mapping. <em>Methods in Molecular Biology</em>. 819:13-27.
<a href='http://structure.bu.edu/sites/default/files/mapping_mmb.pdf'>pdf</a>
</p>

</div>
<?php
$page->footer();
