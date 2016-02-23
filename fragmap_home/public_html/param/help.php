<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'sessionstart.php';


$page = new App_Page( array('name' => 'Help', 'js' => $js) );

$page->header();

?>

<h3>Help</h3>

<p align="justify">The user needs to input the <b>Job Name</b>, <b>PDB ID</b>, and <b>Chain IDs</b> that are to be used in computational mapping. The user can also upload his or her own PDB file. If the Chains field is not entered, all PDB chains are used in computational mapping.</p>
<br>
<p align="justify"> Select <b> Original </b> to use the standard set of 16 small molecules, or select <b> Upload Small Molecules </b> to upload small molecules.</p>
<br>
<li> Copy and past up to 10 SMILES strings into a box </li>
<img src="image/SMILES_2.png" width="350" height="450" alt="SMIlES_1"/>
<br>

<br>
<br>
<hr noshade size=3>
<br>
<br>
<h1>Mapping Results</h1>
<br>
<p align="justify"> Under the <b>Examples Tab</b>, see the "<b>Results</b>" section for a simple analysis of the mapping results. The user can inspect the results online (interactive session can be found here: <a href='http://ftmap.bu.edu/param/model/1704f5922a7e3fad/' target="_blank">Sample Results</a>). </p>
<br>
<p align="justify"> The interactive session shows the submitted protein and its <b>Consensus Sites/ Hot Spots</b>. The top hot spot (<b>Site 1</b>) has the highest number of small molecule clusters based only on the original 16 types of small molecules, followed by the second hot spot (<b>Site 2</b>), so on and so forth.</p>
<br>
<img src="image/Help_map_1.png" width="500" height="425" alt="Help_map_1"/>
<br>
<br>
<br>
<br>
<br>
<br>
<p align="justify"> If the user submitted his or her own small molecules, the interactive session also shows, for each small molecule, the lowest energy cluster representative of each of its conformers. E.g. If the user's submitted Probe 1 has 3 conformers, then, the lowest energy cluster representative of each of Probe 1's 3 conformers are displayed in the interactive session.</p>
<br>
<p align="justify"> <b>NOTE 1</b>: The user can also download the PyMol session; the session shows the top 6 lowest energy cluster representatives for each of the conformers for each of the submitted small molecules.</p>
<br> 
<p align="justify"> <b>NOTE 2</b>: Only cluster representatives in a distance of 4 angstroms to at least one consensus site or hot spot are retained in the final results  <b>i.e. if small molecules submitted by the user bind in locations other than the hot spots identified by the 16 original small molecules, <u>the user's molecules are not shown</u></b></p> 
<br>
<img src="image/Help_map_2.png" width="500" height="425" alt="Help_map_2"/>
<br>


<hr noshade size=3>
<br>
<br>
<h1>Parameterization</h1>
<br>
<p align="justify">Small molecule parameterization can be run independently of computational mapping, and the field accepts a list of small molecules in the format:
<br>
<br>
<p>
<kbd>
#&nbsp&nbsp&nbsp&nbsp&nbspSMILES<br>
#&nbsp&nbsp&nbsp&nbsp&nbspSMILES<br>
#&nbsp&nbsp&nbsp&nbsp&nbspSMILES<br></kbd>
</p>
<br>
<p align="justify"> The sign <b>#</b> is the formal charge of the molecule and the <b>SMILES</b> string must be a correct SMILES string representing a molecule to be used in computational mapping. Enter this text into the field and click <b> Submit</b>. These instructions are reprinted on the <b>Parameterization Tab</b> for the user's reference. Click the <b>Run Example</b> button to run parameterization on one example SMILES string.
<br>
<br>
The server can only generate parameters for a number of small molecules. First, each submission is <b> limited to 10 SMILES strings</b>.  Second, <b> each molecule submitted is limited to 99 conformations</b>.  Should a single molecule exceed 99 conformers, that molecule and its conformers are <b> omitted </b> in computational mapping.
<br>
<br>
The generation of parameters can be run independently of computational mapping by selecting the <b> Parameterization Tab</b> on the FTMap homepage. Users can access and use this function to generate force field parameters to be used in compatible applications. 
</p>
<br>
<br>

<h3><a name="param_results">Parameterization Results</a></h3>
<br>
<p align="Justify">Once the user clicks <b> Submit</b>, parameterization of the small molecules (and conformers, if present) begins. The user is automatically forwarded to the results page when the calculations are done.</b>
<br>
<br>
<p align="Justify">On the results page, the user can download the parameters for inspection or editing. The user can then TAR the edited files and upload the TAR file on the FTMap homepage for use in computational mapping.</p> 
<br>
<br>
<br>
<br>
<p align="justify"> *** *** <br>The downloaded TAR file unpacks into a "<b> Probes </b>" folder:</p>
<br>
<div align="justify" style="margin-left:auto; margin-right:25px">
<ol>
<li>Each small molecule is given an alphanumeric name (e.g. A11), and its own folder containing its topology file (<b> pdbamino_new.rtf </b>) and its parameters file (<b> parm_new.prm </b>).</li>
<br>
<li>The user can upload at most 10 small molecules. The first molecule is assigned a name beginning with "<b>A</b>", and the second beginning with "<b>B</b>", so on and so forth up to and including "<b>J</b>".</li>
<br>
<li>Each small molecule can generate at most 99 conformers. The first conformer for the first small molecule is named "<b>A01</b>, the second conformer "<b>A02</b>", so on and so forth. For the second small molecule, the first conformer is named "<b>B01</b>", the second conformer is named "<b>B02</b>", so on and so forth.</li>
<br>
<li>The following is the <b>pdbamino_new.rtf</b> topology file for the conformer "<b>A01</b>" circled in red, and its atom types are also circled in red. The user can change these parameters in the <b>pdbamino_new.rtf</b> file, as well as the parameters corresponding to the circled atoms types in the <b>parm_new.prm</b> file. It is the user's responsiblity to ensure physical sanity of the changes.</li>
<br>
<br>
<img src=image/Parameters_c.png width="385" height="720" alt="Parameters_c">

</ol>
</div>

</div>
<?php
   $page->footer();
