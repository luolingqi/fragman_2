<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'sessionstart.php';
?>
<a name="contents"></a>
<?php
$page = new App_Page( array('name' => 'Help' ) );
$page->header();
?>
<h3>FTMap Help</h3>

<div class='serverhelp'>
<p>Welcome to FTMap! This page is support for each facet of teh FTMap server. Feel free to use one of the test cases from
the <a href="examples.php">Examples Page</a> or use a protein of your own. The images below use the Angiotensin Converting Enzyme example.</p>

<h3>Contents</h3>
<p><ul>
<li><a href="#login">Login</a></li>
<li><a href="#signup">Sign Up</a></li>
<li><a href="#home">Map</a></li>
<li><a href="#advanced">Advanced Options</a></li>
<li><a href="#queue">Queue</a></li>
<li><a href="#status">Status</a></li>
<li><a href="#status_abbr">Status Abbreviations</a></li>
<li><a href="#results">Results List</a></li>
<li><a href="#models">Results Detail</a></li>
<li><a href="#probes">Probes</a></li>
<li><a href="#preferences">Preferences</a></li>
<li><a href="#publications">Publications</a></li>
<li><a href="#contact">Contact</a></li>
</ul></p>
<br><br>

<a name="login"><h3>Login</h3></a>
<p>When you arrive at the FTMap page, the first screen you see is the login screen, shown here:</p>
<img src="image/login_crop.png" alt="Login Page" />
<p><ul><li>If you already have a username and password, fill in the boxes and click <b>Login</b>. Proceed to <a href="#home">Map</a>.</li><br>
<li>If you prefer not to use an account, click the option below Login to use the server without an ID. Know that your job submissions will be publicly accessible.</li><br>
<li>If you would like to sign up for an account, click the <b>Sign Up</b> option at the bottom of the screen. Proceed to <a href="#signup">Sign Up</a>.</li><br>
<li>If you are a returning user, but have forgotten your password, click <b>Reset Password</b> to have a new password emailed to you</li></ul></p>
<br><br><p><small><a href="#contents">Return to Top</a></small></p>

<a name="signup"><h3>Sign Up</h3></a>
<p>Use the form on this page to create a new user account. Provide your <i>First Name</i>, <i>Last Name</i>, a unique <i>Username</i>, your <i>Affiliation</i>, 
and your <i>Email</i>. Enter the text in the image to <i>Word Verification</i> and click the <i>Agreement</i> checkbox, and then click <b>Create Account</b>. 
If you have entered everything correctly, an email will be sent to you with a password to <a href="#login">Login</a>. 
You can change this password later in <a href="#preferences">Prefereences</a>.</p>
<img src="image/signup_crop.png" alt="Signup Page" />
<br><br><p><small><a href="#contents">Return to Top</a></small></p>

<a name="home"><h3>Map</h3></a>
<p>This is the server home screen. From this page you can submit a new job. See below for descriptions of the information to enter into each box. 
Jump to <a href="advanced">Advanced Options</a> to find further information.</p>
<img src="image/home_crop.png" alt="Home Page" />
<p><ol><li><i>Job Name</i>: [Optional] Provide a job name for this submission. If you choose to leave this blank, a unique ID will be created for this field.</li><br>
<li><i>PDB</i> [Choose one]<ul><li><i>PDB ID</i>: Enter the four digit ID of a PDB or click <i>Upload PDB</i>.</li>
<li><i>Upload PDB</i>: Select <i>Browse</i> to upload a file containing a structure in PDB format or click <i>Use PDB ID</i>.</li></ul></li><br>
<li><i>Chains</i>: Designate specific chains to be used by the server with whitespace separated chain IDs or leave blank for all chains.</li><br>
<li><i>Advanced Options</i>: Click here to see detailed information about <a href="#advanced">Advanced Options</a>.</li><br>
<li><b>Map</b>: Click the Map button when all the fields are filled to your satisfaction to submit your job. Your job will enter the <a href="#queue">Queue</a></li></ol></p>
<br><br><p><small><a href="#contents">Return to Top</a></small></p>

<a name="advanced"><h3>Advanced Options</h3></a>
<p>Click the <i>Advanced Options</i> label on the <a href="#home">Map</a> page to see more options.</p>
<img src="image/advanced_crop.png" alt="Advanced Options" />
<p><ul><li><i>Protein Mask</i>: Upload a file to mask specified atoms as non-interacting. See <a href="http://cluspro.bu.edu/tut_att_rep.php#mask">this guide</a>
for how to create a mask file.</li><br>
<li><i>PPI Mode</i>: Detect binding hot spots for protein-protein interactions using an alternative set of parameters.</li><br>
<li><i>Nucleic Acids</i>: PDB includes nucleic acids. This option uses a specific set of parameters and must be selected for any file containing nucleic acids.</li></ul></p>
<br><br><p><small><a href="#contents">Return to Top</a></small></p>

<a name="queue"><h3>Queue</h3></a>
<p>Use this page to see the status of your submitted jobs. Jobs are run in order of there submission. Runtime will vary by job. 
Click on a job <i>ID</i> to see a detailed <a href="#status">Status Page</a>. See the <a href="#status_abbr">Status Abbreviation</a> table 
for a description of the various job statuses and their meanings.</p>
<img src="image/queue_crop.png" alt="Queue Page" />
<br><br><p><small><a href="#contents">Return to Top</a></small></p>

<a name="status"><h3>Status</h3></a>
<p>This page shows the detailed status of a job that is currently running.</p>
<img src="image/status_crop.png" alt="Status Page" />
<p><ul><li><i>Job ID</i>: Unique job ID number and job name.</li><br>
<li><i>Job Info</i>: Information about this job including status, submission time stamp, and PDB ID. 
Definitions of statuses can be found in <a href="#status_abbr">this table</a>.</li><br>
<li><i>Input</i>: Pictoral representations of the uploaded and processed inputs.</li><br>
<li><i>Probe Status</i>: Status for each small molecule probe. See the <a href="#probe_status_abbr">Probe Status Abbreviation</a> table.</li></ul></p>
<br><br><p><small><a href="#contents">Return to Top</a></small></p>

<a name="status_abbr"><h3>Status Abbreviations</h3></a>
<p>This table shows a list of job statuses and their descriptions that are seen on the <a href="#queue">Queue</a> and <a href="#status">Status</a> pages.</p>
<table border="1" style="width: 100%">
<tr><th>Status</th><th>Description</th></tr>
<tr>
<td>processing pdb files</td>
<td>downloading pdbs from the PDB web site, processing chain information, extracting the chains the user specified</td>
</tr><tr>
<td>pre-docking minimization</td>
<td>running charmm to add missing atoms and polar hydrogens, minimizing the added atoms in the presence of the protein</td>
</tr><tr>
<td>calculating PB potential</td>
<td>using charmm to calculate the Poisson-Boltzmann potential around the protein</td>
</tr><tr>
<td>copying to supercomputer</td>
<td>copying the pdbs and PB potential to the cluster where FTMap will run</td>
</tr><tr>
<td>held on supercomputer</td>
<td>files are on cluster, but job is not yet submitted</td>
</tr><tr>
<td>in queue on supercomputer</td>
<td>jobs for all probes have been submitted on the cluster, but they have not started running</td>
</tr><tr>
<td>running on supercomputer</td>
<td>jobs have begun running on the cluster</td>
</tr><tr>
<td>finished on supercomputer</td>
<td>all probes have run FFT and minimization on the cluster</td>
</tr><tr>
<td>coping to local computer</td>
<td>results are being copied from the cluster back to the FTMap server</td>
</tr><tr>
<td>clustering and minimization</td>
<td>individual probes are clustered, then consensus sites are generated by clustering across probes</td>
</tr><tr>
<td>calculating interactions</td>
<td>calculate nonbonded and hydrogen bonded interactions between probes and the protein using HBPlus</td>
</tr><tr>
<td>finished</td>
<td>everything is complete</td>
</tr>
</table>
<a name="probe_status_abbr"></a>
<br><p>This table shows a list of probe statuses that can be found on the <a href="#status">Status Page</a>.</p>
<table border="1" style="width: 100%">
<tr><th>Status</th><th>Description</th></tr>
<tr>
<td>Queued</td>
<td>waiting in the cluster queue to run FFT</td>
</tr><tr>
<td>Running FFT</td>
<td>probe is running FFT on the cluster</td>
</tr><tr>
<td>Clustering</td>
<td>post-FFT clustering of probes</td>
</tr><tr>
<td>Queued for minimization</td>
<td>waiting in the cluster queue to minimize probes</td>
</tr><tr>
<td>Running Minimization</td>
<td>probe is being minimized on the cluster</td>
</tr><tr>
<td>Finished</td>
<td>probe has finished minimization on the cluster</td>
</tr>
</table>
<br><br><p><small><a href="#contents">Return to Top</a></small></p>

<a name="results"><h3>Results List</h3></a>
<p>See a list of your completed jobs. Click a job <i>ID</i> to see <a href="#models">Result Details</a>. 
Click <i>next-&#62</i> or <i>&#60-prev</i> to see more completed jobs.</p>
<img src="image/results_crop.png" alt="Results Page" />
<br><br><p><small><a href="#contents">Return to Top</a></small></p>

<a name="models"><h3>Results Detail</h3></a>
<p>See detailed information about a completed job.</p>
<img src="image/results_detail_crop.png" alt="Results Detail" />
<p><ul><li><i>Job Details</i>: Job name.</li><br>
<li><i>Download</i>: Download your results in PDB format or as a PyMol Session.</li><br>
<li><a href="#visualization"><i>Visualization</i></a>: Click the image of your results to interact with you uploaded structure using JSmol. 
Allow several seconds for JSmol to load.</li><br>
<li><a href="#graphs"><i>Contact Graphs</i></a>: Bar graphs show contact frequency per residue from FTMap results.</li><br>
<li><i>Probe Summary</i>: Download a text file summary of clusters and their probe compositions.</li></ul></p>
<br>
<a name="visualization"><h4>Visualization</h4></a>
<p>To manipulate your structure, click the image and allow a few seconds for JSmol to load.</p>
<img src="image/visualization_crop.png" alt="Visualization" />
<p>Use the left mouse button and drag to rotate your structure. Use your mouse wheel to zoom in and out.</p>
<p>Use the checkboxes along the bottom to select/deselect specific FTMap clusters. Use the color options to change the structure coloring. 
Execute JSmol scripts from the text field. See the <a href="http://chemapps.stolaf.edu/jmol/jsmol/jsmol.htm">JSmol webpage</a> for further details on JSmol.</p>
<br><br>
<a name="graphs"><h4>Contact Graphs</h4></a>
<img src="image/graphs_crop.png" alt="Contact Graphs" />
<p>The contact graphs display the contact rate of structure residues as a percent of total contacts. One graph is for nonbonded residue interactions, 
the second for hydrogen bond interactions. These results can be downloaded as a tab separated file with exact residue contact counts by click the link 
located directly beneath each respective graph.</p>
<br><br><p><small><a href="#contents">Return to Top</a></small></p>

<a name="probes"><h3>Probes</h3></a>
<table border="1" style="width: 100%">
<tr><th>ID  </th><th>Name</th></tr>
<tr><td>1acd</td><td>acetamide</td></tr>
<tr><td>1acn</td><td>acetonitrile</td></tr>
<tr><td>1act</td><td>acetone</td></tr>
<tr><td>1ady</td><td>acetaldehyde</td></tr>
<tr><td>1amn</td><td>methylamine</td></tr>
<tr><td>1bdy</td><td>benzaldehyde</td></tr>
<tr><td>1ben</td><td>benzene</td></tr>
<tr><td>1but</td><td>isobutanol</td></tr>
<tr><td>1chx</td><td>cyclohexane</td></tr>
<tr><td>1dfo</td><td>N,N-dimethylformamide</td></tr>
<tr><td>1dme</td><td>dimethyl ether</td></tr>
<tr><td>1eol</td><td>ethanol</td></tr>
<tr><td>1eth</td><td>ethane</td></tr>
<tr><td>1phn</td><td>phenol</td></tr>
<tr><td>1ths</td><td>isopropanol</td></tr>
<tr><td>1ure</td><td>urea</td></tr>
</table>
<br><br><p><small><a href="#contents">Return to Top</a></small></p>

<a name="preferences"><h3>Preferences</h3></a>
<p>Set options for your FTMap account.</p>
<img src="image/preferences_crop.png" alt="Preferences Page" />
<p><ul><li>Enter a new password for your account.</li><br>
<li>Select or deselect email notification for completed jobs.</li></ul></p>
<br><br><p><small><a href="#contents">Return to Top</a></small></p>

<a name="publications"><h3>Publications</h3></a>
<p>Please cite these papers in any publications that make use of FTMap.</p>
<img src="image/publications_crop.png" alt="Publications Page" />
<br><br><p><small><a href="#contents">Return to Top</a></small></p>

<a name="contact"><h3>Contact</h3></a>
<p>Feel free to contact us with any additional questions.</p>
<img src="image/contact_crop.png" alt="Contact Page" />
<br><br><p><small><a href="#contents">Return to Top</a></small></p>

</div>
<?php
$page->footer();
