<html>
<head>
<title>FTMap Examples</title>
<meta http-equiv="content-type" content="text/html; charset=utf-8" />
<link rel='stylesheet' type='text/css' href='css/style.css' />
<link rel='stylesheet' type='text/css' href='css/jobsform.css' />
<link rel='stylesheet' type='text/css' href='css/goodform.css' />
<link rel="stylesheet" type="text/css" href="css/grids-min.css" /> 
<link rel="shortcut icon" href="favicon.png" type="image/png" />
<style type='text/css' >
p {width:700px; margin:auto}
</style>
</head>
<body>
<div id="doc" style="text-align:center">
<div id="hd">
<ul id="tabs-menu">
          <li><a id='tabContact' href='contact.php'>Contact</a></li>
	  <li><a id='tabExamples' href='examples.php'>Examples</a></li>
	  <li><a id='tabHelp' href='help.php'>Help</a></li>
          <li><a id='tabPapers' href='publications.php'>Papers</a></li>
	  <li><a id="tabParam" href="parameterization_home.php">Parameterization</a></li>
	  <li><a id="tabFTMap" href="home.php">Map</a></li>
</ul>
<img src="image/banner_ftmap.png" width="750" height="160" alt="FTMap"/>
</div>
<div id="bd">
<br>
<h1>Sample Run</h1>
<hr noshade size=3>
<br>
<br>
<p><h4> Follow Steps 1 to 3, and 7 - 8 to submit a test job using the original 16 small molecules</p> 
<br>
<p>or</p>
<br>
<p> Continue through Step 6, and 7 - 8 to input a new small molecule to be used in mapping</h4></p>
<br>
<br>

<hr noshade size=3>
<br>
<br>

<ol>
<li>Open <a href="home.php" target="_blank">new FTMap session</a>; <b> Original</b> 16 small molecules chosen by default</li>
<img src="image/Help_new_1c.png" width="600" height="450" alt="help_1"/>
<br>
<br>
<br>
<br>

<li>Enter any name in the <b>Job Name</b> field. In this case, the name <b>Test</b> is entered</li>
<br>Enter <b>1HXF</b> into the <b>PDB ID</b> field
<br>Enter &nbsp<b>H  L</b>&nbsp&nbsp into the <b>Chains</b> field
<br>
<br>
<p>*The protein PDB ID 1HXF (Chains H and L) has been chosen for this analysis</p>
<p> but these fields can be replaced by other PDB ID and Chain ID as appropriate</p>
<br>
<img src="image/Help_new_2c.png" width="600" height="450" alt="help_7"/>
<br>
<br>
<br>
<br>

<li> <font color="red"> <b> SKIP TO STEP 4 IF UPLOADING NEW MOLECULES FOR MAPPING&nbsp! </font></b>
<br>
<br>
<p>or</p>
<br>
Click <b> Map </b> to submit the job. The user is then forwarded to the <b>Status Page</b>  </li>
<img src="image/Help_new_3c.png" width="600" height="450" alt="help_7"/>
<br>
<br>
<br>
<br>

<li>Select button <b>Upload Small Molecules</b> to upload molecules in SMILES format</li>
<img src="image/Help_new_4c.png" width="600" height="540" alt="help_2"/>
<br>
<br>
<br>
<br>

<li>Insert Formal Charge: <b>0</b> and SMILES string: <b>NCc1cccc(Cl)c1</b> into the fields</li>
<img src="image/Help_new_5c.png" width="600" height="540" alt="help_3"/>
<br>
<br>
<br>
<br>

<li>Click <b> Map </b> to submit the job. The user is then forwarded to the <b>Status Page</b>  </li>
<img src="image/Help_new_6c.png" width="600" height="540" alt="help_7"/>
<br>
<br>
<br>
<br>

<li><b>Bookmark</b> URL for status updates and results. Status refreshes itself <b>every 30s</b></li>
<br>
<img src="image/help_8_2.png" width="600" height="20" alt="help_8_2"/>
<br>
<img src="image/help_8_1.png" width="600" height="500" alt="help_8_1"/>
<br>
<br>
<br>
<br>

<li>The job is done (<b>Status: finished</b>). Click <b>View Map</b> to view/download results</li>
<br>
<img src="image/help_9.png" width="600" height="480" alt="help_9"/>
<br>

</ol></p>

<br>
<br>
<br>
<br>
<br>
<br>


<hr noshade size=3>
<br>
<h1> Sample Results</h1>
<ol>
<li> <b> Sample Run </b> computational mapping results can be found here : <a href='http://ftmap.bu.edu/param/model/1704f5922a7e3fad/' target="_blank">Sample Run</a> <br><br> The interactive page shows submitted protein and consensus sites / hot spots </li>
<img src="image/Results_new_1c.png" width="450" height="450" alt="Results_new_1c"/>
<br>
<br>
<br>

<li> FTMap generated multiple conformers for the user submitted small molecule</li>
<br>
<p> Only the lowest energy representative of &nbsp"<b>Conformer 1</b>"&nbsp is on display here </p></li>
<img src="image/Results_new_2c.png" width="450" height="450" alt="Results_new_2c"/> 
<br>
<br>
<div align="justify" style="margin-left:auto; margin-right:25px">
<li> Downloadable PyMol session shows the top consensus site (crosscluster.000, in cyan) identifying the most important hot spot, and the lowest energy cluster representative of the first conformer of the user submitted molecule (0101cluster, in purple). The 6 lowest energy cluster representatives are made available for each conformer in the PyMol session.</li>
<br>
<img src="image/Results_new_3_1.png" width="685" height="520" alt="Results_2"/>
<br>
<br>
<img src="image/Results_new_3_2.png" width="685" height="520" alt="Results_3"/>
<br>
<br>
<br>

<li> In the top consensus site, the lowest energy representative of the first conformer of the user submitted molecule (shown as purple sticks) has an almost identical pose to the bound ligand from PDB 2C8Z (shown as cyan sticks). Note that this test case is based on the mapping of the unbound PDB 1HXF.</li>
<img src="image/Results_new_4.png" width="685" height="520" alt="Results_4"/>

</ol>
</div>
<br>
<br>
<br>

</div>
<div id="ft">
  <a href='http://structure.bu.edu/'>Structural Bioinformatics Lab</a>
  <br/>
  <a href='http://www.bu.edu/'>Boston University</a>
</div>
</div>
</body>
</html>
