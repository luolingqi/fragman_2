<?php
require 'env/appvars.php';
require 'sessionstart.php';

$page = new App_Page( array('require_login' => true) );
$page->header();
?>
<h3>Description of coeffs.prm</h3>
<p>This file defines the energy weights and thresholds used
in choosing a docking result from each rotation. Each
line corresponds to a different set of coefficients.</p>

<p>Each line consists of the following:</p>

<ul>
   <li>coefficient set index</li>

   <li>repulsive vdW weight</li>

   <li>attractive vdW weight</li>

   <li>coulombic electrostatics weight</li>

   <li>generalized born approximation
      electrostatics weight</li>

   <li>pairwise potential weight</li>

   <li>repulsive vdW threshold (if repulsive vdW is greater
      than this value, the result at that
      translation will not be considered).</li>

   <li>attractive vdW threshold (if attractive vdW is less
      than this value, the result at that
      translation will not be considered).</li>
</ul>

<p>Example coeffficients file</p>

<?php
$page->footer();
?>
