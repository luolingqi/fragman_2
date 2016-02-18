<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'sessionstart.php';


$page = new App_Page( array('name' => 'Contact', 'js' => $js) );

$page->header();

?>
<h3>Conctact Us</h3>

<p style="width:70%;margin:auto;text-align:left;">Please use the <a href='http://structure.bu.edu/contact'>contact form</a> on our lab page and select "FTMap" as the category.<br/>
If contacting us about a particular job, please include the job number in your message.
</p>
<?php

$page->footer();
