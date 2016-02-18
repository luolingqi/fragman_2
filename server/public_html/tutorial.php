<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'sessionstart.php';

$page = new App_Page( array('name' => 'Tutorial' ) );
$page->header();
?>
<h3>Tutorial</h3>

<iframe src="http://docs.google.com/gview?url=http://structure.bu.edu/sites/default/files/mapping_mmb.pdf&embedded=true" style="width:100%; height:1000px;" frameborder="0"></iframe>

</div>
<?php
$page->footer();
