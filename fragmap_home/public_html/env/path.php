<?php
require_once 'appvars.php';

$path=getenv('PATH');
putenv("PATH=$path:$rootdir/bin:$rootdir/bin/phplibbin:/usr/bin:$rootdir/bin/atlas_parameterization-1.0.2/bin");
$perllib=getenv('PERLLIB');
putenv("PERLLIB=$perllib:$rootdir/bin");
