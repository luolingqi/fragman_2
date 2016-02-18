<?php
require_once 'appvars.php';

$path=getenv('PATH');
putenv("PATH=$path:$rootdir/bin:$rootdir/bin/phplibbin");
$perllib=getenv('PERLLIB');
putenv("PERLLIB=$perllib:$rootdir/bin");
