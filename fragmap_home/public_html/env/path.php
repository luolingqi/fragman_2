<?php
require_once 'appvars.php';

$path=getenv('PATH');
putenv("PATH=$path:$rootdir/bin:$rootdir/bin/phplibbin:/usr/bin");
$perllib=getenv('PERLLIB');
putenv("PERLLIB=$perllib:$rootdir/bin");
