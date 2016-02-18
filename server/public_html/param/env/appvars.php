<?php
$servername='FTMap';
define('SERVER_NAME', $servername);
define('SERVER_TITLE', SERVER_NAME.': A Small Molecule Mapping Server');
define('WEB_ROOT'   , 'http://ftmap.bu.edu/');
define('SERVER_EMAIL', 'ftmap@gmail.com');

# directories
$rootdir='/home/ftmap';
define('ROOT_DIR', $rootdir);
define('CLASS_DIR', "$rootdir/classes");
$storagedir = "/data/ftmap";
define('STORAGE_DIR', $storagedir);
$prmsdir = "$rootdir/prms";
define('PRMS_DIR', $prmsdir);
define('BIN_DIR', "$rootdir/bin");
define('SCV_USER', 'gychuang');
define('SCV_PASS', '223CluaP');
define('SCV_STORAGE', '/projectnb/mhcpep/ftmap/jobs');
define('SCV_BIN', '/projectnb/mhcpep/ftmap/bin');

define('DEFAULT_ATOMS' , 'atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec');
define('DEFAULT_COEFFS', 'coeffs.0.0.6');
define('PPI_COEFFS'    , 'ppi_coeffs.0.0.6');
define('DEFAULT_ROTS'  , 'rot70k.0.0.4.prm');

# filenames
# original user submitted files
$recorigpre = 'protorig';
$recorig = $recorigpre . '.pdb';

# processed files
$rec='prot.pdb';

# form fields
$recupfile='prot';
$recpdb='protpdb';

# parameters
$nrots=500;
$chopres=2;
#$nrots=1500;

# website viewing
$jobsperpage=10;

function __autoload($class)
{  
   $path = str_replace('_', DIRECTORY_SEPARATOR, $class) . '.php';
   if ( is_readable(CLASS_DIR . "/$path") )
   {
      require CLASS_DIR . "/$path";
   }
   elseif ( is_readable(CLASS_DIR . "/phplibclasses/$path") )
   {
      require CLASS_DIR . "/phplibclasses/$path";
   }
}
