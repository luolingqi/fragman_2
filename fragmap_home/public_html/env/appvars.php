<?php
$servername='FragMap';
define('SERVER_NAME', $servername);
define('SERVER_TITLE', SERVER_NAME.': A Small Fragment Molecule Mapping Server');
define('WEB_ROOT'   , 'http://127.0.0.1/~fragmap_home/');
define('SERVER_EMAIL', 'bufragmap@gmail.com');

# directories
$rootdir='/vagrant/fragmap_home';
define('ROOT_DIR', $rootdir);
define('CLASS_DIR', "$rootdir/classes");
$storagedir = "/home/vagrant/fragmap_data";
define('STORAGE_DIR', $storagedir);
$prmsdir = "$rootdir/prms";
define('PRMS_DIR', $prmsdir);
define('BIN_DIR', "$rootdir/bin");
define('SCV_USER', 'gychuang');
define('SCV_PASS', '223CluaP');
define('SCV_STORAGE', '/projectnb/docking/marsluo/fragmap/jobs');
define('SCV_BIN', '/projectnb/docking/marsluo/atlas/atlas-1.5.4-alpha.2/atlas_package/python_bin');

define('DEFAULT_ATOMS' , 'atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec');
define('DEFAULT_COEFFS', 'coeffs.0.0.6');
define('PPI_COEFFS'    , 'ppi_coeffs.0.0.6');
define('DEFAULT_ROTS'  , 'rot70k.0.0.4.prm');


# filenames
# original user submitted files
$recorigpre = 'ftmaporig';
//$recorig = $recorigpre; // temporarily just name without extension
$fragorigpre = 'fragorig';
//$fragorig = $fragorigpre; // temporarily just name without extension


# processed files
$rec='prot.pdb';

# form fields
$recupfile='ftmap';
$fragmentupfile='fraglist';
#$recpdb='protpdb';

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
