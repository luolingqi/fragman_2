<?php
class App_PDB {

    public static function charmm($dir, $file, $chain) {
        PDB::charmm($dir, $file, $chain, null, PRMS_DIR.'/pdbamino.rtf', PRMS_DIR.'/parm.prm', 'charmm28');
    }

    public static function psf($dir, $file) {
        //PDB::charmm($dir, $file, '?', null, PRMS_DIR.'/pdbamino.rtf', PRMS_DIR.'/parm.prm', 'charmm28', 1, 0);
        PDB::charmm($dir, $file, '?', null, 'pdbamino_new.rtf', 'parm_new.prm', 'charmm28', 1, 0);
    }

    public static function nmin($dir, $file, $chain, $smod='R') {
        PDB::nmin($dir, $file, $chain, $smod, PRMS_DIR.'/pdbamino.rtf', PRMS_DIR.'/parm.prm', true);
    }
}
