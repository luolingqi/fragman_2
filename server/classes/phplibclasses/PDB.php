<?php
class PDB {

    public static function get($pdbcode, $ofile='')
    {
        // $name should be /\w{4}/
        $pdbcode = escapeshellarg($pdbcode);

        $pdbgetcmd = "pdbget.pl $pdbcode $ofile";
        echo "pdbgetcmd: $pdbgetcmd\n";
        exec($pdbgetcmd, $output, $status);
        return (! $status); //status==0 indicates success
    }
    public static function get2($pdbcode, $ofile='')
    {
        // $name should be /\w{4}/
        return file_put_contents($ofile, file_get_contents("http://www.rcsb.org/pdb/files/".$pdbcode.".pdb"));
        //$pdbgetcmd = "pdbget.pl $pdbcode $ofile";
        //echo "pdbgetcmd: $pdbgetcmd\n";
        //exec($pdbgetcmd, $output, $status);
        //return (! $status); //status==0 indicates success
    }

    public static function prep($dir, $file)
    {
        $file = escapeshellarg($file);
        $pwd = getcwd();
        chdir($dir) or die("couldn't chdir $dir\n");

        //mac os<9 line endings (\r) cause problems, windows
        //appears to be handled seamlessly
        $cmd = "mac2unix -q $file";
        echo "make all unixy: $cmd\n";
        exec($cmd);

        $pdbprepcmd = "pdbprep.pl $file";
        echo "pdbprepcmd: $pdbprepcmd\n";
        exec($pdbprepcmd, $output, $status);
        chdir($pwd);
        return (! $status); //status==0 indicates success
    }

    public static function charmm($dir, $file, $chain, $smod=null, $rtf=null, $prm=null, $charmm=null, $clean=0, $nsteps=null)
    {
        $pwd = getcwd();
        chdir($dir) or die("couldn't chdir $dir\n");

        // chain should have been checked for format
        $chain = trim($chain);
        if (preg_match('/^\s*$/', $chain))
        {
            $chain = '?';
        }

        $chains = explode(' ', $chain);

        $cmd = 'pdbchm.pl ';
        if (! is_null($smod))
        {
            $smod = escapeshellarg($smod);
            $cmd .= "--smod=$smod ";
        }
        if ($clean === 0)
        {
            $cmd .= '--dont-clean ';
        }
        if (! is_null($rtf))
        {
            $rtf  = escapeshellarg($rtf);
            $cmd .= "--rtf=$rtf ";
        }
        if (! is_null($prm))
        {
            $prm  = escapeshellarg($prm);
            $cmd .= "--prm=$prm ";
        }
        if (! is_null($charmm))
        {
            $charmm = escapeshellarg($charmm);
            $cmd   .= "--charmm=$charmm ";
        }
        if (! is_null($nsteps))
        {
            $nsteps = escapeshellarg($nsteps);
            $cmd   .= "--nsteps=$nsteps ";
        }

        $file   = escapeshellarg($file);
        $chains = array_map('escapeshellarg', $chains);
        $chain  = implode(' ', $chains);
        $cmd   .= "$file $chain";
        echo "dir: $dir, pdbchmcmd: $cmd\n";
        exec($cmd, $output, $status);
        chdir($pwd);
        return (! $status); //status==0 indicates success
    }

    public static function glyseealpha($file)
    {
        $file = escapeshellarg($file);
        $pdbgcacmd = "pdbglyseealpha.pl $file $file";
        echo "pdbgcacmd: $pdbgcacmd\n";
        exec($pdbgcacmd, $output, $status);
        return (! $status); //status==0 indicates success
    }

    public static function resiTypeSuffix($file, $suffix)
    {
        $file = escapeshellarg($file);
        $suffix = escapeshellarg($suffix);
        $cmd = "pdb_resitypesuffix.pl $file $suffix $file";
        echo "resitypesuffix: $cmd\n";
        exec($cmd, $output, $status);
        return (! $status); //status==0 indicates success
    }

    public static function get_chain($pdb)
    {
        $client = new SoapClient("http://www.pdb.org/pdb/services/pdbws?wsdl");
        return $client->getChains($pdb);
    }

    public static function get_id_status($pdb)
    {
        $client = new SoapClient("http://www.pdb.org/pdb/services/pdbws?wsdl");
        return $client->getIdStatus($pdb);
    }

    public static function ngen($dir, $file, $chain, $smod=null, $rtf=null, $prm=null, $xplor_psf=false, $clean=0)
    {
        $pwd = getcwd();
        chdir($dir) or die("couldn't chdir $dir\n");

        // chain should have been checked for format
        $chain = trim($chain);
        if (preg_match('/^\s*$/', $chain))
        {
            $chain = '?';
        }

        $chains = explode(' ', $chain);

        $cmd = 'pdbnmd.pl --dont-minimize ';
        if (! is_null($smod))
        {
            $smod = escapeshellarg($smod);
            $cmd .= "--smod=$smod ";
        }
        if (! is_null($rtf))
        {
            $rtf  = escapeshellarg($rtf);
            $cmd .= "--rtf=$rtf ";
        }
        if (! is_null($prm))
        {
            $prm  = escapeshellarg($prm);
            $cmd .= "--prm=$prm ";
        }
        if ($xplor_psf === true)
        {
            $cmd .= '--xplor-psf ';
        }
        if ($clean === 0)
        {
            $cmd .= '--dont-clean ';
        }

        $file   = escapeshellarg($file);
        $chains = array_map('escapeshellarg', $chains);
        $chain  = implode(' ', $chains);
        $cmd   .= "$file $chain";
        echo "dir: $dir, pdbnmdcmd: $cmd\n";
        exec($cmd, $output, $status);
        chdir($pwd);
        return (! $status); //status==0 indicates success
    }

    public static function nmin($dir, $file, $chain, $smod=null, $rtf=null, $prm=null, $xplor_psf=false, $clean=0, $ph=null)
    {
        $pwd = getcwd();
        chdir($dir) or die("couldn't chdir $dir\n");

        // chain should have been checked for format
        $chain = trim($chain);
        if (preg_match('/^\s*$/', $chain))
        {
            $chain = '?';
        }

        $chains = explode(' ', $chain);

        $cmd = 'pdbnmd.pl ';
        if (! is_null($smod))
        {
            $smod = escapeshellarg($smod);
            $cmd .= "--smod=$smod ";
        }
        if (! is_null($rtf))
        {
            $rtf  = escapeshellarg($rtf);
            $cmd .= "--rtf=$rtf ";
        }
        if (! is_null($prm))
        {
            $prm  = escapeshellarg($prm);
            $cmd .= "--prm=$prm ";
        }
        if (! is_null($ph))
        {
            $ph   = escapeshellarg($ph);
            $cmd .= "--ph=$ph ";
        }
        if ($xplor_psf === true)
        {
            $cmd .= '--xplor-psf ';
        }
        if ($clean === 0)
        {
            $cmd .= '--dont-clean ';
        }

        $file   = escapeshellarg($file);
        $chains = array_map('escapeshellarg', $chains);
        $chain  = implode(' ', $chains);
        $cmd   .= "$file $chain";
        echo "dir: $dir, pdbnmdcmd: $cmd\n";
        exec($cmd, $output, $status);
        chdir($pwd);
        return (! $status); //status==0 indicates success
    }
}
