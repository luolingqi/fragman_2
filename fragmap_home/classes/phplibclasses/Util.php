<?php
class Util {

   //random string generator for things like passwords
   /**
    * generates a random string, with selected letters and numbers
    * for readability
    * @param int $nchars
    * @return string of length $nchars
    */
   public static function randString($nchars)
   {
      $possible = '23456789bcdfghjkmnpqrstvwxyzBCDFGHJKMNPQRSTVWXYZ';
      $string = '';
      $i = 0;
      for ($i =1; $i <= $nchars; $i++)
      { 
      	$string .= substr($possible, mt_rand(0, strlen($possible)-1), 1);
      }
      return $string;
   }

   /**
    * copies a directory from a local source to the remote destination
    * using ssh2 extension
    * @param ssh2_session $c, string $src, string $dest
    */
   public static function scp_send_dir($c, $src, $dest)
   {
      $src  = rtrim($src , '/');
      $dest = rtrim($dest, '/');
   
      ssh2_exec($c, "mkdir -p $dest");

      $cwd = getcwd();

      chdir($src);

      $tar = tempnam('/tmp', 'FT_TAR_');
      $tarbase = basename($tar);

      $cmd = "tar cJf $tar *";
      echo $cmd, "\n";
      system($cmd);
   
      echo "copying $tar to $dest/$tarbase\n";
      #ssh2_scp_send($c, "$tar", "$dest/$tarbase");
      $sftp = ssh2_sftp($c);
      #copy("$tar", "ssh2.sftp://$sftp/scratch/$tarbase");
      copy("$tar", "ssh2.sftp://$sftp/$dest/$tarbase");
   
      echo "untarring $tarbase\n"; 
      #ssh2_exec($c, "cd $dest;tar xf /scratch/$tarbase");
      ssh2_exec($c, "cd $dest;tar xJf /$dest/$tarbase");

      #$cmd = "rm -f /scratch/$tarbase";
      $cmd = "rm -f /$dest/$tarbase";
      echo $cmd, "\n";
      ssh2_exec($c, $cmd);

      unlink($tar);

      chdir($cwd);
   }

   /**
    * copies a directory from a remote source to the local destination
    * using ssh2 extension
    * @param ssh2_session $c, string $src, string $dest
    */
   public static function scp_recv_dir($c, $src, $dest)
   {
      $src  = rtrim($src , '/');
      $dest = rtrim($dest, '/');
   
      system("mkdir -p $dest");

      $cwd = getcwd();

      chdir($dest);

      $tar = 'FT_RECV_TAR_'.self::randString(5);

      $cmd  = "cd $src;";
      $cmd .= "tar czf /scratch/$tar *";
      echo $cmd, "\n";
      ssh2_exec($c, $cmd);
   
      echo "copying $tar to $dest/$tar\n";
      #ssh2_scp_recv($c, "$src/$tar", "/tmp/$tar");
      $sftp = ssh2_sftp($c);
      copy("ssh2.sftp://$sftp/scratch/$tar", "/tmp/$tar");
   
      echo "untarring $tar\n"; 
      chdir($dest);
      system("tar xzf /tmp/$tar");

      echo "deleting $tar\n";
      unlink("/tmp/$tar");

      $cmd = "rm -f /scratch/$tar";
      echo $cmd, "\n";
      ssh2_exec($c, $cmd);

      chdir($cwd);
   }
}
