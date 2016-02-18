<?php
abstract class File {
   
   abstract public function fullpath();

   abstract public function content_type();
   
   abstract public function gzipped();
   
   abstract public function http_query($sep = '&amp;');
      
   public function readable()
   {
      return is_readable( $this->fullpath() );
   }

   public function send()
   {
      if (! $this->readable() ) {
         header('HTTP/1.0 403 Forbidden');
         echo 'File Not Readable';
         return false;
      }
      
      
      header('Content-Type: ' . $this->content_type());
      header('Content-Length: ' . filesize( $this->fullpath() ) );
      header('Content-Disposition: attachment; filename=' . basename( $this->fullpath(), '.gz' ) );
      if ($this->gzipped()) {
         header('Content-Encoding: gzip');
      }

      if ( filesize( $this->fullpath() ) < 50000000  ) {
           readfile($this->fullpath());
      } else {
           $out = fopen($this->fullpath(), 'rb');
           while ( !feof($out) ) {
               echo fread($out, 50000000);
           }
           fclose($out);
      }
      return true;
   }
}
