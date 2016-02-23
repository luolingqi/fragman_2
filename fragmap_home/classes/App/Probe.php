<?php
class App_Probe  {

   public $job = 0;
   public $probe = '';
   protected $properties = array();

   public function __construct($jobid, $probe)
   {
      $this->job = new App_Job($jobid);
      $this->probe = $probe;

      //check for existence in database
      $e_jid   = pg_escape_string($this->job->jobid);
      $e_probe = pg_escape_string($this->probe);
      $query = "SELECT COUNT(*) FROM jobsprobes WHERE jobid = $e_jid AND probe = '$e_probe'";
      $result = pg_query($query) or $this->job->le('Could not fetch probe');;
      list($exists) = pg_fetch_row($result);
      if (0 == $exists)
      {
         $query = "INSERT INTO jobsprobes (jobid, probe) VALUES ('$e_jid', '$e_probe')";
         pg_query($query) or $this->job->le('Could not create probe');
      }
   }

   public function __set($name, $value)
   {
      $e_value = pg_escape_string($value);
      $e_jid   = pg_escape_string($this->job->jobid);
      $e_probe = pg_escape_string($this->probe);
      $query = "UPDATE jobsprobes SET $name = '$e_value' WHERE jobid = '$e_jid' AND probe = '$e_probe'";
      pg_query($query) or die( 'Query failed: ' . pg_last_error() );
      $this->properties[$name] = $value;
   }

   public function __get($property)
   {
      if (! array_key_exists($property, $this->properties) )
      {
         $e_jid = pg_escape_string($this->job->jobid);
         $e_probe = pg_escape_string($this->probe);
         $query = "SELECT $property FROM jobsprobes WHERE jobid = $e_jid AND probe = '$e_probe'";

         $result = pg_query($query) or die( 'Query failed: ' . pg_last_error() );
         
         list($this->properties[$property]) = pg_fetch_row( $result );
      }
      return $this->properties[$property];
   }
   
   public function touch()
   {
      $e_jid = pg_escape_string($this->job->jobid);
      $e_probe = pg_escape_string($this->probe);
      $query = "UPDATE jobsprobes SET touched=NOW()  WHERE jobid = '$e_jid' AND probe = '$e_probe'";
      pg_query($query) or die( 'Query failed: ' . pg_last_error() );

      $this->job->touch();
   }

}
