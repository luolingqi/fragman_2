<?php
abstract class User {

   public  $userid     = 0;
   protected $properties = array();

   public function __construct($id)
   {
      $this->userid = $id;
   }
   
   public function __set($name, $value)
   {
		global $dbh;
      $query = "UPDATE users SET $name=:value WHERE userid=:userid";
		$sth = $dbh->prepare($query);
		$sth->execute(array(':value' => $value, ':userid' => $this->userid));
      $this->properties[$name] = $value;
   }

   public function __get($property)
   {
      if (! array_key_exists($property, $this->properties) )
      {
			global $dbh;
         $query = "SELECT $property FROM users WHERE userid =:userid";
			$sth = $dbh->prepare($query);
			$sth->execute(array(':userid' => $this->userid));
         
         $this->properties[$property] = $sth->fetchColumn();
      }
      return $this->properties[$property];
   }
   
   public function ownsJob($jobid)
   {
      if ($this->isSuperUser() )
      {
         return true;
      } 
      
		global $dbh;
      $query = "SELECT COUNT(*) FROM jobs WHERE id =:jobid AND userid =:userid";
		$sth = $dbh->prepare($query);
		$sth->execute(array(':jobid' => $jobid, ':userid' => $this->userid));
		
      
      $owned = $sth->fetchColumn();
      
      return $owned;
   }

   public function isSuperUser()
   {
      return $this->privilege >= 2;
   }
   
   public function isPrivileged()
   {
      return $this->privilege >= 1;
   }
}
