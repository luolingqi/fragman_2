<?php
abstract class Job {

   public  $jobid      = 0;
   protected $properties = array();
   
   abstract public function owner();

   public function __construct($id)
   {
      $this->jobid = (int) $id;
   }

   public function __set($name, $value)
   {
		global $dbh;
      $query = "UPDATE jobs SET $name = :value WHERE id = :jobid";
		$sth = $dbh->prepare($query);
		$sth->execute(array(':value' => $value, ':jobid' => $this->jobid));
      $this->properties[$name] = $value;
   }

   public function __get($property)
   {
      if (! array_key_exists($property, $this->properties) )
      {
			global $dbh;
         $query = "SELECT $property FROM jobs WHERE id =:jobid";
			$sth = $dbh->prepare($query);
			$sth->execute(array(':jobid' => $this->jobid));
         
         $this->properties[$property] = $sth->fetchColumn();
      }
      return $this->properties[$property];
   }

   //this function allows multiple properties to quickly be loaded
   //to reduce the number of sql calls 
   public function quickPropertyLoad()
   {
      if (func_num_args() === 0)
      {
         exit;
      }
      
		global $dbh;
      $props = func_get_args();
      $fields = implode(', ', $props);
      $query = "SELECT $fields FROM jobs WHERE id=:jobid";
		$sth = $dbh->prepare($query);
		$sth->execute(array(':jobid' => $this->jobid));
      $values = $sth->fetch( PDO::FETCH_ASSOC );
      $this->properties = array_merge($this->properties, $values);
   }

   public function id_dir()
   {
      return STORAGE_DIR."/$this->jobid";
   }

   public function prms_dir()
   {
      return $this->id_dir().'/prms';
   }
   
   //This function enters a local error for a job
   public function le($errstring)
   {
      $this->errstring = $errstring;
      $this->status = 'l.e';
      //email user if notification is on
      if ($this->owner()->emailJobComplete) {
          $mail = App_Mail::instance();
          $mail->addTo($this->owner()->email,'');
          $mail->setSubject(SERVER_NAME." job failed");
          $message  = "Your job '$this->jobname' has failed to run on ".SERVER_NAME.".\n";
          $message .= "It failed with an error of: $errstring\n\n";
          $message .= "-".SERVER_NAME." Team\n";
          $message .= "You can turn off e-mail notification at ".WEB_ROOT."preferences.php";
          $mail->setBodyText($message);
          $mail->setFrom(SERVER_EMAIL, $servername);
          $mail->send();
      }
      die("$errstring\n");
   }

   //This function enters a remote error for a job
   public function re($errstring)
   {
      $this->errstring = $errstring;
      $this->status = 'r.e';
      //email user if notification is on
      if ($this->owner()->emailJobComplete) {
          $mail = App_Mail::instance();
          $mail->addTo($this->owner()->email,'');
          $mail->setSubject(SERVER_NAME." job failed");
          $message  = "Your job '$this->jobname' has failed to run on ".SERVER_NAME.".\n";
          $message .= "It failed with an error of: $errstring\n\n";
          $message .= "-".SERVER_NAME." Team\n";
          $message .= "You can turn off e-mail notification at ".WEB_ROOT."preferences.php";
          $mail->setBodyText($message);
          $mail->setFrom(SERVER_EMAIL, $servername);
          $mail->send();
      }
      die("$errstring\n");
   }
   
   public function touch()
   {
		global $dbh;
      $query = "UPDATE jobs SET touched=NOW()  WHERE id =:jobid";
		$sth = $dbh->prepare($query);
		$sth->execute(array(':jobid' => $this->jobid));
   }

   /**
    * generates medium length status from 3 letter status
    * @param string $status
    * @return long status
    */
   public static function convertMedStatus($status)
   {
		switch ($status) {
		   case 'r.e':
			   return "error on supercomputer";
		      break;
		   default:
			   return self::convertLongStatus($status);
		}
   }

   /**
    * generates long status from 3 letter status
    * @param string $status
    * @return long status
    */
   public static function convertLongStatus($status)
   {
		switch ($status) {
            case 'l.p':
                return "fragment parameterization";
		   case 'l.e':
			   return "error on local system";
		      break;
		   case 'l.d':
			   return "deleted";
		      break;
		   case 'r.e':
			   return "error on supercomputer";
		      break;
		   case 'l.g':
		   	return "processing pdb files";
		      break;
		   case 'l.m':
		   	return "pre-docking minimization";
		      break;
		   case 'l.f':
		   	return "finished";
		      break;
		   case 'r.i':
		   	return "held on supercomputer";
		      break;
		   case 'r.q':
		   	return "in queue on supercomputer";
		      break;
		   case 'r.r':
		   	return "running on supercomputer";
		      break;
		   case 'l.r':
		   	return "copying to supercomputer";
		      break;
		   case 'r.o':
		   	return "copying to local computer";
		      break;
		   case 'l.c':
		   case 'l.i':
		   case 'r.c':
		   	return "clustering and minimization";
		      break;
		   default:
			   return $status;
		}
   }
}
