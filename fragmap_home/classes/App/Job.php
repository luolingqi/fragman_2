<?php
class App_Job extends Job {

   public  $jobid      = 0;
   protected $properties = array();

   public function __construct($id)
   {
      parent::__construct($id);
   }

   public function le($errstring)
   {
		if (strlen($this->callbackurl) > 0) {
			$request = new HttpRequest($this->callbackurl, HTTP_METH_POST);
			$request->addPostFields(array('id'=>$this->jobid, 'status' => 'l.e', 'errstring' => $errstring));
			$request->send();
		}
		parent::le($errstring);
   }

   public function owner()
   {
      return new App_User($this->userid);
   }
   
   public function bzip_html()
   {
      $file = new App_File($this, 'model_bz2');
      if ( $file->readable() )
      {
         return "<a href='file.php?" . $file->http_query() . "'>Download all Models for all Coefficients</a>";
      }
   }
   
   public function model_file_html($model, $coeffi)
   {
      $model  = (int) $model;
      $coeffi = (int) $coeffi;
      $file = new App_File($this, 'model_file', $model, $coeffi);
      if ( $file->readable() )
      {
         return "<a href='file.php?" . $file->http_query() . "'>$model</a>";
      }
      return $model;
   }
   
   public function model_img_html($model, $coeffi)
   {
      $model  = (int) $model;
      $coeffi = (int) $coeffi;
      $file = new App_File($this, 'model_img', $model, $coeffi);
      if ( $file->readable() )
      {
         return "<img src='file.php?" . $file->http_query() . "' alt='$model' />";
      }
      
      return '(image not found)';
   }

   public function user_file_html($prot)
   {
      $filetype = "user_${prot}_file";
      $file = new App_File($this, $filetype);
      return "<a href='file.php?" . $file->http_query() . "'>pdb</a>";
   }
   
   public function proc_file_html($prot)
   {
      $filetype = "proc_${prot}_file";
      $file = new App_File($this, $filetype);
      if ( $file->readable() )
      {
         return "<a href='file.php?" . $file->http_query() . "'>pdb</a>";
      }

      return $prot;
   }
   
   public function pdb_img($state, $prot)
   {
      $filetype = "${state}_${prot}_img";
      $file = new App_File($this, $filetype);
      if ( $file->readable() )
      {
         return "<img src='file.php?" . $file->http_query() . "' alt='$state $prot' />";
      }
      return '(image not found)';
   }

	public function probes_w_status($status)
	{
		$query  = "SELECT probe FROM jobsprobes WHERE jobid= '$this->jobid' ";
		$query .= "AND status = '$status'";
		$result = pg_query($query);
		return pg_fetch_all_columns($result);
	}

   public function probes_running()
   {
      $query  = "SELECT COUNT(*) FROM jobsprobes WHERE jobid='$this->jobid' ";
      $query .= "AND status <> 'f.q' AND status IS NOT NULL";
      $result = pg_query($query);
      list($count) = pg_fetch_row($result);
      
      //if there are any probes in the result, then it is running remotely
      return ($count);
   }

   public function probes_errored()
   {
      $query  = "SELECT COUNT(*) FROM jobsprobes WHERE jobid='$this->jobid' ";
      $query .= "AND substring(status from 3 for 1) = 'e'";
      $result = pg_query($query);
      list($count) = pg_fetch_row($result);
      
      //if there are any probes in the result, then it is running remotely
      return ($count);
   }


   public function probes_finished()
   {
      $query  = "SELECT COUNT(*) FROM jobsprobes WHERE jobid='$this->jobid' ";
      $query .= "AND status <> 'm.f' AND substring(status from 3 for 1) <> 'e'";
      $result = pg_query($query);
      list($count) = pg_fetch_row($result);
     
      //only when this count is 0, are the jobs finished 
      return (! $count);
   }

   public function probes_table_html()
   {
      $table = "<table class='nice'><tr><th>Probe</th><th>Status</th></tr>";
      
      $query  = "SELECT p.name, jp.probe, jp.status FROM jobsprobes AS jp LEFT OUTER JOIN probes AS p ";
      $query .= "ON jp.probe = p.probe WHERE jp.jobid ='$this->jobid' ORDER BY p.name, probe";
      $result = pg_query($query);
      while ($probe = pg_fetch_assoc($result))
      {
	      $probe['name'] = is_null($probe['name']) ? $probe['probe'] : $probe['name'];
	      $table .= $this->probe_row_html($probe['name'], $probe['status']);
      }
      $table .= "</table>";
      
      return $table;
   }
   
   public function probe_row_html($name, $status)
   {
      $row = "<tr><td>$name</td><td>";
      switch ($status) {
         case 'f.e':
         case 'c.e':
            $row .= 'Error';
            break;
         case 'f.r':
            $row .= 'Running FFT';
            break;
         case 'f.f':
            $row .= 'Clustering';
            break;
         case 'm.q':
            $row .= 'Queued for Minimization';
            break;
         case 'm.r':
            $row .= 'Running Minimization';
            break;
         case 'm.f':
            $row .= 'Finished';
            break;
         case 'm.e':
            $row .= 'Minimization Error';
            break;
         case 'f.q':
         default:
            $row .= 'Queued';
      }
      $row .= '</td></tr>';
      
      return $row;
   }
}
