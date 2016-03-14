<?php
class App_File extends File {

   private $job;
   private $filetype;
   private $model;
   private $coeffi;

   public function __construct(Job $job, $filetype, $model = null, $coeffi = null)
   {
      $this->job = $job;
      $this->filetype = $filetype;
      $this->model  = (int) $model;
      $this->coeffi = (int) $coeffi;
   }

   public function fullpath()
   {
      switch ($this->filetype) {
         case 'user_prot_file':
            $path = $this->job->jobid . '/user/' . $this->job->protname;
            break;
         case '1sidems_file':
            $path = $this->job->jobid . '/probes/rec/1rec.1sidehp3.ms';
            break;
         case 'prot_psf':
            $path = $this->job->jobid . '/protorig_cmin.psf';
            break;
         case 'user_prot_img':
            $path = $this->job->jobid . '/user/prot.png';
            break;
         case 'proc_prot_file':
            $path = $this->job->jobid . '/prot.pdb';
            break;
         case 'proc_prot_img':
            $path = $this->job->jobid . '/prot.png';
            break;
         case 'model_bz2':
            $path = $this->job->jobid . '/cluspro.' . $this->job->jobid . '.tar.bz2';
            break;
         case 'model_file':
            $path = $this->job->jobid . '/fftmap.' . $this->job->jobid . '.pdb';
            break;
         case 'model_pse':
            $path = $this->job->jobid . '/fftmap.' . $this->job->jobid . '.pse';
            break;
         case 'modelb64_pse':
            $path = $this->job->jobid . '/fftmap.' . $this->job->jobid . '.b64.pse';
            break;
         case 'nb_file':
            $path = $this->job->jobid . '/nonbonded.' . $this->job->jobid . '.rawextract';
            break;
         case 'hb_file':
            $path = $this->job->jobid . '/hbonded.' . $this->job->jobid . '.rawextract';
            break;
         case 'summary_file':
            $path = $this->job->jobid . '/probes/cluster/crossclustersummary';
            break;
            /*
         case 'ftsite_zip':
            $path = $this->job->jobid . '/ftsite.zip';
				if (!file_exists(STORAGE_DIR.'/'.$path)) { //this file is not created by default
					$zip = new ZipArchive;
					$zip->open(STORAGE_DIR.'/'.$path, ZIPARCHIVE::OVERWRITE);
					$zip->addFile(STORAGE_DIR.'/'.$this->job->jobid.'/probes/rec/1rec.pdb', 'rec/1rec.pdb');
					$zip->addFile(STORAGE_DIR.'/'.$this->job->jobid.'/probes/rec/1rec.psf', 'rec/1rec.psf');
					$user_rec_pdb = gzfile(STORAGE_DIR.'/'.$this->job->jobid.'/user/'.$this->job->protname.'.gz');
					$user_rec_pdb_string = implode($user_rec_pdb);
					$zip->addFromString('rec/user_rec.pdb', $user_rec_pdb_string);
					exec('ls '.STORAGE_DIR.'/'. $this->job->jobid.'/probes/cluster/1???.min.pdb', $content);
					foreach ($content as $line )
					{
						$zip->addFile($line, 'cluster/'.basename($line));
					}
					exec('ls '.STORAGE_DIR.'/'. $this->job->jobid.'/probes/cluster/crosscluster.*.pdb', $content);
					foreach ($content as $line )
					{
						$zip->addFile($line, 'cluster/'.basename($line));
					}
					$zip->close();
				}
            break;
         case 'ftflex_zip':
            $path = $this->job->jobid . '/ftflex.zip';
				if (!file_exists(STORAGE_DIR.'/'.$path)) { //this file is not created by default
					$zip = new ZipArchive;
					$zip->open(STORAGE_DIR.'/'.$path, ZIPARCHIVE::OVERWRITE);
					$zip->addFile(STORAGE_DIR.'/'.$this->job->jobid.'/probes/rec/1rec.pdb', 'rec/1rec.pdb');
					$zip->addFile(STORAGE_DIR.'/'.$this->job->jobid.'/probes/rec/1rec.1sidehp3.ms', 'rec/1rec.1sidehp3.ms');
					$zip->addFile(STORAGE_DIR.'/'.$this->job->jobid.'/probes/rec/1rec.psf', 'rec/1rec.psf');
					$user_rec_pdb = gzfile(STORAGE_DIR.'/'.$this->job->jobid.'/user/'.$this->job->protname.'.gz');
					$user_rec_pdb_string = implode($user_rec_pdb);
					$zip->addFromString('rec/user_rec.pdb', $user_rec_pdb_string);
					$zip->close();
				}
            break;
         case 'contacts_zip':
             $path = $this->job->jobid . '/contacts.zip';
             if (!file_exists(STORAGE_DIR.'/'.$path)) { //this file is not created by default
                 $zip = new ZipArchive;
                 $zip->open(STORAGE_DIR.'/'.$path, ZIPARCHIVE::OVERWRITE);
                 $zip->addFile(STORAGE_DIR.'/'.$this->job->jobid.'/hbonded.'.$this->job->jobid.'.rawextract', 'hbond.rawextract');
                 $zip->addFile(STORAGE_DIR.'/'.$this->job->jobid.'/nonbonded.'.$this->job->jobid.'.rawextract', 'nonbonded.rawextract');
                 $zip->close();
             }
             break;
         case 'ftdyn_zip':
             $path = $this->job->jobid . '/ftdyn.zip';
             if (!file_exists(STORAGE_DIR.'/'.$path)) { //this file is not created by default
                 $zip = new ZipArchive;
                 $zip->open(STORAGE_DIR.'/'.$path, ZIPARCHIVE::OVERWRITE);
                 $zip->addFile(STORAGE_DIR.'/'.$this->job->jobid.'/probes/rec/1rec.pdb', '1rec.pdb');
                 $zip->addFile(STORAGE_DIR.'/'.$this->job->jobid.'/hbonded.'.$this->job->jobid.'.rawextract', 'hbond.rawextract');
                 $zip->addFile(STORAGE_DIR.'/'.$this->job->jobid.'/nonbonded.'.$this->job->jobid.'.rawextract', 'nonbonded.rawextract');
                 $zip->addFile(STORAGE_DIR.'/'.$this->job->jobid.'/fftmap.'.$this->job->jobid.'.pse', 'ftmap.pse');
                 $zip->addFile(STORAGE_DIR.'/'.$this->job->jobid.'/fftmap.'.$this->job->jobid.'.b64.pse', 'ftmap.b64.pse');
                 $model_pdb = gzfile(STORAGE_DIR.'/'.$this->job->jobid.'/fftmap.'.$this->job->jobid.'.pdb.gz');
                 $model_pdb_string = implode($model_pdb);
                 $zip->addFromString('ftmap.pdb', $model_pdb_string);
                 $zip->close();
             }
             break;*/
         case 'model_img':
            $path = $this->job->jobid . '/model.' . sprintf('%03d', $this->coeffi) . '.' . sprintf('%02d', $this->model) . '.png';
            break;
         case 'model_img_old':
            $path = $this->job->jobid . '/lig' . $this->model . '.png';
            break;
         case 'result_img':
            $path = $this->job->jobid . '/result.png';
            break;
         default:
            return null;
            exit;
      }

      if ($this->gzipped()) {
         //$path .= '.gz';
      }


      return STORAGE_DIR . '/'. $path;
   }

   public function content_type()
   {
      list($null, $ext) = preg_split('/_/', $this->filetype);
      switch ($ext) {
         case 'file':
            $content_type = 'chemical/x-pdb';
            break;
         case 'img':
            $content_type = 'image/png';
            break;
         case 'bz2':
            $content_type = 'application/x-bzip';
            break;
         case 'zip':
            $content_type = 'application/zip';
            break;
         case 'pse':
            $content_type = 'text/plain; charset=x-user-defined';
            break;
         default:
            return null;
            exit;
         }
      return $content_type;
   }

   public function http_query($sep = '&amp;')
   {
      $data = array(
               'jobid'  => $this->job->jobid,
               'coeffi' => $this->coeffi,
               'model'  => $this->model,
               'filetype' => $this->filetype
               );
      return http_build_query($data, '', $sep);
   }

   public function gzipped() {
      switch ($this->filetype) {
         case 'user_prot_file':
         case 'proc_prot_file':
         case 'model_file':
            return true;
            break;
         default:
            return false;
      }
   }
}
