#!/usr/bin/php -q
<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'env/path.php';

if ($argc < 2)
{
   echo "usage: $argv[0] jobid\n";
   exit;
}

$jobid = (int) $argv[1];
$job = new App_Job($jobid);
$job->touch();

$remote_dir = SCV_STORAGE . "/$jobid";
$local_dir  = STORAGE_DIR . "/$jobid";
$prms_dir   = "$remote_dir/prms";

//most diffs between servers here
//but there are other diffs further
//down
if ('scc2' == $job->server)
{
   $r_server = 'scc2.bu.edu';
   $r_port = 22;
   $mpirun = '~tanggis/bin/mympirun.scc';
   $np = 16;
   $piper = '~tanggis/bin/piper.0.0.6-0.mpi';
   $status_column = 5;
}
if ('katana' == $job->server)
{
   $r_server = 'katana';
   $r_port = 22;
   $mpirun = '~/bin/mympirun.katana.long';
   $np = 16;
   $piper = '~/ft.bin/piper.katana.mpi.20090601';
   $status_column = 5;
}
elseif ('lee' == $job->server)
{
   $r_server = 'lee';
   $r_port = 22;
   $mpirun = '~/bin/mympirun.bgl.wrapper';
   $np = 32;
   $piper = '~/ft.bin/piper.cavity.bgl.20080327';
	$status_column = 8;
}


//kill the script and break out of infinite loop if ssh disconnects
//we'll pick it back up with a script
function disconnecting()
{
   exit('ssh2 disconnected');
}

$c = ssh2_connect($r_server, $r_port, array(), array('disconnect' => 'disconnecting'));

ssh2_auth_password($c, SCV_USER, SCV_PASS) or die("no connection\n");

while(1)
{
   echo "status: $job->status\n";
   $job->touch();
   switch ($job->status) {
      case 'l.b':
	#check how many things are already running pbpdb.pl
	$cmd = "pgrep -c pbpdb";
	echo $cmd, "\n";
	$output = array();
	exec($cmd, $output);
	if ($output[0] > 3)
	{
		break;
	}
         $cmd  = "pbpdb.pl " . PRMS_DIR . "/pdbamino.rtf ";
         $cmd .= PRMS_DIR . "/parm.prm ";
         $cmd .= "$local_dir/ protorig_cmin.pdb $job->pbmode";
         echo $cmd, "\n";
         system($cmd);

         $cwd = getcwd();

         chdir($local_dir);
         if (!file_exists("protorig_cmin_accs_$job->pbmode.pot") || filesize("protorig_cmin_accs_$job->pbmode.pot") === 0) {
		$job->le('PB Potential not created');
	}
         rename("protorig_cmin_accs_$job->pbmode.pot", "pb.pot");
         unlink("protorig_cmin_accs_$job->pbmode.kap");
         rename("protorig_cmin_accs_$job->pbmode.pdb", "1rec.pdb");
         copy("1rec.pdb", "1rec.temp.pdb");

#         $cmd = "pdb.renumberResidues.pl 1rec.temp.pdb";
#         echo $cmd, "\n";
#         system($cmd);

         if ($job->keep_metals == 't') {
	    $cmd = copy("1rec.temp.pdb", "1rec.charmm.pdb");
         } else {
            $cmd  = "grep -Ev '(ZN)|(HE(C|M))|(CO)|(MN)' 1rec.temp.pdb | "; //remove metals
            $cmd .= "cut -c -67 > 1rec.charmm.pdb"; //and only keep first 67 columns for charmm
         }
         echo $cmd, "\n";
         system($cmd);

         $cmd = "echo 'END' >> 1rec.charmm.pdb"; //append END to make charmm happy
         echo $cmd, "\n";
         system($cmd);

         $cmd = 'pdb.renumberAtoms.pl 1rec.charmm.pdb';
         echo $cmd, "\n";
         system($cmd);

         $cmd  = "1sidehphobe_markms 1rec.pdb 1rec.1sidehp3.ms ";
         $cmd .= PRMS_DIR . "/atoms.0.0.4.prm.ms.3cap+0.5ace_hphobe3.Hr0rec";
         echo $cmd, "\n";
         system($cmd);

         //give us our probes
         $cmd  = "tar xf " . PRMS_DIR . "/probes-$job->probeset.tar";
         echo $cmd, "\n";
         system($cmd);

         rename("1rec.pdb", "probes/rec/1rec.pdb");
         rename("protorig_cmin.psf", "probes/rec/1rec.psf");
         rename("1rec.charmm.pdb", "probes/rec/1rec.charmm.pdb");
         rename("pb.pot", "probes/rec/pb.pot");
         rename("1rec.1sidehp3.ms", "probes/rec/1rec.1sidehp3.ms");
         if ($job->protmask != '')
         {
            rename("mask.pdb", "probes/rec/mask.pdb");

            $cmd = "fake_msur.pl probes/rec/mask.pdb probes/rec/mask.ms";
            echo $cmd, "\n";
            system($cmd);
         }

         $cmd = "touch complexs";
         echo $cmd, "\n";
         system($cmd);

	 chdir('probes');
         exec('ls -d */', $content);

         foreach ($content as $line )
         {
            $probeid = rtrim($line, '/');
            if ($probeid !== 'rec')
            {
	        chdir($local_dir);
		$probe = new App_Probe($jobid, $probeid);

                $cmd = "echo $probeid >> complexs";
                echo $cmd, "\n";
                system($cmd);

		chdir("$local_dir/probes/$probeid");

                $cmd = 'cat 1rec.charmm.pdb 1lig.ms > 1rec1lig.pdb';
                echo $cmd, "\n";
                system($cmd);

                $cmd = 'pdb.renumberAtoms.pl 1rec1lig.pdb';
                echo $cmd, "\n";
                system($cmd);

                $cmd = 'pdbprep.pl 1rec1lig.pdb';
                echo $cmd, "\n";
                system($cmd);

		App_PDB::psf(getcwd(), '1rec1lig', '?');

		$cmd = 'tail -n `wc -l 1lig.ms | awk \'{print $1}\'` 1rec1lig.pdb > 1lig.pdb';
                echo $cmd, "\n";
                system($cmd);

            }
         }

         //move back into directory we started in
         chdir($cwd);

         $job->status = 'l.r';
      //job is ready to transmit files
      case 'l.r':
         $job->touch();

         Util::scp_send_dir($c, $local_dir, $remote_dir);

         //prevent us from overwriting ourselves on the way back
         $cmd = "rm -f $remote_dir/jobrunner.log";
         echo $cmd, "\n";
         ssh2_exec($c, $cmd);

         $job->touch();

         $job->status = 'r.i';
      case 'r.i':
         $job->touch();

         echo "submit job\n";
         $job->status = 'r.q';
      case 'r.q':
      case 'r.r':
         $job->touch();

			$unsubmitted_fft = $job->probes_w_status('f.i');
			if ( sizeof($unsubmitted_fft) > 0 )
			{
         	$mpicmd  = "$mpirun -np $np ";
         	$mpicmd .= "$piper ";
         	$mpicmd .= "-p $job->atomsfile ";
         	$mpicmd .= "-f $prms_dir/$job->coeffsfile ";
         	$mpicmd .= "-r $prms_dir/$job->rotsfile ";
         	$mpicmd .= "-c0.8 ";
         	$mpicmd .= "-k6 ";
         	$mpicmd .= "-t$job->kperrot ";
         	$mpicmd .= "-d3.0 ";
         	$mpicmd .= "-R$job->nrots ";
         	$mpicmd .= "--water_sigma=10.0 ";
         	$mpicmd .= "--pb=pb.pot ";
         	$mpicmd .= "--pbt=15.0 ";
         	$mpicmd .= "--rvdw_sa_scale=1.0 ";
         	$mpicmd .= "--avdw_sa_scale=1.0 ";
         	$mpicmd .= "--avdw_sa_inc=3.5 ";
		$mpicmd .= "--avdw_nsa_scale=1.0 ";
         	$mpicmd .= "--avdw_nsa_inc=3.5 ";
         	$mpicmd .= "--ff-psf-rec=1rec.psf --ff-psf-lig=1lig.psf ";
         	$mpicmd .= "-vv ";
         	if ($job->protmask != '')
         	{
            	$mpicmd .= "--maskr=1.0 ";
            	#$mpicmd .= "--mask=mask.ms ";
            	$mpicmd .= "--maskrec=mask.ms ";
         	}
         	$mpicmd .= "1rec.1sidehp3.ms 1lig.ms";

		foreach ( $unsubmitted_fft AS $probeid)
				{
					$probe = new App_Probe($jobid, $probeid);
					$cmd = "cd $remote_dir/probes/$probeid;$mpicmd";
					echo "cmd: $cmd\n";

					$res = ssh2_exec($c, $cmd);
					stream_set_blocking($res, true);
         		//stream contents like
         		if ('scc2' == $job->server)
			{
				$parts = explode(' ', stream_get_contents($res));
				$probe->fftqueueid = $parts[2];

			}
         		if ('katana' == $job->server)
         		{
            		$parts = explode(' ', stream_get_contents($res));
            		$probe->fftqueueid = $parts[2];
         		}
         		elseif ('lee' == $job->server)
         		{
            		//llsubmit: The job "fe2.bgl.bu.edu.103287" has been submitted.
            		$parts = explode('"', stream_get_contents($res));
            		$probe->fftqueueid = $parts[1];
			$cmd = "llprio -1 $parts[1]";
			echo "cmd: $cmd\n";
			ssh2_exec($c, $cmd);
         		}
         		fclose($res);
         		$probe->status = 'f.q';
				}
			}

			$queued_fft = $job->probes_w_status('f.q');
			$running_fft = $job->probes_w_status('f.r');
			$submitted_fft = array_merge($queued_fft, $running_fft);

			if (count($submitted_fft) > 0)
			{
				$cmd = "qstat";
				echo "cmd: $cmd\n";
				$res = ssh2_exec($c, $cmd);
				stream_set_blocking($res, true);
				$qstat = explode("\n", stream_get_contents($res));
				fclose($res);


				//limit qstat output to the user running these jobs
				$qstat = preg_grep('/'.SCV_USER.'/', $qstat);
			}

			foreach ( $submitted_fft AS $probeid )
			{
				$probe = new App_Probe($jobid, $probeid);

				$probeline = array_values(preg_grep("/$probe->fftqueueid/", $qstat));

				if (count($probeline) > 0)
				{
					$probeparts = preg_split("/\s+/", trim($probeline[0]));
					$status = $probeparts[$status_column-1];
				} else {
					$status = '';
				}

				echo "$probeid status: $status\n";

				if (strtoupper($status) === 'R')
				{
					$probe->status = 'f.r';
				}
				if ($status === '')
				{
					$probe->status = 'f.f';
				}
			}

			$finished_fft = $job->probes_w_status('f.f');
			foreach ($finished_fft AS $probeid)
			{
				$probe = new App_Probe($jobid, $probeid);
				$probe_dir  = "$remote_dir/probes/$probeid";

				$k_per_rot = $job->kperrot;
				$top_energies = $job->topenergies;
				for ($k = 0; $k < $k_per_rot; $k++)
				{
					$k = sprintf('%02d', $k);
					$cmd = "wc -l $probe_dir/ft.000.$k";
					echo $cmd, "\n";
					$wc = ssh2_exec($c, $cmd);
					stream_set_blocking($wc, true);
					$lines = (int) stream_get_contents($wc);
					fclose($wc);
					if ($lines != $job->nrots)
					{
						$probe->status = 'f.e';
					}
				}

				//sometimes lee is stupid and doesn't boot
				//checking for a blank ll.out allows us to see if that happened
				if ('f.e' === $probe->status)
				{
					if ('scc2' == $job->server)
					{
						$outfile = "sge.jcf.po$probe->fftqueueid";
						$errfile = "sge.jcf.o$probe->fftqueueid";
					}
					if ('katana' == $job->server)
					{
						$outfile = "sge.jcf.o$probe->fftqueueid";
						$errfile = "sge.jcf.e$probe->fftqueueid";
					}
					elseif ('lee' == $job->server)
					{
            					$parts = explode('.', $probe->fftqueueid);
						$queueid = $parts[4];
						$outfile = "ll.$queueid.out";
						$errfile = "ll.$queueid.err";
					}
					$cmd = "wc -l $probe_dir/$outfile";
					echo $cmd, "\n";
					$wc = ssh2_exec($c, $cmd);
					stream_set_blocking($wc, true);
					$lines = (int) stream_get_contents($wc);
					fclose($wc);
					if ($lines == 0)
					{
						$probe->status = 'f.i';
						break;
					}
					else
					{
						$cmd = "grep -c deallocate $probe_dir/$errfile";
						echo $cmd, "\n";
						$wc = ssh2_exec($c, $cmd);
						stream_set_blocking($wc, true);
						$dealloc = (int) stream_get_contents($wc);
						fclose($wc);
						if ($dealloc != 0)
						{
							$probe->status = 'f.i';
							break;
						} else {
						$job->re("Not enough lines in output for $probeid");
						}
					}
				}


				//create basename for clusters
				$file_k = sprintf('%02d', $k_per_rot-1);
				$cluster_ft = "ft.000.00-$file_k";

				//concatenate files into single ft file
				$file_k = $k_per_rot-1;
				$cmd = "cat $probe_dir/ft.000.0[0-$file_k] | sort -snk5 | head -n $top_energies > $probe_dir/$cluster_ft.h";
				echo $cmd, "\n";
				ssh2_exec($c, $cmd);

				//generate clusters
				$cmd  = "cd $probe_dir;";
				$cmd .= SCV_BIN . "/gen_ftclusterswithneighborsbysize.0.0.4.pl $cluster_ft.h 10 4";
				echo $cmd, "\n";
				ssh2_exec($c, $cmd);

				$top_clusters = $job->topclusters;
				$file_tc = $top_clusters-1;
				$cmd = "cat $probe_dir/$cluster_ft.h.clusters.[0-$file_tc].nghbrs > $probe_dir/ftfinal";
				echo $cmd, "\n";
				ssh2_exec($c, $cmd);

				$probe->status = 'c.f';

				if ( $job->skipcharmm == 't' ) {
					$probe->status = 'm.f';
				}
				$job->touch();
			}

			$finished_cluster = $job->probes_w_status('c.f');
			foreach ($finished_cluster AS $probeid)
			{
				$probe = new App_Probe($jobid, $probeid);
				$probe_dir  = "$remote_dir/probes/$probeid";

				$cmd  = "cd $probe_dir;";
				$cmd .= "$mpirun -np 16 ";
				$cmd .=	SCV_BIN . "/ace_min.mpi 1rec1lig_cmin.psf parm_new.prm pdbamino_new.rtf /dev/null 1rec1lig.pdb 1rec.charmm.pdb 0.5 $prms_dir/$job->rotsfile 1lig.pdb ftfinal ~tanggis/src/ace_minmpi/trunk/allh.tl 1lig.pdb";
				echo "cmd: $cmd\n";

				$res = ssh2_exec($c, $cmd);
				stream_set_blocking($res, true);

				if ('scc2' == $job->server)
				{
					$parts = explode(' ', stream_get_contents($res));
					$probe->charmmqueueid = $parts[2];

				}
				elseif ('lee' == $job->server)
				{
					$parts = explode('"', stream_get_contents($res));
					$probe->charmmqueueid = $parts[2];
				}
				fclose($res);


				$probe->status = 'm.q';
			}

			$queued_min = $job->probes_w_status('m.q');
			$running_min = $job->probes_w_status('m.r');
			$submitted_min = array_merge($queued_min, $running_min);

			if ( count( $submitted_min ) > 0)
			{
				$cmd = "qstat";
				echo "cmd: $cmd\n";
				$res = ssh2_exec($c, $cmd);
				stream_set_blocking($res, true);
	 			$qstat = explode("\n", stream_get_contents($res));
				fclose($res);

				//limit qstat output to the user running these jobs
				$qstat = preg_grep('/'.SCV_USER.'/', $qstat);
			}

			foreach ($submitted_min AS $probeid)
			{
				$probe = new App_Probe($jobid, $probeid);

				$probeline = array_values(preg_grep("/$probe->charmmqueueid/", $qstat));

				if (count($probeline) > 0)
				{
					$probeparts = preg_split("/\s+/", trim($probeline[0]));
					$status = $probeparts[$status_column-1];
				} else {
					$status = '';
				}

				echo "$probeid status: $status\n";

				if (strtoupper($status) === 'R')
				{
					$probe->status = 'm.r';
				}

				if ($status === '')
				{
					$probe->status = 'm.f';
					$probe_dir  = "$remote_dir/probes/$probeid";

					//check if output file is empty => lee failed to allocate partition
					if ('scc2' == $job->server)
					{
						$outfile = "sge.jcf.po$probe->charmmqueueid";
						$errfile = "sge.jcf.o$probe->charmmqueueid";
					}

					$cmd = "wc -l $probe_dir/$outfile";
					echo $cmd, "\n";
					$wc = ssh2_exec($c, $cmd);
					stream_set_blocking($wc, true);
					$lines = (int) stream_get_contents($wc);
					fclose($wc);
					if ($lines == 0)
					{
						$probe->status = 'c.f';
						break;
					}

					$cmd = "grep -c deallocate $probe_dir/$errfile";
					echo $cmd, "\n";
					$wc = ssh2_exec($c, $cmd);
					stream_set_blocking($wc, true);
					$dealloc = (int) stream_get_contents($wc);
					fclose($wc);
					if ($dealloc != 0)
					{
						$probe->status = 'c.f';
						break;
					}

					//check if all conformers are minimized
					$cmd = "wc -l $probe_dir/ftfinal";
					echo $cmd, "\n";
					$wc = ssh2_exec($c, $cmd);
					stream_set_blocking($wc, true);
					$lines = (int) stream_get_contents($wc);
					echo $lines, "\n";
					fclose($wc);

					$cmd = "ls $probe_dir/min_out????.pdb|wc -l";
					echo $cmd, "\n";
					$wc = ssh2_exec($c, $cmd);
					stream_set_blocking($wc, true);
					$lines1 = (int) stream_get_contents($wc);
					echo $lines1, "\n";
					fclose($wc);

					if ($lines != $lines1)
					{
						$probe->status = 'm.e';
						break;
					}

					$cmd  = SCV_BIN . "/process_acemin_results_1probe.sh $probe_dir";
					echo $cmd, "\n";
					ssh2_exec($c, $cmd);

				}

			}

			$job->touch();

         if ($job->probes_running())
         {
            $job->status = 'r.r';
         }
         if ($job->probes_finished())
         {
            $job->status = 'r.f';
         }

         if ($job->status !== 'r.f')
         {
            break;
         }
      case 'r.f':
         $job->touch();
         //do any verifying of results
         //check if any of the probes had errors
         if ($job->probes_errored())
         {
            $job->status = 'r.e';
            die;
         }

	 $cmd  = SCV_BIN . "/process_acemin_results_modi.sh $remote_dir/probes ";
	 $cmd .= "$remote_dir/complexs ";
	 $cmd .= "cluster";
	 echo $cmd, "\n";
	 ssh2_exec($c, $cmd);

         $job->status = 'r.o';
      case 'r.o':
         $job->touch();
         //time to drop off the files locally
         Util::scp_recv_dir($c, $remote_dir, $local_dir);

         //clean up data on server
         ssh2_exec($c, "sh -c 'rm -rf $remote_dir &'");

         $job->status = 'l.i';

			if ( $job->skipcharmm == 't' ) {
				$job->status = 'c.f';
			}

      case 'l.i':
      case 'l.c':
         $job->touch();

         chdir("$local_dir/probes/cluster");

         $complexs = fopen("$local_dir/complexs", 'r');
         while (!feof($complexs)) {
            $pdb = rtrim(fgets($complexs));
            if (strlen($pdb) > 0) {
                $cmd = "ln -s ../$pdb/1lig.psf $pdb.psf";
		#$cmd = 'ln -s ' . ROOT_DIR . "/pdbpsf/${pdb}_cmin.psf $pdb.psf";
                echo "$cmd\n";
                system($cmd);
                $cmd = "ln -s ../$pdb/1lig.ms $pdb.pdb";
                echo "$cmd\n";
                system($cmd);
            }
         }
         fclose($complexs);

         $top_clusters = $job->topclusters;
	 $top_energies = $job->topenergies;
	 if ($job->probeset == 'customize') {
		 $cmd = "smcluspro -n $top_energies -u $top_clusters -r 4 -R 4 -s 10 -S 2 $local_dir/probes/probes-ext;";
		 $cmd .= "rm crosscluster*pdb;";
		 $cmd .= "smcluspro -n $top_energies -u $top_clusters -r 4 -R 4 -s 10 -S 2 $local_dir/probes/probes";
	 }else {
		 $cmd = "smcluspro -n $top_energies -u $top_clusters -r 4 -R 4 -s 10 -S 2 $local_dir/complexs";
	 }
         echo $cmd, "\n";
         system($cmd);

	 $cmd = "ls crosscluster*pdb | wc -l";
	 echo $cmd, "\n";
	 $output = array();
	 exec($cmd, $output);
	 if ($output[0] == 0) {
		 $job->le('Failed during clustering');
	 }

	 if ($job->probeset == 'customize') {
		 $cmd = "select_clusters.pl 4 ../probes-ext";
		 echo $cmd, "\n";
		 system($cmd);
	 }

         chdir($local_dir);

         echo "copying protein 1rec.pdb\n";
         copy('probes/rec/1rec.pdb', '1rec.pdb');

         $cmd = "addchain.pl 1rec.pdb";
         echo "$cmd\n";
         system($cmd);

	 //make pymol session
	 $cmd = "create_pymol_session.sh fftmap.$jobid.pse";
	 echo "$cmd\n";
	 system($cmd);

         $cmd = "echo 'HEADER protein' >fftmap.$jobid.pdb";
         echo "$cmd\n";
         system($cmd);

	//use protorig_join for proteins with over 1000 residues
	if (file_exists('protorig_join.pdb')) {
        	 $cmd = "cat protorig_join.pdb >> fftmap.$jobid.pdb";
        	 echo "$cmd\n";
        	 system($cmd);
	} else {
        	 $cmd = "cat 1rec.pdb >> fftmap.$jobid.pdb";
        	 echo "$cmd\n";
        	 system($cmd);
	}

         echo "unlinking copy of protein 1rec.pdb\n";
         unlink('1rec.pdb');

         chdir("$local_dir/probes/cluster");
         $output = array();
         $cmd = "ls crosscluster.0*.pdb";
         echo "$cmd\n";
         exec($cmd, $output);
         foreach ($output AS $crosscluster)
         {
		 $cmd = "pdb.renumberAtoms.pl $crosscluster";
		 echo "$cmd\n";
		 system($cmd);
	 }

         chdir($local_dir);
         foreach ($output AS $crosscluster)
         {
            $cmd = "echo 'HEADER $crosscluster' >> fftmap.$jobid.pdb";
            echo "$cmd\n";
            system($cmd);

            $cmd = "cat probes/cluster/$crosscluster >> fftmap.$jobid.pdb";
            echo "$cmd\n";
            system($cmd);
         }

         $chain = $job->protchains;
         $cmd = "addprobechaintag21.pl fftmap.$jobid.pdb '$chain'";
         echo "$cmd\n";
         system($cmd);

         $cmd = "pdb.renumberAtoms.pl fftmap.$jobid.pdb";
         echo "$cmd\n";
         system($cmd);

         $job->status = 'c.f';
      case 'c.f':
         $job->touch();
         chdir($local_dir);

         mkdir("$local_dir/workdir");

         copy("probes/rec/1rec.pdb", "workdir/1rec.pdb");

         $cmd = "grep -h ATOM probes/cluster/*.min.pdb > workdir/probes.allpdb";
         echo "$cmd\n";
         system($cmd);

         chdir("$local_dir/workdir");

         $cmd = "rmchain.pl probes.allpdb";
         echo "$cmd\n";
         system($cmd);

         $cmd = "ls *.allpdb >listallpdb";
         echo "$cmd\n";
         system($cmd);

         $cmd  = "splitallpdb.pl listallpdb ../complexs";
         echo "$cmd\n";
         system($cmd);

         $cmd = "splitintokfileslig.pl splitlist 10";
         echo "$cmd\n";
         system($cmd);

         $cmd = "addchain.pl 1rec.pdb";
         echo "$cmd\n";
         system($cmd);

         $file = fopen('splitlist', 'r');
         while ( !feof($file) )
         {
            $line = trim(fgets($file));
            if ($line === '') {
                continue;
            }
            $cmd = 'seq 0 9|xargs -i bash -c "cat 1rec.pdb '.$line.'{} > '.$line.'{}rec.pdb"';
            echo "$cmd\n";
            system($cmd);
         }
         fclose($file);

         $cmd = "ls *rec.pdb> listfinal";
         echo "$cmd\n";
         system($cmd);

         $cmd = "Mapping_FixHET.sh";
         echo "$cmd\n";
         system($cmd);

         $cmd = "convert_probe_atoms.pl listfinal";
         echo "$cmd\n";
         system($cmd);

         $file = fopen('listfinal', 'r');
         while ( !feof($file) )
         {
            $job->touch();
            $line = trim(fgets($file));
            $cmd  = "/usr/bin/python ";
            $cmd .= BIN_DIR . "/calc_interactions.py ";
            $cmd .= BIN_DIR . "/hbplus $line";
            echo "$cmd\n";
            system($cmd);
         }
         fclose($file);

         mkdir('nnbout');

         $cmd = "cat *.nnb.out > nnbout/result.all.nnb.out";
         echo "$cmd\n";
         system($cmd);

         $cmd = "cat *.hhb.out > nnbout/result.all.hhb.out";
         echo "$cmd\n";
         system($cmd);

         chdir('nnbout');

         $cmd = "extractinteractionnumber.pl result.all.nnb.out";
         echo "$cmd\n";
         system($cmd);

         $cmd = "extractinteractionnumber.pl result.all.hhb.out";
         echo "$cmd\n";
         system($cmd);

         chdir($local_dir);

         $cmd = "mv workdir/nnbout/*.rawextract ./";
         echo "$cmd\n";
         system($cmd);

         $cmd = "sort -k2,2 -k1,1n result.all.nnb.out.rawextract > nonbonded.$jobid.rawextract";
         echo "$cmd\n";
         system($cmd);

         $cmd = "sort -k2,2 -k1,1n result.all.hhb.out.rawextract > hbonded.$jobid.rawextract";
         echo "$cmd\n";
         system($cmd);

	//clean everything up
	//$cmd = "rm -r workdir probes/*/pbg.dx probes/*/r_CA_nothphil_cavity_sigma_10.0.dx";
	$cmd = "find . -name min_out*pdb|xargs rm";
	echo "$cmd\n";
	system($cmd);
	$cmd = "rm -r probes/rec/pb.pot workdir";
	echo "$cmd\n";
	system($cmd);
	system("gzip --best $local_dir/*.pdb");
	system("gzip --best $local_dir/user/*.pdb");

        $job->status = 'l.f';

        if (strlen($job->callbackurl) > 0) {
           $request = new HttpRequest($job->callbackurl, HTTP_METH_POST);
           $request->addPostFields(array('id'=>$job->jobid, 'status' => 'l.f'));
           $request->send();
        }


      default:
         exit;
    }
	//close the db connection while we're asleep
	pg_close($dbconn);
	$cmd = "pgrep -fc ".realpath($_SERVER['SCRIPT_FILENAME']);
	echo $cmd, "\n";
	$output = array();
	exec($cmd, $output);
	echo 'sleeping for ', 60 * $output[0], " seconds\n";
	sleep(60 * $output[0]);
	require 'env/dbvars.php';
}
