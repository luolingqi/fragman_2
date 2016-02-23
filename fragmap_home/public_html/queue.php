<?php
require 'env/appvars.php';
require 'sessionstart.php';
require 'env/dbvars.php';
require 'env/path.php';

//use the queue page to look for jobs who haven't been touched
//for at least 100 minutes
$query  = "SELECT id FROM jobs WHERE touched < NOW()-interval '100minutes' AND status NOT IN ('l.e', 'r.e', 'l.f', 'l.d') FOR UPDATE;";
$query .= "UPDATE jobs SET touched=NOW() WHERE touched < NOW()-interval '300minutes' AND status NOT IN ('l.e', 'r.e', 'l.f','l.d')";

pg_send_query($dbconn, $query);
$result = pg_get_result($dbconn);
$sleep = 0;
while ( $line = pg_fetch_assoc($result) )
{
   $id = $line['id'];
   $local_dir = STORAGE_DIR . "/$id";
   $jobcmd = "jobrunner.php $id";
   exec("sleep $sleep; $jobcmd 1>>$local_dir/jobrunner.log &");
   $sleep += 1;
}
//clear out results from update part of query
$result = pg_get_result($dbconn);

$page = new App_Page( array( 'name' => 'Queue', 'require_login' => true, 'autorefresh' => true) );

$page->header();

?>
    <h3>Running Jobs</h3>
	 <table class='nice'>
      <tr>
      	<th>Id</th>
      	<th>Name</th>
      	<th>User</th>
      	<th>Status</th>
   	</tr>

<?php
$query = "
	select jobs.id, jobs.jobname, users.username, jobs.status
	from jobs, users
	where status NOT IN ('l.e', 'r.e', 'l.f', 'l.d')
	and jobs.userid = users.userid
	order by id desc
	";

$result = pg_query($query) or die('Query failed: ' . pg_last_error());

echo "<p><b>", pg_num_rows($result), "</b> jobs running for all users</p>";

while ( $line = pg_fetch_assoc($result) )
{
	echo '<tr>';

		$val = $line['id'];
      if ( $liuser->ownsJob($line['id']) )
      {
         $val = "<a href='jobdetail.php?job=$val'>$val</a>";
      }
		echo "<td>$val</td>";

		echo "<td>$line[jobname]</td>";

		echo "<td>$line[username]</td>";

		echo "<td>".Job::convertMedStatus($line['status'])."</td>";

	echo '</tr>';
}

?>

	</table>
<?php
$page->footer();
