<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'sessionstart.php';

$page = new App_Page( array( 'name' => 'Results', 'require_login' => true) );
$page->header();
?>

   <h3>Finished Jobs</h3>
<p>Note, jobs will be deleted 2 months after being submitted</p>
<?php
if ( $liuser->isSuperUser() )
{
    echo "<a href='?user=$liuser->username'>$liuser->username</a>";
    echo " | <a href='?user='>everyone</a>";
}
?>
	 <table class="nice"> 

	<tr>
	<th>Id</th>
	<th>Name</th>
<?php if ( $liuser->isSuperUser() )
      {
echo "<th>User</th>"; 
      };
?>
	<th>Status</th>
	</tr>

<?php

$user = isset($_GET['user']) ? $_GET['user'] : '';
$user_url = ($user == '') ? '' : "&user=$user";
$offset = isset($_GET['offset']) ? (int) $_GET['offset'] : 0;

$query = "
	select id, jobname, status
	from jobs
	where status IN ('l.e', 'r.e', 'l.f', 'l.d')
   and userid = $_SESSION[userid]
	order by id desc
	limit $jobsperpage+1
	offset $offset*$jobsperpage
	";
if ( $liuser->isSuperUser() )
{
    if ($user == '')
    {
        $query = "
            select jobs.id, jobs.jobname, users.username, jobs.status
            from jobs, users
            where status IN ('l.e', 'r.e', 'l.f', 'l.d')
            and jobs.userid = users.userid
            order by id desc
            limit $jobsperpage+1
            offset $offset*$jobsperpage
            ";
    }
    else
    {
        $query = "
            select jobs.id, jobs.jobname, users.username, jobs.status
            from jobs, users
            where status IN ('l.e', 'r.e', 'l.f', 'l.d')
            and jobs.userid = users.userid
            and users.username = '$user'
            order by id desc
            limit $jobsperpage+1
            offset $offset*$jobsperpage
            ";
    }
}
#echo "$query<br>";
$result = pg_query($query) or die('Query failed: ' . pg_last_error());

for ($i = 0; $i < pg_num_rows ($result) && $i < $jobsperpage; $i++)
{
	$line = pg_fetch_assoc($result);

	echo "<tr>";
		$val = $line['id'];
		echo "<td>";
		if ($line['status'] == 'l.f')
		{
			echo "<a href='models.php?job=$val'>$val</a>";
		}
		else // job error
		{
			echo "<a href='jobdetail.php?job=$val'>$val</a>";
		}
		echo "</td>";

		echo "<td>$line[jobname]</td>";
   
      if ( $liuser->isSuperUser() )
      {
                echo "<td><a href='?user=$line[username]'>$line[username]</a></td>";
      }
		
		echo "<td>".Job::convertMedStatus($line['status'])."</td>";
	echo "</tr>";
}

echo "	</table>";
echo "<br />";

$prev=0; // remember if prev is set
if ($offset > 0)
{
	$prevoffset = $offset-1;
	echo "<a href='?offset=$prevoffset$user_url'><- prev</a>";
	$prev=1;
}
if ( pg_num_rows($result) > $jobsperpage)
{
	if ($prev)
	{
		//echo "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;";
		//echo "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;";
		echo "&nbsp;&nbsp;&nbsp;";
	}
	$nextoffset = $offset+1;
	echo "<a href='?offset=$nextoffset$user_url'>next -></a>";
}

$page->footer();
