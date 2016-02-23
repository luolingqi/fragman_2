<?php
//init error
$errors = array();

if (! isset($_SESSION['key']))
{
	header('Location: signup.php');
	exit;
}

//trim user input
$_POST = array_map('trim', $_POST);

// check user input

// first name
$lenmax = 24;
if ( $_POST['firstname'] === '' )
{
	$errors[] = 'Please enter your first name.';
}
if (strlen($_POST["firstname"]) > $lenmax)
{
	$_POST['firstname'] = substr($_POST['firstname'], 0, $lenmax);	
}
if (! preg_match('/^[\w.]*$/', $_POST['firstname']))
{
	$errors[] = 'First name can contain only letters (a-z), numbers (0-9), periods (.).';
}

// last name
$lenmax = 24;
if ( $_POST['lastname'] === '' )
{
	$errors[] = 'Please enter your last name.';
}
if (strlen($_POST["lastname"]) > $lenmax)
{
   //if their last name is too long, just cut it off
	$_POST['lastname'] = substr($_POST['lastname'], 0, $lenmax);	
}
if (! preg_match('/^[\w.]*$/', $_POST['lastname']))
{
	$errors[] = 'Last name can contain only letters (a-z), numbers (0-9), periods (.).';
}

// desired username
$lenmin = 5; $lenmax = 24;
if (strlen($_POST['username']) < $lenmin ||
	strlen($_POST['username']) > $lenmax
		)
{
	$errors[] = "Your username must be between $lenmin and $lenmax characters long.";
}
if (! preg_match('/^[\w.]*$/', $_POST['username']))
{
	$errors[] = 'Your username can contain only letters (a-z), numbers (0-9), periods (.).';
}
if (preg_match('/^admin$/i', $_POST['username']))
{
	$errors[] = 'Sorry, the username you have chosen is reserved.';
}
if ( preg_match('/'.SERVER_NAME.'/i', $_POST['username']) )
{
	$errors[] = 'Sorry, variants of '.SERVER_NAME.' are not allowed in your username.';
}
if (preg_match('/piper/i', $_POST['username']))
{
	$errors[] = 'Sorry, variants of Piper are not allowed in your username.';
}
else // check that it's unique
{
	$esc_username = pg_escape_string($_POST['username']);
	$query = "
		select COUNT(*)
		from users
		where username = '$esc_username'
		";
	//echo "$query<br>";
	$result = pg_query($query) or die('Query failed: ' . pg_last_error());
   list($taken) = pg_fetch_row($result);

	if ($taken > 0)
	{
		$errors[] = 'Sorry, the username you have chosen is already in use.';
	}
}

// affiliation
$lenmax = 50;
if ( $_POST['affiliation'] === '')
{
	$errors[] = 'Please enter your affiliation.';
}
if (strlen($_POST['affiliation']) > $lenmax)
{
	$_POST['affiliation'] = substr($_POST['affiliation'], 0, $lenmax);	
}

// email
$lenmax = 100;
if ( $_POST['email'] === '')
{
	$errors[] = 'Please enter your email.';
}
if (strlen($_POST['email']) > $lenmax)
{
	$errors[] = "Your email must be fewer than $lenmax characters long.";
}
if (! preg_match('/^[_a-z0-9-]+(\.[_a-z0-9-]+)*@[a-z0-9-]+(\.[a-z0-9-]+)*(\.[a-z]{2,3})$/i',
			$_POST['email'])
		)
{
	$errors[] = 'The provided email address is invalid.';
}

//exclude certain email addresses
$endings = array( '\.aero',
                   '\.asia',
                   '\.biz',
                   '\.cat',
                   '\.com',
                   '\.coop',
                   '\.info',
                   '\.int',
                   '\.jobs',
                   '\.mobi',
                   '\.name',
                   '\.net',
                   '\.org',
                   '\.pro',
                   '\.tel',
                   '\.travel',
                   '\.us' ,
                   '@.*alum.*\.edu',
                   '\.com\.[a-z]{2}',
                   '\.co\.[a-z]{2}',
                   '\.me\.[a-z]{2}',
                   '\.net\.[a-z]{2}',
                   '\.nom\.[a-z]{2}',
                   '\.org\.[a-z]{2}'
                  );

$pattern = '/(' . implode(')|(', $endings) . ')$/i';
if ( preg_match($pattern, $_POST['email']) )
{
	$errors[] = 'Sorry, we are only accepting educational and government email addresses at this time';
}
$esc_email = pg_escape_string($_POST['email']);
$query = "
	select COUNT(*)
	from users
	where email = '$esc_email'
	";
//echo "$query<br>";
$result = pg_query($query) or die('Query failed: ' . pg_last_error());
list($exists) = pg_fetch_row($result);

if ($exists > 0)
{
	$errors[] = 'Sorry, the email address you have chosen is already in use.';
}


// word verification
// note the generated string is always lowercase
if (md5(strtolower($_POST['wordverification'])) != $_SESSION['key'])
{
	$errors[] = 'The word verification code was incorrect. Please try again.';
}

//tos
if (! isset($_POST['tos']))
{
	$errors[] = 'You must agree to the Terms of Service';
}

// if error, leave user on signup page
if ( sizeof($errors) === 0)
{
   // generate password for user
   $nchars = 8;
   $password = Util::randstring($nchars);
   
   // insert user into db
   $e = array_map('pg_escape_string', $_POST);
   $e['password'] = pg_escape_string($password);
/*   $esc_firstname = pg_escape_string($_SESSION['firstname']);
   $esc_lastname = pg_escape_string($_SESSION['lastname']);
   $esc_username = pg_escape_string($_SESSION['username']);
   $esc_affiliation = pg_escape_string($_SESSION['affiliation']);
   $esc_email = pg_escape_string($_SESSION['email']);
 */  
   $query = "
   	insert into users
   	(username, password, email,
   	 firstname, lastname, affiliation,
   	 timecreated)
   	values
   	('$e[username]', MD5('$e[password]'), '$e[email]',
   	 '$e[firstname]', '$e[lastname]', '$e[affiliation]',
   	 'now')
   	 ";
   
   $result = pg_query($query) or die('query failed: ' . pg_last_error());
   
   
   // send email and thank user
   $to = $_POST['email'];
   $subject = "$servername account";
   $message = "Thank you for your interest in $servername.\n";
   $message .= "The password for your account is:\n\n";
   $message .= "$password\n\n";
   $message .= "Please login to automatically";
   $message .= " finalize the account creation process.\n\n";
   $message .= WEB_ROOT."\n";
   $message .= "\n\n\n";
   $message .= "Note: This was an automatically generated message.";
   $message .= " Please do not reply.";
   $from = "$servername <no-reply@bu.edu>";
   $headers = "From: $from";
   mail($to,$subject,$message,$headers);
   
   
   $page = new App_Page();
   $page->header();
   
   echo "<h3>Thank you for your interest in $servername.</h3>";
   echo "<div class='main-sub-aljustify'>";
   echo "<br>";
   echo "An email containing your password has been sent to the address:";
   echo "</div>";
   
   echo "<h4>" . $_POST["email"] . "</h4>";
   
   echo "<div class='main-sub-aljustify'>";
   echo "This email has the subject: \"$subject\".";
   echo " (If it does not arrive momentarily, please check that it was not filtered into";
   echo " your \"spam\" folder.)";
   echo " Also, please do not forget your username:";
   echo "</div>";
   
   echo "<h4>$_POST[username]</h4>";
   
   echo "<div class='main-sub-aljustify'>";
   echo "Your username was not included in the email for security reasons.";
   echo " Once you have received your password, please login to automatically";
   echo " finalize the account creation process.";
   echo "<br>";
   echo "<br>";
   echo "</div>";
   echo "<a href='login.php'>login</a>";
   
   $page->footer();
   
   // session is no longer needed
   session_destroy();
   exit;
}
