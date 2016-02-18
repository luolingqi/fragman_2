<?php
require 'env/appvars.php';
require 'env/dbvars.php';

if(!isset($_GET['k']))
{
   exit;
}

$key = pg_escape_string($_GET['k']);
$query = "SELECT userid, email FROM users
          WHERE forgotpasskey=md5('$key')";
$result = pg_query($query);

if (pg_num_rows($result) > 0)
{
   $line = pg_fetch_assoc($result);
   
   $user = new App_User($line['userid']);
   $password = Util::randstring(8);
   $user->password = md5($password);
   $user->forgotpasskey = '';

   $to = $line['email'];
   $subject  = "New $servername Password";
   $message  = "You have successfully reset your passord for $servername. Your new password is:\n";
   $message .= "$password\n\n";

   $from     = "$servername <no-reply@bu.edu>";
   $headers  = "From: $from";
   mail($to, $subject, $message, $headers);
   
   $page = new App_Page();
   $page->header();
   echo "<p>Your new password has been sent to your email.</p>";
   $page->footer();
}
