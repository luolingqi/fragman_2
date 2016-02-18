<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'sessionstart.php';

if (isset($_POST['action']) )
{
   $email = pg_escape_string($_POST['email']);
   $query = "SELECT userid FROM users WHERE email='$email'";
   $result = pg_query($query);

   if (pg_num_rows($result) > 0)
   {
      $urlkey = Util::randString(35);
      list( $userid ) = pg_fetch_row( $result );
      $user = new App_User($userid);
      $user->forgotpasskey = md5($urlkey);
      
      $to = $_POST['email'];
      $subject  = "$servername Password Request";
      $message  = "Someone (hopefully you) requested your password on $servername be reset. ";
      $message .= "If this was not you, you may ignore this message.\n\n";
      $message .= "Please click the link below to be sent a new password:\n";
      $message .= WEB_ROOT."resetpassword.php?k=$urlkey";
      
      $from     = "$servername <no-reply@bu.edu>";
      $headers  = "From: $from";
      mail($to, $subject, $message, $headers);

      $page = new App_Page();
      $page->header();
      echo "<p>An email with a link to reset your password has been sent to your address</p>.";
      $page->footer();
      exit;
   }
   else
   {
      $errors[] = 'No user found with that email address.';
   }
}

$page = new App_Page();
$page->header();



if ( isset($errors) && sizeof($errors) >0 )
{
   echo '<div id="errorbox">';
   echo '<h2>There were some errors in processing your request:</h2>';
   echo '<ul>';
   foreach ($errors AS $error)
   {
      echo "<li>$error</li>";
   }
   echo '</ul>';
   echo '</div>';
}
?>
<h3>Forgot Password</h3>
<form class="goodform" action="forgotpassword.php" method="post">
   <div>
      <label for="email">Your Email:</label>
      <input type="text" name="email" id="email" />
   </div>
<div id="submit">
    <input name="action" type="submit" value="Send New Password" />
</div>
</form>
<?php
$page->footer();
