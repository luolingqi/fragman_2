<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'sessionstart.php';

if ( isset($_POST['action']) && ('piper' !== $liuser->username) )
{
   if (isset($_POST['newpassword']) && ($_POST['newpassword'] != '')  )
   {
      if ($_POST['newpassword'] === $_POST['newpassword2'])
      {
         $liuser->password = md5($_POST['newpassword']);
         $successes[] = 'Updated Password';
      }
      else
      {
         $errors[] = 'New Passwords Do Not Match';
      }
   }

   //if it's set, then that means set to true, else, set to false
   if ( isset($_POST['emailuser']) && (!$liuser->emailJobComplete) ) {
       //user has turned on email
       $liuser->emailJobComplete = true;

       $successes[] = 'Turned on Email Notification';
   }
   if ( (! isset($_POST['emailuser']) ) && ($liuser->emailJobComplete) ) {
       //user has turned off email
       $liuser->emailJobComplete = 0;

       $successes[] = 'Turned off Email Notification';
   }
}

$page = new App_Page( array('name' => 'Preferences', 'require_login' => true) );

$page->header();

if ( isset($successes) && sizeof($successes) >0 )
{
   echo '<div id="successbox">';
   echo '<h2>The following preferences were updated:</h2>';
   echo '<ul>';
   foreach ($successes AS $success)
   {
      echo "<li>$success</li>";
   }
   echo '</ul>';
   echo '</div>';
}

if ( isset($errors) && sizeof($errors) >0 )
{
   echo '<div id="errorbox">';
   echo '<h2>There were some errors in updating your preferences</h2>';
   echo '<ul>';
   foreach ($errors AS $error)
   {
      echo "<li>$error</li>";
   }
   echo '</ul>';
   echo '</div>';
}

if ($liuser->username == 'piper')
{
	echo "<h3>Preferences</h3><p>Sorry, anonymous user's preferences cannot be changed";
}
else {
?>
<h3>Preferences</h3>

<form class="goodform" action="preferences.php" method="post">
   <div>
      <label for="newpassword">New Password:</label>
      <input type="password" name="newpassword" id="newpassword" />
   </div>
   <div>
      <label for="newpassword2">Confirm New Password:</label>
      <input type="password" name="newpassword2" id="newpassword2" />
   </div>
   <div>
      <input type="checkbox" name="emailuser" id="emailuser" <?php echo ($liuser->emailJobComplete) ? 'checked="checked"' : ''; ?> value="checked" style="float:left;left:150px;position:relative;" />
      <label for"emailuser" style="float:right;textalign:left;">Email when jobs finished</label>
   </div>
<div id="submit">
    <input name="action" type="submit" value="Change Preferences" />
</div>
</form>
<?php
}
$page->footer();
