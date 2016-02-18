<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'sessionstart.php';

if (! $liuser->isSuperUser() )
{
   header('Location: home.php');
   exit;
}

if ( isset($_POST['action']) )
{
   if ($_POST['makeprivilege'] != '')
   {
      $user = new App_User($_POST['makeprivilege']);
      $user->privilege = 1;
      $successes[] = "Gave $user->username Privileges";
   }
   if ($_POST['removeprivilege'] != '')
   {
      $user = new App_User($_POST['removeprivilege']);
      $user->privilege = 0;
      $successes[] = "Removed Privileges from $user->username";
   }
}

$page = new App_Page( array('name' => 'Admin', 'require_login' => true) );

$page->header();

if ( isset($successes) && sizeof($successes) >0 )
{
   echo '<div id="successbox">';
   echo '<h2>The following privileges were updated:</h2>';
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
   echo '<h2>There were some errors:</h2>';
   echo '<ul>';
   foreach ($errors AS $error)
   {
      echo "<li>$error</li>";
   }
   echo '</ul>';
   echo '</div>';
}

?>
<h3>Privileges</h3>

<form class="goodform" action="admin.php" method="post">
   <div>
      <label for="makeprivilege">Make User Privileged:</label>
      <select name="makeprivilege" id="makeprivilege">
         <option value=""></option>
         <?php
         $query = "SELECT userid, username FROM users WHERE privilege='0'";
         $result = pg_query($query);
         while ( $user = pg_fetch_assoc($result) )
         {
            echo "<option value='$user[userid]'>$user[username]</option>";
         }
         ?>
      </select>
   </div>
   <div>
      <label for="removeprivilege">Remove User Privileges:</label>
      <select name="removeprivilege" id="removeprivilege">
         <option value=""></option>
         <?php
         $query = "SELECT userid, username FROM users WHERE privilege='1'";
         $result = pg_query($query);
         while ( $user = pg_fetch_assoc($result) )
         {
            echo "<option value='$user[userid]'>$user[username]</option>";
         }
         ?>
      </select>
   </div>
<div id="submit">
    <input name="action" type="submit" value="Submit" />
</div>
</form>
<?php

$page->footer();
