<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'sessionstart.php';

if ( isset($_POST['action']) )
{
   $e_password = pg_escape_string($_POST['password']);
   $e_username = pg_escape_string($_POST['username']);
   $query = "SELECT userid FROM users
             WHERE username = '$e_username' 
             AND password = md5('$e_password')";
   $result = pg_query($query) or die('Query failed: ' . pg_last_error());
   
   if ( $line = pg_fetch_assoc($result) )
   {
      $_SESSION['userid'] = $line['userid'];
      $liuser = new App_User($_SESSION['userid']);

      $liuser->hasloggedin = 't';
      $liuser->isloggedin = 't';
      //header('Location: home.php');
      header('Location: '.$_POST['redir']);
      exit;
   }
   else
   {
      $error = 'Incorrect Username or Password';
   }
}

$page = new App_Page(array('name' => 'Login') );
$page->header();

$defaultform = array(
         'username' => '',
         'password' => ''
);

$form = array_merge($defaultform, $_POST);

$form['redir'] = isset($_GET['redir']) ? $_GET['redir'] : 'home.php';

?>
<div style="padding-top:1em;width:60%;">
<h2>
Welcome to the new FragMap Server!
</h2>
</div>

<h3>Login</h3>

<?php
if ( isset($error) )
{
	echo '<div id="errorbox">';
	echo "<h2 style='border-bottom-width:0;'>$error</h2>";
	echo '</div>';
}
?>

<form id="loginform" enctype="multipart/form-data" action="login.php" method="post">
<div>
   <label for="username">Username:</label>
   <input type="text" id="username" name="username" value="<?php echo $form['username']; ?>" />
</div>
<div>
   <label for="password">Password:</label>
   <input type="password" id="password" name="password" value="<?php echo $form['password']; ?>" />
</div>
<input type="hidden" id="redir" name="redir" value="<?php echo $form['redir']; ?>" />
<div id="submit">
   <input type="submit" name="action" value="Login" />
</div>
</form>

<br/>
<p><a href="nousername.php">Use FragMap without the benefits of your own account</a></p>

<br />
Need an Account?
<br />
<a href='signup.php'>Sign Up</a>
<br /><br />
Forgot Password?
<br />
<a href='forgotpassword.php'>Reset Password</a>

<?php
$page->footer();
