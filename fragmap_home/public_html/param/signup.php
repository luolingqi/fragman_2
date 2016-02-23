<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'sessionstart.php';

if ( isset($_POST['action']) )
{
   require 'createaccount.php';
}

$js = '
$(document).ready(function(){
   $("#refreshcaptcha").click(function(){
      $("#wordimage").removeAttr("src").attr("src","image/captcha.php?"+Math.random());
   });
});
';

$page = new App_Page( array('name' => 'Signup', 'js' => $js) );
$page->header();

$defaultform = array(
   'firstname'   => '',
   'lastname'    => '',
   'username'    => '',
   'affiliation' => '',
   'email'       => '',
   'tos'         => ''
);

$form = array_merge($defaultform, $_POST);

?>
<h3>Create Account</h3>

<?php
if (isset($errors) && sizeof($errors) > 0)
{
	echo '<div id="errorbox">';
	echo '<h2>Sorry, your registration contained errors:</h2>';
   echo '<ul>';
   foreach ($errors AS $error)
   {
      echo "<li>$error</li>";
   }
   echo '</ul>';
	echo '</div>';
}

?>

<form id="signupform" enctype='multipart/form-data' action='signup.php' method='post'>

<div>
   <label for="firstname">First Name:</label>
   <input type="text" name="firstname" id="firstname" value="<?php echo $form['firstname']; ?>" />
</div>

<div>
   <label for="lastname">Last Name:</label>
   <input type="text" name="lastname" id="lastname" value="<?php echo $form['lastname']; ?>" />
</div>

<div>
   <label for="username">Username:</label>
   <input type="text" name="username" id="username" value="<?php echo $form['username']; ?>" />
</div>

<div>
   Examples: JSmith, John.Smith
</div>

<div>
   <label for="affiliation">Affiliation:</label>
   <input type="text" name="affiliation" id="affiliation" value="<?php echo $form['affiliation']; ?>" />
</div>

<div>
   Example: Boston University
</div>

<div>
   <label for="email">Email:</label>
   <input type="text" name="email" id="email" value="<?php echo $form['email']; ?>" />
</div>

<div id="emailrequirements">
 Your email must be an educational or government address.
 Your login password will be sent to this address
 (you can change this password after you've logged in).
</div>

<div>
   <label for="wordverification">Word Verification:</label>
   <input type="text" name="wordverification" id="wordverification" />
</div>
<div>
   <span style="float:left;width:50%;">Type the characters seen in the picture </span>
   <img src="image/captcha.php" alt="" id="wordimage" /><br /><a id="refreshcaptcha">Refresh Image</a>
</div>
<div style="width:150%;">
   <input type="checkbox" name="tos" id="tos" <?php echo ($form['tos']) ? 'checked="checked"' : ''; ?> value="checked" />
   <label for="tos">I agree to use <?php echo SERVER_NAME; ?> only for noncommercial purposes.</label>
</div>
<div>
   
</div>
<div id="submit">
   <input type="submit" name="action" value="Create Account" />
</div>
</form>

<br />
Already Have an Account?
<br />
<a href="login.php">Login</a>


<?php
$page->footer();
?>
