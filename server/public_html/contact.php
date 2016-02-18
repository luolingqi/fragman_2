<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'sessionstart.php';

if ( isset($_POST['action']) )
{
   if (! preg_match('/^[_a-z0-9-]+(\.[_a-z0-9-]+)*@[a-z0-9-]+(\.[a-z0-9-]+)*(\.[a-z]{2,3})$/i',
			$_POST['email'])
		)
   {
	   $errors[] = 'The provided email address is invalid.';
   }
   // word verification
   // note the generated string is always lowercase
   if (md5(strtolower($_POST['wordverification'])) != $_SESSION['key'])
   {
	   $errors[] = 'The word verification code was incorrect. Please try again.';
   }
   if (! isset($errors) || sizeof($errors) === 0 )
   {
   
      //send email to ourselves
      $to = SERVER_EMAIL;
      $subject = 'Contact Form';
      $message = "$_POST[message]";
      $headers = "From: $_POST[email]";
      mail($to, $subject, $message, $headers);
      
      //send email to them
      $to = $_POST['email'];
      $subject = "$servername Contact Form";
      $message = "Thank you your sending a message to $servername.  We appreciate your feedback.";
      $from = "$servername <no-reply@bu.edu>";
      $headers = "From: $from";
      mail($to, $subject, $message, $headers);
   
      $successes[] = 'yes';
      $_POST = array();
   }
}

$js = '
$(document).ready(function(){
   $("#refreshcaptcha").click(function(){
      $("#wordimage").removeAttr("src").attr("src","image/captcha.php?"+Math.random());
   });
});
';

$page = new App_Page( array('name' => 'Contact', 'js' => $js) );

$page->header();

if ( isset($successes) && sizeof($successes) >0 )
{
   echo '<div id="successbox" class="jobsform">';
   echo '<h2 style="border-bottom:0;">Message Sent</h2>';
   echo '</div>';
}

if ( isset($errors) && sizeof($errors) >0 )
{
   echo '<div id="errorbox">';
   echo '<h2>There were some errors in processing your message:</h2>';
   echo '<ul>';
   foreach ($errors AS $error)
   {
      echo "<li>$error</li>";
   }
   echo '</ul>';
   echo '</div>';
}

$defaultform = array(
   'email'   => '',
   'message' => '',
);

$form = array_merge($defaultform, $_POST);
   
?>
<h3>Conctact Us</h3>

<p style="width:70%;margin:auto;text-align:left;">Please use the <a href='http://structure.bu.edu/contact'>contact form</a> on our lab page and select "FTMap" as the category.<br/>
If contacting us about a particular job, please include the job number in your message.
</p>
<?php

$page->footer();
