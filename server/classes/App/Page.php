<?php
class App_Page extends Page {

   protected $name          = null;
   protected $require_login = false;
   protected $autorefresh   = false;
   protected $js            = null;

   public function __construct($params = array())
   {
      parent::__construct($params);
   }

   public function header()
   {
      ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN" "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <title><?php echo SERVER_TITLE; ?></title>
    <?php
     if ($this->autorefresh) {
        echo '<meta http-equiv="refresh" content="300"/>';
     } ?>
     <meta http-equiv="content-type" content="text/html; charset=utf-8" />
     <link rel='stylesheet' type='text/css' href='css/style.css' />
     <link rel='stylesheet' type='text/css' href='css/loginform.css' />
     <link rel='stylesheet' type='text/css' href='css/signupform.css' />
     <link rel='stylesheet' type='text/css' href='css/contactform.css' />
     <link rel='stylesheet' type='text/css' href='css/jobsform.css' />
     <link rel='stylesheet' type='text/css' href='css/goodform.css' />
     <link rel="stylesheet" type="text/css" href="css/grids-min.css" />
     <?php
         if (! is_null($this->name) )
         {
            echo '<style type="text/css">';
            echo "#tab$this->name { font-weight:bold; }";
            echo '</style>';
         }
     ?>
     <link rel="shortcut icon" href="favicon.png" type="image/png" />
     <script src="http://code.jquery.com/jquery-1.11.0.min.js"></script>
     <?php
         if (! is_null($this->js) )
         {
            echo <<<JS
               <script type="text/javascript">
                  $this->js
               </script>
JS;
         }
         ?>

  </head>

  <body>
  <div id="doc">
    <div id="hd">
        <ul id='tabs-menu'>
          <li><a id='tabContact' href='contact.php'>Contact</a></li>
          <li><a id='tabPapers' href='publications.php'>Papers</a></li>
          <!--<li><a id='tabTutorial' href='tutorial.php'>Tutorial</a></li>-->
          <li><a id='tabHelp' href='serverhelp.php'>Help</a></li>
          <li><a id='tabExamples' href='examples.php'>Examples</a></li>
<?php if ( !isset($_SESSION['userid']) ) { ?>
          <li><a id='tabSignup' href='signup.php'>Sign Up</a></li>
          <li><a id='tabLogin' href='login.php'>Login</a></li>
<?php }
      if ( $this->require_login || isset($_SESSION['userid']) ) {
         $user = new App_User($_SESSION['userid']);
?>
<?php
         if ( $user->isSuperUser() )
         {
            echo "<li><a id='tabAdmin' href='admin.php'>Admin</a></li>";
         }
       ?>
<?php if ( $user->username != 'piper' ) { ?>
          <li><a id='tabPreferences' href='preferences.php'>Preferences</a></li>
<?php } ?>
          <li><a id='tabResults' href='results.php'>Results</a></li>
          <li><a id='tabQueue' href='queue.php'>Queue</a></li>
          <li><a id='tabMap' href='home.php'>Map</a></li>
<?php } ?>
        </ul>
      <img src='image/banner_ftmap.png' width='750' height='160' alt=''/>
    </div>
      <div id="bd">
      <!--<p style="margin: 1em auto;padding:0; width;70%">
      <strong><font color='red'>Warning: </font></strong> There is a <a href="http://www.bu.edu/tech/support/research/whats-happening/updates/">scheduled maintenance</a> on Monday, August 11. FTMap will be unavailable from Sunday evening to Tuesday morning.
      </p> -->
     <?php
     if ( $this->require_login || isset($_SESSION['userid']) ) {
       echo "
        <div id='main-header-right'>
          <a href='logout.php'>sign out</a>
        </div>
       ";
     }
   }

   public function footer()
   { ?>
</div>
<div id="ft">
  <a href='http://structure.bu.edu/'>Structural Bioinformatics Lab</a>
  <br/>
  <a href='http://www.bu.edu/'>Boston University</a>
</div>
</div>
</body>
</html>
<?php
   }
}
