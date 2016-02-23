<?php
//session_set_cookie_params(0, '/', 'ftmap.bu.edu', false, true);
session_start();

if (isset($_SESSION['userid']) )
{
   $liuser = new App_User($_SESSION['userid']);
}
