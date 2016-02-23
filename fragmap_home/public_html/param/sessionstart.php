<?php
session_start();

if (isset($_SESSION['userid']) )
{
   $liuser = new App_User($_SESSION['userid']);
}
