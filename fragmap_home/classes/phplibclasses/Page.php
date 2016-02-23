<?php
abstract class Page {

   protected $name          = null;
   protected $require_login = false;
   protected $autorefresh   = false;
   protected $js            = null;
   protected $css            = null;

   public function __construct($params = array())
   {
      foreach ($params AS $key=>$value) {
         $this->$key = $value;
      }
      if ($this->require_login) {
         if (! isset($_SESSION['userid']) ) {
            if ($_SERVER['REQUEST_URI'] === '/') {
            	header('Location: login.php');
            } else {
                header('Location: login.php?redir='.htmlentities($_SERVER['REQUEST_URI']));
            }
            exit;
         }
      }
   }

   abstract public function header();

   abstract public function footer();
}
