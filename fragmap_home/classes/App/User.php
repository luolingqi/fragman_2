<?php
class App_User extends User {

   public  $userid     = 0;
   protected $properties = array();

   public function __construct($id)
   {
      parent::__construct($id);
   }
}
