$(function() {
    $('input:checkbox').button();
$('input[id^="cc"]').change(function() {
  "use strict";
  var cc_no = $(this).attr('id').substr(2);
  if($(this).is(':checked')) {
                document.av.execute('display lines on molecule cc'+cc_no+';');
  }else
  {
                document.av.execute('display lines off molecule cc'+cc_no+';');
  }

});
$('input[id^="lig"]').change(function() {
  "use strict";
  var lig_no = $(this).attr('id').substr(3);
  if($(this).is(':checked')) {
                document.av.execute('display sticks on molecule lig'+lig_no+';');
  }else
  {
                document.av.execute('display sticks off molecule lig'+lig_no+';');
  }

});
});
