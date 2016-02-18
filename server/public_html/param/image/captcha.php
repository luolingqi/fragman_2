<?php
session_start();


$x = 120;
$y = 40;
$nchars = 6;
$fontid = imageloadfont('backlash.gdf');

$captchaimg = imagecreate($x, $y) or die("error imagecreate\n");
$background_color = imagecolorallocate($captchaimg, 255, 255, 255); // (first call to imagecolorallocate sets background color)
//$background_color = imagecolorallocate($captchaimg, 220, 220, 220); // (first call to imagecolorallocate sets background color)
$text_color = imagecolorallocate($captchaimg, 0, 0, 0);

$possible = '23456789bcdfghjkmnpqrstvwxyz';
$string = '';
$i = 0;

while ($i < $nchars)
{ 
	$string .= substr($possible, mt_rand(0, strlen($possible)-1), 1);
	$i++;
}

imagestring($captchaimg, $fontid, 0, 0, $string, $text_color);
$angle = mt_rand(-15,15);
$captchaimg = imagerotate($captchaimg, $angle, $background_color);
$y = imagesy($captchaimg);
$x = imagesx($captchaimg);


//$noise_color = imagecolorallocate($captchaimg, 190, 190, 190);
//$noise_color = imagecolorallocate($captchaimg, 90, 90, 90);
$noise_color = imagecolorallocate($captchaimg, 0, 0, 0);


// generate random dots in background 
//for( $i=0; $i<($x*$y)/3; $i++ ) {
//for( $i=0; $i<($x*$y)/5; $i++ ) {
for ( $i=0; $i<($x*$y)/5; $i++ ) {
	imagefilledellipse($captchaimg, mt_rand(0,$x), mt_rand(0,$y), 1, 1, $noise_color);
}
// generate random lines in background 
//for( $i=0; $i<($x*$y)/150; $i++ ) {
for ( $i=0; $i<($x*$y)/500; $i++ ) {
	imageline($captchaimg, mt_rand(0,$x), mt_rand(0,$y), mt_rand(0,$x), mt_rand(0,$y), $noise_color);
}

/*
imagefilter ($captchaimg, IMG_FILTER_EMBOSS);
imagefilter ($captchaimg, IMG_FILTER_GAUSSIAN_BLUR);
imagefilter ($captchaimg, IMG_FILTER_GAUSSIAN_BLUR);
imagefilter ($captchaimg, IMG_FILTER_GAUSSIAN_BLUR);
*/


//   Encrypt and store the key inside of a session
$_SESSION['key'] = md5($string);
//$_SESSION['key'] = $string;

header('Content-Type: image/png');
imagepng($captchaimg);

?>
