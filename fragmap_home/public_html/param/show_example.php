<?php
require 'env/path.php';

$example = $_GET['example'];
$EXAMPLES = array(
    "ace",
    "ca",
    "neuraminidase",
    "phospholipase",
    "urokinase",
    "thrombin"
);

if ( ! isset($example) or !in_array($example, $EXAMPLES)) {
    header('Location: examples.php');
    exit();
}

$EXAMPLE_NAMES = [
    "ace" => "Angiotensin Converting Enzyme",
    "ca" => "Carbonic Anhydrase I",
    "neuraminidase" => "Neuraminidase N2",
    "phospholipase" => "Phospholipase A2",
    "urokinase" => "Urokinase Type Plasminogen Activator",
    "thrombin" => "Thrombin"
];

$EXAMPLE_DETAILS = [
    # example => [unbound_pdb, molecule_name, molecule_link, bound_pdb]
    "ace" => ['1o8a', 'Lisinopril', 'http://en.wikipedia.org/wiki/Lisinopril', '1o86'],
    "ca" => ['2cab', 'Dorzolamide', 'http://en.wikipedia.org/wiki/Dorzolamide', '1cil'],
    "neuraminidase" => ['1ivg', 'sialic acid', 'http://en.wikipedia.org/wiki/Sialic_acid', '2bat'],
    "phospholipase" => ['1bbc', 'Niflumic acid', 'http://en.wikipedia.org/wiki/Niflumic_acid', '1td7'],
    "urokinase" => ['2o8t', '8 nM inhibitor', 'http://dx.doi.org/10.1016/S1074-5521(01)00084-9', '1gjc'],
    "thrombin" => ['1hxf', '1-(3-Chlorophenyl)Methanamine', 'http://www.rcsb.org/pdb/ligand/ligandsummary.do?hetId=C2A&sid=2C8Z', '2c8z']
];

$DEFAULT_CLUSTERS = [
    "ace" => [0, 1, 4, 5],
    "ca" => [0, 2, 3],
    "neuraminidase" => [0, 4, 7],
    "phospholipase" => [0, 2, 5, 6],
    "urokinase" => [0, 1, 3],
    "thrombin" => [0]
];

$pse = "examples/$example/$example"."_pse.b64";
$example_dir = "/home/ftmap/public_html/param/examples/".$example;
$example_unbound = $EXAMPLE_DETAILS[$example][0];
$example_lig_name = $EXAMPLE_DETAILS[$example][1];
$example_lig_href = $EXAMPLE_DETAILS[$example][2];
$example_bound = $EXAMPLE_DETAILS[$example][3];

$page = new App_Page( array( 'require_login' => false, 'name' => 'Examples' ) );

$page->header();
?>
<h2><?php echo $EXAMPLE_NAMES[$example]; ?></h2>

<p>
  Mapping of unbound
  <a href="http://www.pdb.org/pdb/explore/explore.do?structureId=<?php echo $example_unbound; ?>" target="_blank">
    PDB <?php echo $example_unbound; ?>
  </a> with
  <a href="<?php echo $example_lig_href; ?>" target="_blank"><?php echo $example_lig_name; ?></a> from
  <a href="http://www.pdb.org/pdb/explore/explore.do?structureId=<?php echo $example_bound; ?>" target="_blank">
    PDB <?php echo $example_bound; ?>
  </a> superimposed.
</p>
<a href='#jsmol_app' id='load_jsmol'>Click to view results using JSmol</a>
<div id="jsmol_app"></div>
<div id="jsmol_control">
<script src="../js/jsmol/JSmol.min.nojq.js"></script>
<script src="../js/jsmol/js/JSmolThree.js"></script>
<script src="../js/jsmol/js/JSmolGLmol.js"></script>

<?php
// read crossclusters for JSMol
chdir($example_dir."/clusters");
$output = array();
$cmd = "ls crosscluster.*.pdb";
exec($cmd, $output);
$clusters=array();
foreach ($output AS $crosscluster)
{
    $parts = explode('.', $crosscluster);
    $clusters[]=$parts[2];
}
?>

<!-- Load pse javascript-->
<script type='text/javascript'>
$(document).ready(function() {
    Jmol.setDocument(0);

    $("#jsmol_control").hide();
    var info = {
        background: '#000000',
            height: 600,      // pixels (but it may be in percent)
            width: '100%',
            use: 'WebGL HTML5',         // 'HTML5' or 'Java' (case-insensitive)
            j2sPath: '../js/jsmol/j2s',           // only used in the HTML5 modality
            script: 'load <?php echo $pse; ?>; background white;',
            addSelectionOptions: false,
            debug: false,
            coverImage: 'examples/<?php echo $example; ?>/result.png',
            coverTitle: 'Click to view model using JSmol',
            deferApplet: true,  // wait to load applet until click
            deferUncover: true, // wait to uncover applet until script completed
            readyFunction: function() {
                $("#load_jsmol").hide();
                $("#jsmol_control").show();
            }
    };

    $('#load_jsmol').click(function(e) {
        info.deferApplet=false;
        $('#jsmol_app').html(Jmol.getAppletHtml('jsmol', info));
    });
    $('#jsmol_app').html(Jmol.getAppletHtml('jsmol', info));

});
</script>

<script language="javascript" type="text/javascript">
document.write('<b>Show/hide clusters:</b> ');
Jmol.jmolCheckbox('jsmol','', '', 'All', true, 'cc_all');
<?php
$jsarg="'cc_all'";
foreach ($clusters as $i => $size) {
    $cc=sprintf('__crosscluster_%03d_%03d', $i, $size);
    $size=sprintf('%d',$size);
    echo "Jmol.jmolCheckbox('jsmol','select {$cc};wireframe 0.3', 'select {$cc};wireframe off', '$i ($size)', true, 'cc$i');\n";
    $jsarg.=sprintf(",'cc%s'",$i);
}
echo "Jmol.jmolSetCheckboxGroup($jsarg);\n";
?>
document.write('<br/><b>Color protein by:</b> ');
Jmol.jmolRadioGroup('jsmol',[
    ['select protein;color group', 'Direction'],
    ['select protein;color structure', 'Secondary structures'],
    ['select protein;color temperature', 'b-factor'],
    //['select all;color chain', 'Chain']
    ] );
document.write('<br/>');
Jmol.jmolCommandInput('jsmol', 'Execute Jmol script', 30);
document.write('<br/>');
</script>
<br/>
<b>Rotate: </b> Left drag |
<b>Zoom: </b> Scroll wheel
</div>

<hr/>

<?php
//prepare data array for google chart
//nonbonded
$lines=file($example_dir."/nonbonded.rawextract", FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);
$js_data_nb='["Sequence","non bonded %", {role: "tooltip"}]';
$cmd = "awk -F'\t' '{sum+=$4} END {print sum}' $example_dir/nonbonded.rawextract";
$output = array();
exec($cmd, $output);
$total_hits=$output[0];
foreach ($lines AS $l) {
    $p = explode("\t", $l);
    $percent=sprintf('%.2f',$p[3]/$total_hits*100);
    $data=array( "$p[1]$p[0]", (float) $percent, "$p[1]/$p[0]/$p[2]: $percent%" );
    $js_data_nb.=",\n".json_encode($data);
}

$lines=file("$example_dir/hbonded.rawextract", FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);
$js_data_hb='["Sequence","Hbond %", {role: "tooltip"}]';
$cmd = "awk -F'\t' '{sum+=$4} END {print sum}' $example_dir/nonbonded.rawextract";
$output = array();
exec($cmd, $output);
$total_hits=$output[0];
foreach ($lines AS $l) {
    $p = explode("\t", $l);
    $percent=sprintf('%.2f',$p[3]/$total_hits*100);
    $data=array( "$p[1]$p[0]", (float) $percent, "$p[1]/$p[0]/$p[2]: $percent%" );
    $js_data_hb.=",\n".json_encode($data);
}
?>
<!-- Google Charts API for data plotting -->
<script type="text/javascript" src="https://www.google.com/jsapi"></script>
<script type="text/javascript">
google.load("visualization", "1", {packages:["corechart"]});
google.setOnLoadCallback(drawChart);
function drawChart() {
    var data_nb = google.visualization.arrayToDataTable([ <?php echo "$js_data_nb\n"; ?> ]);
    var options_nb = {
        hAxis: {title: 'Sequence', textStyle: {fontSize: 8}}
    };
    var chart_nb = new google.visualization.ColumnChart(document.getElementById('chart_div_nb'));
    chart_nb.draw(data_nb, options_nb);
<?php if (isset($js_data_hb)) { ?>
    var data_hb = google.visualization.arrayToDataTable([ <?php echo "$js_data_hb\n"; ?> ]);
    var options_hb = {
        hAxis: {title: 'Sequence', textStyle: {fontSize: 8}}
    };
    var chart_hb = new google.visualization.ColumnChart(document.getElementById('chart_div_hb'));
    chart_hb.draw(data_hb, options_hb);
<?php } ?>
}
</script>

<div id='chart_div_nb'></div>
<a id='nblink' href='examples/<?php echo $example; ?>/nonbonded.rawextract'>
  Nonbonded Interactions
</a>
<br/>

<hr/>
<div id='chart_div_hb'></div>
<a id='hblink' href='examples/<?php echo $example; ?>/hbonded.rawextract'>
  Hbond Interactions
</a>
<br/>

<hr/>
<a id='summarylink' href='examples/<?php echo $example; ?>/crossclustersummary'>
  Probes summary
</a>
<br/>

<?php
$page->footer();
