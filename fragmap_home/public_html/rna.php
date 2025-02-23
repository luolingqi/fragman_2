<?php
require 'env/appvars.php';
require 'env/dbvars.php';
require 'sessionstart.php';

$page = new App_Page( array('name' => 'RNA' ) );
$page->header();
?>
<br>

<h3>Supported RNA Residue and Atom Names</h3>
<br>

adenine:<br>
<pre>
ATOM    195  P   A      19      38.769 -32.636 -69.811  1.00 20.00
ATOM    196  OP1 A      19      37.589 -31.806 -70.153  1.00 20.00
ATOM    197  OP2 A      19      39.783 -32.084 -68.881  1.00 20.00
ATOM    198  O5' A      19      38.277 -34.029 -69.261  1.00 20.00
ATOM    199  C5' A      19      37.580 -34.939 -70.120  1.00 20.00
ATOM    200  H5' A      19      38.178 -35.173 -71.031  1.00  0.00
ATOM    201 H5'' A      19      36.624 -34.466 -70.440  1.00  0.00
ATOM    202  C4' A      19      37.264 -36.221 -69.364  1.00 20.00
ATOM    203  H4' A      19      36.371 -36.722 -69.794  1.00  0.00
ATOM    204  O4' A      19      38.352 -37.174 -69.486  1.00 20.00
ATOM    205  C1' A      19      38.520 -37.878 -68.282  1.00 20.00
ATOM    206  H1' A      19      38.426 -38.957 -68.524  1.00  0.00
ATOM    207  N9  A      19      39.861 -37.623 -67.758  1.00 20.00
ATOM    208  C5  A      19      41.867 -37.778 -66.862  1.00 20.00
ATOM    209  N7  A      19      41.732 -36.491 -67.299  1.00 20.00
ATOM    210  C8  A      19      40.538 -36.431 -67.826  1.00 20.00
ATOM    211  H8  A      19      40.096 -35.548 -68.296  1.00  0.00
ATOM    212  N1  A      19      42.734 -39.755 -65.905  1.00 20.00
ATOM    213  C2  A      19      41.573 -40.313 -66.237  1.00 20.00
ATOM    214  H2  A      19      41.489 -41.361 -65.964  1.00  0.00
ATOM    215  N3  A      19      40.499 -39.803 -66.837  1.00 20.00
ATOM    216  C4  A      19      40.715 -38.503 -67.139  1.00 20.00
ATOM    217  C6  A      19      42.929 -38.438 -66.214  1.00 20.00
ATOM    218  N6  A      19      44.084 -37.823 -65.906  1.00 20.00
ATOM    219  H61 A      19      44.238 -36.833 -66.081  1.00  0.00
ATOM    220  H62 A      19      44.948 -38.319 -65.727  1.00  0.00
ATOM    221  C2' A      19      37.421 -37.435 -67.341  1.00 20.00
ATOM    222 H2'' A      19      37.831 -37.436 -66.350  1.00  0.00
ATOM    223  O2' A      19      36.313 -38.302 -67.364  1.00 20.00
ATOM    224  H2' A      19      35.637 -37.725 -66.980  1.00  0.00
ATOM    225  C3' A      19      37.066 -36.053 -67.861  1.00 20.00
ATOM    226  H3' A      19      37.817 -35.331 -67.461  1.00  0.00
ATOM    227  O3' A      19      35.717 -35.644 -67.598  1.00 20.00
</pre>

<br>
cytosine:<br>
<pre>
ATOM    164  P   C      23      36.435 -32.568 -75.197  1.00 20.00
ATOM    165  OP1 C      23      35.759 -33.788 -74.695  1.00 20.00
ATOM    166  OP2 C      23      35.665 -31.303 -75.246  1.00 20.00
ATOM    167  O5' C      23      37.750 -32.326 -74.367  1.00 20.00
ATOM    168  C5' C      23      38.716 -33.376 -74.221  1.00 20.00
ATOM    169  H5' C      23      38.953 -33.907 -75.158  1.00  0.00
ATOM    170 H5'' C      23      38.263 -34.152 -73.556  1.00  0.00
ATOM    171  C4' C      23      39.949 -32.862 -73.494  1.00 20.00
ATOM    172  H4' C      23      40.662 -33.686 -73.260  1.00  0.00
ATOM    173  O4' C      23      40.667 -31.918 -74.307  1.00 20.00
ATOM    174  C1' C      23      41.229 -30.926 -73.480  1.00 20.00
ATOM    175  H1' C      23      42.327 -31.114 -73.552  1.00  0.00
ATOM    176  N1  C      23      40.909 -29.586 -74.002  1.00 20.00
ATOM    177  C6  C      23      39.695 -29.335 -74.596  1.00 20.00
ATOM    178  H6  C      23      38.994 -30.172 -74.605  1.00  0.00
ATOM    179  C5  C      23      39.410 -28.129 -75.107  1.00 20.00
ATOM    180  H5  C      23      38.450 -27.924 -75.580  1.00  0.00
ATOM    181  C2  C      23      41.854 -28.582 -73.887  1.00 20.00
ATOM    182  O2  C      23      42.923 -28.838 -73.324  1.00 20.00
ATOM    183  N3  C      23      41.586 -27.358 -74.392  1.00 20.00
ATOM    184  C4  C      23      40.408 -27.110 -74.994  1.00 20.00
ATOM    185  N4  C      23      40.188 -25.880 -75.492  1.00 20.00
ATOM    186  H41 C      23      40.933 -25.202 -75.393  1.00  0.00
ATOM    187  H42 C      23      39.567 -25.748 -76.264  1.00  0.00
ATOM    188  C2' C      23      40.781 -31.156 -72.042  1.00 20.00
ATOM    189 H2'' C      23      40.452 -30.226 -71.531  1.00  0.00
ATOM    190  O2' C      23      41.825 -31.750 -71.305  1.00 20.00
ATOM    191  H2' C      23      41.504 -31.840 -70.393  1.00  0.00
ATOM    192  C3' C      23      39.609 -32.106 -72.215  1.00 20.00
ATOM    193  H3' C      23      38.668 -31.527 -72.356  1.00  0.00
ATOM    194  O3' C      23      39.494 -33.049 -71.153  1.00 20.00
</pre>

<br>
guanine:<br>
<pre>
ATOM    130  P   G      97      34.258 -29.771 -79.476  1.00 20.00
ATOM    131  OP1 G      97      33.165 -30.596 -78.896  1.00 20.00
ATOM    132  OP2 G      97      34.317 -28.335 -79.120  1.00 20.00
ATOM    133  O5' G      97      35.657 -30.416 -79.128  1.00 20.00
ATOM    134  C5' G      97      35.936 -31.782 -79.451  1.00 20.00
ATOM    135  H5' G      97      35.804 -32.016 -80.514  1.00  0.00
ATOM    136 H5'' G      97      35.206 -32.417 -78.895  1.00  0.00
ATOM    137  C4' G      97      37.314 -32.171 -78.940  1.00 20.00
ATOM    138  H4' G      97      37.525 -33.251 -79.095  1.00  0.00
ATOM    139  O4' G      97      38.354 -31.444 -79.633  1.00 20.00
ATOM    140  C1' G      97      39.371 -31.080 -78.719  1.00 20.00
ATOM    141  H1' G      97      40.319 -31.530 -79.094  1.00  0.00
ATOM    142  N9  G      97      39.520 -29.620 -78.727  1.00 20.00
ATOM    143  C4  G      97      40.633 -28.887 -78.401  1.00 20.00
ATOM    144  N2  G      97      43.940 -28.724 -77.281  1.00 20.00
ATOM    145  H21 G      97      44.653 -28.026 -77.126  1.00  0.00
ATOM    146  H22 G      97      44.238 -29.681 -77.264  1.00  0.00
ATOM    147  N3  G      97      41.794 -29.378 -77.949  1.00 20.00
ATOM    148  C2  G      97      42.718 -28.422 -77.748  1.00 20.00
ATOM    149  N1  G      97      42.466 -27.083 -77.993  1.00 20.00
ATOM    150  H1  G      97      43.214 -26.417 -77.827  1.00  0.00
ATOM    151  C6  G      97      41.274 -26.551 -78.475  1.00 20.00
ATOM    152  O6  G      97      41.145 -25.336 -78.665  1.00 20.00
ATOM    153  C5  G      97      40.299 -27.571 -78.664  1.00 20.00
ATOM    154  N7  G      97      38.996 -27.478 -79.115  1.00 20.00
ATOM    155  C8  G      97      38.566 -28.713 -79.132  1.00 20.00
ATOM    156  H8  G      97      37.568 -29.033 -79.436  1.00  0.00
ATOM    157  C2' G      97      39.025 -31.672 -77.362  1.00 20.00
ATOM    158 H2'' G      97      39.309 -31.005 -76.524  1.00  0.00
ATOM    159  O2' G      97      39.674 -32.915 -77.204  1.00 20.00
ATOM    160  H2' G      97      39.438 -33.210 -76.317  1.00  0.00
ATOM    161  C3' G      97      37.516 -31.836 -77.459  1.00 20.00
ATOM    162  H3' G      97      37.025 -30.859 -77.240  1.00  0.00
ATOM    163  O3' G      97      36.984 -32.896 -76.655  1.00 20.00
</pre>

<br>
uracil:<br>
<pre>
ATOM     32  P   U     694      37.267 -15.812 -89.135  1.00 20.00
ATOM     33  OP1 U     694      36.245 -14.839 -88.690  1.00 20.00
ATOM     34  OP2 U     694      38.445 -16.062 -88.272  1.00 20.00
ATOM     35  O5' U     694      36.572 -17.187 -89.443  1.00 20.00
ATOM     36  C5' U     694      35.493 -17.289 -90.375  1.00 20.00
ATOM     37  H5' U     694      35.639 -16.668 -91.279  1.00  0.00
ATOM     38 H5'' U     694      34.580 -16.892 -89.870  1.00  0.00
ATOM     39  C4' U     694      35.214 -18.752 -90.708  1.00 20.00
ATOM     40  H4' U     694      34.395 -18.858 -91.452  1.00  0.00
ATOM     41  O4' U     694      36.384 -19.374 -91.283  1.00 20.00
ATOM     42  C1' U     694      36.503 -20.704 -90.811  1.00 20.00
ATOM     43  H1' U     694      36.490 -21.354 -91.714  1.00  0.00
ATOM     44  N1  U     694      37.780 -20.878 -90.097  1.00 20.00
ATOM     45  C6  U     694      38.457 -19.786 -89.600  1.00 20.00
ATOM     46  H6  U     694      37.904 -18.864 -89.603  1.00  0.00
ATOM     47  C2  U     694      38.275 -22.155 -89.918  1.00 20.00
ATOM     48  O2  U     694      37.702 -23.154 -90.341  1.00 20.00
ATOM     49  N3  U     694      39.478 -22.214 -89.251  1.00 20.00
ATOM     50  H3  U     694      39.812 -23.140 -88.980  1.00  0.00
ATOM     51  C4  U     694      40.211 -21.160 -88.743  1.00 20.00
ATOM     52  O4  U     694      41.282 -21.363 -88.178  1.00 20.00
ATOM     53  C5  U     694      39.629 -19.869 -88.967  1.00 20.00
ATOM     54  H5  U     694      40.101 -18.983 -88.545  1.00  0.00
ATOM     55  C2' U     694      35.296 -21.001 -89.940  1.00 20.00
ATOM     56 H2'' U     694      35.548 -21.631 -89.056  1.00  0.00
ATOM     57  O2' U     694      34.288 -21.646 -90.692  1.00 20.00
ATOM     58  H2' U     694      33.632 -21.859 -90.009  1.00  0.00
ATOM     59  C3' U     694      34.862 -19.611 -89.497  1.00 20.00
ATOM     60  H3' U     694      35.469 -19.266 -88.628  1.00  0.00
ATOM     61  O3' U     694      33.459 -19.503 -89.233  1.00 20.00
</pre>






<br>
<br>
<br>
<br>
<?php
$page->footer();
