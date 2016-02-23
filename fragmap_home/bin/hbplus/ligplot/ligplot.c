/***********************************************************************

  Ligplot.c - Program to produce a schematic 2D plot of a given protein
              ligand, showing the protein residues to which it is
              hydrogen-bonded and those with which it makes hydrophobic
              contact

************************************************************************

Version:-      v.4.0  -  16 November 1997

Written by:-   Andrew Wallace & Roman Laskowski

               Department of Biochemistry and Molecular Biology
               University College London
               Gower Street, London WC1E 6BT, UK.

Additional
routines:-     qikfit and matfit - least-squares fitting routines by
(v.3.0)        Andrew Martin, Dept. Biochem. & Mol. Biol., University
	       College London.

               Eigen value calculation routines by Ian Tickle & David
               Moss, Dept. Crystallography, Birkbeck College.

Reference:-    Wallace A C, Laskowski R A & Thornton J M .
               LIGPLOT: A program to generate schematic diagrams
               of protein-ligand interactions. Prot. Eng., 8,
               127-134.

Further info:- http://www.biochem.ucl.ac.uk/bsm/ligplot/ligplot.html

Requirements:- The list of hydrogen bonds and non-bonded contacts
               required by the program can either be user-generated
               in the format described in the Operating Manual, or
               produced by the HBPLUS program, written by Ian McDonald.
               To obtain a copy of the latter program, see

               http://sjh.bi.umist.ac.uk/naccess.html

               The ligand atoms can also be shaded (coloured) according
               to their solvent accessibility. Accessibilities can be
               calculated using the NACCESS program, written by
               Simon Hubbard. To obtain a copy, see

               http://www.biochem.ucl.ac.uk/bsm/naccess/naccess.html

Amendments:-   Amendments after v.3.0 are labelled by version number
               v.m.n--> and v.m.n<-- where m.n is the version
               number corresponding to the change

v.1.0          Original release version.                 September 1994
                                     (Andrew Wallace & Roman Laskowski)

v.2.0          Revised input parameters, including many new parameters,
               particularly those defining the colours and sizes of the
               different items and labels plotted. Portrait/landscape
               option added. Name of input parameter file changed from
               ligplot.prm to ligplot.par to prevent confusion with
               former version.
               Changes to PostScript output routines to make the output
               PostScript file simpler to understand and edit.
               Several bugs in v.1.0 found and fixed.
	       Double- and triple-bond option added (these have to be
               manually defined, though).
               Read-in and use of CONECT records in input PDB file.
               Output of ligplot.frm file giving original coordinates,
               in PDB format, of all the residues shown on the LIGPLOT
               diagram.
                                                         Jan 1996 (RAL)

v.2.0.1        Addition of routine to check for hydrogen atoms so that
               they can be discarded.                  2 Feb 1996 (RAL)

               Replacement of "connect" array by the bond structure
               to improve program's storage efficiency. New routines:
               update_bond, get_bond_type.
               Extra option (only used when parameter Include_Hydrogens
               set to TRUE and program recompiled) to include hydrogen
               atoms on the plot (Note: corresponding H-bonds have to be
               defined to H-atoms in the .hhb file).
                                                      26 Feb 1996 (RAL)

v.2.1          Bug-fix to incorrect conversion of atom-names from ligplot-
               hbond file format to .hhb hbond file format.
               Addition of read_conec_records and verify_conec_records
               routines to pick up all the non-ligand residues that
               are covalently bonded to the ligand (as per the CONECT
               records), and to store these bonds as pseudo H-honds.
                                                   11-14 Nov 1996 (RAL)
               Routines to store all ligand residue-names and numbers,
               as encountered in get_ligand_start_end, to compare against
               when reading in the H-bonds.
                                                      21 Nov 1996 (RAL)
v.3.0          Major revisions, including change of ligplot.par back
               to ligplot.prm.
                                             29 Jan - 31 Mar 1997 (RAL)
v.3.0.1        Bug-fix to routine determining whether line read in from
               HBPLUS file is valid. Lines with blanks at the start were
               ignored, but HBPLUS fills line-starts with blanks if it
               cannot identify a PDB 4-letter code in the name of the
               input PDB file.
                                                       8 Apr 1997 (RAL)
v.3.1          Extension of combine_objects routine to ligand objects,
               specifically to handle phosphotyrosines in the ligand. In
               the process, also fixed a bug in the join_residues function
               which goes wrong if block of residues or atoms to be moved
               are already adjacent in the list.
                                                      11 Apr 1997 (RAL)
               Amendments to hbadd.c program.
                                                      13 Apr 1997 (RAL)
               Improvements to some of the most computationally
               intensive routines to speed up whole program.
                                                      14 Apr 1997 (RAL)
               Check for change in residue name as well as number in
               read_pdb_file function, otherwise H-bond info not picked
               up for HET group attachements to residues
               (eg phosphotyrosines)
                                                      15 Apr 1997 (RAL)
               Improved program's ability to cope with duplicate
               residue numbering in PDB files.
                                                      16 Apr 1997 (RAL)
v.3.1.1        Bug-fix for certain cases where cyclic peptides identified
               (eg covalently-linked sugars in 1gya). New function to
               check for such cases prior to the split_object function
               so that the different units are split into separate
               objects for ease of flattening.
               Check for ENDMDL added to function get_ligand_start_end
               so that warning messages not repeated unnecessarily.
                                                      20 Apr 1997 (RAL)
v.3.1.2        Amendment primarily for PDBsum. Where ligand residues
               are duplicated and last "apparent" ligand is consists of
               a single atom only, then take the "true" ligand to be the
               penultimate one, rather than the last.
                                                      22 Apr 1997 (RAL)
               Speeded up processing of CONECT records for NMR
               structures having large numbers of them.
               Additional flag to specify that ligand is a metal.
                                                      23 Apr 1997 (RAL)
v.3.2          Bug-fix in stretch_mainchain function causing crash when
               no carbonyl oxygen present.
                                                      28 Apr 1997 (RAL)
               Addition of graph-matching routines to hbadd.c.
                                                       1 May 1997 (RAL)
               Bug-fix. When two residues were combined into a hydrophobic
               group, and then one of the residues' atoms all deleted as
               unwanted, gave problems in the plot with the residue name
               displaced from the "eye-lash" representing it.
               "Hydrophobic contacts" extended to include non-bonded
               contact between any pair of atoms where at least one of
               the pair is a carbon or sulphur.
                                                       2 May 1997 (RAL)
               Addition of two new parameters to the ligplot.prm file,
               one dealing with the definition of non-bonded contacts and
               the other with a final rotation of the picture.
               Removal of one of the warning messages.
                                                       8 May 1997 (RAL)
               Amendment to distance-energy for non-bonded contacts.
               Amendment to read-in of CONECT records to check for early
               line-end (as in output from the Babel program).
                                                      12 May 1997 (RAL)
               Amendment to initially place all hydrophobic interactions
               a long way from the ligand to give the H-bonded groups a
               better chance of finding better positions.
               New option to show covalent bonds to external groups as
               solid lines rather than as dotted lines.
                                                      13 May 1997 (RAL)
               Bug-fix. Residues belonging to deleted objects now marked
               as deleted to prevent them being written out to the
               ligplot.frm file.
               Bug-fix. Problems in transformation routines when ligand
               consisting of a single atom only.
                                                      15 May 1997 (RAL)
               New output file, ligplot.res, listing all the residues
               included in the plot (for use by the raslig program in
               the PDBsum LIGPLOT pages).
               Improvements to function that interprets the command-line
               parameters.
                                                      16 May 1997 (RAL)
               Bug-fix to lig_search routines which were failing to find
               and display the first "cofactor".
                                                      19 May 1997 (RAL)
               Bug-fix to remove any bonds between hydrophobic groups 
               and any non-ligand residues.
                                                      25 May 1997 (RAL)
               Amendment to cope with real problem cases where squashing
               flat of the ligand is required (eg as in 2cmm). Function
               align_object added to orientate the ligand as best as
               possible before it is brutally forced flat.
               Added check for Hg and Ho in hydrogen-check.
                                                      28 May 1997 (RAL)
v.3.2.1        Addition of marker to identify the ligand residues in the
               output ligplot.res file.
                                                      19 Jun 1997 (RAL)
v.4.0          Addition of routines to read in interactions details across
               a dimer or domain interface, as generated by Alice
               Wuensche's program, dimer.c
                                             28 Jun - 10 Aug 1997 (RAL)
               Addition of output ligplot.rcm, which holds the residue
               centres of mass for use in generating similar LIGPLOT
               layouts for similar input structures. Routines to read in
               the ligplot.rcm file and use its data for the minimization
               procedure.
                                                      10 Aug 1997 (RAL)
               Amendments to hbadd (rewrite of clique-detection routine).
                                                   16-18 Sep 1997 (RAL)
               Addition of check that maximum accessibility non-zero,
               otherwise led to divide by zero. Spotted by Dave Renouf.
                                                      21 Sep 1997 (RAL)
               Amendment of .rcm file routines to have centres of mass
               of residues, rather than object, and to allow restraint
               of individual atoms also.
                                                      16 Oct 1997 (RAL)
               Increase number of atoms stored from 1200 to 5000 as was
               not enough for some structures (eg NAD 336 of 1gd1).
               Fix to exclude CONECT records giving covalent bonds
               between non-ligand residues (as in 1gd1).
               Minor fixes to ligplot.res file (only used by PDBsum).
                                                      10 Nov 1997 (RAL)
               Additional command-line input of residue names, as well as
               residue numbers, for defining the range of the ligand's
               residues. Particularly useful where PDB file contains
               multiple residues having exactly the same residue number.
               Amendment to prevent protein residues being incorporated
               into the ligand if they have the same residue number and
               chain ID as a residue in the ligand.
                                                      11 Nov 1997 (RAL)
               Addition of list of covalent bonds between ligand and
               protein residues to .res file.
                                                      12 Nov 1997 (RAL)
               Addition of option to parameters to allow water atoms to
               be shown as spheres, or not.
                                                      16 Nov 1997 (RAL)
               Storage of maximum accessibility value in the ligplot.pdb
               file so that shading preserved when plot regenerated.
                                                      18 Nov 1997 (RAL)
v.4.0.1        Bug-fix. Ligand atoms were being printed even if option
               not to print them had been selected.
                                                       3 Dec 1997 (RAL)
v.4.0.2        Bug-fix to deal with cases where ligand residues are
               numbered backwards (as in 9lpr).
               Bug-fix to increase maximum number of command-line
               parameters that can be stored (MAXTOKENS) as giving
               problems with certain structures (eg MG 435 -> HAD 438
               in 1gin).
                                                       8 Dec 1997 (RAL)


----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name       File pointer    Description
-----------       ------------    -----------
<filename>.pdb    fil_pdb         Input PDB file holding coords of
                                  protein and ligand
<filename>.hhb    fil_hbplus      Input files produced by HBPLUS (or
<filename>.nnb                    otherwise). The .hhb file gives all
				  the hydrogen bonds in the structure,
				  while the .nnb file gives all the
                                  non-bonded contacts. These files can
                                  either be in LIGPLOT format or in
                                  HBPLUS format.
ligplot.bonds     ligplot_bonds   Output file listing of bonds and bond
                  (fil_bonds)     types, for the molecules shown in the
                                  plot.                                 [*]
ligplot.frm       ligplot_frm     Output file in PDB format of the
                                  molecules shown in the plot, prior to
                                  flattening.
ligplot.hhb       ligplot_hhb     Output file containing just the
                                  hydrogen bonds from the original .hhb
                                  file actually used in the plot.       [*]
ligplot.nnb       ligplot_nnb     Output file containing just the
                                  non-bonded contacts from the original
                                  .nnb file actually used in the plot.  [*]
ligplot.prm       ligplot_par     Input parameter file governing
                                  appearance of plot produced.
ligplot.pdb       ligplot_pdb     Output file in PDB format of the
                                  flattened molecules shown in the
                                  plot.                                 [*]
ligplot.ps        ps_file         Output PostScript file of the plot.
ligplot.rcm       ligplot_rcm     Output file listing all the residue
                                  names and numbers in the plot, plus their
                                  centre of mass position on the plot. For
                                  use when generating a plot that is to
                                  have a similar residue layout to a plot
                                  for a related protein/complex. File is
                                  used if LIGPLOT detects <filename>.rcm
                                  file, where <filename> is the name of
                                  the PDB file being plotted.
ligplot.res       ligplot_res     Output file listing all the residue
                                  names and numbers in the plot. For use
                                  by the raslig program in PDBsum LIGPLOT
                                  pages. (File only generated with -r
                                  command-line option).

[*] The asterisked output files are taken as input when the program
    is run with ligplot.pdb as the input PDB file (in which case the
    Print_as_is flag is set to TRUE).

*/

/* Function calling-tree
   ---------------------

main
   -> initialise_parameters
   -> read_in_parameters
         -> update_current_parameter
   -> get_random_number
   -> define_object_sizes
   -> define_text_sizes
   -> get_command_arguments
         -> check_for_dimplot
         -> lig_search
               -> chk_cofactor
               -> first_res
               -> print_cofactor
               -> lig_size
               -> print_chains
               -> ligand_input
         -> get_residue_number
         -> get_filename
               -> getnam
   -> get_ligand_start_end
         -> check_if_in_ligand
   -> read_conec_records
   -> verify_conec_records
   -> read_hbplus_file
         -> get_hbplus_info
         -> check_for_metal
         -> check_if_in_ligand
   -> special_resinc
         -> get_hbplus_info
   -> read_pdb_file
         -> check_if_in_ligand
         -> is_hydrogen_atom
         -> update_bond
   -> define_objects
   -> update_all_boundaries
         -> update_residue_boundaries
   -> open_pdb_output
   -> open_bonds_output
   -> covalent_connectivity
         -> update_bond
   -> combine_objects
         -> check_for_join
         -> check_standard_aa
         -> join_residues
               -> get_last
               -> get_pointer_to_residue
               -> get_pointer_to_atom
   -> delete_false_ligands
         -> check_for_metal
   -> add_hbplus_bonds
         -> check_if_in_ligand
         -> update_bond
   -> bond_linkage
   -> check_connections
         -> get_test_bond
         -> remove_unconnected_atoms
               -> initialise_atoms
               -> mark_downstream_atoms
                    -> initialise_bonds
         -> delete_atoms_and_bonds
               -> delete_unwanted_atoms
               -> delete_unwanted_bonds
               -> delete_unwanted_bond_links
               -> delete_unwanted_bond_connections
   -> check_for_cyclics
   -> split_objects
         -> check_for_attachment
   -> assign_bonds_to_objects
         -> delete_object
   -> delete_unwanted_bond_connections
   -> mark_reachable_atoms
         -> initialise_atoms
   -> delete_unreachable_residues
         -> delete_object
   -> show_deleted_residues_message
   -> get_atom_sizes
   -> atom_linkage
   -> identify_rotatable_bonds
         -> ring_check
               -> initialise_bonds
   -> write_frm_file
   -> write_lig_prot_bonds
   -> pare_down_hydrophobics
         -> initialise_atoms
   -> delete_atoms_and_bonds
         -> delete_unwanted atoms
         -> delete_unwanted_bonds
         -> delete_unwanted_bond_links
         -> delete_unwanted_bond_connections
   -> calc_fit_coordinates
         -> centre_of_mass
         -> adjust_stored_coords
         -> principal_components
               -> calc_eigen_values
               -> eigen_sort
               -> orien3
                     -> vecprd
                     -> dot3
         -> apply_rotation
   -> flatten_all_objects
         -> create_simple_hgroup
               -> initialise_atoms
         -> delete_atoms_and_bonds
               -> delete_unwanted_atoms
               -> delete_unwanted_bonds
               -> delete_unwanted_bond_links
               -> delete_unwanted_bond_connections
         -> remove_mainchains
               -> initialise_atoms
         -> update_residue_boundaries
         -> place_group_at_origin
         -> flatten_all_rings
               -> get_other_ring_bonds
                     -> initialise_bonds
               -> flatten_ring
                     -> get_ring_atom_coords
                           -> initialise_atoms
                     -> centre_of_mass
                     -> move_object
                     -> adjust_stored_coords
                     -> principal_components
                           -> calc_eigen_values
                           -> eigen_sort
                           -> orien3
                                 -> vecprd
                                 -> dot3
                     -> rotate_object
                     -> squash_z
         -> correct_bond_lengths
               -> line_transformation
                     -> get_angle
                     -> calculate_matrix
                     -> apply_transformation
               -> transform_object
                     -> apply_transformation
               -> stretch_bond
                     -> initialise_bonds
                     -> initialise_atoms
         -> flatten_ring_offshoots
               -> offshoot_splay
                     -> plane_transformation
                           -> get_angle
                           -> calculate_matrix
                           -> apply_transformation
                     -> transform_object
                           -> apply_transformation
                     -> line_transformation
                           -> get_angle
                           -> calculate_matrix
                     -> apply_transformation
                     -> transform_downstream_atoms
                           -> initialise_bonds
                           -> initialise_atoms
                           -> apply_transformation
         -> unroll_structure
               -> initialise_atoms
               -> flatten_object
                     -> line_transformation
                           -> get_angle
                           -> calculate_matrix
                     -> transform_object
                     -> flatten_bonds
                           -> plane_transformation
                                 -> get_angle
                                 -> calculate_matrix
                                 -> apply_transformation
                           -> transform_downstream_atoms
                                 -> initialise_bonds
                                 -> initialise_atoms
                                 -> apply_transformation
                           -> line_transformation
                                 -> get_angle
                                 -> calculate_matrix
                           -> apply_transformation
                     -> move_object
                     -> flip_object
               -> move_object
         -> align_object
               -> centre_of_mass
               -> adjust_to_origin
               -> adjust_stored_coords
               -> principal_components
                     -> calc_eigen_values
                     -> eigen_sort
                     -> orien3
                           -> vecprd
                           -> dot3
               -> rotate_object
         -> force_flat
         -> create_simple_ligand
               -> stretch_mainchain
                     -> line_transformation
                     -> transform_object
                     -> transform_downlist_atoms
                           -> transform_downstream_atoms
                                 -> initialise_bonds
                                 -> initialise_atoms
                                 -> apply_transformation
                     -> move_object
                     -> apply_transformation
               -> mark_sidechain_atoms
         -> adjust_for_schematic
               -> initialise_atoms
         -> delete_atoms_and_bonds
               -> delete_unwanted atoms
	       -> delete_unwanted_bonds
	       -> delete_unwanted_bond_links
               -> delete_unwanted_bond_connections
         -> untangle_object
               -> get_random_number
               -> save_restore_coords
               -> calc_internal_energy
                     -> update_residue_boundaries
                     -> fit_object
                           -> matfit
                                 -> qikfit
                           -> get_object_distance_match
                           -> move_object
                           -> flip_object
                           -> rotate_object
               -> line_transformation
                     -> get_angle
                     -> calculate_matrix
               -> transform_object
               -> flip_about_bond
                     -> initialise_bonds
               -> flip_back
   -> score_objects
   -> sort_objects
   -> calc_residue_means
   -> read_rcm_file
         -> get_filename
   -> lay_interface_objects_out
         -> get_object_means
         -> centre_of_mass
         -> adjust_stored_coords
         -> principal_components
         -> rotate_centres_of_mass
         -> place_interacting_residues
   -> lay_objects_out
         -> fit_object
               -> matfit
                     -> qikfit
               -> get_object_distance_match
               -> move_object
               -> flip_object
               -> rotate_object
         -> update_residue_boundaries
   -> energy_minimize
         -> update_all_boundaries
         -> get_random_number
         -> get_total_energy
               -> calculate_energy
                     -> get_atom_clash_energy
                     -> line_transformation
                           -> get_angle
                           -> calculate_matrix
                     -> get_bond_clash_energy
                           -> apply_transformation
                           -> get_distance_energy
                     -> get_bond_overlap_energy
                           -> apply_transformation
                     -> get_anchor_energy
                     -> get_rel_pstn_energy
                     -> get_boundary_energy
         -> save_restore_all_objects
               -> save_restore_coords
         -> save_restore_coords
         -> calculate_energy
               -> get_atom_clash_energy
               -> line_transformation
                     -> get_angle
                     -> calculate_matrix
               -> get_bond_clash_energy
                     -> apply_transformation
                     -> get_distance_energy
               -> get_bond_overlap_energy
                     -> apply_transformation
               -> get_anchor_energy
               -> get_rel_pstn_energy
               -> get_boundary_energy
         -> update_residue_boundaries
         -> make_random_move
               -> get_random_number
               -> move_object
         -> make_random_rotation
               -> get_random_number
               -> move_object
               -> rotate_object
         -> make_random_flip
               -> get_random_number
               -> move_object
               -> line_transformation
                     -> get_angle
                     -> calculate_matrix
               -> transform_object
               -> flip_about_bond
               -> calc_internal_energy
                     -> update_residue_boundaries
                     -> fit_object
                           -> matfit
                                 -> qikfit
                           -> move_object
                           -> flip_object
                           -> rotate_object
                -> flip_attached_groups
                      -> move_object
                      -> transform_object
                      -> flip_object
         -> make_random_atom_flip
               -> get_random_number
               -> move_object
               -> flip_object
         -> get_total_energy
               -> calculate_energy
                     -> get_atom_clash_energy
                     -> line_transformation
                           -> get_angle
                           -> calculate_matrix
                     -> get_bond_clash_energy
                           -> apply_transformation
                           -> get_distance_energy
                     -> get_bond_overlap_energy
                           -> apply_transformation
                     -> get_anchor_energy
                     -> get_rel_pstn_energy
                     -> get_boundary_energy
   -> ligand_orientate
         -> centre_of_mass
         -> adjust_to_origin
         -> adjust_stored_coords
         -> principal_components
               -> calc_eigen_values
               -> eigen_sort
               -> orien3
                     -> vecprd
                     -> dot3
         -> rotate_whole_structure
         -> swap_x_and_y
         -> check_ligand_direction
               -> flip_about_x
               -> flip_about_y
         -> update_residue_boundaries
   -> calc_rotation_matrix
   -> rotate_whole_structure
   -> update_all_boundaries
   -> write_pdb_file
   -> write_bonds
   -> write_rcm_file
   -> psopen_
         -> match_colour_names
         -> define_object_colours
         -> define_text_colours
         -> define_object_sizes
         -> define_text_sizes
         -> pscomm_
         -> pshade_
         -> psubox_
   -> adjust_for_landscape
   -> psrot_
   -> print_title
         -> pscomm_
         -> pscolb_
         -> psctxt_
   -> print_symbol_key
         -> pscomm_
         -> pssave_
         -> pscolb_
         -> pslwid_
         -> print_key_text
               -> pstext_
         -> psline_
         -> psrest_
         -> psphcl_
         -> pspher_
         -> psdash_
         -> psctxt_
         -> plot_contact_group
               -> pssave_
               -> pscolb_
               -> pslwid_
               -> psarc_
               -> plot_spokes
                     -> psline_
               -> psrest_
         -> accessibility_circle
               -> pscrgb_
               -> psucir_
         -> psccol_
         -> psucir_
   -> extract_minmax
   -> scale
   -> centre_point
   -> define_object_sizes
   -> define_text_sizes
   -> plot_postscript_picture
         -> shade_accessibilities
               -> pscomm_
               -> initialise_atoms
               -> sort_accessibilities
               -> accessibility_cirle
                     -> pscrgb_
                     -> psucir_
               -> psccol_
               -> psucir_
         -> plot_hbond_lines
               -> pscomm_
               -> pssave_
               -> pscolb_
               -> pslwid_
               -> psdash_
               -> psline_
               -> print_hbond_lengths
                     -> pssave_
                     -> pscolb_
                     -> pslwid_
                     -> psubox_
                     -> psctxt_
                     -> psrest_
               -> psrest_
         -> plot_ext_bond_lines
               -> pscomm_
               -> pssave_
               -> pscolb_
               -> pslwid_
               -> psdash_
               -> psline_
               -> psrest_
         -> plot_bond_lines
               -> pscomm_
               -> pssave_
               -> pscolb_
               -> pslwid_
               -> double_bond_split
               -> get_atom_colour
               -> psrest_
         -> plot_atoms
               -> pscomm_
               -> pssave_
               -> initialise_atoms
               -> print_current_atom
                     -> pscolb_
                     -> get_atom_colour
                     -> psphcl_
                     -> pspher_
               -> print_atom_name
                     -> position_names
                     -> get_atom_text_colour
                     -> pscolb_
                     -> psctxt_
               -> psctxt_
               -> psrest_
         -> plot_hydrophobics
               -> pscomm_
               -> pssave_
               -> pslwid_
               -> print_name
                     -> convert_residue_name
                     -> pssave_
                     -> pslwid_
                     -> pscolb_
                     -> pshade_
                     -> psbbox_
                     -> psubox_
                     -> psctxt_
               -> plot_contact_group
                     -> pssave_
                     -> pscolb_
                     -> pslwid_
                     -> psarc_
                     -> plot_spokes
                           -> psline_
                     -> psrest_
               -> psrest_
         -> plot_simple_ligand_residues
               -> pscomm_
               -> pssave_
               -> pslwid_
               -> psccol_
               -> pscolb_
               -> pscirc_
               -> print_name
                     -> convert_residue_name
                     -> pssave_
                     -> pslwid_
                     -> pscolb_
                     -> pshade_
                     -> psbbox_
                     -> psubox_
                     -> psctxt_
               -> psrest_
         -> plot_simple_hgroups
               -> pscomm_
               -> pssave_
               -> pslwid_
               -> print_name
                     -> convert_residue_name
                     -> pssave_
                     -> pslwid_
                     -> pscolb_
                     -> pshade_
                     -> psbbox_
                     -> psubox_
                     -> psctxt_
               -> psrest_
         -> print_residue_names
               -> pssave_
               -> pscolb_
               -> convert_residue_name
                     -> find_ca_position
                     -> coords_resname
               -> resname_grid_search
                     -> initialise_atoms
                     -> extract_search_atoms
                           -> initialise_atoms
                     -> perform_search
                           -> label_dist
               -> psrest_
   -> psclos_

*/

#include "ligplot.h"


/***********************************************************************

initialise_parameters  -  Initialise the plot parameters

***********************************************************************/

void initialise_parameters(void)
{
  int icol, icolour, iline, isize;

  /* Initialise counts */
  n_labels = 0;

  /* Initialise global variables */
  Page_Min_x = BBOXX1;
  Page_Max_x = BBOXX2;
  Page_Min_y = BBOXY1;
  Page_Max_y = BBOXY2;
  Margin_y = (float)Page_Min_y;
  Split_Colour_Ligand_Bonds = FALSE;
  Split_Colour_Nonligand_Bonds = FALSE;
/* v.4.0--> */
  Have_Anchors = FALSE;
/* <--v.4.0 */
  
  /* Initialise pointers */
  first_bond_ptr = NULL;
  first_hhb_info_ptr = NULL;
  last_hhb_info_ptr = NULL;
  
  /* Initialise picture options */
  Picture = (struct Picture_Options *)
    malloc(sizeof(struct Picture_Options));
  Picture->In_Colour = FALSE;
  Picture->Portrait = TRUE;

  /* Title of plot */
  root_name[0] = '\0';
  Print_Title[0] = '\0';
  
  /* Special flags and values */
  xsite_file = FALSE;
  Max_atom_radius = 0.0;
  Maximum_accessibility = 100.0;
  Nobjects = 0;
  Plot_Centre_x = (Page_Min_x + Page_Max_x) / 2.0;
  Plot_Centre_y = (Page_Min_y + Page_Max_y) / 2.0;
  
  /* Parameters determining what's to be included in the plot */
  Include = (struct Include_Parameters *)
    malloc(sizeof(struct Include_Parameters));
  Include->Waters = FALSE;
  Include->Mainchain_Atoms = TRUE;
  Include->Hydrophobics = TRUE;
  Include->Hbonds = TRUE;
  Include->External_Bonds = TRUE;
  Include->Hydrophobic_Bonds = FALSE;
  Include->Internal_Hbonds = FALSE;
  Include->Simple_Ligand_Residues = FALSE;
  Include->Simple_Nonligand_Residues = FALSE;
  Include->Accessibilities = FALSE;
  Include->Ligand_Atoms = TRUE;
  Include->Nonligand_Atoms = TRUE;
/* v.4.0--> */
  Include->Water_Atoms = TRUE;
  Include->Ligand_Accessibilities_Only = TRUE;
/* <--v.4.0 */
  Include->Atom_Names = TRUE;
  Include->Residue_Names = TRUE;
  Include->Hbond_Lengths = TRUE;
  Include->Linked_Residues = FALSE;
  Include->Key = TRUE;
  Include->Double_Bonds = FALSE;
  Include->Filename_for_Title = TRUE;
/* v.3.2--> */
  Include->External_Bonds_Solid = FALSE;
  Include->Contact_Type = 0;
/* <--v.3.2 */
  Include_Hydrogens = FALSE;
  
  /* Residues H-bonded to non-ligand residues */
  for (iline = 0; iline < MAXSPECIAL_RES; iline++)
    strcpy(special_res[iline],"       ");
  
  /* Sizes of objects on plot */
  for (isize = 0; isize < OBJECT_SIZES; isize++)
    object_size[isize] = 0.0;
  Size = (struct Sizes *) malloc(sizeof(struct Sizes));
  Size_Val = (struct Size_Vals *) malloc(sizeof(struct Size_Vals));
  strcpy(Size->Ligand_Atoms,           "Ligatom_radius ");
  strcpy(Size->Nonligand_Atoms,        "Nligatom_radius");
  strcpy(Size->Waters,                 "Water_radius   ");
  strcpy(Size->Hydrophobics,           "Hphobic_radius ");
  strcpy(Size->Simple_Residues,        "Simple_radius  ");
  strcpy(Size->Ligand_Bonds,           "Ligbond_width  ");
  strcpy(Size->Nonligand_Bonds,        "Nligbond_width ");
  strcpy(Size->Hydrogen_Bonds,         "Hbond_width    ");
  strcpy(Size->External_Bonds,         "Ext_bond_width ");

  /* Text sizes on plot */
  for (isize = 0; isize < TEXT_SIZES; isize++)
    text_size[isize] = 0.0;
  Text_Size = (struct Text_Sizes *) malloc(sizeof(struct Text_Sizes));
  Text_Size_Val 
    = (struct Text_Size_Vals *) malloc(sizeof(struct Text_Size_Vals));
  strcpy(Text_Size->Ligand_Residue_Names,   "Ligresname_size");
  strcpy(Text_Size->Nonligand_Residue_Names,"Nligresnam_size");
  strcpy(Text_Size->Water_Names,            "Watername_size ");
  strcpy(Text_Size->Hydrophobic_Names,      "Hydrophnam_size");
  strcpy(Text_Size->Simple_Residue_Names,   "Simpletext_size");
  strcpy(Text_Size->Ligand_Atom_Names,      "Ligatmname_size");
  strcpy(Text_Size->Nonligand_Atom_Names,   "Nligatmnam_size");
  strcpy(Text_Size->Hbond_Lengths,          "HBtext_size    ");
  
  /* Object colours */
  for (icolour = 0; icolour < OBJECT_COLOURS; icolour++)
    strcpy(object_colour[icolour],"Grey00");
  strcpy(object_colour[0],"Grey10");
  strcpy(object_colour[6],"Grey10");
  strcpy(object_colour[8],"Grey08");
  strcpy(object_colour[9],"Grey04");
  strcpy(object_colour[11],"Grey09");
  strcpy(object_colour[12],"Grey02");
  strcpy(object_colour[13],"Grey04");
  strcpy(object_colour[14],"Grey04");
  strcpy(object_colour[15],"Grey04");
  Colour = (struct Colours *) malloc(sizeof(struct Colours));
  strcpy(Colour->Background,       "Background_col ");
  strcpy(Colour->Ligand_Bonds,     "Ligbond_colour ");
  strcpy(Colour->Nonligand_Bonds,  "Nligbond_colour");
  strcpy(Colour->Hydrogen_Bonds,   "Hbond_colour   ");
  strcpy(Colour->External_Bonds,   "Ext_bond_colour");
  strcpy(Colour->Hydrophobics,     "Hydrophobic_col");
  strcpy(Colour->Accessibility_Min,"Accessib_min   ");
  strcpy(Colour->Accessibility_Max,"Accessib_max   ");
  strcpy(Colour->Nitrogen,         "Nitrogen_colour");
  strcpy(Colour->Oxygen,           "Oxygen_colour  ");
  strcpy(Colour->Carbon,           "Carbon_colour  ");
  strcpy(Colour->Sulphur,          "Sulphur_colour ");
  strcpy(Colour->Water,            "Water_colour   ");
  strcpy(Colour->Phosphorus,       "Phosphorus_col ");
  strcpy(Colour->Iron,             "Iron_colour    ");
  strcpy(Colour->Other,            "Other_colour   ");
  strcpy(Colour->Atom_Edges,       "Atom_edge_col  ");
  strcpy(Colour->Simple_Residues,  "Simpleres_col  ");
  
  /* Accessibility maximum and minimum */
  for (icol = 0; icol < 3; icol++)
    {
      Accessibility_Max[icol] = 1.0;
      Accessibility_Min[icol] = 0.2;
    }
  
  /* Text colours */
  for (icolour = 0; icolour < TEXT_COLOURS; icolour++)
    strcpy(text_colour[icolour],"Grey00");
  Text_Colour = (struct Text_Colours *) malloc(sizeof(struct Text_Colours));
  strcpy(Text_Colour->Title,                  "Title_colour   ");
  strcpy(Text_Colour->Key_Text,               "Key_Text_colour");
  strcpy(Text_Colour->Ligand_Residue_Names,   "Ligresname_col ");
  strcpy(Text_Colour->Nonligand_Residue_Names,"Nligresname_col");
  strcpy(Text_Colour->Water_Names,            "Watername_col  ");
  strcpy(Text_Colour->Hydrophobic_Names,      "Hydrophname_col");
  strcpy(Text_Colour->Ligand_Atom_Names,      "Ligatmname_col ");
  strcpy(Text_Colour->Nonligand_Atom_Names,   "Nligatmname_col");
  strcpy(Text_Colour->Hbond_Lengths,          "HBlength_col   ");
  
  /* Colour definitions */
  for (icolour = 0; icolour < MAX_COLOURS; icolour++)
    {
      strcpy(Colour_Table_Name[icolour],"BLACK          ");
      Colour_Table[icolour][0] = 0.0000;
      Colour_Table[icolour][1] = 0.0000;
      Colour_Table[icolour][2] = 0.0000;
    }
  
  /* Minimization parameters */
  Atom_Atom_Clash = 10.0;
  Bond_Atom_Clash = 0.20;
  Overlap_Score = 50.0;
  HB_Weight = 0.5;
  Nonbond_Weight = 0.2;
  Internal_Energy_Weight = 10.0;
  Max_HBdist = 15.0;
  Bond_Stretch = 1.0;
  Energy_Drop_Maximum = 0.0;
  Max_Loops = 1000;
  Random_Start = FALSE;
/* v.4.0--> */
  Anchor_Weight = 2.0;
  Boundary_Energy_Weight = 1.0;
  Min_Boundary_Dist = 2.0;
  Rel_Pstn_Energy_Weight = -1.0;
/* <--v.4.0 */
}
/***********************************************************************

update_current_parameter  -  Update the current plot parameter with
                             the value just read in

***********************************************************************/

void update_current_parameter(int parameter_group, int parameter,
			      int true_false, float size,
			      char colour[COL_NAME_LEN + 1],
			      float red, float green, float blue,
			      char residue_pair[8],
			      char special_res[MAXSPECIAL_RES][15])
{
  /* Process according to the parameter-group */
  switch (parameter_group)
    {
      /* Colour option */
    case 0 :
      
      /* Update the appropriate parameter */
      switch (parameter)
	{
	case 0 : Picture->In_Colour = true_false; break;
	case 1 : Picture->Portrait = true_false; break;
/* v.3.2--> */
	case 2 : Picture->Rotation_Angle = size; break;
/* <--v.3.2 */
	default :
	  printf("*** Program error. Invalid Picture parameter: %d\n",
		 parameter);
	  exit(1);
	}
      break;

      /* Include options */
    case 1 :

      /* Update the appropriate parameter */
      switch (parameter)
	{
	case  0 : Include->Hydrophobics = true_false; break;
	case  1 : Include->Waters = true_false; break;
	case  2 : Include->Mainchain_Atoms = true_false; break;
	case  3 : Include->Linked_Residues = true_false; break;
	case  4 : Include->Hbonds = true_false; break;
	case  5 : Include->Internal_Hbonds = true_false; break;
	case  6 : Include->External_Bonds = true_false; break;
	case  7 : Include->Hydrophobic_Bonds = true_false; break;
	case  8 : Include->Simple_Ligand_Residues = true_false; break;
	case  9 : Include->Simple_Nonligand_Residues = true_false; break;
	case 10 : Include->Accessibilities = true_false; break;
	case 11 : Include->Ligand_Atoms = true_false; break;
	case 12 : Include->Nonligand_Atoms = true_false; break;
	case 13 : Include->Double_Bonds = true_false; break;
	case 14 : Include->Key = true_false; break;
	case 15 : Include->Residue_Names = true_false; break;
	case 16 : Include->Atom_Names = true_false; break;
	case 17 : Include->Hbond_Lengths = true_false; break;
	case 18 : Include->Filename_for_Title = true_false; break;
/* v.3.2--> */
	case 19 : Include->External_Bonds_Solid = true_false; break;
	case 20 : Include->Contact_Type = size; break;
/* <--v.3.2 */
/* v.4.0--> */
	case 21 : Include->Water_Atoms = true_false; break;
	case 22 : Include->Ligand_Accessibilities_Only = true_false; break;
/* <--v.4.0 */
	default :
	  printf("*** Program error. Invalid Include parameter: %d\n",
		 parameter);
	  exit(1);
	}
      break;

      /* Linked residues */
    case 2 :
      
      /* Update the appropriate parameter */
      strncpy(special_res[parameter],residue_pair,8);
      special_res[parameter][8] = '\0';
      break;

      /* Sizes */
    case 3 :
      
      /* Update the appropriate parameter */
      object_size[parameter] = size;
      break;

      /* Text sizes */
    case 4 :
      
      /* Update the appropriate parameter */
      text_size[parameter] = size;
      break;

      /* Object colours */
    case 5 :

      /* Update the appropriate object colour */
      if (Picture->In_Colour == TRUE)
	strcpy(object_colour[parameter],colour);

      /* If black-and-white, update only if this is an ATOM option
	 on ligand or non-ligand bonds */
      else if ((parameter == 1 || parameter == 2) &&
	       !strncmp(colour,"ATOM           ",15))
	strcpy(object_colour[parameter],colour);
      break;

      /* Text colours */
    case 6 :
      
      /* Update the appropriate text colour */
      if (Picture->In_Colour == TRUE)
	strcpy(text_colour[parameter],colour);
      break;

      /* Colour definitions */
    case 7 :

      /* Update the appropriate parameter */
      Colour_Table[parameter][0] = red;
      Colour_Table[parameter][1] = green;
      Colour_Table[parameter][2] = blue;
      strcpy(Colour_Table_Name[parameter],colour);
      break;

      /* Minimization parameters */
    case 8 :

      /* Update the appropriate parameter */
      switch (parameter)
	{
	case  0 : Atom_Atom_Clash = size; break;
	case  1 : Bond_Atom_Clash = size; break;
	case  2 : Overlap_Score = size; break;
	case  3 : HB_Weight = size; break;
	case  4 : Nonbond_Weight = size; break;
	case  5 : Internal_Energy_Weight = size; break;
	case  6 : Max_HBdist = size; break;
	case  7 : Bond_Stretch = size; break;
	case  8 : Max_Loops = size; break;
	case  9 : Energy_Drop_Maximum = size; break;
	case 10 : Random_Start = true_false; break;
/* v.4.0--> */
	case 11 : Anchor_Weight = size; break;
	case 12 : Boundary_Energy_Weight = size; break;
	case 13 : Min_Boundary_Dist = size; break;
	case 14 : Rel_Pstn_Energy_Weight = size; break;
/* <--v.4.0 */
	default :
	  printf("*** Program error. \n");
	  printf("*** Invalid Minimization parameter: %d\n",
		 parameter);
	  exit(1);
	}
      break;

      /* Error message */
    default :
      printf("*** Program error. Invalid parameter group: %d\n",
	     parameter_group);
      exit(1);
    }
}
/***********************************************************************

read_in_parameters  -  Read in the plot parameters from the
                       ligplot.par parameter file

***********************************************************************/

void read_in_parameters(char special_res[MAXSPECIAL_RES][15])
{
  char line[LINELEN + 1];
  char colour_name[COL_NAME_LEN + 1], residue_pair[8];
  char show_string[14];
/* v.3.2--> */
/*  char ligplot_version[] = "LIGPLOT v.3.0"; */
/* v.4.0--> */
/*  char ligplot_version[] = "LIGPLOT v.3.2"; */
/* v.4.0.2--> */
/*  char ligplot_version[] = "LIGPLOT v.4.0"; */
  char ligplot_version[] = "LIGPLOT v.4.0.2";
/* <--v.4.0.2 */
/* <--v.4.0 */
  char old_version[16];
/* <--v.3.2 */
  char parameter_type, parameter_types[] = {'P', 'Y', 'L', 'S', 'S',
					      'C', 'C', 'D', 'S'};
  int iline, last_group, parameter, parameter_group, parameter_lines;
  int number_of_parameters[] = { PICTURE_OPTIONS, INCLUDE_PARAMETERS,
				   MAXSPECIAL_RES, OBJECT_SIZES,
				   TEXT_SIZES, OBJECT_COLOURS,
				   TEXT_COLOURS, MAX_COLOURS,
				   MINIMIZATION_PARAMETERS };
  int true_false;
  float blue, green, red, size;

/* v.4.0--> */
  /* Print version number */
  printf("\n");
  printf("Running %s ...\n",ligplot_version);
  printf("\n");
/* <--v.4.0 */

  /* Open the parameter file */
  if ((ligplot_par = fopen("ligplot.prm","r")) == NULL)
    {
      printf("\n*** Unable to open parameter file: ligplot.prm\n\n");
      exit(1);
    }

  /* Initialise flags and counts */
  strcpy(colour_name,"WHITE");
  iline = 0;
  last_group = -1;
  parameter = 0;
  parameter_group = -1;
  parameter_lines = -1;
  size = 0.0;
  true_false = FALSE;

  /* Loop while reading records from the parameter file */
  while (fgets(line,LINELEN,ligplot_par)!= NULL)
    {
      /* If this is the first line, then check that the parameter
	 file corresponds to the current program version number */
      if (iline == 0)
	{
	  if (strncmp(line,ligplot_version,13))
	    {
/* v.3.2--> */
	      /* Check for a version which may work OK */
/* v.4.0--> */
/*	      if (!strncmp(line,ligplot_version,11)) */
	      if (!strncmp(line,"LIGPLOT v.3",11))
/* <--v.4.0 */
		{
		  printf("*** Warning. Your ligplot.prm file is for ");
		  strncpy(old_version,line,15);
		  old_version[15] = '\0';
		  printf("%s\n",old_version);
		  printf("***          You will get better results if you");
		  printf(" delete the file and re-run the\n");
/* v.4.0--> */
/*		  printf("***          program which is %s\n",
			 ligplot_version); */
		  printf("***          program to generate a new ");
		  printf("parameter file for %s\n",ligplot_version);
		  printf("\n");
/* <--v.4.0 */
		  Nwarnings++;
		}
/* <--v.3.2 */
	      /* If this is a ligplot file for any other version, then
		 unlikely to work OK */
/* v.3.2--> */
/* 	      if (!strncmp(line,"LIGPLOT",7)) */
	      else if (!strncmp(line,"LIGPLOT",7))
/* <--v.3.2 */
		{
		  printf("*** ERROR. Incorrect version of parameter file");
		  printf(", ligplot.prm\n");
		  strncpy(show_string,line,13);
		  show_string[13] = '\0';
		  printf("           File is for %s, whereas program ",
			 show_string);
		  printf("requires %s\n",ligplot_version);
/* v.3.2--> */
		  exit(1);
/* <--v.3.2 */
		}

	      /* Otherwise, if not a LIGPLOT parameter file at all, then
		 forget it */
	      else
		{
		  printf("*** ERROR. Unrecognised first line in parameter");
		  printf(" file, ligplot.prm\n");
		  printf("%s\n",line);
/* v.3.2--> */
		  exit(1);
/* <--v.3.2 */
		}
/* v.3.2--> */
/*	      exit(1); */
/* <--v.3.2 */
	    }
	}

      /* Check for one of the keyword lines */
      if (!strncmp(line,"PRINT OPTIONS",13))
	parameter_group = 0;
      else if (!strncmp(line,"PLOT PARAMETERS",15))
	parameter_group = 1;
      else if (!strncmp(line,"LINKED RESIDUES",15))
	parameter_group = 2;
      else if (!strncmp(line,"SIZES",5))
	parameter_group = 3;
      else if (!strncmp(line,"TEXT SIZES",10))
	parameter_group = 4;
      else if (!strncmp(line,"COLOURS",7))
	parameter_group = 5;
      else if (!strncmp(line,"TEXT COLOURS",12))
	parameter_group = 6;
      else if (!strncmp(line,"COLOUR DEFINITIONS",18))
	parameter_group = 7;
      else if (!strncmp(line,"MINIMIZATION PARAMETERS",23))
	parameter_group = 8;

      /* If this is a new parameter-group, then initialise
	 line-counts */
      if (parameter_group != last_group)
        {
	  parameter = -2;
	  last_group = parameter_group;
	  parameter_lines = number_of_parameters[parameter_group];
	  parameter_type = parameter_types[parameter_group];
        }

      /* If this line is one of the current group of parameters,
	 then extract the values and store */
      if (parameter > -1 && parameter < parameter_lines)
        {
	  /* Check for a Yes/No response */
	  if (line[0] == 'Y' || line[0] == 'y')
	    true_false = TRUE;
	  else
	    true_false = FALSE;

	  /* If this is a yes/no response, then get its value */
	  if (parameter_type == 'Y')
            {
	      if (line[0] == 'Y' || line[0] == 'y')
		true_false = TRUE;
	      else
		true_false = FALSE;
/* v.3.2--> */
	      if (line[0] != 'Y' && line[0] != 'y' &&
		  line[0] != 'N' && line[0] != 'n')
		sscanf(line,"%f ",&size);
/* <--v.3.2 */
            }

	  /* If this is one of the print options, interpret it */
	  else if (parameter_type == 'P')
            {
	      if (parameter == 0 && (line[0] == 'Y' || line[0] == 'y'))
		true_false = TRUE;
	      else if (parameter == 1 && (line[0] == 'P' || line[0] == 'p'))
		true_false = TRUE;
	      else
		true_false = FALSE;
/* v.3.2--> */
	      if (parameter == 2)
		sscanf(line,"%f ",&size);
/* <--v.3.2 */
            }

	  /* If this is a colour, then store its name */
	  else if (parameter_type == 'C')
            {
	      strncpy(colour_name,line,COL_NAME_LEN);
	      colour_name[COL_NAME_LEN] = '\0';
            }

	  /* If this is a size, then retrieve the value */
	  else if (parameter_type == 'S')
            {
	      /* Read in the value of the size */
	      sscanf(line,"%f ",&size);
            }

	  /* If this is a link-residue pair, retrieve the residue names */
	  else if (parameter_type == 'L')
            {
	      strncpy(residue_pair,line,7);
	      residue_pair[7] = '\0';
            }

	  /* If this is a colour definition line, read in the RGB
	     values and the colour name */
	  else if (parameter_type == 'D')
            {
	      /* Read in the RGB values */
	      sscanf(line,"%f %f %f ",&red,&green,&blue);
	      
	      /* Read in the colour name */
	      strncpy(colour_name,line+22,COL_NAME_LEN);
	      colour_name[COL_NAME_LEN] = '\0';
            }

	  /* Update the current parameter */
	  update_current_parameter(parameter_group,parameter,true_false,
				   size,colour_name,red,green,blue,
				   residue_pair,special_res);
        }
      /* Increment the parameter line-count */
      parameter++;
      iline++;
    }

  /* Correct any incompatible parameters */
/* v.4.0--> */
/*  if (Include->Simple_Ligand_Residues == TRUE ||
      Include->Simple_Nonligand_Residues == TRUE)
    {
      Include->Accessibilities = FALSE;
      Include->External_Bonds = FALSE;
    } */
  if (Include->Simple_Ligand_Residues == TRUE)
    {
      Include->Accessibilities = FALSE;
      Include->External_Bonds = FALSE;
    }
  if (Include->Simple_Nonligand_Residues == TRUE)
    {
      Include->External_Bonds = FALSE;
      Include->Ligand_Accessibilities_Only = TRUE;
    }
/* <--v.4.0 */

  /* For interface plots, need to ensure that certain options aren't set */
  if (Interface_Plot == TRUE)
    {
      Include->Simple_Ligand_Residues = FALSE;
      Include->Simple_Nonligand_Residues = FALSE;

      /* Initialise weight for boundary energy calcs */
      if (Rel_Pstn_Energy_Weight < 0.0)
	Rel_Pstn_Energy_Weight = 0.3;
    }
}
/***********************************************************************

getnam  -  Peel off the directory path and extension from the full name
           of the .pdb file

***********************************************************************/

int getnam(char pdbfil[FILENAME_LEN],int *istart,int *iend)
{
    char pchar;
    int  finish, gotdot, ipos, istate;


    /* Initialise variables */
    finish = FALSE;
    *iend = -1;
    *istart = 0;
    istate = 1;
    ipos = strlen(pdbfil) - 1;
    gotdot = FALSE;

    /* Check through the filename from right to left */
    while (finish == FALSE && ipos > -1)
    {

        /* Pick off next character */
        pchar = pdbfil[ipos];


        /* State 1: Searching for first non-blank character */
        if (istate == 1)
        {
            if (pchar == '/' || pchar == '\\' || pchar == ']')
            {
                printf("*** ERROR in supplied name of file: [%s\n",
                       pdbfil);
                return(-1);
            }

            if (pchar != ' ' && pchar != '.')
            {
                *iend = ipos;
                istate = 2;
            }
        }

        /* State 2: Searching for end of extension, or end of
           directory path */
        else if (istate == 2)
        {
            /* If character is a dot, and is the first dot, then
               note position */
            if (pchar == '.' && gotdot == FALSE)
            {
                *iend = ipos - 1;
                gotdot = TRUE;
            }

            /* If character signifies the end of a directory path,
               note pstn */
            else if (pchar == '/' || pchar == '\\' || pchar == ']')
            {
                *istart = ipos + 1;
                finish = TRUE;
            }
        }

        /* Step back a character */
        ipos = ipos - 1;
    }

    /* Check whether file name is sensible */
    if (*istart > *iend)
    {
        printf("*** ERROR in supplied name of file: %s\n", pdbfil);
        printf("    No name found\n");
        return(-1);
    }
    return(0);
}
/***********************************************************************

get_filename  -  Strip off directory path of PDB filename
                 and append required extention to root

***********************************************************************/

void get_filename(char pdb_name[FILENAME_LEN], char file_name[FILENAME_LEN],
                  char extension[10])
{
    int      iend, ierror, inew, ipos, istart;

    /* Initialise returned filename */
    file_name[0] = '\0';

    /* Check that PDB filename is valid */
    if (pdb_name[0] != '\0' && pdb_name[0] != '\n')
    {
        /* Peel off the directory path and extension, returning just
           the root file name */
        ierror = getnam(pdb_name,&istart,&iend);
        if (ierror != -1)
        {
            /* Transfer just the root name into new file name */
            for (ipos = istart, inew = 0; ipos < iend + 1;
                 ipos++, inew++)
                file_name[inew] = pdb_name[ipos];
            file_name[inew] = '\0';

            /* Store the root name */
            strcpy(root_name,file_name);

            /* Append the extension */
            strcat(file_name,extension);
        }
    }
}
/***********************************************************************

get_residue_number  -  Interpret the residue-number entered by the
                       user

***********************************************************************/

/* v.4.0--> */
/* void get_residue_number(char tmpstr[LINELEN],char res_num[6]) */
int get_residue_number(char tmpstr[LINELEN],char res_num[6],char res_name[4])
/* <--v.4.0 */
{
  char insertion_code, res_number[6], new_number[6];
/* v.4.0--> */
  char temp_string[LINELEN];
/* <--v.4.0 */
  int cpos, got_number, idigit, ipos, ndigits;
/* v.4.0--> */
  int at_end, residue_name, residue_number;
  int len, nchars;
/* <--v.4.0 */

  /* Initialise position of insertion code */
  cpos = 0;
  insertion_code = ' ';
  got_number = FALSE;
  ndigits = 0;
  for (idigit = 0; idigit < 5; idigit++)
    new_number[idigit] = ' ';
/* v.4.0--> */
  at_end = FALSE;
  nchars = 0;
  residue_name = FALSE;
  residue_number = FALSE;

  /* Check whether this has a -n prefix indicating a residue name */
  if (!strncmp(tmpstr,"-n",2))
    {
      residue_name = TRUE;

      /* Excise the -n */
      strcpy(temp_string,tmpstr+2);
      strcpy(tmpstr,temp_string);
    }
/* <--v.4.0 */

  /* Loop through the characters in the entered residue-number
     searching for end of number or insertion-code */
/* v.4.0--> */
/*  for (ipos = 0; ipos < 5 && cpos == 0; ipos++) */
  for (ipos = 0; ipos < 5 && at_end == FALSE; ipos++)
/* <--v.4.0 */
    {
      /* If this is the end of the string, store end-position */
      if (tmpstr[ipos] == '\0')
/* v.4.0--> */
	{
/* <--v.4.0 */
	  cpos = ipos;
/* v.4.0--> */
	  at_end = TRUE;
	}
/* <--v.4.0 */

      /* If this is a digit, then have part of a residue number */
/* v.4.0--> */
/*      else if (tmpstr[ipos] >= '0' && tmpstr[ipos] <= '9') */
      else if ((tmpstr[ipos] >= '0' && tmpstr[ipos] <= '9') ||
	       tmpstr[ipos] == '-')
/* <--v.4.0 */
	{
	  res_number[ndigits] = tmpstr[ipos];
	  ndigits++;
	  got_number = TRUE;

/* v.4.0--> */
	  /* If already have a character, then this must be a residue name */
	  if (nchars > 0)
	    {
	      got_number = FALSE;
	      residue_name = TRUE;
	    }
/* <--v.4.0 */
	}

      /* If this is a character it might be a prefix to the
	 residue number, or the insertion code */
      else if ((tmpstr[ipos] >= 'A' && tmpstr[ipos] <= 'Z') ||
	       tmpstr[ipos] == '\'')
	{
	  /* If we already have a number, this must be an
	     insertion code */
	  if (got_number == TRUE)
	    {
	      cpos = ipos;
	      insertion_code = tmpstr[ipos];
	    }

	  /* Otherwise it must be a residue name */
	  else
	    {
/* v.4.0--> */
/*	      res_number[ndigits] = tmpstr[ipos];
	      ndigits++; */
	      residue_name = TRUE;
	      got_number = FALSE;
/* <--v.4.0 */
	    }
/* v.4.0--> */
	  /* Increment count of characters */
	  nchars++;
/* <--v.4.0 */
	}

      /* If this is a space, have reached the end of the number */
      else if (tmpstr[ipos] == ' ')
        {
	  if (got_number == TRUE)
	    cpos = ipos;
        }

      /* Otherwise, have an invalid character in the residue
	 number */
      else
	cpos = -1;
    }

  /* If have too many digits, then last one is an insertion code */
  if (ndigits > 4)
    {
      insertion_code = res_number[4];
      ndigits--;
    }

  /* If entry doesn't resemble a residue number, then see whether it
     might be a residue name instead */
/* v.4.0--> */
/*  if (cpos < 0 || got_number == FALSE) */
  if (residue_name == TRUE || cpos < 0 || got_number == FALSE ||
      nchars > 1)
/* <--v.4.0 */
    {
/* v.4.0--> */
/*        printf("*** Invalid residue-number: [%s]\n",tmpstr);
        exit(1);
      } */
      /* Get length of residue name */
      len = strlen(tmpstr);

      /* If too long, then have an error */
      if (len > 3)
	{
	  printf("*** Invalid residue-number or name:");
	  printf(" [%s]\n",tmpstr);
	  exit(1);
	}

      /* Store the name */
      strcpy(res_name,tmpstr);

      residue_number = FALSE;
    }

  /* Otherwise, reformat the number */
  else
    {
/* <--v.4.0 */
      if (cpos == 0)
	cpos = 4;

      /* Now re-form the number such the first four characters
	 contain the digits and the 5th contains the insertion
	 code (if any) */
      for (idigit = 0; idigit < ndigits; idigit++)
	{
	  ipos = idigit + (4 - ndigits);
	  new_number[ipos] = res_number[idigit];
	}

      /* Add the insertion code */
      new_number[4] = insertion_code;
      new_number[5] = '\0';

      /* Transfer the new number across */
      strncpy(res_num,new_number,5);
      res_num[5] = '\0';
/* v.4.0--> */

      /* Set flag to say that residue number is OK */
      residue_number = TRUE;
    }

  /* Return flag indicating whether residue number identified */
  return(residue_number);
/* <--v.4.0 */

}
/***********************************************************************

chk_cofactor  -  Check for a cofactor

***********************************************************************/

int chk_cofactor(int a,char data[MAXCYCLE][4],char pdb_line[81])
{
    int 
	i,
	cofactor;
/* v.3.2--> */
/*    cofactor = OFF; */
    cofactor = -1;
/* <--v.3.2 */
    for (i = 0;i < a;i++) 
	if (!strncmp(data[i],pdb_line+17,3))
	    cofactor = i;
    return cofactor;
}
/***********************************************************************

print_cofactor  -  Print out the cofactor information

***********************************************************************/

void print_cofactor(char pdb_line[81],char name[4],char number[6],
		    int cofactor,char infoc[MAXCYCLE][41])
{
    infoc[cofactor][40] = '\0';
    printf("\t%s %s %c\t%s\n",name,number,
	   pdb_line[21],infoc[cofactor]);
}
/***********************************************************************

first_res  -  Save the first residues details

***********************************************************************/

void first_res(char pdb_line[81],char num_start[6],char *chain_start,
	       char name_start[3])
{
    strncpy(num_start,pdb_line+22,5);
    num_start[5] = '\0';
    *chain_start = pdb_line[21];
    strncpy(name_start,pdb_line+17,3);
    name_start[3] = '\0';
}
/***********************************************************************

print_chains  -  Print out the information

***********************************************************************/

void print_chains(char num_start[6],
	     char chain_start,  
	     char name_start[3],
	     char res_num[6],
	     char chain_id,
	     char res_name[4])
{
    num_start[5] = '\0';
    name_start[3] = '\0';
    res_num[5] = '\0';
    res_name[3] = '\0';

    if ((atoi(res_num) - atoi(num_start)) > 10)
    {
	printf("\t%s %s %c to",
	       name_start,num_start,chain_start);
	printf(" %s %s %c \n",res_name,
	       res_num,chain_id);
    }
}
/***********************************************************************

ligand_input  -  For the ligand input string 

***********************************************************************/

void ligand_input(char character,
	      char ligand2[6],
	      char ligand3[6],
	      char ligand4[6],
	      char ligand5[6])
{
    static int
	num = 2,
	i = 0, 
	copied = OFF;

    if (character != ' ')
    {
	switch (num)
	{
	case(2):
	{
	    ligand2[i] = character;
	    i++;
	    ligand2[i] = '\0';
	    copied = ON;
	}
	    break;
	case(3):
	{
	    ligand3[i] = character;
	    i++;
	    ligand3[i] = '\0';
	}
	    break;
	case(4):
	{
	    ligand4[i] = character;
	    i++;
	    ligand4[i] = '\0';
	}
	    break;
	case(5):
	{
	    ligand5[i] = character;
	    i++;
	    ligand5[i] = '\0';
	}
	    break;
	}
    }
    else if (character == ' ' &&
	    (num != 2 || copied == ON))
    {
	i = 0;
	num++;
    }
    if (!strncmp(ligand4,"-h",2))
    {
	strncpy(ligand5,"-h",2);
	ligand5[2] = '\0';
	ligand4[0] = '\0';
    }
    if (!strncmp(ligand4,"-w",2))
    {
	strncpy(ligand5,"-w",2);
	ligand5[2] = '\0';
	ligand4[0] = '\0';
    }
}

/***********************************************************************

lig_size  -  See how big a ligand chain is

***********************************************************************/

int lig_size(FILE *fil_pdb,char num[5],char chain)
{
    char
	num_temp[6],
	pdb_line[81];
    int
	count;
    static int
	once = 1;
    if (once == 1)
    {
	once++;
	count = 0;
	strncpy(num_temp,"#####",5);
	while (fgets (pdb_line,80,fil_pdb)!= NULL)    
	{
	    if (strncmp(num_temp,pdb_line+22,5))
	    {
		strncpy(num_temp,pdb_line+22,5);
		if (pdb_line[21] == chain && 
		   strncmp(num,pdb_line+22,5))
		    count++;
		else
		    count = 0;
	    }
	    if (count == 12)
		return OFF;
	}
	return ON;
    }
    return OFF;
}
/***********************************************************************
  
lig_search  -  Search for the ligand information in PDB file
  
***********************************************************************/

/* v.3.2--> */
/* void lig_search(char pdb_name[FILENAME_LEN],char ligand2[6],char ligand3[6],
		char ligand4[6],char ligand5[6]) */
void lig_search(char pdb_name[FILENAME_LEN],
		char token[MAXTOKENS][FILENAME_LEN],int *ntoken)
/* <--v.3.2 */
{
  FILE
    *fil_pdb;
  char
    chain,
    num[6],
    nme[4],
    name[4],
    number[6],
    num_temp0[6],
    num_temp1[6],
    pdb_line[81],
    res_old[6],
    name_old[4],
    chain_old,
    num_start[6],
    chain_start,
    name_start[4],
    params[81],
    infoc[MAXCYCLE][41],
    data[MAXCYCLE][4];
/* v.3.2--> */
  char ligand2[6], ligand3[6], ligand4[6], ligand5[6];
/* <--v.3.2 */
  int
    print_lig,
    position,
    stop,
    printed_cofactor,
    a,i,
    cofactor,
    save_first_res;

  /* Initialise */
  a = 0;
  chain = '\0';
  save_first_res = ON;
  strncpy(name,"###",3);
  strncpy(number,"#####",5);
  strncpy(nme,"###",3);
  strncpy(num,"#####",5);
  printed_cofactor = OFF;
  stop = OFF;
/* v.3.2--> */
  cofactor = -1;
/* <--v.3.2 */

  printf("\nThe ligands or peptide chains that occur in your protein,\n");
  printf("%s, are as follows:\n\n",pdb_name);
    
  /* Open the input PDB file */
  fil_pdb = fopen(pdb_name,"r");
  if (fil_pdb == NULL)
    {
      printf("WARNING - Unable"); 
      printf(" to open %s file - EXITING PROGRAM \n",pdb_name);
      exit(2);
    }

  /* Loop while reading in the PDB file */
  while (fgets(pdb_line,80,fil_pdb) != NULL)
    {
      if (!strncmp("ENDMDL",pdb_line,6))
	stop = ON;
      if (stop == OFF)
	{
	  if (!strncmp(pdb_line,"HET   ",6))
	    {
	      strncpy(infoc[a],pdb_line+30,50);
	      infoc[a][50] = '\0';
	      strncpy(data[a],pdb_line+7,3);
	      data[a][3]  = '\0';
	      a++;
	      if (a > MAXCYCLE)
		{
		  printf("Too many ligands in protein -- EXITING PROGRAM\n");
		  exit(2);
		}
	    }

	  /* Get the protein details */
	  if ((!strncmp(pdb_line,"ATOM",4) ||
	       !strncmp(pdb_line,"HETATM",6) ||
	       !strncmp(pdb_line,"TER",3)) &&
	      strncmp("HOH",pdb_line+17,3))
	    {
	      cofactor = chk_cofactor(a,data,pdb_line);
	      if (save_first_res == ON && 
/* v.3.2--> */
/*		  cofactor == OFF) */
		  cofactor == -1)
/* <--v.3.2 */
		{
		  strncpy(res_old,pdb_line+22,5);
		  res_old[5] = '\0';
		  chain_old = pdb_line[21];
		  strncpy(name_old,pdb_line+17,3);

		  /* Get the first main chain residue details */
		  first_res(pdb_line,num_start,&chain_start,
			    name_start);
		  save_first_res = OFF;
		}

	      /* Check to see whether it's a cofactor*/
/* v.3.2--> */
/*	      if (cofactor != OFF) */
	      if (cofactor > -1)
/* <--v.3.2 */
		{
		  if (strncmp(name,pdb_line+17,3) ||
		      strncmp(number,pdb_line+22,5) ||
		      chain != pdb_line[21])
		    {
		      strncpy(name,pdb_line+17,3);
		      name[3] = '\0';
		      strncpy(number,pdb_line+22,5);
		      number[5] = '\0';
		      chain = pdb_line[21];
		      print_cofactor(pdb_line,
				     name,
				     number,
				     cofactor,
				     infoc);
		      printed_cofactor = ON;
		    }
		}
		if (printed_cofactor == ON && 
		   !strncmp(pdb_line,"ATOM ",5) &&
		   (strncmp(pdb_line+17,nme,3) ||
		    strncmp(pdb_line+22,num,5)) &&
		    pdb_line[21] == chain)
		{
		    strncpy(nme,pdb_line+17,3);
		    strncpy(num,pdb_line+22,5);
		    nme[3] = '\0';
		    num[5] = '\0';
		    position =  ftell(fil_pdb);
		    print_lig = lig_size(fil_pdb,num,chain);
		    fseek(fil_pdb, position, SEEK_SET);
		    if (print_lig == ON)
			printf("\t%s %s %c\n",nme,num,pdb_line[21]);
		}
		strncpy(num_temp0,pdb_line+22,5);
		num_temp0[5] = '\0';
		if (((atoi(num_temp1) - atoi(num_temp0)) > 10 ||
		    !strncmp(pdb_line,"TER",3)) &&
/* v.3.2--> */
/*		    cofactor == OFF ) */
		    cofactor == -1)
/* <--v.3.2 */
		{
		    printed_cofactor = OFF;
		    strncpy(res_old,pdb_line+22,5);
		    res_old[5] = '\0';
		    chain_old = pdb_line[21];
		    strncpy(name_old,pdb_line+17,3);
		    save_first_res = ON;
		    print_chains(num_start,
				 chain_start,
				 name_start,
				 res_old,
				 chain_old,
				 name_old);
		}
		strncpy(num_temp1,pdb_line+22,5);
		num_temp1[5] = '\0';
	    }
	}
    }
    fclose(fil_pdb);
	
    printf("\n\nNow enter one of these ligand ranges.\n");
    printf("For example :'1  4  I', leaving a <space> after \n");
    printf("each parameter.  If want a heading for your plot\n");
    printf("enter a '-h' as your last parameter\n");
    gets(params);
    i = 0;
    ligand2[0] = '\0';
    ligand3[0] = '\0';
    ligand4[0] = '\0';
    ligand5[0] = '\0';
    while (params[i])
    {
	ligand_input(params[i],
		     ligand2,
		     ligand3,
		     ligand4,
		     ligand5);
	i++;
    }
    if (!ligand3[0])
 	strcpy(ligand3,ligand2);
    if (!strncmp(ligand3,"-h",2))
    {
	strncpy(ligand4,ligand3,2);
	ligand4[2] = '\0';
	strcpy(ligand3,ligand2);
    }
    if (isalpha(ligand3[0]))
    {
	ligand4[0] = ligand3[0];
 	strcpy(ligand3,ligand2);
    }
/* v.3.2--> */
  /* Store the tokens entered */
  strcpy(token[1],ligand2);
  strcpy(token[2],ligand3);
  strcpy(token[3],ligand4);
  strcpy(token[4],ligand5);
  *ntoken = 5;
  if (token[4][0] == '\0')
    {
      *ntoken = 4;
      if (token[3][0] == '\0')
	{
	  *ntoken = 3;
	  if (token[2][0] == '\0')
	    {
	      *ntoken = 2;
	      if (token[1][0] == '\0')
		*ntoken = 1;
	    }
	}
    }
/* <--v.3.2 */
  
}
/* v.4.0--> */
/***********************************************************************

check_for_dimplot  -  Check the ligplot.pdb file to test whether it
                      contains a flattened DIMPLOT from a previous run

***********************************************************************/

int check_for_dimplot(void)
{
  char line[LINELEN + 1];
  int dimplot, keep_reading;

  FILE *fil_pdb;

  /* Open the ligplot.pdb file */
  if ((fil_pdb = fopen("ligplot.pdb","r")) == NULL)
    {
      printf("\n*** Unable to open the ligplot.pdb file\n");
      exit(1);
    }

  /* Initialise variables */
  dimplot = FALSE;
  keep_reading = TRUE;

  while (fgets(line,LINELEN,fil_pdb) != NULL && keep_reading == TRUE)
    {
      /* If this is a DIMPLOT record indicating that this is a flattened
	 DIMPLOT, then can stop reading */
      if (!strncmp(line,"DIMPLOT",7))
	{
	  /* Set the flags */
	  dimplot = TRUE;
	  keep_reading = FALSE;
	}

      /* Otherwise, if have hit the first ATOM or HETATM record, stop
	 reading */
      else if (!strncmp(line,"ATOM",4) || !strncmp(line,"HETATM",6))
	keep_reading = FALSE;
    }

  /* Close the ligplot.pdb file */
  fclose(fil_pdb);

  /* Return whether this is a DIMPLOT */
  return(dimplot);
}
/* <--v.4.0 */
/***********************************************************************

get_command_arguments  -  Interpret the command line arguments

***********************************************************************/

void get_command_arguments(char *strptrs[],int num,
/* v.3.2--> */
			   char pdb_name[FILENAME_LEN],
/* <--v.3.2 */
/* v.4.0--> */
			   char res_name1[4],
/* <--v.4.0 */
			   char res_num1[6],
/* v.4.0--> */
			   char res_name2[4],
/* <--v.4.0 */
			   char res_num2[6],
			   char *chain_identi,
			   char hhb_name[FILENAME_LEN],
			   char nnb_name[FILENAME_LEN],
			   char bonds_name[FILENAME_LEN],
			   char title[TITLE_LEN])
{
/* v.4.0-->
  int
    lig_name;
<--v.4.0 */
/* v.3.2--> */
  int irange, itoken;
/* v.4.0--> */
  int res_num;
/* <--v.4.0 */
  char token[MAXTOKENS][FILENAME_LEN];
/* <--v.3.2 */
  char
/* v.4.0--> */
/*    ligand2[5],
    ligand3[5],
    ligand4[5],
    ligand5[5], */
/* <--v.4.0 */
    tmpstr1[LINELEN];


  /* Initialise */
  Print_as_is = FALSE;
  title[0] = '\0';
/* v.3.2--> */
  pdb_name[0] = '\0';
  res_num1[0] = '\0';
  res_num2[0] = '\0';
/* v.4.0--> */
  res_name1[0] = '\0';
  res_name2[0] = '\0';
/* <--v.4.0 */
  *chain_identi = ' ';
  for (itoken = 0; itoken < MAXTOKENS; itoken++)
    token[itoken][0] = '\0';
/* <--v.3.2 */

/* v.3.2--> (Replacement for previous code now in ligout.c) */
  /* If no parameters present, not even filename, then stop */
  if (num == 0)
    {
      printf("*** ERROR. No command-line arguments entered!\n");
      exit(-1);
    }

  /* If only a single parameter entered, then take this to be the filename
     and prompt user for remaining residue information */
  else if (num == 1)
    {
/* v.4.0--> */
/*      lig_name = ON; */
/* <--v.4.0 */

/* v.3.2--> */
/*      lig_search(strptrs[1],ligand2,ligand3,ligand4,ligand5); */
      strcpy(pdb_name,strptrs[1]);
      strcpy(token[0],strptrs[1]);

/* v.4.0--> */
      /* Check whether the file is the output file from dimer, giving
	 the residue interactions across a dimer or domain interface */
      if (!strcmp(pdb_name,"dimplot.pdb"))
	{
	  /* Have dimplot file, so set appropriate flag */
	  Interface_Plot = TRUE;
	}

      /* If this is a ligplot.pdb file, then check whether it contains
	 a flattened DIMPLOT from a previous run */
      if (!strcmp(pdb_name,"ligplot.pdb"))
	Interface_Plot = check_for_dimplot();

      /* If this is not a DIMPLOT, ask user for the ligand details */
      if (Interface_Plot == FALSE)
/* <--v.4.0 */
	lig_search(pdb_name,token,&num);
/* <--v.3.2 */
    }

  /* If two or more parameters, transfer the command-line parameters
     into the array to be processed */
  else
    {
      for (itoken = 0; itoken < num; itoken++)
	strcpy(token[itoken],strptrs[itoken + 1]);
    }

  /* Process the entered tokens */

  /* Initialise count of non-flag tokens */
  irange = 0;

  /* Loop over all the command-line arguments */
  for (itoken = 0; itoken < num; itoken++)
    {
      /* Check for the -h option */
      if (!strncmp(token[itoken],"-h",2))
	{
	  /* Title for plot to be entered by user - accept it */
	  printf("\nEnter title of your plot\n");
	  gets(title);

	  /* Blank out this token */
	  token[itoken][0] = '\0';
	}

      /* Check for the -m option signifying that ligand is a metal */
      else if (!strncmp(token[itoken],"-m",2))
	Metal_as_Ligand = TRUE;

      /* Check for the -w option signifying that ligand is a water */
      else if (!strncmp(token[itoken],"-w",2))
	Water_as_Ligand = TRUE;

      /* Check for the -r option signifying that ligplot.res file
	 to be generated */
      else if (!strncmp(token[itoken],"-r",2))
	Write_Res_File = TRUE;

      /* Check for the -x option signifying that ligand to be plotted
	 exactly as it is, without any flattening */
      else if (!strncmp(token[itoken],"-x",2))
	Print_as_is = TRUE;

      /* Otherwise, take the token to be the filename or part of
	 the residue-range, according to order of entry */
      else if (token[itoken][0] != '\0')
	{
	  /* If this is the first non-flag parameter, then it must
	     be the file name */
	  if (irange == 0)
	    {
	      strcpy(pdb_name,token[itoken]);

/* v.4.0--> */
	      /* Check whether the file is dimplot.pdb */
	      if (!strcmp(pdb_name,"dimplot.pdb"))
		Interface_Plot = TRUE;

	      /* If this is a ligplot.pdb file, then check whether it
		 contains a flattened DIMPLOT from a previous run */
	      if (!strcmp(pdb_name,"ligplot.pdb"))
		Interface_Plot = check_for_dimplot();
/* <--v.4.0 */
	    }

/* v.4.0--> */
	  /* If not an interface plot, check for residue range and
	     chain-id information */
	  if (Interface_Plot == FALSE)
	    {
/* <--v.4.0 */
	  /* If this is the first of the residue-range parameters,
	     then take it to be the first residue in the range */
/* v.4.0--> */
/*	  else if (irange == 1) */
	      if (irange == 1)
/* <--v.4.0 */
		{
		  /* Interpret the entered residue number */
		  strcpy(tmpstr1,token[itoken]);
/* v.4.0--> */
/*		  get_residue_number(tmpstr1,res_num1); */
		  res_num = get_residue_number(tmpstr1,res_num1,res_name1);

		  /* If not a residue number, then decrement range as
		     still need to pick up the residue number */
		  if (res_num == FALSE)
		    irange--;
/* <--v.4.0 */
		}

	      /* If this is the second of the residue-range parameters,
		 then take it to be the second residue in the range */
	      else if (irange == 2)
		{
		  /* Interpret the entered residue number */
		  strcpy(tmpstr1,token[itoken]);
/* v.4.0--> */
/*		  get_residue_number(tmpstr1,res_num2); */
		  res_num = get_residue_number(tmpstr1,res_num2,res_name2);

		  /* If not a residue number, then decrement range as
		     still need to pick up the residue number */
		  if (res_num == FALSE)
		    irange--;
/* <--v.4.0 */
		}

	      /* The third residue-range parameter will be the chain-id */
	      else if (irange == 3)
		*chain_identi = token[itoken][0];
/* v.4.0--> */
	    }
/* <--v.4.0 */

	  /* Increment count of residue-range parameters */
	  irange++;
	}
    }

  /* If only one residue entered for the range, then make the end-residue
     number the same */
  if (res_num2[0] == '\0')
/* v.4.0--> */
    {
/* <--v.4.0 */
      strcpy(res_num2,res_num1);
/* v.4.0--> */
      strcpy(res_name2,res_name1);
    }
/* <--v.4.0 */

  /* If no PDB filename entered, then abort */
  if (pdb_name[0] == '\0')
    {
      printf("*** ERROR. No filename entered!\n");
      exit(-1);
    }

  /* Determine what the names of the corresponding .hhb, .nnb and
     .bonds files are from the entered PDB filename */
  printf("Input file names:-\n");
  printf("   PDB file              %s\n",pdb_name);
  get_filename(pdb_name,hhb_name,".hhb");
  printf("   H-bonds file          %s\n",hhb_name);
  get_filename(pdb_name,nnb_name,".nnb");
  printf("   Non-bonded contacts   %s\n",nnb_name);
  get_filename(pdb_name,bonds_name,".bonds");
  printf("\n");

  /* The .bonds option only currently works where the .bonds file
     has been output as ligplot.bonds */
  if (strncmp(bonds_name,"ligplot.bonds",13))
    bonds_name[0] = '\0';

  /* Print confirmation of entered details */
/* v.4.0--> */
/*      printf("Residue range:  %s to %s   ",res_num1,res_num2); */
  if (Interface_Plot == FALSE)
    {
      printf("Residue range:  ");
      if (res_name1[0] != '\0')
	printf("%s ",res_name1);
      printf("%s to ",res_num1);
      if (res_name2[0] != '\0')
	printf("%s ",res_name2);
      printf("%s   ",res_num2);
/* <--v.4.0 */
      if (*chain_identi != '#')
	printf("Chain [%c]\n",*chain_identi);
      else
	printf("\n");
/* v.4.0--> */
    }
/* <--v.4.0 */

  /* If the input file is called 'ligplot.pdb' then want
     to leave the structure as it is */
  if (!strncmp(pdb_name,"ligplot.pdb",11))
      Print_as_is = TRUE;
/* <--v.3.2 */

   /* If no title entered, and option to use filename as the default
      title, then store the filename as the title */
   if (Include->Filename_for_Title == TRUE && Print_Title[0] == '\0' &&
/* v.4.0--> */
/*       Print_as_is == FALSE) */
       Print_as_is == FALSE && Interface_Plot == FALSE)
/* <--v.4.0 */
        strcpy(Print_Title,root_name);

  /* If the ligand is a water molecule, then make sure that water
     molecules are allowed */
  if (Water_as_Ligand == TRUE)
    Include->Waters = TRUE;

}
/***********************************************************************

check_if_in_ligand  -  Determine whether the given residue is a ligand
                       residue

***********************************************************************/

/* v.3.1--> */
/*int check_if_in_ligand(char chain, char res_name[4], char res_num[6]) */
int check_if_in_ligand(char chain, char res_name[4], char res_num[6],
		       int partial_match,char duplicate_res_name[4])
/* <--v.3.1 */
{
  int ires, inligand;

  /* Initialise answer */
  inligand = FALSE;
  ires = 0;
  duplicate_res_name[0] = '\0';

  /* Loop through the stored names of the residues making up the ligand
     to see if given residue is amongst them */
  while (ires < ligand_residues && inligand == FALSE)
    {
      /* First check if have match on residue number and chain-id */
      if (chain == lig_chain_store[ires] &&
	  !strncmp(res_name,lig_res_name_store[ires],3) &&
	  !strncmp(res_num,lig_res_num_store[ires],5))
	inligand  = TRUE;

/* v.3.1--> */
      /* If only a partial match required, check for match on residue
	 number and chain-id only */
      if (partial_match == TRUE)
	if (chain == lig_chain_store[ires] &&
	    !strncmp(res_num,lig_res_num_store[ires],5))
	  {
	    inligand  = TRUE;
	    strncpy(duplicate_res_name,lig_res_name_store[ires],3);
	    duplicate_res_name[3] = '\0';
	  }
/* <--v.3.1 */

      /* Increment ligand count */
      ires++;
    }

  /* Return answer */
  return(inligand);
}
/* v.3.1--> */
/***********************************************************************

get_ligand_start_end  -  Read through PDB file to get the start
                         and end residue positions of the ligand

***********************************************************************/

/* v.4.0--> */
/* void get_ligand_start_end(char pdb_name[FILENAME_LEN],char res_num1[6],
			  char res_num2[6],
			  char chain_id) */
void get_ligand_start_end(char pdb_name[FILENAME_LEN],char res_name1[4],
			  char res_num1[6],char res_name2[4],char res_num2[6],
			  char chain_id)
/* <--v.4.0 */
{
  char line[LINELEN + 1];
/* v.3.1.1--> */
/*  char atmnum[6], duplicate_resnam[4], resnam[4]; */
  char duplicate_resnam[4], resnam[4];
/* <--v.3.1.1 */
  char chain, last_chain, last_resnum[6], resnum[6];
/* v.3.1.2--> */
  char last_resnam[4];
/* <--v.3.1.2 */
  char icode, icode_from, icode_to, temp_num[5];

  int resno, resno_from, resno_to;
  int inligand, residue_count;
  int duplicate;
  int fstpos, ipos, ires;
/* v.3.1.1--> */
  int keep_reading;
/* <--v.3.1.1 */
/* v.4.0--> */
  int first_residue, got_end, got_ligand, last_residue;
/* <--v.4.0 */
/* v.4.0.2--> */
  int backwards;
/* <--v.4.0.2 */

  /* Open the pdb file */
  if ((fil_pdb = fopen(pdb_name,"r")) == NULL)
    {
      printf("\n*** Unable to open your PDB file [%s]\n",pdb_name);
      exit(1);
    }

  /* Initialise variables */
/* v.3.1.1--> */
  keep_reading = TRUE;
/* <--v.3.1.1 */
/* v.4.0--> */
  got_end = FALSE;
  got_ligand = FALSE;
/* <--v.4.0 */
  inligand = FALSE;
  ligand_residues = 0;
  residue_count = 0;

  /* Read in the data from the PDB file */
  printf("\nSearching for ligand residues in PDB file ...\n");

  /* Initialise variables */
/* v.3.1.2--> */
  last_resnam[0] = '\0';
/* <--v.3.1.2 */
  last_resnum[0] = '\0';
  last_chain = '\0';

  /* Split the residue numbers giving the range of residues making up the
     ligand into their residue number and insertion code components */

  /* First residue of range */
  strncpy(temp_num,res_num1,4);
  temp_num[4] = '\0';
  resno_from = atoi(temp_num);
  icode_from = res_num1[4];
    
  /* Second residue of range */
  strncpy(temp_num,res_num2,4);
  temp_num[4] = '\0';
  resno_to = atoi(temp_num);
  icode_to = res_num2[4];
    
/* v.4.0.2--> */
  /* Check whether residues are numbered backwards (!) */
  backwards = FALSE;
  if (resno_from > resno_to ||
      (resno_from ==resno_to && icode_from > icode_to))
    backwards = TRUE;
/* <--v.4.0.2 */

  /* Search through the PDB file to find residue-range defining the
     ligand */
/* v.3.1.1--> */
/*  while (fgets(line,LINELEN,fil_pdb) != NULL) */
  while (fgets(line,LINELEN,fil_pdb) != NULL && keep_reading == TRUE)
/* <--v.3.1.1 */
    {
      /* If this is an ATOM or HETATM record, process it */
      if (!strncmp(line,"ATOM",4) || !strncmp(line,"HETATM",6))
	{
	  /* Get the residue name, number and chain-ID */
	  strncpy(resnam,line+17,3);
	  resnam[3] = '\0';
	  strncpy(resnum,line+22,5);
	  resnum[5] = '\0';
	  chain = line[21];

	  /* Check whether this is a new residue */
/* v.3.1.2--> */
/*	  if (chain != last_chain || strncmp(resnum,last_resnum,5)) */
	  if (chain != last_chain || strncmp(resnum,last_resnum,5) ||
	      strncmp(resnam,last_resnam,3))
/* <--v.3.1.2 */
	    {
	      /* Split the residue number in an integer and insertion code */
	      strncpy(temp_num,resnum,4);
	      temp_num[4] = '\0';
	      resno = atoi(temp_num);
	      icode = resnum[4];

	      /* Determine whether this residue falls in the range defining
		 the ligand */
	      inligand = FALSE;
/* v.4.0--> */
	      first_residue = FALSE;
	      last_residue = FALSE;
/* <--v.4.0 */
/* v.4.0.2--> */
/*	      if (chain == chain_id && resno >= resno_from &&
		  resno <= resno_to) */
	      if (chain == chain_id &&
		  ((backwards == FALSE && resno >= resno_from &&
		    resno <= resno_to) ||
		  (backwards == TRUE && resno <= resno_from &&
		    resno >= resno_to)))
/* <--v.4.0.2 */
		{
		  /* Set flag to say residue is in ligand */
		  inligand = TRUE;

		  /* If residue number is start of range, check the
		     insertion code is within the range, too */
/* v.4.0.2--> */
/*		  if (resno == resno_from && icode < icode_from) */
		  if (resno == resno_from &&
		      ((backwards == FALSE && icode < icode_from) ||
		       (backwards == TRUE && icode > icode_from)))
/* <--v.4.0.2 */
		    inligand = FALSE;

/* v.4.0--> */
		  /* Check whether this is the start of the residue range */
		  if (resno == resno_from && icode == icode_from)
		    {
		      first_residue = TRUE;
		      got_ligand = TRUE;
		    }

		  /* Check whether this is the end of the residue range */
		  if (resno == resno_to && icode == icode_to)
		    last_residue = TRUE;
/* <--v.4.0 */

		  /* If residue number is end of range, check the
		     insertion code is within the range, too */
/* v.4.0.2--> */
/*		  if (resno == resno_to && icode > icode_to) */
		  if (resno == resno_to &&
		      ((backwards == FALSE && icode > icode_to) ||
		       (backwards == TRUE && icode < icode_to)))
/* <--v.4.0.2 */
		    inligand = FALSE;
		}
/* v.4.0--> */
	      else
		got_ligand = FALSE;
/* <--v.4.0 */

	      /* If this is a water molecule, then cannot be in the
		 ligand unless -w option selected in command-line
		 parameters */
	      if (!strncmp(resnam,"HOH",3) && Water_as_Ligand == FALSE)
		inligand = FALSE;

/* v.4.0--> */
	      /* If this is the start-residue and have a name for the
		 start-residue, then check whether this name corresponds */
	      if (first_residue == TRUE && res_name1[0] != '\0')
		{
		  /* If the residue name doesn't match, then it is not the
		     first ligand residue */
		  if (strncmp(resnam,res_name1,3))
		    got_ligand = FALSE;
		}

	      /* If haven't encountered a proper start for ligand, or have
                 hit the recognised end, then even if residue is in the
		 right range, then don't want it */
	      if (got_ligand == FALSE || got_end == TRUE)
		inligand = FALSE;

	      /* Check whether this is the predefined residue end */
	      if (last_residue == TRUE && res_name2[0] != '\0')
		{
		  /* If the residue name matches, then have the last
		     ligand residue */
		  if (!strncmp(resnam,res_name2,3))
		    got_end = TRUE;
		}
/* <--v.4.0 */

	      /* If residue is in the ligand, then process */
	      if (inligand == TRUE)
		{
		  /* Check whether a residue with the same residue
		     number has already been stored */
		  duplicate = check_if_in_ligand(chain,resnam,resnum,
						 TRUE,duplicate_resnam);

		  /* If have a duplicate, then print warning message */
		  if (duplicate == TRUE)
		    {
		      printf("\n");
		      printf("*** Warning. Duplicate res. no. in");
		      printf(" PDB file for ligand res: %s %c (%s & %s)\n",
			     resnum,chain,resnam,duplicate_resnam);
		      printf("\n");
		      inligand = TRUE;
		      Nwarnings++;
		    }

		  /* Check that maximum number of ligand residues
		     not exceeded */
		  if (ligand_residues > MAXLIGRES - 1)
		    {
		      printf("*** Warning. Maximum number of ligand ");
		      printf("residues (%d) exceeded!\n",MAXLIGRES);
		      Nwarnings++;
		    }

		  /* Store the residue details */
		  else
		    {
		      strncpy(lig_res_name_store[ligand_residues],resnam,3);
		      lig_res_name_store[ligand_residues][3] = '\0';
		      strncpy(lig_res_num_store[ligand_residues],resnum,5);
		      lig_res_num_store[ligand_residues][5] = '\0';
		      lig_chain_store[ligand_residues] = chain;
		      ligand_residues++;
		    }
		}

	      /* Save the residue details */
/* v.3.1.2--> */
	      strncpy(last_resnam,resnam,3);
	      last_resnam[3] = '\0';
/* <--v.3.1.2 */
	      strncpy(last_resnum,resnum,5);
	      last_resnum[5] = '\0';
	      last_chain = chain;

	      /* Increment the residue-count */
	      residue_count++;
	    }
	}

/* v.3.1.1--> */
      /* If this is an ENDMDL record, stop reading here */
      else if (!strncmp(line,"ENDMDL",6))
	keep_reading = FALSE;
/* <--v.3.1.1 */
    }

  /* Print counts of data read in */
  printf("   Number of residues read in        = %7d\n",residue_count);
  printf("   Number of ligand residues         = %7d\n",ligand_residues);

  if (ligand_residues == 0)
    {
      printf("\n*** No ligand residues found. Nothing to plot\n");
      exit(1);
    }
  else
    {
      /* Print the ligand residues identified */
      printf("\n");
      printf("   Ligand residues: ");
      ipos = 20;
      for (ires = 0; ires < ligand_residues; ires++)
	{
	  fstpos = 0;
	  if (lig_res_num_store[ires][0] == ' ')
	    fstpos = 1;
	  if (lig_res_num_store[ires][1] == ' ')
	    fstpos = 2;
	  if (lig_res_num_store[ires][2] == ' ')
	    fstpos = 3;
	  printf("%s %s",lig_res_name_store[ires],
		 lig_res_num_store[ires]+fstpos);
	  ipos = ipos + 9 - fstpos;
	  if (lig_chain_store[ires] != ' ')
	    {
	      printf("(%c)",lig_chain_store[ires]);
	      ipos = ipos + 3;
	    }
	  if (ires < ligand_residues - 1)
	    {
	      printf(" - ");
	      ipos = ipos + 3;
	    }
	  if (ipos > 60)
	    {
	      printf("\n                    ");
	      ipos = 20;
	    }
	}
      printf("\n\n");
    }
  
  /* Close the PDB file */
  fclose(fil_pdb);
}
/* <--v.3.1 */
/***********************************************************************

read_conec_records  -  Read through CONECT records to pick up any
                       non-ligand residues that are covalently
		       bonded to the ligand, and to store these.

***********************************************************************/

/* v.3.1--> */
/*void read_conec_records(int lig_atom_start,int lig_atom_end,
			char pdb_name[FILENAME_LEN],int *nlinks) */
void read_conec_records(char pdb_name[FILENAME_LEN],int *nlinks)
/* <--v.3.1 */
{
  char line[101];
  char atmnum[6];

  int done, have_link;
/* v.3.1--> */
/*  int iatom, icon, ipos, jatom, ncon, ntoligand;
  int inlig_atom[8], store_atom[8]; */
  int iatom, icon, ipos, jatom, ncon;
  int store_atom[8];
/* <--v.3.1 */
/* v.3.1.2--> */
  int atom_number, last_atom_number, nconnect;
  int keep_reading;
/* <--v.3.1.2 */

  struct hhb_info *hhb_info_ptr;

  /* Open the pdb file */
  if ((fil_pdb = fopen(pdb_name,"r")) == NULL)
    {
      printf("\n*** Unable to open your PDB file [%s]\n",pdb_name);
      exit(1);
    }
    
  /* Initialise variables */
  *nlinks = 0;
/* v.3.1--> */
/*  ntoligand = 0; */
/* <--v.3.1 */
/* v.3.1.2--> */
  last_atom_number = 0;
  nconnect = 0;
  keep_reading = TRUE;
/* <--v.3.1.2 */
   
  /* Read in the data from the PDB file */
  printf("\nChecking CONECT records in the PDB file ...\n");

  /* Search through the PDB file to scan the CONECT records */
  while (fgets(line,100,fil_pdb) != NULL)
    {
/* v.3.1.2--> */
/*      if (!strncmp(line,"CONECT",6))
	{ */

      /* If this is an ATOM or HETATM record, store the atom number */
      if ((!strncmp(line,"ATOM",4) || !strncmp(line,"HETATM",6)) &&
	       keep_reading == TRUE)
	{
	  /* Store atom-number of last atom read in (for use in whittling
	     out which CONECT records are to be processed) */
	  strncpy(atmnum,line+6,5);
	  atmnum[5] = '\0';
	  last_atom_number = atoi(atmnum);
	}

      /* If this is an ENDMDL record, then not interested in any
         of the following ATOM and HETATM records */
      else if (!strncmp("ENDMDL",line,6))
	keep_reading = FALSE;

      /* If this is a CONECT record, process it */
      else if (!strncmp(line,"CONECT",6))
	{
	  /* Get the atom-number of the first number in the CONECT record */
	  strncpy(atmnum,line+6,5);
	  atmnum[5] = '\0';
	  atom_number = atoi(atmnum);

	  /* If the atom number is within the range of atoms read in,
	     then process */
	  if (atom_number <= last_atom_number)
	    {
/* <--v.3.1.2 */
	      /* Extract all the atom numbers in the CONECT record */
	      ipos = 6;
	      ncon = 0;
	      done = FALSE;
	      while (done == FALSE && ncon < 5)
		{
		  /* Extract this atom number and store */
		  strncpy(atmnum,line+ipos,5);
/* v.3.2--> */
		  if (line[ipos] == '\0' || line[ipos] == '\n')
		    strncpy(atmnum,"     ",5);
/* <--v.3.2 */
		  atmnum[5] = '\0';

		  /* If have blank, then end of connections reached */
		  if (!strncmp(atmnum,"     ",5))
		    done = TRUE;

		  /* Otherwise, store the number */
		  else
		    {
		      store_atom[ncon] = atoi(atmnum);

		      /* Check whether the atom belongs to the ligand */
/* v.3.1--> */
/*		  if (store_atom[ncon] >= lig_atom_start &&
		      store_atom[ncon] <= lig_atom_end)
		    inlig_atom[ncon] = TRUE;
		  else
		    inlig_atom[ncon] = FALSE; */
/* <--v.3.1 */
		  
		      /* Increment variables */
		      ncon++;
		      ipos = ipos + 5;
		    }
		}

	      /* Consider each connection in turn */
	      iatom = store_atom[0];
	      for (icon = 1; icon < ncon; icon++)
		{
		  jatom = store_atom[icon];

		  /* Check whether we already have this link */
		  have_link = FALSE;

		  /* Set pointer to the first of the links */
		  hhb_info_ptr = first_hhb_info_ptr;

		  /* Search through all the stored links */
		  while (hhb_info_ptr != NULL)
		    {
		      /* Check if we have this link already */
		      if ((hhb_info_ptr->atom_number1 == iatom &&
			   hhb_info_ptr->atom_number2 == jatom) ||
			  (hhb_info_ptr->atom_number1 == jatom &&
			   hhb_info_ptr->atom_number2 == iatom))
			have_link = TRUE;

		      /* Get the pointer to the next link */
		      hhb_info_ptr = hhb_info_ptr->next_hhb_info_ptr;
		    }
		      
		  /* Add if this is a new link */
		  if (have_link == FALSE)
		    {
		      /* Allocate memory for structure to hold
			 current interaction */
		      hhb_info_ptr = (struct hhb_info *)
			malloc(sizeof(struct hhb_info));

		      /* If this is the first interaction stored, then
			 update pointer to head of linked list */
		      if (first_hhb_info_ptr == NULL)
			first_hhb_info_ptr = hhb_info_ptr;
	      
		      /* Otherwise, make previous interaction point to
			 the current one */
		      else
			last_hhb_info_ptr->next_hhb_info_ptr
			  = hhb_info_ptr;
	      
		      /* Save the current interaction's pointer */
		      last_hhb_info_ptr = hhb_info_ptr;

		      /* If either atom is in the ligand store this
			 connection as a proper CONECT */
/* v.3.1--> */
/*		  if (inlig_atom[0] == TRUE || inlig_atom[icon] == TRUE)
		    {
		      hhb_info_ptr->source = CONECT;
		      ntoligand++;
		    } */
/* <--v.3.1 */

		      /* Otherwise, store it as an EXTRA, possibly
			 for later use */
/* v.3.1--> */
/*		  else
		    hhb_info_ptr->source = EXTRA; */
		      hhb_info_ptr->source = CONECT;
/* <--v.3.1 */

		      /* Store the data */
		      hhb_info_ptr->atom_number1 = iatom;
		      hhb_info_ptr->atom_number2 = jatom;

		      /* Initialise the other fields in the
			 structure */

		      /* First atom's details */
		      strncpy(hhb_info_ptr->atom_type1,"    ",4);
		      hhb_info_ptr->atom_type1[4] = '\0';
		      strncpy(hhb_info_ptr->res_name1,"   ",3);
		      hhb_info_ptr->res_name1[3] = '\0';
		      strncpy(hhb_info_ptr->res_num1,"     ",5);
		      hhb_info_ptr->res_num1[5] = '\0';
		      hhb_info_ptr->x1 = 0.0;
		      hhb_info_ptr->y1 = 0.0;
		      hhb_info_ptr->z1 = 0.0;

		      /* Second atom's details */
		      strncpy(hhb_info_ptr->atom_type2,"    ",4);
		      hhb_info_ptr->atom_type2[4] = '\0';
		      strncpy(hhb_info_ptr->res_name2,"   ",3);
		      hhb_info_ptr->res_name2[3] = '\0';
		      strncpy(hhb_info_ptr->res_num2,"     ",5);
		      hhb_info_ptr->res_num2[5] = '\0';
		      hhb_info_ptr->chain1 = ' ';
		      hhb_info_ptr->chain2 = ' ';
		      hhb_info_ptr->x2 = 0.0;
		      hhb_info_ptr->y2 = 0.0;
		      hhb_info_ptr->z2 = 0.0;

		      /* Additional data fields */
		      hhb_info_ptr->bond_length = 0.0;
		      hhb_info_ptr->next_hhb_info_ptr = NULL;

		      /* Increment number of connections stored */
		      (*nlinks)++;
		    }
		}
/* v.3.1.2--> */
	    }
/* <--v.3.1.2 */
	}
    }

  /* Write out the number of pseudo H-bonds */
  printf("   Total number of CONECT records read in = %7d\n",*nlinks);
/* v.3.1--> */
/*  printf("   Number of covalent bonds to ligand     = %7d\n",ntoligand); */
/* <--v.3.1 */

  /* Close the PDB file */
  fclose(fil_pdb);
}
/***********************************************************************

verify_conec_records  -  Verify the stored CONECT links by calculating
                         the appropriate distances

***********************************************************************/

void verify_conec_records(char pdb_name[FILENAME_LEN])
{
  char line[101];
  char atmnum[6];

  int iatom;
  int delete_link, wanted;

  float dist, x, x1, x2, y, y1, y2, z, z1, z2;

  struct hhb_info *hhb_info_ptr, *previous_hhb_info_ptr;

  /* Initialise variables */
  previous_hhb_info_ptr = NULL;

  /* Open the pdb file */
  if ((fil_pdb = fopen(pdb_name,"r")) == NULL)
    {
      printf("\n*** Unable to open your PDB file [%s]\n",pdb_name);
      exit(1);
    }
    
  /* Read in the data from the PDB file */
  printf("\nVerifying covalent bonds to ligand ...\n");

  /* Search through the PDB file to pick up the atoms involved in the
     CONECT records */
  while (fgets(line,100,fil_pdb) != NULL)
    {
      /* If this is an ATOM or HETATM record, process it */
      if ((!strncmp(line,"ATOM",4) || !strncmp(line,"HETATM",6)))
	{
	  /* Extract the atom number */
	  strncpy(atmnum,line+6,5);
	  atmnum[5] = '\0';
	  iatom = atoi(atmnum);

	  /* Set pointer to the first of the links */
	  hhb_info_ptr = first_hhb_info_ptr;

	  /* Search through all the stored links to see if this atom
	     is amongst them */
	  while (hhb_info_ptr != NULL)
	    {
	      /* Check if the atom belongs to the current link */
	      if (iatom == hhb_info_ptr->atom_number1 ||
		  iatom == hhb_info_ptr->atom_number2)
		{
		  /* Retrieve the atomic coordinates */
		  sscanf(line+30," %f %f %f ",&x,&y,&z);

		  /* Determine whether this is the first or second
		     atom of the pair */
		  if (iatom == hhb_info_ptr->atom_number1)
		    {
		      /* Transfer the appropriate details from the PDB
			 record to the H-bond info record */
		      strncpy(hhb_info_ptr->atom_type1,line+12,4);
		      hhb_info_ptr->atom_type1[4] = '\0';
		      strncpy(hhb_info_ptr->res_name1,line+17,3);
		      hhb_info_ptr->res_name1[3] = '\0';
		      strncpy(hhb_info_ptr->res_num1,line+22,5);
		      hhb_info_ptr->res_num1[5] = '\0';
		      hhb_info_ptr->chain1 = line[21];

		      /* Store this atom's coordinates */
		      hhb_info_ptr->x1 = x;
		      hhb_info_ptr->y1 = y;
		      hhb_info_ptr->z1 = z;
		    }

		  /* If this is the second atom of the link, store the
		     details */
		  else
		    {
		      /* Transfer the appropriate details from the PDB
			 record to the H-bond info record */
		      strncpy(hhb_info_ptr->atom_type2,line+12,4);
		      hhb_info_ptr->atom_type2[4] = '\0';
		      strncpy(hhb_info_ptr->res_name2,line+17,3);
		      hhb_info_ptr->res_name2[3] = '\0';
		      strncpy(hhb_info_ptr->res_num2,line+22,5);
		      hhb_info_ptr->res_num2[5] = '\0';
		      hhb_info_ptr->chain2 = line[21];

		      /* Store this atom's coordinates */
		      hhb_info_ptr->x2 = x;
		      hhb_info_ptr->y2 = y;
		      hhb_info_ptr->z2 = z;
		    }
		}

	      /* Get the pointer to the next link */
	      hhb_info_ptr = hhb_info_ptr->next_hhb_info_ptr;
	    }
	}
    }

  /* Loop through the stored links to calculate distances and discard
     any links with missing atoms or which are clearly wrong */

  /* Set pointer to the first of the links */
  hhb_info_ptr = first_hhb_info_ptr;

  /* Loop through all the stored links */
  while (hhb_info_ptr != NULL)
    {
      /* Initialise variables */
      delete_link = FALSE;

      /* If have located both atoms, then use stored coordinates to
	 compute the distance between them */
      if (strncmp(hhb_info_ptr->atom_type1,"    ",4) &&
	  strncmp(hhb_info_ptr->atom_type2,"    ",4))
	{
	  /* Retrieve both atoms' coordinates */
	  x1 = hhb_info_ptr->x1;
	  y1 = hhb_info_ptr->y1;
	  z1 = hhb_info_ptr->z1;
	  x2 = hhb_info_ptr->x2;
	  y2 = hhb_info_ptr->y2;
	  z2 = hhb_info_ptr->z2;

	  /* Calculate the distance between them */
	  dist = (x1 - x2) * (x1 - x2)
	    + (y1 - y2) * (y1 - y2)
	      + (z1 - z2) * (z1 - z2);
	  dist = sqrt((double) dist);

	  /* Check the bond-length */
	  wanted = TRUE;
	  if (dist > 3.0)
	    {
	      printf("*** Warning. Long bond from CONECT ");
	      printf("[%5d ->%5d] Length %7.2f",
		     hhb_info_ptr->atom_number1,
		     hhb_info_ptr->atom_number2,dist);
	      if (dist > 4.0)
		{
		  printf(" - discarded\n");
		  wanted = FALSE;
		}
	      else
		printf("\n");
	      Nwarnings++;
	    }

	  /* Store the distance if not too long */
	  if (wanted == TRUE)
	    {
	      hhb_info_ptr->bond_length = dist;
	    }
	  else
	    delete_link = TRUE;
	}

      /* Otherwise, if one or both atoms not located, show warning message */
      else
	{
	  /* Print warning message */
	  printf("* Warning. Error in CONECT record joining atoms ");
	  printf("%d and %d. \n",hhb_info_ptr->atom_number1,
		 hhb_info_ptr->atom_number2);
	  printf("*               ");
	  if (!strncmp(hhb_info_ptr->atom_type1,"    ",4))
	    {
	      printf("Atom %d not found. ",hhb_info_ptr->atom_number1);
	    }
	  if (!strncmp(hhb_info_ptr->atom_type2,"    ",4))
	    {
	      printf("Atom %d not found. ",hhb_info_ptr->atom_number2);
	    }
	  printf("\n");
	  Nwarnings++;

	  /* Mark for deletion */
	  delete_link = TRUE;
	}

      /* If link to be deleted, do so by making pointers by-pass it */
      if (delete_link == TRUE)
	{
	  /* If this is the first link, then set the next link as the
	     first */
	  if (previous_hhb_info_ptr == NULL)
	    first_hhb_info_ptr = hhb_info_ptr->next_hhb_info_ptr;

	  /* Otherwise, make pointers bypass the current link */
	  else
	    previous_hhb_info_ptr->next_hhb_info_ptr
	      = hhb_info_ptr->next_hhb_info_ptr;

	  /* If this was the last in the list, then set the previous
	     pointer as the new terminus */
	  if (hhb_info_ptr->next_hhb_info_ptr == NULL)
	    last_hhb_info_ptr = previous_hhb_info_ptr;
	}

      /* Save the current pointer */
      previous_hhb_info_ptr = hhb_info_ptr;

      /* Get the pointer to the next link */
      hhb_info_ptr = hhb_info_ptr->next_hhb_info_ptr;
    }

  /* Close the PDB file */
  fclose(fil_pdb);
}
/***********************************************************************

get_hbplus_info  -  Extract the atom details from the interaction just
                    read in from the .hhb or .nnb file

***********************************************************************/

int get_hbplus_info(char line[LINELEN + 1],char atom_type1[5],
		    char atom_type2[5],char *chain1,char *chain2,
		    char res_name1[4],char res_name2[4],
		    char res_num1[6],char res_num2[6],float *bond_length,
		    int file_type,int ligplot_format)
{
  char bond_lenstr[5];

  int ipos, nonzero, valid;

  /* Initialise variables */
  valid = TRUE;

  /* Check for obviously unwanted records */
/* v.3.0.1--> */
  /* if (!strncmp(line,"ligplot",7) || !strncmp(line,"   ",3))  */
/* v.4.0--> */
/*  if (!strncmp(line,"ligplot",7) || !strncmp(line + 10,"   ",10)) */
  if ((int) strlen(line) < 10 ||
      (ligplot_format == TRUE && !strncmp(line + 21,"   ",3)) ||
      (ligplot_format == FALSE && !strncmp(line + 10,"   ",3)))
/* <--v.4.0 */
/* <--v.3.0.1 */
    valid = FALSE;

  /* Otherwise, extract the data according to the file format */
  else
    {
      /* LIGPLOT-format of the .hhb and .nnb files */
      if (ligplot_format == TRUE)
	{
	  /* Atom name */
	  strncpy(atom_type1,line+12,4);
	  atom_type1[4] = '\0';
	  strncpy(atom_type2,line+33,4);
	  atom_type2[4] = '\0';
		
	  /* Chain-id */
	  *chain1 = line[4];
	  *chain2 = line[25];
	  if ((*chain1) == '-')
	    *chain1 = ' ';
	  if ((*chain2) == '-')
	    *chain2 = ' ';
	  
	  /* Residue details */
	  strncpy(res_num1,line+6,5);
	  res_num1[5] = '\0';
	  strncpy(res_num2,line+27,5);
	  res_num2[5] = '\0';
	  strncpy(res_name1,line,3);
	  res_name1[3] = '\0';
	  strncpy(res_name2,line+21,3);
	  res_name2[3] = '\0';

	  /* Bond length */
	  strncpy(bond_lenstr,line+41,4);
	  bond_lenstr[4] = '\0';
	  *bond_length = atof(bond_lenstr);
	}

      /* HBPLUS-format hydrogen bonds file */
      else if (file_type == HHB_FILE || file_type == NNB_FILE)
	{
	  /* Atom name */
	  strncpy(atom_type1,line+19,4);
	  atom_type1[4] = '\0';
	  strncpy(atom_type2,line+39,4);
	  atom_type2[4] = '\0';
		
	  /* Chain-id */
	  *chain1 = line[10];
	  *chain2 = line[30];
	  if ((*chain1) == '-')
	    *chain1 = ' ';
	  if ((*chain2) == '-')
	    *chain2 = ' ';
	  
	  /* Residue details */
	  strncpy(res_num1,line+11,5);
	  res_num1[5] = '\0';
	  strncpy(res_num2,line+31,5);
	  res_num2[5] = '\0';
	  strncpy(res_name1,line+16,3);
	  res_name1[3] = '\0';
	  strncpy(res_name2,line+36,3);
	  res_name2[3] = '\0';

	  /* Bond length */
	  strncpy(bond_lenstr,line+45,4);
	  bond_lenstr[4] = '\0';
	  *bond_length = atof(bond_lenstr);
	}

      /* For the residue numbers, may need to replace any leading
	 zeros by spaces */

      /* First residue */
      nonzero = FALSE;
      for (ipos = 0; ipos < 3 && nonzero == FALSE; ipos ++)
	{
	  if (res_num1[ipos] == '0')
	    res_num1[ipos] = ' ';
	  else
	    nonzero = TRUE;
	}

      /* Replace hyphen in insertion code with a space */
      if (res_num1[4] == '-')
	res_num1[4] = ' ';
		
      /* Repeat for second residue */
      nonzero = FALSE;
      for (ipos = 0; ipos < 3 && nonzero == FALSE; ipos ++)
	{
	  if (res_num2[ipos] == '0')
	    res_num2[ipos] = ' ';
	  else
	    nonzero = TRUE;
	}

      /* Replace hyphen in insertion code with a space */
      if (res_num2[4] == '-')
	res_num2[4] = ' ';
    }

  /* Return status */
  return(valid);
}
/* v.3.1.2--> */
/***********************************************************************

check_for_metal  -  Check whether the given atom is a metal ion

***********************************************************************/

int check_for_metal(char atom_name[5])
{
  int imetal;
  int metal;

  static int nmetals = 25;

  static char *metal_name[] = {
    "CA  ", "FE  ", "ZN  ", " U  ", "CU  ",
    "MN  ", "CL  ", "MG  ", "NA  ", "HG  ",
    "CO  ", " K  ", "CD  ", "AU  ", "AG  ",
    "IN  ", "HO  ", " F  ", "PB  ", "AL  ",
    "YB  ", "BR  ", "....", "....", "...."
  };

  /* Initialise variables */
  metal = FALSE;

  /* Loop through all the metals */
  for (imetal = 0; imetal < nmetals && metal == FALSE; imetal++)
    {
      /* If name matches, then have a metal */
      if (!strncmp(atom_name,metal_name[imetal],2))
	metal = TRUE;
    }

  /* Return the answer */
  return(metal);
}
/* <--v.3.1.2 */
/***********************************************************************

read_hbplus_file  -  Read in the required H-bond and nonbonded contacts
                     data from the appropriate HBPLUS output file

***********************************************************************/

void read_hbplus_file(char file_name[FILENAME_LEN],int file_type,
		      int *nbonds)
{
  char file_desc[20], line[LINELEN + 1];
  char atom_type1[5], atom_type2[5], chain1, chain2, res_name1[4],
  res_name2[4], res_num1[6], res_num2[6];
/* v.3.1--> */
  char dummy[4];
/* <--v.3.1 */

  int have_file, inrange1, inrange2, ligplot_format, valid, wanted;

  float bond_length;

  struct hhb_info *hhb_info_ptr;

  /* Initialise variables */
  have_file = TRUE;
  *nbonds = 0;
  ligplot_format = FALSE;
  if (file_type == HHB_FILE)
    strcpy(file_desc,"hydrogen bonds");
  else
    strcpy(file_desc,"non-bonded contacts");

  /* Open the HBPLUS input file */
  if ((fil_hbplus = fopen(file_name,"r")) == NULL)
    {
      printf("\n*** Unable to open your %s file [%s]\n",file_desc,
	     file_name);
      printf("*** No %s will be shown on the plot\n",file_desc);
      have_file = FALSE;
      Nwarnings++;
    }

  /* If input file exists, read in the interactions */
  if (have_file == TRUE)
    {
      printf("\nReading in %s from file %s ...\n",file_desc,file_name);

      /* Loop while reading through the file */
      while (fgets(line,LINELEN,fil_hbplus) != NULL)
	{
	  /* Check if this line identifies the file as being in
	     LIGPLOT format */
	  if (!strncmp(line,"ligplot",7))
	    ligplot_format = TRUE;

	  /* Extract the relevant information from the line just
	     read in */
	  valid = get_hbplus_info(line,atom_type1,atom_type2,&chain1,
				  &chain2,res_name1,res_name2,
				  res_num1,res_num2,&bond_length,
				  file_type,ligplot_format);

	  /* If this is a valid interaction check whether it is
	     wanted */
	  if (valid == TRUE)
	    {
	      /* Initialise flag which determines whether this
		 interaction is required */
	      wanted = TRUE;

	      /* For non-bonded contacts, check whether this is a desired
                 contact */
	      if (file_type == NNB_FILE)
		{
/* v.3.2--> */
/*		  if ((atom_type1[1] != 'C' && atom_type1[1] != 'S') ||
		      (atom_type2[1] != 'C' && atom_type2[1] != 'S'))
		    wanted = FALSE; */

		  /* If hydrophobics only (as in the old version) check
		     that both atoms are either a carbon or a sulphur */
		  if (Include->Contact_Type == HYDROPHOBIC_ONLY)
		    {
		      if ((atom_type1[1] != 'C' && atom_type1[1] != 'S') ||
			  (atom_type2[1] != 'C' && atom_type2[1] != 'S'))
			wanted = FALSE;
		    }

		  /* If at least one of the atoms has to be a hydrophobic
		     then check for this */
		  else if (Include->Contact_Type == HYDROPHOBIC_ANY)
		    {
		      if (atom_type1[1] == 'C' || atom_type1[1] == 'S' ||
			  atom_type2[1] == 'C' || atom_type2[1] == 'S')
			wanted = TRUE;
		      else
			wanted = FALSE;
		    }

		  /* Otherwise, accept any contact */
		  else
		    wanted = TRUE;

		  /* If either atom is a water, then don't want a
		     non-bonded contact */
		  if (!strncmp(res_name1,"HOH",3) ||
		      !strncmp(res_name2,"HOH",3))
		    wanted = FALSE;

		  /* If either is a metal, then again don't want this
		     contact */
		  if (wanted == TRUE)
		    {
		      if (check_for_metal(atom_type1) == TRUE)
			wanted = FALSE;
		      else if (check_for_metal(atom_type2) == TRUE)
			wanted = FALSE;
		    }
/* <--v.3.2 */
		}

	      /* If waters not required and one of the atoms is a
		 water, then exclude */
	      if (Include->Waters == FALSE &&
		  (!strncmp("HOH",res_name1,3) ||
		   !strncmp("HOH",res_name2,3)))
		wanted = FALSE;

	      /* Check if either residue is in the ligand */
	      else if (wanted == TRUE)
		{
/* v.4.0--> */
		  /* If plotting interface interactions, then check whether
		     this interaction is wanted */
		  if (Interface_Plot == TRUE)
		    {
		      /* Set interaction as wanted */
		      wanted = TRUE;
		    }

		  /* Otherwise check whether one of the residues belongs
		     to a ligand */
		  else
		    {
/* <--v.4.0 */
		      /* Check if first residue belongs to the ligand */
/* v.3.1--> */
/*		  inrange1 = check_if_in_ligand(chain1,res_name1,res_num1); */
		      inrange1 = check_if_in_ligand(chain1,res_name1,
						    res_num1,FALSE,dummy);
/* <--v.3.1 */

		      /* Repeat for second residue */
/* v.3.1--> */
/*		  inrange2 = check_if_in_ligand(chain2,res_name2,res_num2); */
		      inrange2 = check_if_in_ligand(chain2,res_name2,
						    res_num2,FALSE,dummy);
/* <--v.3.1 */
		    
		      /* If either residue falls within the range defining the
			 ligand, then want to save this interaction */
		      if (inrange1 == FALSE && inrange2 == FALSE)
			wanted = FALSE;
		    }
		}

	      /* If this interaction is still wanted, then store it */
	      if (wanted == TRUE)
		{
		  /* Allocate memory for structure to hold current 
		     interaction */
		  hhb_info_ptr = (struct hhb_info *)
		    malloc(sizeof(struct hhb_info));

		  /* If this is the first interaction stored, then
		      update pointerto head of linked list */
		  if (first_hhb_info_ptr == NULL)
		    first_hhb_info_ptr = hhb_info_ptr;
	      
		  /* Otherwise, make previous interaction point to the
		     current one */
		  else
		    last_hhb_info_ptr->next_hhb_info_ptr = hhb_info_ptr;
	      
		  /* Save the current interaction's pointer */
		  last_hhb_info_ptr = hhb_info_ptr;

		  /* Store the data */
		  hhb_info_ptr->source = file_type;
		  hhb_info_ptr->atom_number1 = 0;
		  hhb_info_ptr->atom_number2 = 0;

		  /* Initialise the other fields in the structure */

		  /* First atom's details */
		  strncpy(hhb_info_ptr->atom_type1,atom_type1,4);
		  hhb_info_ptr->atom_type1[4] = '\0';
		  strncpy(hhb_info_ptr->res_name1,res_name1,3);
		  hhb_info_ptr->res_name1[3] = '\0';
		  strncpy(hhb_info_ptr->res_num1,res_num1,5);
		  hhb_info_ptr->res_num1[5] = '\0';
		  hhb_info_ptr->chain1 = chain1;
		  hhb_info_ptr->x1 = 0.0;
		  hhb_info_ptr->y1 = 0.0;
		  hhb_info_ptr->z1 = 0.0;

		  /* Second atom's details */
		  strncpy(hhb_info_ptr->atom_type2,atom_type2,4);
		  hhb_info_ptr->atom_type2[4] = '\0';
		  strncpy(hhb_info_ptr->res_name2,res_name2,3);
		  hhb_info_ptr->res_name2[3] = '\0';
		  strncpy(hhb_info_ptr->res_num2,res_num2,5);
		  hhb_info_ptr->res_num2[5] = '\0';
		  hhb_info_ptr->chain2 = chain2;
		  hhb_info_ptr->x2 = 0.0;
		  hhb_info_ptr->y2 = 0.0;
		  hhb_info_ptr->z2 = 0.0;

		  /* Distance between atoms */
		  hhb_info_ptr->bond_length = bond_length;

		  /* Additional data fields */
		  hhb_info_ptr->next_hhb_info_ptr = NULL;

		  /* Increment count of bonds stored */
		  (*nbonds)++;
		}
	    }
	}

      /* Close the input file */
      fclose(fil_hbplus);
    }
}
/***********************************************************************

special_resinc  -  Include specially defined residues in the parameter
                   file - ie  residues that are connected to other
		   residues connected to the ligand

***********************************************************************/

void special_resinc(char file_name[FILENAME_LEN],int file_type,
		    char special_res[MAXSPECIAL_RES][15])
{
  char file_desc[20], line[LINELEN + 1];
  char atom_type1[5], atom_type2[5], chain1, chain2, res_name1[4],
  res_name2[4], res_num1[6], res_num2[6];
  char check_chain, check_res_name[4], check_res_num[6], other_res_name[4];

  int have_file, ligplot_format, valid;
  int match1, match2, res_to_ligand, nspecial_res, sp_res;

  float bond_length;

  struct hhb_info *hhb_info_ptr;

  /* Initialise variables */
  have_file = TRUE;
  ligplot_format = FALSE;
  if (file_type == HHB_FILE)
    strcpy(file_desc,"hydrogen bonds");
  else
    strcpy(file_desc,"non-bonded contacts");

  /* Open the HBPLUS input file */
  if ((fil_hbplus = fopen(file_name,"r")) == NULL)
    {
      printf("\n*** Unable to open your %s file [%s]\n",file_desc,
	     file_name);
      printf("*** No %s will be shown on the plot\n",file_desc);
      have_file = FALSE;
      Nwarnings++;
    }

  /* If input file exists, read in the interactions */
  if (have_file == TRUE)
    {
      printf("\nReading in %s from file %s ...\n",file_desc,file_name);

      /* Determine how many special-residue records there are */
      nspecial_res = 0;
      for (sp_res = 0; sp_res < MAXSPECIAL_RES; sp_res++)
	{
	  /* If entry non-blank, then save its position */
	  if (strncmp(special_res[sp_res],"       ",7))
	    nspecial_res = sp_res + 1;
	}

      /* Loop while reading through the file */
      while (fgets(line,LINELEN,fil_hbplus) != NULL)
	{
	  /* Check if this line identifies the file as being in
	     LIGPLOT format */
	  if (!strncmp(line,"ligplot",7))
	    ligplot_format = TRUE;

	  /* Extract the relevant information from the line just
	     read in */
	  valid = get_hbplus_info(line,atom_type1,atom_type2,&chain1,
				  &chain2,res_name1,res_name2,
				  res_num1,res_num2,&bond_length,
				  file_type,ligplot_format);

	  /* If this is a valid interaction check whether it is
	     wanted */
	  if (valid == TRUE)
	    {
	      /* Check if this interaction matches any of the selection
		 criteria */
	      for (sp_res = 0; sp_res < nspecial_res; sp_res++)
		{
		  /* For the case where two residues that are non-ligand 
		     and hydrogen bonded to the ligand and the user 
		     wishes to H-bond these together */
		  match1 = FALSE;
		  if ((!strncmp(special_res[sp_res],res_name1,3) ||
		       !strncmp(special_res[sp_res],"***",3)) &&
		      (!strncmp(special_res[sp_res]+4,res_name2,3) ||
		       !strncmp(special_res[sp_res]+4,"***",3)))
		    match1 = TRUE;
		    
		  match2 = FALSE;
		  if ((!strncmp(special_res[sp_res],res_name2,3) ||
		       !strncmp(special_res[sp_res],"***",3)) &&
		      (!strncmp(special_res[sp_res]+4,res_name1,3) ||
		       !strncmp(special_res[sp_res]+4,"***",3)))
		    match2 = TRUE;
		  
		  /* If have a possible match, need to check if either
		     of the residues is H-bonded to the ligand */
		  if (match1 == TRUE || match2 == TRUE)
		    {
		      /* Need to loop through the stored H-bond links
			 to check if either residue is H-bonded to the
			 ligand */
		      res_to_ligand = FALSE;

		      /* Store which of the two residues needs to be
			 H-bonded to the ligand */
		      if (match1 == TRUE)
			{
			  strncpy(check_res_name,res_name1,3);
			  check_res_name[3] = '\0';
			  strncpy(check_res_num,res_num1,5);
			  check_res_num[5] = '\0';
			  check_chain = chain1;
			  strncpy(other_res_name,res_name2,3);
			  other_res_name[3] = '\0';
			}
		      else
			{
			  strncpy(check_res_name,res_name2,3);
			  check_res_name[3] = '\0';
			  strncpy(check_res_num,res_num2,5);
			  check_res_num[5] = '\0';
			  check_chain = chain2;
			  strncpy(other_res_name,res_name1,3);
			  other_res_name[3] = '\0';
			}

		      /* Set pointer to the first of the H-bond links
			 stored so far */
		      hhb_info_ptr = first_hhb_info_ptr;
		      
		      /* Search through all the stored links */
		      while (hhb_info_ptr != NULL && res_to_ligand == FALSE)
			{
			  /* Check if first residue is H-bonded to
			     the ligand */
			  if (hhb_info_ptr->source == HHB_FILE &&
			      ((!strncmp(hhb_info_ptr->res_name1,
					 check_res_name,3) &&
				!strncmp(hhb_info_ptr->res_num1,
					 check_res_num,5) &&
				hhb_info_ptr->chain1 == check_chain) ||
			       (!strncmp(hhb_info_ptr->res_name2,
					 check_res_name,3) &&
				!strncmp(hhb_info_ptr->res_num2,
					 check_res_num,5) &&
				hhb_info_ptr->chain2 == check_chain)))
			    res_to_ligand = TRUE;
			  
			  /* Get the pointer to the next link */
			  hhb_info_ptr = hhb_info_ptr->next_hhb_info_ptr;
			}

		      /* If the outer residue is a water, then only
			 include if waters to be shown on plot */
		      if (!strncmp(other_res_name,"HOH",3) &&
			  Include->Waters == FALSE)
			res_to_ligand = FALSE;

		      /* If have found an interaction between the
			 residue and the ligand, want to add the current
			 HBPLUS interaction to the list of hydrogen bonds */
		      if (res_to_ligand == TRUE)
			{
			  /* Allocate memory for structure to hold current 
			     interaction */
			  hhb_info_ptr = (struct hhb_info *)
			    malloc(sizeof(struct hhb_info));

			  /* If this is the first interaction stored, then
			     update pointer to head of linked list */
			  if (first_hhb_info_ptr == NULL)
			    first_hhb_info_ptr = hhb_info_ptr;

			  /* Otherwise, make previous interaction point
			     to the current one */
			  else
			    last_hhb_info_ptr->next_hhb_info_ptr
			      = hhb_info_ptr;
	      
			  /* Save the current interaction's pointer */
			  last_hhb_info_ptr = hhb_info_ptr;

			  /* Store the data */
			  hhb_info_ptr->source = file_type;
			  hhb_info_ptr->atom_number1 = 0;
			  hhb_info_ptr->atom_number2 = 0;

			  /* Initialise the other fields in the structure */

			  /* First atom's details */
			  strncpy(hhb_info_ptr->atom_type1,atom_type1,4);
			  hhb_info_ptr->atom_type1[4] = '\0';
			  strncpy(hhb_info_ptr->res_name1,res_name1,3);
			  hhb_info_ptr->res_name1[3] = '\0';
			  strncpy(hhb_info_ptr->res_num1,res_num1,5);
			  hhb_info_ptr->res_num1[5] = '\0';
			  hhb_info_ptr->chain1 = chain1;
			  hhb_info_ptr->x1 = 0.0;
			  hhb_info_ptr->y1 = 0.0;
			  hhb_info_ptr->z1 = 0.0;

			  /* Second atom's details */
			  strncpy(hhb_info_ptr->atom_type2,atom_type2,4);
			  hhb_info_ptr->atom_type2[4] = '\0';
			  strncpy(hhb_info_ptr->res_name2,res_name2,3);
			  hhb_info_ptr->res_name2[3] = '\0';
			  strncpy(hhb_info_ptr->res_num2,res_num2,5);
			  hhb_info_ptr->res_num2[5] = '\0';
			  hhb_info_ptr->chain2 = chain2;
			  hhb_info_ptr->x2 = 0.0;
			  hhb_info_ptr->y2 = 0.0;
			  hhb_info_ptr->z2 = 0.0;

			  /* Distance between atoms */
			  hhb_info_ptr->bond_length = bond_length;

			  /* Additional data fields */
			  hhb_info_ptr->next_hhb_info_ptr = NULL;
			}
		    }
		}
	    }
	}

      /* Close the input file */
      fclose(fil_hbplus);
    }
}
/***********************************************************************

is_hydrogen_atom  -  Check whether this is a hydrogen or a pseudo atom
                     from any of the many possible representations of
                     such in the PDB

***********************************************************************/

int is_hydrogen_atom(char atname[5])
{
  char atom_name[5];

  int is_hydrogen;

  /* Initialise */
  is_hydrogen = FALSE;
  strncpy(atom_name,atname,4);
  atom_name[4] = '\0';

  /* Check the atom name for any indicators to suggest this is either
     a hydrogen or a pseudo atom */
  if (atom_name[0] == 'H' || atom_name[1] == 'H' || atom_name[1] == 'Q' ||
      (atom_name[1] == 'D' && (atom_name[0] == ' ' || atom_name[0] == '1' ||
			       atom_name[0] == '2' || atom_name[0] == '3')))
    is_hydrogen = TRUE;

/* v.3.2--> */
  /* Check for one of the possible metals */
  if (!strncmp(atom_name,"HG  ",4) || !strncmp(atom_name,"HO  ",4))
    is_hydrogen = FALSE;
/* <--v.3.2 */

  return is_hydrogen;
}
/***********************************************************************

update_bond  -  Add/remove bond to/from the bond-list, or update its
                     bond-type

***********************************************************************/

void update_bond(struct coordinate *first_atom_ptr, 
		      struct coordinate *second_atom_ptr,int bond_type,
		      int bond_source,float bond_length)
{
  int got_it_already;

  struct bond *bond_ptr;

  /* Check that bond not already stored */

  /* Loop through the stored bond-list to see if already have bond */
  bond_ptr = first_bond_ptr;
  got_it_already = FALSE;
  while (bond_ptr != NULL && got_it_already == FALSE)
    {
      /* If this bond contains both atoms, then have it already */
      if ((bond_ptr->first_atom_ptr == first_atom_ptr &&
	   bond_ptr->second_atom_ptr == second_atom_ptr) ||
	  (bond_ptr->first_atom_ptr == second_atom_ptr &&
	   bond_ptr->second_atom_ptr == first_atom_ptr))
	{
	  got_it_already = TRUE;

	  /* If bond to be deleted, mark it so */
	  if (bond_type == DELETED)
	    bond_ptr->bond_type = DELETED;
	}

      /* Go to the next bond in the linked list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* If bond is a new bond to be added, then do so */
  if (bond_type != DELETED && got_it_already == FALSE)
    {
      /* Get pointer for memory next entry in bond-list */
      bond_ptr = (struct bond *)
	malloc(sizeof(struct bond));

      /* Store the bond details */
      bond_ptr->first_atom_ptr = first_atom_ptr;
      bond_ptr->second_atom_ptr = second_atom_ptr;
      bond_ptr->bond_type = bond_type;
      if (bond_type == COVALENT)
	bond_ptr->elastic = FALSE;
      else
	bond_ptr->elastic = TRUE;
      bond_ptr->rotatable_bond = FALSE;
      bond_ptr->ring_bond = FALSE;
      bond_ptr->end_bond = FALSE;
      bond_ptr->bond_order = 1;
      bond_ptr->bond_source = bond_source;
      bond_ptr->bond_length = bond_length;
      bond_ptr->checked = FALSE;
      bond_ptr->flattened = FALSE;
      bond_ptr->nfirst_atom_links = 0;
      bond_ptr->nsecond_atom_links = 0;
      bond_ptr->nbonds_from_first_atom = 0;
      bond_ptr->nbonds_from_second_atom = 0;
      bond_ptr->first_bond_link_ptr = NULL;
      bond_ptr->next_stack_ptr = NULL;
      bond_ptr->next_link_ptr = NULL;
      bond_ptr->next_bond_ptr = first_bond_ptr;

      /* Store this pointer as the last-encountered bond */
      first_bond_ptr = bond_ptr;
    }
}
/***********************************************************************

read_pdb_file  -  Read in the required atom coordinates from the
                  PDB file

***********************************************************************/

int read_pdb_file(char pdb_name[FILENAME_LEN])
{
  FILE *fa;
  char line[LINELEN + 1];
/* v.4.0--> */
/*  char accstring[8], res_name[4], res_num[6]; */
  char accstring[10], res_name[4], res_num[6];
/* <--v.4.0 */
/* v.3.1.2--> */
  char atmnum[6];
/* <--v.3.1.2 */
/* v.3.1--> */
  char dummy[4], last_res_name[4], occup[4];
/* <--v.3.1 */
  char last_res_num[6];
  char atom_partner[6];
  char atom_type[5], chain, last_chain;

  int asa_file, atom_counter, residue_counter, counter;
  int cpos, ipos, inligand, new_residue, found, keep_reading,
  want_atom, residue_in_ligand, residue_in_group, natoms_in_ligand,
  wanted;
  int iamino, ncopies, residues_stored;
  int hydrogen;
/* v.3.1.2--> */
  int atom_number, last_atom_number, nconnect;
/* <--v.3.1.2 */
/* v.4.0--> */
  int in_lig1, in_lig2, ligplot_pdb_file, surface1;
  int len;
/* <--v.4.0 */

  float accessibility, length, length2, x1, x2, y1, y2, z1, z2;

  struct coordinate *atom_ptr, *atom1_ptr, *atom2_ptr, *last_atom_ptr,
  *other_atom_ptr;
  struct hhb_info *hhb_info_ptr;
  struct residue *residue_ptr, *last_residue_ptr, *other_residue_ptr;

  static char *amino[] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
    "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
    "PRO", "SER", "THR", "TRP", "TYR", "VAL"
  };
    
  /* Open PDB file */
  if ((fa = fopen(pdb_name,"r")) ==  NULL)
    {
      printf("\n*** Unable to open your PDB file [%s]\n",pdb_name);
      exit(1);
    }

  /* Initialise variables */
  atom_counter = 0;
  residue_counter = 0;
  residues_stored = 0;
  natoms_in_ligand = 0;
  Maximum_accessibility = 0.0;
  first_atom_ptr = NULL;
  last_atom_ptr = NULL;
  first_residue_ptr = NULL;
  last_residue_ptr = NULL;
  residue_ptr = NULL;

  /* Read in the data from the PDB file */
  printf("\nReading in atom coordinates from PDB file ...\n");

  /* Initialise variables */
  asa_file = FALSE;
/* v.4.0--> */
  ligplot_pdb_file = FALSE;
/* <--v.4.0 */
  xsite_file = FALSE;
  xsite_probe[0] = '\0';
  counter = 0;
  keep_reading = TRUE;
/* v.3.1.2--> */
  last_atom_number = 0;
/* <--v.3.1.2 */
  last_res_num[0] = '\0';
/* v.3.1--> */
  last_res_name[0] = '\0';
/* <--v.3.1 */
  last_chain = '\0';
  inligand = FALSE;
/* v.3.1.2--> */
  nconnect = 0;
/* <--v.3.1.2 */
  new_residue = TRUE;
  residue_in_ligand = FALSE;
  residue_in_group = FALSE;
/* v.4.0--> */
  surface1 = TRUE;
/* <--v.4.0 */
  want_atom = FALSE;

/* v.4.0--> */
  /* Check whether the file is a .asa file output by NACCESS */
  len = strlen(pdb_name);
  if (!strncmp(pdb_name+(len - 4),".asa",4))
    asa_file = TRUE;
  if (!strncmp(pdb_name,"ligplot.pdb",11))
    ligplot_pdb_file = TRUE;
/* <--v.4.0 */

  /* Loop through the PDB file, picking up all the atoms that
     are required */
  while (fgets(line,LINELEN,fa) != NULL)
    {
/* v.4.0--> */
      /* Check for key words in a dimplot file */
      if (Interface_Plot == TRUE)
	{
	  /* If this is the first surface, then treat as ligand */
	  if (!strncmp(line,"Title: ",7) &&
	      Include->Filename_for_Title == TRUE)
	    strcpy(Print_Title,line+7);
	  else if (!strncmp(line,"SURFACE 1",9))
	    surface1 = TRUE;
	  else if (!strncmp(line,"SURFACE 2",9))
	    surface1 = FALSE;
	}
/* <--v.4.0 */

      /* If this is a HEADER record, check whether the PDB file is
	 a .asa file, output by ACCESS, containing all the atom
	 accessibilities */
      if (!strncmp(line,"HEADER",6))
	{
	  if (!strncmp(line+13,"ATOMIC ASAs",11))
	    asa_file = TRUE;
	}
	
      /* Check for the end of an NMR structure model */
      else if (!strncmp("ENDMDL",line,6))
	keep_reading = FALSE;

      /* If this is the TITLE record in a ligplot.pdb file, then
	 store it */
      else if (!strncmp("TITLE ",line,6) && Print_as_is == TRUE)
	{
	  if (Include->Filename_for_Title == TRUE
	      && Print_Title[0] == '\0')
	    {
	      /* Get the title and chop off at last non-blank character */
	      ipos = 0;
	      for (cpos = 0; cpos < TITLE_LEN && line[cpos + 7] != '\n';
		   cpos++)
		{
		  if (line[cpos + 7] != ' ')
		    {
		      ipos = cpos + 8;
		    }
		}
	      line[ipos] = '\0';
	      strcpy(Print_Title,line+7);
	    }
	}
      
      /* Check whether this is an XSITE file */
      else if (!strncmp(line,"XSITE",5))
	{
	  xsite_file = TRUE;
	  strncpy(xsite_probe,line+5,6);
	  xsite_probe[6] = '\0';
	}

      /* If this is an ATOM or HETATM record, process it */
      else if ((!strncmp(line,"ATOM",4) || !strncmp(line,"HETATM",6)) &&
	       keep_reading == TRUE)
	{
	  /* Get the residue name, number, chain-ID and atom type */
	  strncpy(res_name,line+17,3);
	  res_name[3] = '\0';
	  strncpy(res_num,line+22,5);
	  res_num[5] = '\0';
	  chain = line[21];
	  strncpy(atom_type,line+12,4);
	  atom_type[4] = '\0';

/* v.3.1.2--> */
	  /* Store atom-number of last atom read in (for use in whittling
	     out which CONECT records are to be stored */
	  strncpy(atmnum,line+6,5);
	  atmnum[5] = '\0';
	  last_atom_number = atoi(atmnum);
/* <--v.3.1.2 */

	  /* Initialise variables for this atom */
	  ncopies = 0;

	  /* Get the accessibility to update the maximum value */
	  
/* v.4.0--> */
	  /* If this is a ligplot.pdb file, then accessibility is beyond
	     the B-value field */
	  if (ligplot_pdb_file == TRUE)
	    {
	      strncpy(accstring,line+66,8);
	      accstring[8] = '\0';
	    }
/* <--v.4.0 */

	  /* If this is an ordinary PDB file, store the B-value
	     as the atom's accessibility */
/* v.4.0--> */
/*	  if (asa_file == FALSE) */
          else if (asa_file == FALSE)
/* <--v.4.0 */
	    {
	      strncpy(accstring,line+60,6);
	      accstring[6] = '\0';
	    }
	    
	  /* Otherwise, if this is a .asa file output by Simon Hubbard's
	     NACCESS program, then pick up the accessibility from the
	     appropriate location */
	  else
	    {
	      strncpy(accstring,line+65,7);
	      accstring[7] = '\0';
	    }
	  accessibility = atof(accstring);
	  if (accessibility > Maximum_accessibility)
	    Maximum_accessibility = accessibility;

	  /* Check whether atom belongs to a new residue */
/* v.3.1--> */
/*	  if (chain != last_chain || strncmp(res_num,last_res_num,5)) */
	  if (chain != last_chain || strncmp(res_num,last_res_num,5) ||
	      strncmp(res_name,last_res_name,3))
/* <--v.3.1 */
	    {
	      /* If residue is a new one, need to check whether it
		 either belongs to the ligand or to one of the
		 residues H-bonded to the ligand */
	      new_residue = TRUE;
	      residue_in_ligand = FALSE;
	      residue_in_group = FALSE;

	      /* Check whether the residue falls within the range
		 defining the ligand */
	      inligand = TRUE;

/* v.4.0--> */
	      /* For interface plot, check whether this residue belongs
		 to the first surface */
	      if (Interface_Plot == TRUE)
		{
		  inligand = surface1;

		  /* Check whether residue from the other surface
		     or water */
		  if (surface1 == FALSE)
		    residue_in_group = TRUE;
		}

	      /* Otherwise, check whether residue falls within the
		 sequential range defining the ligand */
	      else
/* <--v.4.0 */

/* v.3.1--> */
/*	      inligand = check_if_in_ligand(chain,res_name,res_num); */
		inligand = check_if_in_ligand(chain,res_name,res_num,
					      FALSE,dummy);
/* <--v.3.1 */

	      /* If the residue is a ligand residue, then want to
		 store all its atom coords */
	      if (inligand == TRUE)
		residue_in_ligand = TRUE;
		  
	      /* If residue is not one of the ligand residues, check
		 whether it is one of the residues H-bonded to the
		 ligand, or involved in non-bonded contacts */
/* v.4.0--> */
/*	      else */
	      else if (Interface_Plot == FALSE)
/* <--v.4.0 */
		{
		  inligand = FALSE;

		  /* If printing structure as it stands (ie reading
		     in from the ligplot.pdb file), then want
		     to read in all the atom records present */ 
		  if (Print_as_is == TRUE)
		    residue_in_group = TRUE;

		  /* Otherwise, check the HBPLUS information */
		  else
		    {
		      /* Set pointer to the first of the links */
		      hhb_info_ptr = first_hhb_info_ptr;

		      /* Search through all the stored links */
		      while (hhb_info_ptr != NULL)
			{
			  /* Initialise flag */
			  wanted = FALSE;

			  /* Check if this residue is in the current
			     interaction */
			  if (chain == hhb_info_ptr->chain1 &&
			      !strncmp(res_name,hhb_info_ptr->res_name1,3) &&
			      !strncmp(res_num,hhb_info_ptr->res_num1,5))
			    wanted = TRUE;
			  else if (chain == hhb_info_ptr->chain2 &&
				   !strncmp(res_name,
					    hhb_info_ptr->res_name2,3) &&
				   !strncmp(res_num,
					    hhb_info_ptr->res_num2,5))
			    wanted = TRUE;

/* v.4.0--> */
			  /* If this is a covalent bond, then check
			     that is is between the ligand and a non-ligand
			     residue */
			  if (wanted == TRUE &&
			      hhb_info_ptr->source == CONECT)
			    {
			      /* Check if first residue is in the ligand */
			      in_lig1
				= check_if_in_ligand(hhb_info_ptr->chain1,
						     hhb_info_ptr->res_name1,
						     hhb_info_ptr->res_num1,
						     FALSE,dummy);

			      /* Check if first residue is in the ligand */
			      in_lig2
				= check_if_in_ligand(hhb_info_ptr->chain2,
						     hhb_info_ptr->res_name2,
						     hhb_info_ptr->res_num2,
						     FALSE,dummy);

			      /* If neither residue is in the ligand,
				 discount this bond */
			      if (in_lig1 == FALSE && in_lig2 == FALSE)
				wanted = FALSE;
			    }
/* <--v.4.0 */

			  /* If residue is included, then need to store
			     all its coordinates */
			  if (wanted == TRUE)
			    residue_in_group = TRUE;

			  /* Get the pointer to the next link */
			  hhb_info_ptr = hhb_info_ptr->next_hhb_info_ptr;
			}
		    }
		}

	      /* Save the residue details */
	      last_chain = chain;
	      strncpy(last_res_num,res_num,5);
	      last_res_num[5] = '\0';
/* v.3.1--> */
	      strncpy(last_res_name,res_name,3);
	      last_res_name[3] = '\0';
/* <--v.3.1 */
	      
	      /* Increment residue-count */
	      residue_counter++;
	    }

	  /* Determine whether current atom is wanted for a ligand
	     or H-group */
	  if (residue_in_ligand == TRUE || residue_in_group == TRUE)
	    want_atom = TRUE;
	  else
	    want_atom = FALSE;
	  
	  if (residue_in_ligand == TRUE)
	    inligand = TRUE;
	  
	  /* NMR structures in the pdb have hydrogens and pseudo-atoms, 
	     so we need to remove these */
	  hydrogen = is_hydrogen_atom(line+12);
	  if (hydrogen == TRUE && Include_Hydrogens != TRUE)
	    want_atom = FALSE;

	  /* If this is an alternate-occupancy atom, check whether
	     we already have it */
/* v.3.2--> */
/*	  if (want_atom == TRUE && line[16] != ' ') */
	  strncpy(occup,line+56,3);
	  occup[3] = '\0';
	  if (want_atom == TRUE && (line[16] != ' ' || strncmp(occup,"1.0",3)))
/* <--v.3.2 */
	    {
	      /* Initialise pointer to start of linked list of atom coords */
	      other_atom_ptr = first_atom_ptr;

	      /* Loop through all the atoms to check whether we
		 already have this one */
	      while (other_atom_ptr != NULL && want_atom == TRUE)
		{
		  /* Get pointer to the atom's residue */
		  other_residue_ptr = other_atom_ptr->residue_ptr;

		  /* If we already have this atom, then don't want
		     the second copy */
		  if (!strncmp(other_residue_ptr->res_name,res_name,3) &&
		      !strncmp(other_residue_ptr->res_num,res_num,5) &&
		      other_residue_ptr->chain == chain &&
		      !strncmp(other_atom_ptr->atom_type,atom_type,4))
		    want_atom = FALSE;

		  /* Get pointer to next atom in linked-list */
		  other_atom_ptr = other_atom_ptr->next;
		}	       
	    }

	  /* If the current atom is required, save its details */
	  if (want_atom == TRUE)
	    {	    	    
	      /* If the current atom is wanted, store it */

	      /* Allocate memory for structure to hold current atom's
		 details */
	      atom_ptr
		= (struct coordinate*)malloc(sizeof(struct coordinate));
	      if (atom_ptr == NULL)
		{
		  printf("*** Can't allocate memory for struct coordinate\n");
		  exit (1);
		}

	      /* If this is the first atom stored, then update pointer
		 to head of linked list */
	      if (first_atom_ptr == NULL)
		first_atom_ptr = atom_ptr;
	      
	      /* Otherwise, make previous atom point to the current
		 one */
	      else
		last_atom_ptr->next = atom_ptr;
	      
	      /* Save the current atom's pointer */
	      last_atom_ptr = atom_ptr;

	      /* Store the atom name and number */
	      strncpy(atom_ptr->atom_type,line+12,4);
	      atom_ptr->atom_type[4] = '\0';
	      strncpy(atom_ptr->print_name,atom_ptr->atom_type,4);
	      atom_ptr->print_name[4] = '\0';
	      strncpy(atom_ptr->atom_number,line+6,5);
	      atom_ptr->atom_number[5] = '\0';

	      /* Store the atom's coordinates */
	      sscanf(line+30," %f %f %f ",&atom_ptr->x,&atom_ptr->y,
		     &atom_ptr->z);
	      atom_ptr->original_x = atom_ptr->x;
	      atom_ptr->original_y = atom_ptr->y;
	      atom_ptr->original_z = atom_ptr->z;
	      atom_ptr->fit_x = atom_ptr->x;
	      atom_ptr->fit_y = atom_ptr->y;
	      atom_ptr->fit_z = atom_ptr->z;
	      atom_ptr->save_x = atom_ptr->x;
	      atom_ptr->save_y = atom_ptr->y;
	      
	      /* Store occupancy, B-balue and accessibility, where
		 present */
	      strncpy(atom_ptr->occupancy,line+54,6);
	      atom_ptr->occupancy[6] = '\0';
	      strncpy(atom_ptr->bvalue,line+60,6);
	      atom_ptr->bvalue[6] = '\0';
	      atom_ptr->accessibility = accessibility;
	      atom_ptr->atom_size = 0.0;
	      atom_ptr->deleted = FALSE;

	      /* If atom in ligand, increment count of ligand atoms */
	      if (inligand == TRUE)
		natoms_in_ligand++;
		
	      /* Mark all residue side chain atoms*/
	      if (!strncmp(atom_ptr->atom_type," C ",3)||
		  !strncmp(atom_ptr->atom_type," O ",3)||
		  !strncmp(atom_ptr->atom_type," CA ",4)||
		  !strncmp(atom_ptr->atom_type," N ",3)||
		  !strncmp(atom_ptr->atom_type," P ",3)||
		  !strncmp(atom_ptr->atom_type," ON ",4)||
		  !strncmp(atom_ptr->atom_type," OXT",4))
		atom_ptr->side_chain = FALSE;
	      else
		atom_ptr->side_chain = TRUE;

	      /* Initialise other fields for this atom */
	      atom_ptr->checked = FALSE;
	      atom_ptr->natom_links = 0;
	      atom_ptr->plot_atom = TRUE;
/* v.4.0--> */
	      atom_ptr->have_anchor = FALSE;
	      atom_ptr->anchor_pstn_x = 0.0;
	      atom_ptr->anchor_pstn_y = 0.0;
/* <--v.4.0 */
	      atom_ptr->next_stack_ptr = NULL;
	      atom_ptr->first_atom_link_ptr = NULL;
	      atom_ptr->next = NULL;

	      /* Check that this is a standard amino acid */
	      found = FALSE;
	      for (iamino = 0; iamino < 20 && found == FALSE; iamino++)
		if (!strncmp(res_name,amino[iamino],3))
		  found = TRUE;

	      /* If this is not a standard amino acid, then
		 mark atom as a main-chain atom */
	      if (Include->Simple_Ligand_Residues == TRUE && found == FALSE)
		atom_ptr->side_chain = FALSE;

	      /* If this is an XSITE file, and the print-field is
		 asterisked, read in the character-string to be
		 printed in place of the atom name on the final plot */
	      if (xsite_file == TRUE && line[67] == '*')
		{
		  strncpy(atom_ptr->print_name,line+68,4);
		  atom_ptr->print_name[4] = '\0';
		}

	      /* Increment counts */
	      counter++;
	      if (atom_counter >= MAXATOMS)
		{
		  printf("*** Maximum number of atoms");
		  printf(" (%d) exceeded. Program aborted.\n",
			 MAXATOMS);
		  exit(1);
		}
	      ncopies++;
	      atom_counter++;

	      /* If this is a new residue, create a residue record,
		 store residue data, and increment the residue-count */
	      if (new_residue == TRUE)
		{
		  /* Store this residue's details */
		  
		  /* Allocate memory for structure to hold current residue's
		     details */
		  residue_ptr
		    = (struct residue*)malloc(sizeof(struct residue));
		  if (residue_ptr == NULL)
		    {
		      printf("*** ERROR. Unable to allocate memory for");
		      printf(" struct residue\n");
		      printf("***        Program ligplot terminated with");
		      printf(" error.\n");
		      exit (1);
		    }

		  /* Store the residue details */
		  strncpy(residue_ptr->res_name,res_name,3);
		  residue_ptr->res_name[3] = '\0';
		  strncpy(residue_ptr->res_num,res_num,5);
		  residue_ptr->res_num[5] = '\0';
		  residue_ptr->chain = chain;
		  residue_ptr->natoms = 0;
		  residue_ptr->residue_type = STANDARD;
		  if (inligand == TRUE &&
		      Include->Simple_Ligand_Residues == TRUE)
		    residue_ptr->residue_type = SIMPLE_LIGAND;
		  if (inligand == FALSE &&
		      Include->Simple_Nonligand_Residues == TRUE)
		    residue_ptr->residue_type = SIMPLE_HGROUP;
		  residue_ptr->inligand = inligand;
		  residue_ptr->minx = 0.0;
		  residue_ptr->maxx = 0.0;
		  residue_ptr->miny = 0.0;
		  residue_ptr->maxy = 0.0;
		  residue_ptr->max_atom_size = 0.0;
/* v.3.2--> */
		  residue_ptr->deleted = FALSE;
/* <--v.3.2 */
/* v.4.0--> */
		  residue_ptr->original_mean_x = 0.0;
		  residue_ptr->original_mean_y = 0.0;
		  residue_ptr->original_mean_z = 0.0;
		  residue_ptr->flattened_mean_x = 0.0;
		  residue_ptr->flattened_mean_y = 0.0;
		  residue_ptr->have_anchor = FALSE;
		  residue_ptr->anchor_pstn_x = 0.0;
		  residue_ptr->anchor_pstn_y = 0.0;
		  residue_ptr->nanchored_atoms = 0;
/* <--v.4.0 */
		  residue_ptr->first_atom_ptr = atom_ptr;
		  residue_ptr->object_ptr = NULL;
		  residue_ptr->next_residue_ptr = NULL;

		  /* If this is the first residue stored, then update
		     pointer to head of linked list */
		  if (first_residue_ptr == NULL)
		    first_residue_ptr = residue_ptr;
		    
		  /* Otherwise, make previous residue point to the current
		     one */
		  else
		    last_residue_ptr->next_residue_ptr = residue_ptr;
		    
		  /* Save the current residue's pointer */
		  last_residue_ptr = residue_ptr;

		  residues_stored++;
		  new_residue = FALSE;
		}

	      /* Increment count of atoms in this residue */
	      residue_ptr->natoms++;

	      /* Get atom record to point to its residue record */
	      atom_ptr->residue_ptr = residue_ptr;
	    }
	}
      
      /* If this is a CONECT record, store the covalent bonds it reports */
      else if (!strncmp(line,"CONECT",6))
	{
/* v.3.1.2--> */
	  /* Get the atom-number of the first number in the CONECT record */
	  strncpy(atmnum,line+6,5);
	  atmnum[5] = '\0';
	  atom_number = atoi(atmnum);

	  /* If the atom number is within the range of atoms read in,
	     then process */
	  if (atom_number <= last_atom_number)
	    {
/* <--v.3.1.2 */

	      /* Search through all the stored atoms for a match with the
		 first atom-number in the CONECT record */
	      found = FALSE;

	      /* Initialise pointer to start of linked list of atom coords */
	      atom1_ptr = first_atom_ptr;

	      /* Loop through all the atoms in the linked list */
	      while (atom1_ptr != NULL && found == FALSE)
		{
		  /* Check if this atom's number matches */
		  if (!strncmp(atom1_ptr->atom_number,line+6,5))
		    {
		      /* If match found, store the atom's coordinates */
		      x1 = atom1_ptr->x;
		      y1 = atom1_ptr->y;
		      z1 = atom1_ptr->z;
		      found = TRUE;
		    }

		  /* If have found a match for the first atom, locate all
		     the other atoms it is reported as being bonded to */
		  if (found == TRUE)
		    {
		      /* Get the next atom-number */
		      ipos = 0;
		      strncpy(atom_partner,line+(5 * ipos + 11),5);
		      atom_partner[5] = '\0';
		  
		      /* Loop through all possible bonding partners */
		      while (strncmp(atom_partner,"     ",5) && ipos < 4)
			{
			  /* Search for this atom-partner in the atom
			     records stored */
			  found = FALSE;

			  /* Initialise pointer to start of linked list
			     of atom coords */
			  atom2_ptr = first_atom_ptr;

			  /* Loop through all the atoms in the linked list */
			  while (atom2_ptr != NULL && found == FALSE)
			    {
			      /* Initialise flag indicating whether bond
				 is required */
			      wanted = FALSE;

			      /* Check if this atom's number matches */
			      if (!strncmp(atom2_ptr->atom_number,
					   atom_partner,5))
				{
				  found = TRUE;

				  /* Want this covalent bond only if it
				     belongs entirely to the ligand (any
				     ligand-nonligand bonds will already have
				     been picked up in a previous pass) */
				  if (atom1_ptr->residue_ptr->inligand
				      == TRUE &&
				      atom2_ptr->residue_ptr->inligand == TRUE)
				    wanted = TRUE;
				}

			      /* If desired match found, store the atom's 
				 coordinates */
			      if (wanted == TRUE)
				{
				  x2 = atom2_ptr->x;
				  y2 = atom2_ptr->y;
				  z2 = atom2_ptr->z;

				  /* Check the bond-length */
				  length2 = (x2 - x1) * (x2 - x1)
				    + (y2 - y1) * (y2 - y1)
				      + (z2 - z1) * (z2 - z1);
				  if (length2 > 9.0)
				    {
				      printf("*** Warning. Long bond from ");
				      printf("CONECT [%5s -> %5s] %7.2f",
					     atom1_ptr->atom_number,
					     atom_partner,
					     sqrt((double) length2));
				      if (length2 > 16.0)
					{
					  printf(" -- discarded\n");
					  found = FALSE;
					}
				      else
					printf("\n");
				      Nwarnings++;
				    }
			      
				  /* Store this covalent bond */
				  if (found == TRUE)
				    {
				      length = sqrt((double) length2);
				      update_bond(atom1_ptr,atom2_ptr,
						  COVALENT,CONECT,length);
/* v.3.1.2--> */
				      nconnect++;
/* <--v.3.1.2 */
				    }
				}

			      /* Get pointer to next atom in linked-list */
			      atom2_ptr = atom2_ptr->next;
			    }

			  /* Get the next atom-number */
			  ipos++;
			  strncpy(atom_partner,line+(5 * ipos + 11),5);
			  atom_partner[5] = '\0';
			}
		    }

		  /* Get pointer to next atom in linked-list */
		  atom1_ptr = atom1_ptr->next;
		}
	    }
	}

/* v.4.0--> */
      /* If this is an ENDLIN record in a ligplot.pdb file, pick up the
	 maximum accessibility value */
      else if (ligplot_pdb_file == TRUE && !strncmp(line,"ENDLIN",6))
	{
	  /* Get the accessibility value */
	  strncpy(accstring,line+66,8);
	  accstring[8] = '\0';
	  
	  /* Save as the maximum value, if indeed higher than current
	     maximum */
	  accessibility = atof(accstring);
	  if (accessibility > Maximum_accessibility)
	    Maximum_accessibility = accessibility;
	}
/* <--v.4.0 */
    }
    
  /* Print counts of data read in */
  printf("   Number of residues stored         = %7d\n",residues_stored);
  printf("   Number of atoms stored            = %7d\n",atom_counter);
  printf("   Number of ligand atoms            = %7d\n",natoms_in_ligand);
/* v.3.1.2--> */
  printf("   Number of CONECT bonds stored     = %7d\n",nconnect);
/* <--v.3.1.2 */
  if (atom_counter == 0)
    {
      printf("\n*** No atoms found. Nothing to plot\n");
      exit(1);
    }
  if (natoms_in_ligand == 0)
    {
      printf("\n*** No ligand atoms found. Nothing to plot\n");
      exit(1);
    }

  /* Close the PDB file */
  fclose(fa);

/* v.4.0--> */
  /* Check that maximum accessibility is non-zero */
  if (Maximum_accessibility == 0.0)
    Maximum_accessibility = 100.0;
/* <--v.4.0 */

  /* Return number of atoms read in */
  return(atom_counter);
}
/***********************************************************************

print_object  -  Print current object (for debugging purposes)

***********************************************************************/

void print_object(struct object *object_ptr)
{
  int atom_no, iatom, iresid, natoms, nresid;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

/* v.4.0--> */
  /* Print object type */
  printf("Object type %d, number of residues = %d\n",
	 object_ptr->object_type,object_ptr->nresidues);
/* <--v.4.0 */
  
  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;
  atom_no = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, printing each one */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* Print atom details */
/* v.4.0--> */
/*	  printf("ATOM  %5d %s %s %c%s   %8.3f%8.3f%8.3f  1.00  1.00\n", */
	  printf("ATOM  %5d %s %s %c%s   %8.3f%8.3f%8.3f%s%s %7.3f\n",
/* <--v.4.0 */
		 atom_no + 1,
		 atom_ptr->atom_type,
		 atom_ptr->residue_ptr->res_name,
		 atom_ptr->residue_ptr->chain,
		 atom_ptr->residue_ptr->res_num,
		 atom_ptr->x,
		 atom_ptr->y,
/* v.4.0--> */
/*		 atom_ptr->z); */
		 atom_ptr->z,
		 atom_ptr->occupancy,
		 atom_ptr->bvalue,
		 atom_ptr->accessibility);
/* <--v.4.0 */

	  /* Get pointer to the next atom */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	  atom_no++;
	}

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

define_objects  -  Define the objects that will appear on the final
                   LIGPLOT diagram

***********************************************************************/

/* v.3.1--> */
/* void define_objects(int *nobjects) */
void define_objects(int *nobjects,int *nligand_objects)
/* <--v.3.1 */
{
  int iatom, natoms;
  int got_ligand_object;
/* v.4.0--> */
  int new_object;
/* <--v.4.0 */

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;
  struct object *object_ptr, *last_object_ptr;

  /* Initialise variables */
  got_ligand_object = FALSE;
  first_object_ptr = NULL;
  last_object_ptr = NULL;
  *nobjects = 0;
/* v.3.1--> */
  *nligand_objects = 0;
/* <--v.3.1 */

  /* Initialise pointer to the first stored residue */
  residue_ptr = first_residue_ptr;

  /* Loop through all residues */
  while (residue_ptr != NULL)
    {
/* v.4.0--> */
      /* Determine whether a new object is to be created */
      new_object = FALSE;
      
      /* If not a ligand residue, or the first of the ligand residues,
         need to create a new object record */
      if (residue_ptr->inligand == FALSE || got_ligand_object == FALSE)
	new_object = TRUE;

      /* If plotting a DIMPLOT from a flattened ligplot.pdb file, then
	 each residue needs to be a separate object */
      if (Interface_Plot == TRUE && Print_as_is == TRUE)
	new_object = TRUE;
/* <--v.4.0 */

/* v.4.0--> */
      /* If this is a new object, create a new object record */
/*      if (residue_ptr->inligand == FALSE || got_ligand_object == FALSE) */
      if (new_object == TRUE)
/* <--v.4.0 */
	{
	  /* Store this object's details */
		  
	  /* Allocate memory for structure to hold current object's
	     details */
	  object_ptr = (struct object*)malloc(sizeof(struct object));
	  if (object_ptr == NULL)
	    {
	      printf("*** ERROR. Unable to allocate memory for");
	      printf(" struct object\n");
	      printf("***        Program ligplot terminated with");
	      printf(" error.\n");
	      exit (1);
	    }

	  /* Initialise the object details */
	  object_ptr->nresidues = 0;
	  object_ptr->nbonds = 0;
	  object_ptr->nrot_bonds = 0;
	  object_ptr->nmain_chain = 0;
	  object_ptr->n_ca = 0;
	  object_ptr->weight = 0;
	  object_ptr->minx = 0.0;
	  object_ptr->miny = 0.0;
	  object_ptr->maxx = 0.0;
	  object_ptr->maxy = 0.0;
/* v.4.0--> */
	  object_ptr->interface = 0;
	  object_ptr->have_anchors = FALSE;
/* <--v.4.0 */
	  object_ptr->max_atom_size = 0.0;
	  object_ptr->place_x = 0.0;
	  object_ptr->place_y = 0.0;
	  object_ptr->internal_energy = 0.0;
	  object_ptr->total_energy = 0.0;
	  object_ptr->first_residue_ptr = residue_ptr;
	  object_ptr->first_object_bond_ptr = NULL;
	  object_ptr->next_object_ptr = NULL;

	  /* Define whether object is the ligand or a peripheral
	     residue */
	  if (residue_ptr->inligand == TRUE)
	    {
	      object_ptr->object_type = LIGAND;
	      got_ligand_object = TRUE;
/* v.3.1--> */
	      (*nligand_objects)++;
/* <--v.3.1 */
	    }
	  else if (!strncmp(residue_ptr->res_name,"HOH",3) ||
		   !strncmp(residue_ptr->res_name,"WAT",3))
	    object_ptr->object_type = WATER;
	  else
	    object_ptr->object_type = HGROUP;

/* v.3.1--> */
	  /* If this is a non-ligand residue, reset the got-ligand
	     flag */
	  if (residue_ptr->inligand == FALSE)
	    got_ligand_object = FALSE;
/* <--v.3.1 */

	  /* If this is the first object stored, then update
	     pointer to head of linked list */
	  if (first_object_ptr == NULL)
	    first_object_ptr = object_ptr;
		    
	  /* Otherwise, make previous object point to the current
	     one */
	  else
	    last_object_ptr->next_object_ptr = object_ptr;
		    
	  /* Save the current object's pointer */
	  last_object_ptr = object_ptr;

	  /* Increment count of objects */
	  (*nobjects)++;
	}

      /* Increment count of residues in this object */
      object_ptr->nresidues++;

      /* Store pointer from residue to object containing it */
      residue_ptr->object_ptr = object_ptr;

      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
	      
      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, counting number of
	 mainchain atoms */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* If this is a water, then mark as non-mainchain */
	  if (object_ptr->object_type == WATER)
	    atom_ptr->side_chain = TRUE;

	  /* If not a sidechain atom, increment count of this object's
	     mainchain atoms */
	  if (atom_ptr->side_chain == FALSE)
	    object_ptr->nmain_chain++;

	  /* Get pointer to the next atom */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

      /* Get pointer to the next residue in the list */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Print number of objects created */
  printf("   Number of objects                 = %7d\n",*nobjects);
}
/***********************************************************************

update_residue_boundaries  -  Calculate the maximum and minimum
                              coordinates of the given object's residues

***********************************************************************/

void update_residue_boundaries(struct object *object_ptr)
{
  int iatom, iresid, natoms, nresid;
/* v.3.2--> */
  int first_residue;
/* <--v.3.2 */
/* v.4.0--> */
  int ncoords;
  float flattened_mean_x, flattened_mean_y;
/* <--v.4.0 */

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

/* v.3.2--> */
  /* Initialise variables */
  first_residue = TRUE;
/* <--v.3.2 */

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
/* v.4.0--> */
      /* Initialise mean residue coordinates */
      flattened_mean_x = 0.0;
      flattened_mean_y = 0.0;
      ncoords = 0;
/* <--v.4.0 */

      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise residue limits to the coords of this first atom */
/* v.3.2--> */
/*      residue_ptr->minx = atom_ptr->x;
      residue_ptr->maxx = atom_ptr->x;
      residue_ptr->miny = atom_ptr->y;
      residue_ptr->maxy = atom_ptr->y; */
/* <--v.3.2 */

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

/* v.3.2--> */
      /* Initialise residue limits to the coords of this first atom,
         if there is one(!) */
      if (natoms > 0)
	{
	  residue_ptr->minx = atom_ptr->x;
	  residue_ptr->maxx = atom_ptr->x;
	  residue_ptr->miny = atom_ptr->y;
	  residue_ptr->maxy = atom_ptr->y;
	}
/* <--v.3.2 */

      /* Loop over all the residue's atoms, checking the x- and
	 y-coordinates */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* If this is a hydrophobic group, or simple residue, need
	     to fix all the atom-coords to be equal */
	  if (object_ptr->object_type == HYDROPHOBIC ||
	      object_ptr->object_type == SIMPLE_HGROUP)
	    {
	      atom_ptr->x = residue_ptr->minx;
	      atom_ptr->y = residue_ptr->miny;
/* v.4.0--> */
	      residue_ptr->flattened_mean_x = residue_ptr->minx;
	      residue_ptr->flattened_mean_y = residue_ptr->miny;
/* <--v.4.0 */
	    }

	  /* Otherwise check this atom's coordinates against the
	     maximum and minimum values */
	  else
	    {
	      if (atom_ptr->x < residue_ptr->minx)
		residue_ptr->minx = atom_ptr->x;
	      if (atom_ptr->x > residue_ptr->maxx)
		residue_ptr->maxx = atom_ptr->x;
	      if (atom_ptr->y < residue_ptr->miny)
		residue_ptr->miny = atom_ptr->y;
	      if (atom_ptr->y > residue_ptr->maxy)
		residue_ptr->maxy = atom_ptr->y;

/* v.4.0--> */
	      /* Update mean coords of the current x-y coords of the
		 atom in the residue's flattened state */
	      flattened_mean_x = flattened_mean_x + atom_ptr->x;
	      flattened_mean_y = flattened_mean_y + atom_ptr->y;
	      ncoords++;
/* <--v.4.0 */
	    }

	  /* Get pointer to the next atom */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

      /* If this is the object's first residue, set its boundaries
	 equal to that of the residue */
/* v.3.2--> */
/*      if (residue_ptr == object_ptr->first_residue_ptr) */
      if (natoms > 0)
	{
	  if (first_residue == TRUE)
/* <--v.3.2 */
	    {
	      object_ptr->minx = residue_ptr->minx;
	      object_ptr->maxx = residue_ptr->maxx;
	      object_ptr->miny = residue_ptr->miny;
	      object_ptr->maxy = residue_ptr->maxy;
/* v.3.2--> */
	      first_residue = FALSE;
/* <--v.3.2 */
	    }

	  /* Otherwise, adjust if necessary */
	  else
	    {
	      if (residue_ptr->minx < object_ptr->minx)
		object_ptr->minx = residue_ptr->minx;
	      if (residue_ptr->maxx > object_ptr->maxx)
		object_ptr->maxx = residue_ptr->maxx;
	      if (residue_ptr->miny < object_ptr->miny)
		object_ptr->miny = residue_ptr->miny;
	      if (residue_ptr->maxy > object_ptr->maxy)
		object_ptr->maxy = residue_ptr->maxy;
	    }
/* v.3.2--> */
	}
/* <--v.3.2 */

/* v.4.0--> */
      /* Get this residue's mean coordinates */
      if (ncoords > 0)
	{
	  /* Calculate the flattened-out residue's mean position, and
	     store the coords */
	  flattened_mean_x = flattened_mean_x / ncoords;
	  flattened_mean_y = flattened_mean_y / ncoords;

	  /* Store the object's original and flattened mean coords */
	  residue_ptr->flattened_mean_x = flattened_mean_x;
	  residue_ptr->flattened_mean_y = flattened_mean_y;
	}
/* <--v.4.0 */

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

update_all_boundaries  -  Update all object and residue boundaries

***********************************************************************/

void update_all_boundaries(void)
{
  struct object *object_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Update this object's boundaries */
      update_residue_boundaries(object_ptr);

      /* Get pointer to the next object */
      object_ptr = object_ptr->next_object_ptr;
    }
}
/* v.4.0--> */
/***********************************************************************

get_anchor_energy  -  Calculate the given object's distance from its
                      anchor-position and covert into an anchor energy

***********************************************************************/

float get_anchor_energy(struct object *object_ptr)
{
  int anchors, iatom, iresid, nanchored, natoms, nresid;

  float anchor_x, anchor_y, x, y;
  float dist2;
  float anchor_energy;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  anchor_energy = 0.0;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Process only if this residue has an anchor position */
      if (residue_ptr->have_anchor == TRUE)
	{
	  /* Get the residue's anchor coordinates */
	  anchor_x = residue_ptr->anchor_pstn_x;
	  anchor_y = residue_ptr->anchor_pstn_y;

	  /* Store the residue's coordinates */
	  x = residue_ptr->flattened_mean_x;
	  y = residue_ptr->flattened_mean_y;

	  /* Calculate distance between the residue's current mean
	     position and its anchoring position */
	  dist2 = (x - anchor_x) * (x - anchor_x)
	    + (y - anchor_y) * (y - anchor_y);

	  /* Calculate the equivalent energy */
	  anchor_energy = anchor_energy + Anchor_Weight * dist2;
	}

      /* If any of the residue's atoms have anchor positions, loop through
	 these */
      if (residue_ptr->nanchored_atoms > 0)
	{
	  /* Get pointer to first residue's first atom */
	  atom_ptr = residue_ptr->first_atom_ptr;
	  iatom = 0;
	  natoms = residue_ptr->natoms;
	  nanchored = 0;
	  anchors = residue_ptr->nanchored_atoms;

	  /* Loop over all this residue's atoms */
	  while (iatom < natoms && atom_ptr != NULL && nanchored < anchors)
	    {
	      /* If this is an anchored atom, process it */
	      if (atom_ptr->have_anchor == TRUE)
		{
		  /* Get the atom's anchor coordinates */
		  anchor_x = atom_ptr->anchor_pstn_x;
		  anchor_y = atom_ptr->anchor_pstn_y;

		  /* Store the atom's coordinates */
		  x = atom_ptr->x;
		  y = atom_ptr->y;

		  /* Calculate distance between the atom's current mean
		     position and its anchoring position */
		  dist2 = (x - anchor_x) * (x - anchor_x)
		    + (y - anchor_y) * (y - anchor_y);

		  /* Calculate the equivalent energy */
		  anchor_energy = anchor_energy + Anchor_Weight * dist2;

		  /* Increment count of anchores encountered */
		  nanchored++;
		}

	      /* Get pointer to the next atom */
	      atom_ptr = atom_ptr->next;
	      iatom++;
	    }
	}

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Return this object's energy */
  return(anchor_energy);
}
/***********************************************************************

get_rel_pstn_energy  -  Calculate the relative distances between the
                        given object's residues and all other objects'
                        residues to get an energy from a comparison of the
                        distances in the current flattened layout with
                        the distances in the original coordinates

***********************************************************************/

float get_rel_pstn_energy(struct object *object_ptr)
{
  int iresid, jresid, nother_resid, nresid;

  float other_x, other_y, other_z, x, y, z;
  float other_flattened_x, other_flattened_y, flattened_x, flattened_y;
  float dist2, original_dist2;
  float energy, rel_pstn_energy, weight;

  struct object *other_object_ptr;
  struct residue *residue_ptr, *other_residue_ptr;

  /* Initialise variables */
  rel_pstn_energy = 0.0;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Store the residue's coordinates */
      x = residue_ptr->original_mean_x;
      y = residue_ptr->original_mean_y;
      z = residue_ptr->original_mean_z;
      flattened_x = residue_ptr->flattened_mean_x;
      flattened_y = residue_ptr->flattened_mean_y;

      /* Initialise pointer to the first stored object */
      other_object_ptr = first_object_ptr;

      /* Loop through all objects */
      while (other_object_ptr != NULL)
	{
	  /* Process if not the same as the object of interest */
	  if (other_object_ptr != object_ptr)
	    {
	      /* Set the pointer to the first of this object's residues */
	      other_residue_ptr = other_object_ptr->first_residue_ptr;
	      nother_resid = other_object_ptr->nresidues;
	      jresid = 0;

	      /* Loop over all this object's residues */
	      while (jresid < nother_resid && other_residue_ptr != NULL)
		{
		  /* Get this residue's original and current mean coords */
		  other_x = other_residue_ptr->original_mean_x;
		  other_y = other_residue_ptr->original_mean_y;
		  other_z = other_residue_ptr->original_mean_z;
		  other_flattened_x = other_residue_ptr->flattened_mean_x;
		  other_flattened_y = other_residue_ptr->flattened_mean_y;

		  /* Calculate distance between these two residues in the
		     original structure */
		  original_dist2 = (x - other_x) * (x - other_x)
		    + (y - other_y) * (y - other_y)
		      + (z - other_z) * (z - other_z);

		  /* Calculate the weighting factor (such that nearby residues
		     in the original structure maintain their proximity better
		     than distant residues) */
		  if (original_dist2 > 0.01)
		    weight = Rel_Pstn_Energy_Weight / original_dist2;
		  else
		    weight = 0.0;

		  /* Repeat for distance between them for the current
		     flattened objects */
		  dist2 = (flattened_x - other_flattened_x)
		    * (flattened_x - other_flattened_x)
		      + (flattened_y - other_flattened_y)
			* (flattened_y - other_flattened_y);

		  /* Calculate the difference, and hence the effective
                     energy */
		  energy = weight * fabs(original_dist2 - dist2);
		  rel_pstn_energy = rel_pstn_energy + energy;

		  /* Get pointer to the next residue */
		  other_residue_ptr = other_residue_ptr->next_residue_ptr;
		  jresid++;
		}
	    }

	  /* Get pointer to the next object */
	  other_object_ptr = other_object_ptr->next_object_ptr;
	}

      /* Get pointer to the next residue of the current object */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Return this object's energy */
  return(rel_pstn_energy);
}
/***********************************************************************

get_boundary_energy  -  Calculate the given object's boundary energy
                        according to how far its atoms have strayed into
                        the wrong side of the plot (only relevant for
                        interface plots when objects are placed on either
                        side of a boundary, with any waters lying on,
                        or close to, the boundary)

***********************************************************************/

float get_boundary_energy(struct object *object_ptr)
{
  int iatom, iresid, natoms, nresid, object_type;
  int interface;

  float x;
  float dist, object_size;
  float boundary_energy, energy;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  boundary_energy = 0.0;

  /* Get this object type */
  object_type = object_ptr->object_type;
  interface = object_ptr->interface;

  /* Get the object size for calculating distance from boundary */
  object_size = 0.0;
  if (object_type == HYDROPHOBIC)
    object_size = Size_Val->Hydrophobics;
  else if (object_type == WATER)
    object_size = Size_Val->Waters;
  else if (object_type == LIGAND)
    object_size = Size_Val->Ligand_Atoms;
  else if (object_type == HGROUP)
    object_size = Size_Val->Nonligand_Atoms;
  else if (object_type == SIMPLE_HGROUP)
    object_size = Size_Val->Simple_Residues;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, checking whether
	 they are involved in interactions */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* Get this atom's x-coordinate */
	  x = atom_ptr->x;

	  /* Initialise "bad" distance from boundary */
	  dist = 0.0;

	  /* If object is a water molecule, want it fairly close to the
	     boundary */
	  if (object_type == WATER)
	     dist = fabs(x);

	  /* If residue belongs to first surface, want all its
	     atoms on the left-hand side of the boundary */
	  else if (interface == 1 && x > - (Min_Boundary_Dist + object_size))
	    dist = x - (Min_Boundary_Dist + object_size);

	  /* If residue belongs to second surface, want all its
	     atoms on the right-hand side of the boundary */
	  else if (interface == 2 && x < (Min_Boundary_Dist + object_size))
	    dist = (Min_Boundary_Dist + object_size) - x;

	  /* If distance is non-zero, calculate a corresponding
	     boundary-violation energy for this atom */
	  if (dist > 0.0)
	    {
	      energy = Boundary_Energy_Weight * dist * dist;

	      /* Add to total energy for this object */
	      boundary_energy = boundary_energy + energy;
	    }

	  /* Get pointer to the next atom */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Return this object's energy */
  return(boundary_energy);
}
/* <--v.4.0 */
/***********************************************************************

open_pdb_output  -  Open output PDB files, ligplot.pdb, ligplot.frm and
                    ligplot.res

***********************************************************************/

/* v.3.2--> */
/* void open_pdb_output(void) */
/* v.4.0--> */
/* void open_pdb_output(char res_num1[6],char res_num2[6],char chain) */
void open_pdb_output(char res_name1[4],char res_num1[6],char res_name2[4],
		     char res_num2[6],char chain)
/* <--v.4.0 */
/* <--v.3.2 */
{
/* v.3.2--> */
  int fstpos, ipos, ires;
/* <--v.3.2 */

  /* Open output PDB file containing the final flattened coordinates
     of the residues on the plot */
  if ((ligplot_pdb = fopen("ligplot.pdb","w")) ==  NULL)
    {
      printf("\n*** Unable to open output PDB file, ligplot.pdb\n");
      exit(1);
    }

  /* Open output PDB file for SURNET link */
  if ((ligplot_frm = fopen("ligplot.frm","w")) ==  NULL)
    {
      printf("\n*** Unable to open output PDB file, ligplot.frm\n");
      exit(1);
    }

/* v.4.0--> */
  /* Open residue centres of mass file, ligplot.rcm */
  if ((ligplot_rcm = fopen("ligplot.rcm","w")) ==  NULL)
    {
      printf("\n*** Unable to open output residue centres of mass file,");
      printf(" ligplot.rcm\n");
      exit(1);
    }
/* <--v.4.0 */

/* v.3.2--> */
  /* If writing out the ligplot.res file, then open it */
  if (Write_Res_File == TRUE)
    {
      if ((ligplot_res = fopen("ligplot.res","w")) ==  NULL)
	{
	  printf("\n*** Unable to open output file, ligplot.res\n");
	  exit(1);
	}

      /* Write the header records to the ligplot.res file */
      else
	{
	  /* Print the residue range */
/* v.4.0--> */
/*	  fprintf(ligplot_res,"#Residue range: %s to %s",res_num1,res_num2); */
	  fprintf(ligplot_res,"#Residue range: ");
	  if (res_name1[0] != '\0')
	    fprintf(ligplot_res,"%s ",res_name1);
	  fprintf(ligplot_res,"%s to ",res_num1);
	  if (res_name2[0] != '\0')
	    fprintf(ligplot_res,"%s ",res_name2);
	  fprintf(ligplot_res,"%s",res_num2);
/* <--v.4.0 */
	  if (chain != ' ')
	    fprintf(ligplot_res,", chain %c",chain);
	  fprintf(ligplot_res,"\n");

	  /* Print the ligand residues identified */
	  fprintf(ligplot_res,"#Ligand residues: ");
	  ipos = 20;
	  for (ires = 0; ires < ligand_residues; ires++)
	    {
	      /* If not been blanked out, then write out */
	      if (lig_chain_store[ires] != '\0')
		{
		  fstpos = 0;
		  if (lig_res_num_store[ires][0] == ' ')
		    fstpos = 1;
		  if (lig_res_num_store[ires][1] == ' ')
		    fstpos = 2;
		  if (lig_res_num_store[ires][2] == ' ')
		    fstpos = 3;
		  fprintf(ligplot_res,"%s %s",lig_res_name_store[ires],
			  lig_res_num_store[ires]+fstpos);
		  ipos = ipos + 9 - fstpos;
		  if (lig_chain_store[ires] != ' ')
		    {
		      fprintf(ligplot_res,"(%c)",lig_chain_store[ires]);
		      ipos = ipos + 3;
		    }
		  if (ires < ligand_residues - 1)
		    {
		      fprintf(ligplot_res," - ");
		      ipos = ipos + 3;
		    }
		  if (ipos > 60)
		    {
/* v.4.0--> */
/*		      fprintf(ligplot_res,"\n                    "); */
		      fprintf(ligplot_res,"\n");
		      fprintf(ligplot_res,"#                 ");
/* <--v.4.0 */
		      ipos = 20;
		    }
		}
	    }
	  fprintf(ligplot_res,"\n");
	}
    }
/* <--v.3.2 */
}
/***********************************************************************
  
get_bonds  -  Read in the bond connectivities from the .bonds file

***********************************************************************/

void get_bonds(char bonds_name[FILENAME_LEN])
{
  char line[LINELEN + 1];
  char atom_type[2][5], chain[2], res_name[2][5], res_num[2][6];
  char bond_type_char, ebond_char, rbond_char;
  char bond_no[5], bond_len_str[8];

  int bond_order, bond_source, bond_type;
  int ipos, lpos;
  int end_bond, elastic, rotatable_bond, wanted;

  float bond_length;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr, *bond_atom_ptr[2];
  struct residue *residue_ptr;
    
  /* Initialise variables */
  printf("\nReading in bonds from file  %s ...\n",bonds_name);

  /* Open the ligplot.bonds file */
  if ((fil_bonds = fopen(bonds_name,"r")) == NULL)
    {
      printf("\n*** Unable to open bonds file, %s\n",bonds_name);
      exit(1);
    }

  /* Read through the .bonds file */
  while (fgets(line,100,fil_bonds) != NULL)
    {
      /* Initialise flag */
      wanted = FALSE;

      /* Check this is a bond record */
      if (!strncmp(line,"Bond",4))
	{
	  /* Get the bond number */
	  strncpy(bond_no,line+4,4);
	  bond_no[4] = '\0';

	  /* Assume this bond is going to be OK, unless it turns out
	     otherwise */
	  wanted = TRUE;

	  /* Extract the bond details from the record */
	  lpos = 22;

	  /* Store both atoms */
	  for (ipos = 0; ipos < 2; ipos++)
	    {
	      /* Atom name */
	      strncpy(atom_type[ipos],line+lpos,4);
	      atom_type[ipos][4] = '\0';
	      lpos = lpos + 5;

	      /* Residue name */
	      strncpy(res_name[ipos],line+lpos,3);
	      res_name[ipos][3] = '\0';
	      lpos = lpos + 4;

	      /* Residue number */
	      strncpy(res_num[ipos],line+lpos,5);
	      res_num[ipos][5] = '\0';
	      lpos = lpos + 6;

	      /* Chain */
	      chain[ipos] = line[lpos];
	      lpos = lpos + 12;
	    }

	  /* Get bond type */
	  bond_type_char = line[66];
	  if (bond_type_char == 'H')
	    bond_type = HBOND;
	  else if (bond_type_char == 'c')
	    bond_type = COVALENT;
	  else if (bond_type_char == 'n')
	    bond_type = CONTACT;
	  else if (bond_type_char == 'i')
	    bond_type = INTERNAL;
	  else
	    {
	      printf("*** Warning. Invalid bond-type in ligplot.bonds");
	      printf(" file  [%c]\n",bond_type_char);
	      wanted = FALSE;
	      Nwarnings++;
	    }

	  /* Determine whether a rotatable bond */
	  end_bond = FALSE;
	  rotatable_bond = FALSE;
	  rbond_char = line[67];
	  if (rbond_char == 'r')
	    rotatable_bond = TRUE;
	  else if (rbond_char == 'e')
	    end_bond = TRUE;

	  /* Determine whether an elastic bond */
	  ebond_char = line[68];
	  elastic = FALSE;
	  if (ebond_char == 'E')
	    elastic = TRUE;

	  /* Get the bond length */
	  strncpy(bond_len_str,line+69,7);
	  bond_len_str[7] = '\0';
	  bond_length = atof(bond_len_str);

	  /* Get the bond order */
	  bond_order = 1;
	  if (line[77] == 'd')
	    bond_order = 2;
	  else if (line[77] == 't')
	    bond_order = 3;

	  /* Get the bond source */
	  bond_source = CALCULATED;
	  if (!strncmp(line+79,"HBPLUS",6))
	    bond_source = HBPLUS;
	  else if (!strncmp(line+79,"CONECT",6))
	    bond_source = CONECT;

	  /* Initialise pointers to the two bond atoms */
	  bond_atom_ptr[0] = NULL;
	  bond_atom_ptr[1] = NULL;

	  /* Initialise pointer to start of linked list of atom coords */
	  atom_ptr = first_atom_ptr;
	      
	  /* Loop through all the atoms to find the two making up the
	     current bond */
	  while (atom_ptr != NULL)
	    {
	      /* Get this atom's residue */
	      residue_ptr = atom_ptr->residue_ptr;

	      /* Check if this matches either of the bond atoms */
	      for (ipos = 0; ipos < 2; ipos++)
		{
		  /* If have match, then store this atom's pointer */
		  if (!strncmp(res_name[ipos],residue_ptr->res_name,3) &&
		      !strncmp(res_num[ipos],residue_ptr->res_num,5) &&
		      chain[ipos] == residue_ptr->chain &&
		      !strncmp(atom_type[ipos],atom_ptr->atom_type,4))
		    bond_atom_ptr[ipos] = atom_ptr;
		}

	      /* Get pointer to next atom in linked-list */
	      atom_ptr = atom_ptr->next;
	    }	       

	  /* Check whether matches for both the bond's atoms were found */
	  if (bond_atom_ptr[0] == NULL || bond_atom_ptr[1] == NULL)
	    {
	      /* Print warning message */
	      if (bond_atom_ptr[0] == NULL)
		{
		  printf("*** WARNING. Failed to find atom %s %s %s %c for ",
			 atom_type[0],res_name[0],res_num[0],chain[0]);
		  printf("bond number %s\n",bond_no);
		  Nwarnings++;
		}
	      if (bond_atom_ptr[1] == NULL)
		{
		  printf("*** WARNING. Failed to find atom %s %s %s %c for ",
			 atom_type[1],res_name[1],res_num[1],chain[1]);
		  printf("bond number %s\n",bond_no);
		  Nwarnings++;
		}
	      wanted = FALSE;
	    }
	}

      /* If this bond is wanted, then create a new bond record and store
	 the details */
      if (wanted == TRUE)
	{
	  /* Create a new bond record */
	  bond_ptr = (struct bond *)malloc(sizeof(struct bond));

	  /* Store the bond details */
	  bond_ptr->first_atom_ptr = bond_atom_ptr[0];
	  bond_ptr->second_atom_ptr = bond_atom_ptr[1];
	  bond_ptr->bond_type = bond_type;
	  bond_ptr->elastic = elastic;
	  bond_ptr->rotatable_bond = rotatable_bond;
	  bond_ptr->ring_bond = FALSE;
	  bond_ptr->end_bond = end_bond;
	  bond_ptr->bond_order = bond_order;
	  bond_ptr->bond_source = bond_source;
	  bond_ptr->bond_length = bond_length;
	  bond_ptr->checked = FALSE;
	  bond_ptr->flattened = FALSE;
	  bond_ptr->nfirst_atom_links = 0;
	  bond_ptr->nsecond_atom_links = 0;
	  bond_ptr->nbonds_from_first_atom = 0;
	  bond_ptr->nbonds_from_second_atom = 0;
	  bond_ptr->first_bond_link_ptr = NULL;
	  bond_ptr->next_stack_ptr = NULL;
	  bond_ptr->next_link_ptr = NULL;
	  bond_ptr->next_bond_ptr = first_bond_ptr;

	  /* Store this pointer as the last-encountered bond */
	  first_bond_ptr = bond_ptr;
	}
    }

  /*    printf("   Number of bonds in structure      = %7d\n",nbonds);
    if (bond_stack == 0)
      {
	printf("*** No bonds found in structure. Nothing to plot\n");
	exit(1);
      }
    printf("   No. of hydrogen bonds             = %7d\n",nhbonds);
    printf("   No. of linked bonds on stack      = %7d\n",bond_stack - 1);
    printf("   No. of ligand bonds on stack      = %7d\n",nligbonds);
    return bond_stack;  */
    
}
/***********************************************************************
  
covalent_connectivity  -  Find the covalent connectivity using a simple
                          distance cut-off

***********************************************************************/

void covalent_connectivity(void)
{
  int object_type1, object_type2, wanted1, wanted2;

  float x1, x2, y1, y2, z1, z2;
  float distance_real, distance_x, distance_y, distance_z, dist_sqrd;
  float length;

  struct coordinate *atom1_ptr, *atom2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  printf("\nCalculating bonds and bond-connectivity ...\n");

  /* Initialise pointer to start of linked list of atom coords */
  atom1_ptr = first_atom_ptr;
  length = 0.0;
    
  /* Loop through all the atoms in the linked list */
  while (atom1_ptr != NULL)
    {
      /* Get the atom's coordinates */
      x1 = atom1_ptr->x;
      y1 = atom1_ptr->y;
      z1 = atom1_ptr->z;

      /* Get the residue and object that the atom belongs to */
      residue1_ptr = atom1_ptr->residue_ptr;
      object_type1 = residue1_ptr->object_ptr->object_type;

      /* Determine whether it is worth proceeding with this atom */
      wanted1 = TRUE;
      if (object_type1 == WATER || object_type1 == HYDROPHOBIC)
	wanted1 = FALSE;

      /* Initialise pointer to start of linked list of atom coords */
      atom2_ptr = first_atom_ptr;
    
      /* Loop through all the atoms in the linked list */
      while (atom2_ptr != NULL && wanted1 == TRUE)
	{
	  /* Get the residue and object that the atom belongs to */
	  residue2_ptr = atom2_ptr->residue_ptr;
	  object_type2 = residue2_ptr->object_ptr->object_type;

	  /* Determine whether it is worth proceeding with this atom */
	  wanted2 = TRUE;
	  if (object_type2 == WATER || object_type2 == HYDROPHOBIC)
	    wanted2 = FALSE;

	  /* If atoms both belong to H-bond groups, but are not in
	     the same residue, then don't want any bonds between them */
	  if (object_type1 == HGROUP && object_type2 == HGROUP &&
	      residue1_ptr != residue2_ptr)
	    wanted2 = FALSE;

	  /* If one atom in the ligand, and the other isn't, then only
	     want this bond if external groups parameter is set */
	  if ((object_type1 == HGROUP && object_type2 == LIGAND) ||
	      (object_type1 == LIGAND && object_type2 == HGROUP))
	    {
	      if (Include->External_Bonds == FALSE)
		wanted2 = FALSE;
	    }

	  /* If both atoms are the same atom, then don't want distance
	     calculation */
	  if (atom1_ptr == atom2_ptr)
	    wanted2 = FALSE;

/* v.4.0--> */
	  /* If this is an interface plot and have different residues
	     belonging to the same surface, then not interested in any
	     covalent bonds between them */
	  if (Interface_Plot == TRUE && residue1_ptr != residue2_ptr)
	    wanted2 = FALSE;
/* <--v.4.0 */

	  /* Check for a covalent bond only if both atoms are wanted */
	  if (wanted2 == TRUE)
	    {
	      /* Get the second atom's coordinates */
	      x2 = atom2_ptr->x;
	      y2 = atom2_ptr->y;
	      z2 = atom2_ptr->z;

	      /* Calculate distance between atoms 1 & 2 */
	      distance_x = (x1 - x2) * (x1 - x2);
	      distance_y = (y1 - y2) * (y1 - y2);
	      distance_z = (z1 - z2) * (z1 - z2);
	      distance_real = distance_x + distance_y + distance_z;

	      /* Calculate cut-off distance for covalent connectivity,
		 based on atom types involved */
	      dist_sqrd = 3.6;
	      if (atom1_ptr->atom_type[1] == 'P' ||
		  atom2_ptr->atom_type[1] == 'P')
		dist_sqrd = 3.8;
	      if (!strncmp(atom1_ptr->atom_type,"FE",2) ||
		  !strncmp(atom2_ptr->atom_type,"FE",2))
		dist_sqrd = 5.0;
	      if (atom1_ptr->atom_type[1] == 'S' ||
		  atom2_ptr->atom_type[1] == 'S')
		dist_sqrd = 5.0;

	      /* If distance satisfies the distance criteria for covalent
		 bonding, then atoms are bonded */
	      if (distance_real < dist_sqrd)
		{
		  length = sqrt((double) distance_real);

		  /* If this is a disulphide bond, store it as an H-bond */
		  if (!strncmp(atom1_ptr->atom_type," S  ",4) &&
		      !strncmp(atom2_ptr->atom_type," S  ",4))
		    update_bond(atom1_ptr,atom2_ptr,HBOND,CALCULATED,
				length);

		  /* Otherwise, store as a covalent bond */
		  else
		    update_bond(atom1_ptr,atom2_ptr,COVALENT,CALCULATED,
				length);
		}
	    }
	  atom2_ptr = atom2_ptr->next;
	}
      atom1_ptr = atom1_ptr->next;
    }
}
/***********************************************************************

get_last  -  Get the last residue and atom of the given object

***********************************************************************/

void get_last(struct object *object_ptr,struct residue **last_residue_ptr,
	      struct coordinate **last_atom_ptr)
{
  int iatom, iresid, natoms, nresid;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise pointers */
  *last_atom_ptr = NULL;
  *last_residue_ptr = NULL;

  /* Loop through the first object's residues to get to the last
     residue and last atom */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop through all this object's residues, one by one */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Store this residue pointer */
      *last_residue_ptr = residue_ptr;

      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, counting number of
	 mainchain atoms */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* Store this atom pointer */
	  *last_atom_ptr = atom_ptr;

	  /* Get pointer to the next atom */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

get_pointer_to_residue  -  Get the residue that points to the given
                             residue

***********************************************************************/

void get_pointer_to_residue(struct residue *current_residue_ptr,
			    struct residue **to_residue_ptr)
{
  struct residue *residue_ptr;

  /* Initialise pointer */
  *to_residue_ptr = NULL;

  /* Get pointer to very first residue */
  residue_ptr = first_residue_ptr;
  
  /* Loop through all residues */
  while (residue_ptr != NULL && *to_residue_ptr == NULL)
    {
      /* If this residue points to the given residue then store
	 its pointer */
      if (residue_ptr->next_residue_ptr == current_residue_ptr)
	*to_residue_ptr = residue_ptr;

      /* Get pointer to next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }
}
/***********************************************************************

get_pointer_to_atom  -  Get the atom that points to the given
                             atom

***********************************************************************/

void get_pointer_to_atom(struct coordinate *current_atom_ptr,
			 struct coordinate **to_atom_ptr)
{
  struct coordinate *atom_ptr;

  /* Initialise pointer */
  *to_atom_ptr = NULL;

  /* Get pointer to very first atom */
  atom_ptr = first_atom_ptr;
  
  /* Loop through all atoms */
  while (atom_ptr != NULL && *to_atom_ptr == NULL)
    {
      /* If this atom points to the given atom then store
	 its pointer */
      if (atom_ptr->next == current_atom_ptr)
	*to_atom_ptr = atom_ptr;

      /* Get pointer to next atom */
      atom_ptr = atom_ptr->next;
    }
}
/***********************************************************************

join_residues  -  Combine the given objects' residues and atoms into a
                  single object

***********************************************************************/

/* v.3.1--> */
/* void join_residues(struct object *object1_ptr,struct object *object2_ptr) */
void join_residues(struct object *object1_ptr,struct object *object2_ptr,
		   int inligand)
/* <--v.4.0 */
/* <--v.3.1 */
{
  int iresid, nresid;

  struct coordinate *swap_atom_ptr;
  struct coordinate *last_atom_object1_ptr, *last_atom_object2_ptr;
  struct coordinate *to_first_atom_object1_ptr;
  struct coordinate *to_first_atom_object2_ptr;
  struct object *object_ptr, *to_object2_ptr;
  struct residue *residue_ptr, *swap_residue_ptr;
  struct residue *last_residue_object1_ptr, *last_residue_object2_ptr;
  struct residue *to_first_residue_object1_ptr;
  struct residue *to_first_residue_object2_ptr;

  /* Print that about to connect the residues */
  printf("\n");
  printf("*** Combining residues %s %s %c and %s %s %c",
	 object1_ptr->first_residue_ptr->res_name,
	 object1_ptr->first_residue_ptr->res_num,
	 object1_ptr->first_residue_ptr->chain,
	 object2_ptr->first_residue_ptr->res_name,
	 object2_ptr->first_residue_ptr->res_num,
	 object2_ptr->first_residue_ptr->chain);
  printf(" into a single object\n");
  Nwarnings++;

  /* Get the first object's last residue and last atom */
  get_last(object1_ptr,&last_residue_object1_ptr,&last_atom_object1_ptr);

  /* Repeat for second object */
  get_last(object2_ptr,&last_residue_object2_ptr,&last_atom_object2_ptr);

  /* Get the object that points to the second object */
  to_object2_ptr = NULL;
  object_ptr = first_object_ptr;
  
  /* Loop through all objects */
  while (object_ptr != NULL && to_object2_ptr == NULL)
    {
      /* If this object points to the second object, store its pointer */
      if (object_ptr->next_object_ptr == object2_ptr)
	to_object2_ptr = object_ptr;

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* Get the residue that points to the first object's first residue */
  get_pointer_to_residue(object1_ptr->first_residue_ptr,
			 &to_first_residue_object1_ptr);

  /* Get the residue that points to the second object's first residue */
  get_pointer_to_residue(object2_ptr->first_residue_ptr,
			 &to_first_residue_object2_ptr);

  /* Get the atom that points to the first object's first atom */
  get_pointer_to_atom(object1_ptr->first_residue_ptr->first_atom_ptr,
		      &to_first_atom_object1_ptr);

  /* Get the atom that points to the second object's first atom */
  get_pointer_to_atom(object2_ptr->first_residue_ptr->first_atom_ptr,
		      &to_first_atom_object2_ptr);

/* v.3.2--> */
  /* If adding a residue to a ligand, make sure the first object's
     residues are all identified as ligand residues */
  if (inligand == TRUE)
    {
      residue_ptr = object1_ptr->first_residue_ptr;
      nresid = object1_ptr->nresidues;
      iresid = 0;
 
      /* Loop through all this object's residues, one by one */
      while (iresid < nresid && residue_ptr != NULL)
	{
	  /* Set ligand marker */
	  residue_ptr->inligand = TRUE;

	  /* Get pointer to the next residue */
	  residue_ptr = residue_ptr->next_residue_ptr;
	  iresid++;
	}
    }
/* <--v.3.2 */

  /* Loop through all the residues of the second object and make them
     point to the first object */
  residue_ptr = object2_ptr->first_residue_ptr;
  nresid = object2_ptr->nresidues;
  iresid = 0;
 
  /* Loop through all this object's residues, one by one */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Amend the object pointer */
      residue_ptr->object_ptr = object1_ptr;

/* v.3.1--> */
      /* If in ligand, then mark as such */
      if (inligand == TRUE)
	residue_ptr->inligand = TRUE;
/* <--v.3.1 */

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Cut out the second object altogether */
  if (to_object2_ptr != NULL)
    to_object2_ptr->next_object_ptr = object2_ptr->next_object_ptr;
  else if (first_object_ptr == object2_ptr)
    first_object_ptr = object2_ptr->next_object_ptr;

/* v.3.1--> */
  /* If residues are already adjacent, then don't need to do anything */
  if (last_residue_object1_ptr->next_residue_ptr
      != object2_ptr->first_residue_ptr)
    {
/* <--v.3.1 */
      /* Move the second object's residues to their new location */
      swap_residue_ptr = last_residue_object1_ptr->next_residue_ptr;
      last_residue_object1_ptr->next_residue_ptr
	= object2_ptr->first_residue_ptr;
      if (to_first_residue_object2_ptr != NULL)
	to_first_residue_object2_ptr->next_residue_ptr
	  = last_residue_object2_ptr->next_residue_ptr;
      else if (first_residue_ptr == object2_ptr->first_residue_ptr)
	first_residue_ptr = last_residue_object2_ptr->next_residue_ptr;
      last_residue_object2_ptr->next_residue_ptr = swap_residue_ptr;
/* v.3.1--> */
    }

  /* If atoms are already adjacent, then don't need to do anything */
  if (last_atom_object1_ptr->next
      != object2_ptr->first_residue_ptr->first_atom_ptr)
    {
/* <--v.3.1 */
      /* Move the second object's atoms to their new location */
      swap_atom_ptr = last_atom_object1_ptr->next;
      last_atom_object1_ptr->next
	= object2_ptr->first_residue_ptr->first_atom_ptr;
      if (to_first_atom_object2_ptr != NULL)
	to_first_atom_object2_ptr->next
	  = last_atom_object2_ptr->next;
      else if (first_atom_ptr
	       == object2_ptr->first_residue_ptr->first_atom_ptr)
	first_atom_ptr = last_atom_object2_ptr->next;
      last_atom_object2_ptr->next = swap_atom_ptr;
/* v.3.1--> */
    }
/* <--v.3.1 */

  /* Update the first object's residue-count */
  object1_ptr->nresidues = object1_ptr->nresidues + object2_ptr->nresidues;
/* v.3.1--> */
  object1_ptr->nmain_chain = object1_ptr->nmain_chain
    + object2_ptr->nmain_chain;
  object1_ptr->n_ca = object1_ptr->n_ca + object2_ptr->n_ca;
/* <--v.3.1 */
}
/* v.3.1--> */
/***********************************************************************

check_for_join  -  Check whether the given two residues are connected by
                   a covalent bond

***********************************************************************/

int check_for_join(struct residue *residue1_ptr,
		   struct residue *residue2_ptr)
{
  char chain1, chain2, res_name1[4], res_name2[4], res_num1[6], res_num2[6];

  int connected, match11, match12, match21, match22;

  struct hhb_info *hhb_info_ptr;

  /* Initialise flag */
  connected = FALSE;

  /* Get pointer to first stored bond item in the
     list */
  hhb_info_ptr = first_hhb_info_ptr;

  /* Loop through all the bonds in the linked list */
  while (hhb_info_ptr != NULL)
    {
      /* If this is a covalent bond, check whether
	 it links the two residues */
/* v.3.1--> */
/*      if (hhb_info_ptr->source == CONECT ||
	  hhb_info_ptr->source == EXTRA) */
      if (hhb_info_ptr->source == CONECT)
/* <--v.3.1 */
	{
	  /* Get this bond's two residue details */
	  strncpy(res_name1,hhb_info_ptr->res_name1,3);
	  res_name1[3] = '\0';
	  strncpy(res_num1,hhb_info_ptr->res_num1,5);
	  res_num1[5] = '\0';
	  chain1 = hhb_info_ptr->chain1;
	  strncpy(res_name2,hhb_info_ptr->res_name2,3);
	  res_name2[3] = '\0';
	  strncpy(res_num2,hhb_info_ptr->res_num2,5);
	  res_num2[5] = '\0';
	  chain2 = hhb_info_ptr->chain2;

	  /* Check if the two residues are represented */

	  /* First residue matches res1 of bond */
	  match11 = FALSE;
	  if (!strncmp(residue1_ptr->res_name,
		       res_name1,3) &&
	      !strncmp(residue1_ptr->res_num,
		       res_num1,5) &&
	      residue1_ptr->chain == chain1)
	    match11 = TRUE;

	  /* First residue matches res2 of bond */
	  match12 = FALSE;
	  if (!strncmp(residue1_ptr->res_name,
		       res_name2,3) &&
	      !strncmp(residue1_ptr->res_num,
		       res_num2,5) &&
	      residue1_ptr->chain == chain2)
	    match12 = TRUE;

	  /* Second residue matches res1 of bond */
	  match21 = FALSE;
	  if (!strncmp(residue2_ptr->res_name,
		       res_name1,3) &&
	      !strncmp(residue2_ptr->res_num,
		       res_num1,5) &&
	      residue2_ptr->chain == chain1)
	    match21 = TRUE;

	  /* Second residue matches res2 of bond */
	  match22 = FALSE;
	  if (!strncmp(residue2_ptr->res_name,
		       res_name2,3) &&
	      !strncmp(residue2_ptr->res_num,
		       res_num2,5) &&
	      residue2_ptr->chain == chain2)
	    match22 = TRUE;

	  /* If it spans the two residues in question,
	     then have a link between them */
	  if ((match11 && match22) ||
	      (match12 && match21))
	    {
	      /* Set flag on */
	      connected = TRUE;

	      /* Upgrade this bond so that included
		 later on */
	      hhb_info_ptr->source = CONECT;
	    }
	}

      /* Get the pointer to the next link */
      hhb_info_ptr = hhb_info_ptr->next_hhb_info_ptr;
    }

  /* Return flag indicating whether connection between the two residues
     found or not */
  return(connected);
}
/* v.4.0--> */
/***********************************************************************

check_standard_aa  -  Check whether given residue is a standard amino
                      acid

***********************************************************************/

int check_standard_aa(char res_name[4])
{
  int standard;
  int iamino;

  static char *amino[] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
    "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
    "PRO", "SER", "THR", "TRP", "TYR", "VAL"
  };
    
  /* Initialise return-value */
  standard = FALSE;

  /* Check whether this is a standard amino acid */
  for (iamino = 0; iamino < 20 && standard == FALSE; iamino++)
    if (!strncmp(res_name,amino[iamino],3))
      standard = TRUE;

  /* Return whether given residue is a standard amino acid */
  return(standard);
}
/* <--v.4.0 */
/***********************************************************************

combine_objects  -  Check for any H-group objects that need to be
                    combined (ie one is an attachment for another)

***********************************************************************/

void combine_objects(void)
{
/* v.3.1.1--> */
/*  int connected, inligand, joined, match11, match12, match21, match22; */
  int connected, inligand, joined;
/* <--v.3.1.1 */
  int iresid1, iresid2, nresid1, nresid2;
/* v.4.0--> */
  int standard;
/* <--v.4.0 */

  struct object *object1_ptr, *object2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise pointer to the first stored object */
  object1_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object1_ptr != NULL)
    {
      /* Only process if this is not a water */
      if (object1_ptr->object_type != WATER)
	{
	  /* Get object's first residue and number of residues */
	  residue1_ptr = object1_ptr->first_residue_ptr;
	  nresid1 = object1_ptr->nresidues;
	  iresid1 = 0;
	  joined = FALSE;

	  /* Loop over all this object's residues */
	  while (iresid1 < nresid1 && residue1_ptr != NULL &&
		 joined == FALSE)
	    {
	      /* Loop through all other non-ligand objects except this one */
	      object2_ptr = object1_ptr->next_object_ptr;

	      /* Loop through all other objects */
	      while (object2_ptr != NULL && joined == FALSE)
		{
		  /* Only consider if this is not a water */
		  if (object2_ptr->object_type != WATER)
		    {
		      /* Get object's first residue and number of residues */
		      residue2_ptr = object2_ptr->first_residue_ptr;
		      nresid2 = object2_ptr->nresidues;
		      iresid2 = 0;

		      /* Loop over all this object's residues */
		      while (iresid2 < nresid2 && residue2_ptr != NULL &&
			     joined == FALSE)
			{
			  /* Check whether residue number and chain match */
			  if (!strncmp(residue1_ptr->res_num,
				       residue2_ptr->res_num,5) &&
			      residue1_ptr->chain == residue2_ptr->chain)
			    {
			      /* Check through the list of bonds to confirm
				 that there is indeed a link between these
				 objects */
			      connected = check_for_join(residue1_ptr,
							 residue2_ptr);

/* v.4.0--> */
			      /* If one object is a ligand, and the other
				 isn't, then check that the residue from
				 the non-ligand object is not a standard
				 amino acid */
			      if (object1_ptr->object_type == LIGAND ||
				  object2_ptr->object_type == LIGAND)
				{
				  standard = FALSE;

				  /* For whichever object is the non-ligand
				     one, check whether the residue is a
				     standard amino acid */
				  if (object1_ptr->object_type != LIGAND)
				    standard = 
				      check_standard_aa(residue1_ptr->res_name);
				  else if (object2_ptr->object_type != LIGAND)
				    standard = 
				      check_standard_aa(residue2_ptr->res_name);

				  /* If non-ligand residue is a standard
				     amino acid, then don't want to join
				     the two objects */
				  if (standard == TRUE)
				    connected = FALSE;
				}
/* <--v.4.0 */

			      /* If have a link between the two residues,
				 then may be an attachment */
			      if (connected == TRUE)
				{
				  /* If either object is a ligand, then
				     make sure both are marked as such */
				  if (object1_ptr->object_type == LIGAND ||
				      object2_ptr->object_type == LIGAND)
/* v.4.0--> */
				    {
				      object1_ptr->object_type = LIGAND;
				      object2_ptr->object_type = LIGAND;
/* <--v.4.0 */
				      inligand = TRUE;
/* v.4.0--> */
				    }
/* <--v.4.0 */
				  else
				    inligand = FALSE;

				  /* Join the residues into a single object */
/* v.4.0--> */
/*				  join_residues(object1_ptr,object2_ptr,
						inligand); */
				  join_residues(object1_ptr,object2_ptr,
						inligand);
/* <--v.4.0 */
				  joined = TRUE;
				}
			    }

			  /* Get pointer to the next residue */
			  residue2_ptr = residue2_ptr->next_residue_ptr;
			  iresid2++;
			}
		    }

		  /* Get pointer to next object */
		  object2_ptr = object2_ptr->next_object_ptr;
		}

	      /* Get pointer to the next residue */
	      residue1_ptr = residue1_ptr->next_residue_ptr;
	      iresid1++;
	    }
	}

      /* Get pointer to next object */
      if (joined == FALSE)
	object1_ptr = object1_ptr->next_object_ptr;
    }
}
/* v.3.1--> */
/***********************************************************************

delete_false_ligands  -  Delete false ligands (eg protein residues with
                         the same residue numbers as the ligand residues)

***********************************************************************/

void delete_false_ligands(void)
{
/* v.3.1.1--> */
/*  int inligand; */
/* <--v.3.1.1 */
  int iresid, nligand, nligand_objects, nobject, nresid;
/* v.3.1.2--> */
  int delete_last, delete_ligand, is_metal, metal, water;
  int natoms;
/* <--v.3.1.2 */
/* v.3.2--> */
  int ires;
/* <--v.3.2 */

  struct object *object_ptr;
  struct residue *residue_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  nligand_objects = 0;
/* v.3.1.1--> */
  nobject = 0;
/* <--v.3.1.1 */
/* v.3.1.2--> */
  delete_last = FALSE;
  natoms = 0;
  nresid = 0;
  metal = FALSE;
  water = FALSE;
/* <--v.3.1.2 */

  /* Loop through all objects to count the number of ligand objects
     and make sure that still have more than one */
  while (object_ptr != NULL)
    {
      /* If this is a ligand, then increment count of ligand objects */
      if (object_ptr->object_type == LIGAND)
/* v.3.1.2--> */
/*	nligand_objects++; */
	{
	  /* Get the object's first residue pointer */
	  residue_ptr = object_ptr->first_residue_ptr;

	  /* If object is a water, and are plotting waters, then this
             is OK as the ligand */
	  water = FALSE;
	  if (!strncmp(residue_ptr->res_name,"HOH",3) ||
	      !strncmp(residue_ptr->res_name,"WAT",3))
	    if (Water_as_Ligand == TRUE)
	      water = TRUE;

	  /* If object is a metal, and are plotting metals, then this
             is OK as the ligand */
	  metal = FALSE;
	  if (object_ptr->nresidues == 1 && residue_ptr->natoms == 1)
	    {
	      /* Check if this might be a metal */
	      is_metal
		= check_for_metal(residue_ptr->first_atom_ptr->atom_type);
	      if (Metal_as_Ligand == TRUE && is_metal == TRUE)
		metal = TRUE;
	    }

	  /* Increment number of apparent ligands and get number of
	     residues in this object */
	  nligand_objects++;
	  nresid = object_ptr->nresidues;

	  /* Get the number of atoms in the object's first residue */
	  natoms = residue_ptr->natoms;
	}
/* <--v.3.1.2 */

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* If only have one ligand left, then don't need to delete any objects */
  if (nligand_objects < 2)
    return;

/* v.3.1.2--> */
  /* If last apparent ligand consisted of a single residue containing a
     single atom (and this wasn't a water), then take the penultimate
     ligand as the "true" one, rather than the last one (ie probably just
     a metal ion with the same residue number as a HET group, and probably
     the same number as some residue in the protein!) */
  if (nresid == 1 && natoms == 1 && water == FALSE && metal == FALSE)
    delete_last = TRUE;
/* <--v.3.1.2 */

  /* Initialise variables */
  nligand = 0;
  printf("*** Warning. Number of ligand objects identified: %d\n",
	 nligand_objects);
  printf("***          Retaining last one as ligand only\n");
  Nwarnings++;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* If this is a ligand, then increment count of ligand objects */
      if (object_ptr->object_type == LIGAND)
	nligand++;

      /* If a ligand object, then delete or keep it */
      if (object_ptr->object_type == LIGAND)
	{
	  printf("    Ligand object:  %s %s %c",
		 object_ptr->first_residue_ptr->res_name,
		 object_ptr->first_residue_ptr->res_num,
		 object_ptr->first_residue_ptr->chain);

	  /* If this is not the last ligand object, then delete */
	  nobject++;

/* v.3.1.2--> */
/*	  if (nobject < nligand_objects) */

	  /* Determine whether to delete this particular object */
	  delete_ligand = FALSE;
	  if (delete_last == FALSE && nobject < nligand_objects)
	    delete_ligand = TRUE;
	  else if (delete_last == TRUE && nobject != nligand_objects - 1)
	    delete_ligand = TRUE;

	  if (delete_ligand == TRUE)
/* <--v.3.1.2 */
	    {
	      /* Reclassify this object as a H-bond group */
	      object_ptr->object_type = HGROUP;
	      printf("  *** Marked as non-ligand\n");

	      /* Get pointer to object's first residue */
	      residue_ptr = object_ptr->first_residue_ptr;
	      nresid = object_ptr->nresidues;
	      iresid = 0;
	      
	      /* Loop over all this object's residues to classify them
		 as non-ligand residues */
	      while (iresid < nresid && residue_ptr != NULL)
		{
		  /* Define residue as non-ligand */
		  residue_ptr->inligand = FALSE;

/* v.3.2--> */
	      /* Loop through the list of ligand residues to remove this
	         entry from it */
	      for (ires = 0; ires < ligand_residues; ires++)
		{
		  /* If this is the one just deleted, then remove it
		     from the list */
		  if (lig_chain_store[ires] == residue_ptr->chain &&
		      !strncmp(lig_res_name_store[ires],
			       residue_ptr->res_name,3) &&
		      !strncmp(lig_res_num_store[ires],
			       residue_ptr->res_num,3))
		    {
		      lig_chain_store[ires] = '\0';
		      lig_res_name_store[ires][0] = '\0';
		      lig_res_num_store[ires][0] = '\0';
		    }
		}
/* <--v.3.2 */
		  /* Get pointer to the next residue */
		  residue_ptr = residue_ptr->next_residue_ptr;
		  iresid++;
		}
	    }

	  /* Keep the last ligand object */
	  else
	    printf("      Retained as ligand\n");
	}

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }
}
/* <--v.3.1 */
/***********************************************************************

add_hbplus_bonds  -  Match the atoms listed as H-bonded in the .hhb file
                     (or as involved in non-bonded contacts in the .nnb
                     file) with those read in from the PDB file

***********************************************************************/

void add_hbplus_bonds(void)
{
  int bond_source, bond_type, match, wanted;
/* v.4.0--> */
  int in_lig1, in_lig2;

  char dummy[4];
/* <--v.4.0 */

  float length, x1, x2, y1, y2, z1, z2;

  struct coordinate *atom_ptr, *atom1_ptr, *atom2_ptr;
  struct hhb_info *hhb_info_ptr;
  struct residue *residue_ptr;

  /* Set pointer to the first of the HBPLUS interactions */
  hhb_info_ptr = first_hhb_info_ptr;

  /* Loop through each of the atom-pairs corresponding to the
     hydrogen-bonds read in from the .hhb file and the non-bonded
     contacts read in from the .nnb file */
  while (hhb_info_ptr != NULL)
    {
      /* Determine bond type */
      wanted = TRUE;
      if (hhb_info_ptr->source == CONECT)
	bond_type = COVALENT;
      else if (hhb_info_ptr->source == HHB_FILE)
	bond_type = HBOND;
      else if (hhb_info_ptr->source == NNB_FILE)
	bond_type = CONTACT;
      else
	wanted = FALSE;

      /* Initialise pointer to start of linked list of atom coords */
      atom_ptr = first_atom_ptr;
      atom1_ptr = NULL;
      atom2_ptr = NULL;
      match = FALSE;

      /* Loop through all the atoms in the linked list */
      while (atom_ptr != NULL && match == FALSE && wanted == TRUE)
	{
	  /* Get this atom's residue */
	  residue_ptr = atom_ptr->residue_ptr;

	  /* Check whether atom corresponds to the first of the H-bond
	     atoms */
	  if (residue_ptr->chain == hhb_info_ptr->chain1 &&
	      !strncmp(residue_ptr->res_name,hhb_info_ptr->res_name1,3) &&
	      !strncmp(residue_ptr->res_num,hhb_info_ptr->res_num1,5) &&
	      !strncmp(atom_ptr->atom_type,hhb_info_ptr->atom_type1,4))
	    atom1_ptr = atom_ptr;

	  /* Check whether atom corresponds to the second of the H-bond
	     atoms */
	  else if (residue_ptr->chain == hhb_info_ptr->chain2 &&
	      !strncmp(residue_ptr->res_name,hhb_info_ptr->res_name2,3) &&
	      !strncmp(residue_ptr->res_num,hhb_info_ptr->res_num2,5) &&
	      !strncmp(atom_ptr->atom_type,hhb_info_ptr->atom_type2,4))
	    atom2_ptr = atom_ptr;

/* v.4.0--> */
	  /* If this is a CONECT record, then check that is between
	     ligand and non-ligand residue */
	  if (atom1_ptr != NULL && atom2_ptr != NULL &&
	      hhb_info_ptr->source == CONECT)
	    {
	      /* Check if first residue is in the ligand */
	      in_lig1 = check_if_in_ligand(hhb_info_ptr->chain1,
					   hhb_info_ptr->res_name1,
					   hhb_info_ptr->res_num1,
					   FALSE,dummy);

	      /* Check if first residue is in the ligand */
	      in_lig2 = check_if_in_ligand(hhb_info_ptr->chain2,
					   hhb_info_ptr->res_name2,
					   hhb_info_ptr->res_num2,
					   FALSE,dummy);

	      /* If neither residue is in the ligand, discount
		 this bond */
	      if (in_lig1 == FALSE && in_lig2 == FALSE &&
		  atom1_ptr->residue_ptr->object_ptr !=
		  atom2_ptr->residue_ptr->object_ptr)
		{
		  atom1_ptr = NULL;
		  atom2_ptr = NULL;
		}
	    }
/* <--v.4.0 */

	  /* If have found the two atoms corresponding to both ends
	     of the bond, then add bond to list */
	  if (atom1_ptr != NULL && atom2_ptr != NULL)
	    {
	      /* Set flag to denote that atoms matched */
	      match = TRUE;

	      /* Calculate H-bond distance */
	      x1 = atom1_ptr->x;
	      y1 = atom1_ptr->y;
	      z1 = atom1_ptr->z;
	      x2 = atom2_ptr->x;
	      y2 = atom2_ptr->y;
	      z2 = atom2_ptr->z;
	      length = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)
		+ (z1 - z2) * (z1 - z2);
	      length = sqrt((double) length);

	      /* Determine bond's source */
	      if (hhb_info_ptr->source == CONECT)
		bond_source = CONECT;
	      else
		bond_source = HBPLUS;

	      /* Add bond to list */
	      update_bond(atom1_ptr,atom2_ptr,bond_type,bond_source,length);
	    }

	   /* Get pointer to next atom in linked-list */
	   atom_ptr = atom_ptr->next;
	 }

      /* Get the pointer to the next HBPLUS interaction */
      hhb_info_ptr = hhb_info_ptr->next_hhb_info_ptr;
    }
}
/***********************************************************************

open_bonds_output  -  Open output bonds files, ligplot.hhb, ligplot.nnb
                      and ligplot.bonds

***********************************************************************/

void open_bonds_output(void)
{
  /* Open the ligplot.bonds file for output */
  if ((ligplot_bonds_out = fopen("ligplot.bonds","w")) ==  NULL)
    {
      printf("\n*** Unable to open output bonds file, ligplot.bonds\n");
      exit(1);
    }

  /* Open the ligplot.hhb file for output */
  if ((ligplot_hhb_out = fopen("ligplot.hhb","w")) ==  NULL)
    {
      printf("\n*** Unable to open output H-bonds file, ligplot.hhb\n");
      exit(1);
    }

  /* Write header records to ligplot.hhb */
  fprintf(ligplot_hhb_out,"ligplot.hhb output:\n");    
  fprintf(ligplot_hhb_out,"\n");    
  fprintf(ligplot_hhb_out,"   Donor                  Acceptor     ");
  fprintf(ligplot_hhb_out,"Distance\n");

  /* Open the ligplot.nnb file for output */
  if ((ligplot_nnb_out = fopen("ligplot.nnb","w")) ==  NULL)
    {
      printf("\n*** Unable to open output contacts file, ligplot.nnb\n");
      exit(1);
    }

  /* Write header records to ligplot.nnb */
  fprintf(ligplot_nnb_out,"ligplot.nnb output:\n");    
  fprintf(ligplot_nnb_out,"\n");    
  fprintf(ligplot_nnb_out,"   Atom 1               Atom 2        ");
  fprintf(ligplot_nnb_out,"Distance\n");
}
/***********************************************************************
      
bond_linkage  -  Determine which other bonds are linked to each bond's
                 two ends

***********************************************************************/

void bond_linkage(void)
{
  int link_to_first, link_to_second;

  struct bond *bond_ptr, *other_bond_ptr;
  struct bond_link *bond_link_ptr;

  /* Initialise pointer to start of linked list of bonds */
  bond_ptr = first_bond_ptr;
    
  /* Loop through all the bonds in the linked list, excluding any elastic
     bonds */
  while (bond_ptr != NULL)
    {
      /* Loop through all other bonds to find the ones that have an atom
	 in common with the current bond */
      other_bond_ptr = first_bond_ptr;
      
      /* Loop through all the bonds in the linked list, excluding current
         one */
      while (bond_ptr->elastic == FALSE && other_bond_ptr != NULL)
	{
	  /* Check that this isn't the current bond */
	  if (other_bond_ptr->elastic == FALSE && other_bond_ptr != bond_ptr)
	    {
	      /* Initialise flags identifying links to current bond */
	      link_to_first = FALSE;
	      link_to_second = FALSE;

	      /* Check if other bond is linked to current bond's first
	         atom */
	      if (bond_ptr->first_atom_ptr ==
		  other_bond_ptr->first_atom_ptr ||
		  bond_ptr->first_atom_ptr ==
		  other_bond_ptr->second_atom_ptr)
		link_to_first = TRUE;

	      /* Check if other bond is linked to current bond's second
	         atom */
	      if (bond_ptr->second_atom_ptr ==
		  other_bond_ptr->first_atom_ptr ||
		  bond_ptr->second_atom_ptr ==
		  other_bond_ptr->second_atom_ptr)
		link_to_second = TRUE;

	      /* If the two bonds are linked, then create a link record */
	      if (link_to_first == TRUE || link_to_second == TRUE)
		{
		  /* Allocate memory for structure to hold this bond link */
		  bond_link_ptr
		    = (struct bond_link*)malloc(sizeof(struct bond_link));
		  if (bond_link_ptr == NULL)
		    {
		      printf("*** ERROR. Unable to allocate memory for");
		      printf(" struct bond_link\n");
		      printf("***        Program ligplot terminated with");
		      printf(" error.\n");
		      exit (1);
		    }

		  /* Store the link and appropriate data */
		  bond_link_ptr->bond_ptr = other_bond_ptr;
		  if (link_to_first == TRUE)
		    {
		      bond_link_ptr->bond_end = FIRST;
		      bond_ptr->nfirst_atom_links++;
		    }
		  else
		    {
		      bond_link_ptr->bond_end = SECOND;
		      bond_ptr->nsecond_atom_links++;
		    }
		  bond_link_ptr->next_bond_link_ptr = NULL;
		  
		  /* If this is not the first link for current bond, make
		     new link point to previous link */
		  if (bond_ptr->first_bond_link_ptr != NULL)
		    bond_link_ptr->next_bond_link_ptr
		      = bond_ptr->first_bond_link_ptr;

		  /* Point from current bond to current link */
		  bond_ptr->first_bond_link_ptr = bond_link_ptr;
		}
	    }

	  /* Get pointer to next bond in linked-list */
	  other_bond_ptr = other_bond_ptr->next_bond_ptr;
	}

      /* Get pointer to next bond in linked-list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/***********************************************************************

initialise_atoms  -  Initialise all the atom flags

***********************************************************************/

void initialise_atoms(void)
{
  struct coordinate *atom_ptr;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms to initialise the flag indicating whether
     coords have already been stored */
  while (atom_ptr != NULL)
    {
      /* Initialise the flag and stack-pointer */
      atom_ptr->checked = FALSE;
      atom_ptr->next_stack_ptr = NULL;
      
      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next;
    }	       
}
/***********************************************************************

get_test_bond  -  Pick a representative bond for the given residue from
                  which to check connectivities to all atoms in the
		  residue

***********************************************************************/

struct bond *get_test_bond(struct residue *residue_ptr)
{
  int desirability, stored_desirability;

  struct bond *bond_ptr, *test_bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise variables */
  test_bond_ptr = NULL;
  stored_desirability = 0;

  /* Initialise bond pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds */
  while (bond_ptr != NULL)
    {
      /* Get the residues to which the bond's two atoms belong
	 to */
      residue1_ptr = bond_ptr->first_atom_ptr->residue_ptr;
      residue2_ptr = bond_ptr->second_atom_ptr->residue_ptr;

      /* Process only if this bond is a covalent bond entirely
	 within the current residue */
      if (bond_ptr->bond_type == COVALENT &&
	  residue1_ptr == residue_ptr &&
	  residue2_ptr == residue_ptr)
	{
	  /* Get the bond's two atoms */
	  atom1_ptr = bond_ptr->first_atom_ptr;
	  atom2_ptr = bond_ptr->second_atom_ptr;

	  /* Calculate this bond's desirability as the test
	     bond */

	  /* Best bond is N-CA bond */
	  if ((!strncmp(atom1_ptr->atom_type," N  ",4) &&
	       !strncmp(atom2_ptr->atom_type," CA ",4)) ||
	      (!strncmp(atom2_ptr->atom_type," N  ",4) &&
	       !strncmp(atom1_ptr->atom_type," CA ",4)))
	    desirability = 14;

	  /* Next-best is any mainchain bond */
	  else if (atom1_ptr->side_chain == FALSE &&
		   atom2_ptr->side_chain == FALSE)
	    {
	      /* If rotatable, then better than not */
	      if (bond_ptr->rotatable_bond == TRUE)
		desirability = 12;
	      else
		desirability = 10;
	    }

	  /* Next any bond involving a mainchain atom */
	  else if (atom1_ptr->side_chain == FALSE ||
		   atom2_ptr->side_chain == FALSE)
	    {
	      /* If rotatable, then better than not */
	      if (bond_ptr->rotatable_bond == TRUE)
		desirability = 8;
	      else
		desirability = 6;
	    }

	  /* Finally, any bond will do, preferably a
	     rotatable one */
	  else if (bond_ptr->rotatable_bond == TRUE)
	    desirability = 4;
	  else
	    desirability = 2;

	  /* If this bond is more desirable than the one we
	     already have, then store it */
	  if (desirability > stored_desirability)
	    {
	      /* Store the current bond */
	      test_bond_ptr = bond_ptr;
	      stored_desirability = desirability;
	    }
	}

      /* Get pointer for next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Return the test bond */
  return(test_bond_ptr);
}
/***********************************************************************

initialise_bonds  -  Initialise all the bond flags and stack pointers

***********************************************************************/

void initialise_bonds(void)
{
  struct bond *bond_ptr;

  /* Set pointer to first bond in linked-list */
  bond_ptr = first_bond_ptr;

  /* Loop through the stored bond-list to initialise flags and 
     stack-pointers */
  while (bond_ptr != NULL)
    {
      /* Initialise flag and stack-pointer */
      bond_ptr->checked = FALSE;
      bond_ptr->next_link_ptr = NULL;
      bond_ptr->next_stack_ptr = NULL;

      /* Go to the next bond in the linked list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/***********************************************************************

mark_downstream_atoms  -  Mark all the atoms downstream of the given bond

***********************************************************************/

void mark_downstream_atoms(struct bond *reference_bond_ptr,int end)
{
  int nbonds;

  struct bond *bond_ptr, *next_free_stack_ptr, *other_bond_ptr,
  *next_stack_ptr;
  struct bond_link *bond_link_ptr;

  /* Initialise variables */
  next_stack_ptr = NULL;
  next_free_stack_ptr = NULL;

  /* Initialise all the bond flags and stack pointers */
  initialise_bonds();

  /* Mark both this bond's atoms as already checked */
  reference_bond_ptr->first_atom_ptr->checked = TRUE;
  reference_bond_ptr->second_atom_ptr->checked = TRUE;

  /* Set pointer to the first of the bonds sprouting off the given
     bond */
  bond_link_ptr = reference_bond_ptr->first_bond_link_ptr;

  /* Set the reference-bond as already checked */
  reference_bond_ptr->checked = TRUE;

  /* Initialise count of bonds encountered */
  nbonds = 0;

  /* Initialise stack pointers for bond-search by putting the bonds
     coming off the selected end of the reference-bond onto the stack
     of bonds to be searched */
  while (bond_link_ptr != NULL)
    {
      /* Add this bond to the stack only if it comes off the relevant
	 end */
      if (bond_link_ptr->bond_end == end)
	{
	  /* Get the bond this link is pointing to and add to stack */
	  bond_ptr = bond_link_ptr->bond_ptr;

	  /* Add to stack provided it is not an elastic bond */
	  if (bond_ptr->elastic == FALSE)
	    {
	      /* If this is the first to be added to the stack, set both
		 stack-pointers to point to it */
	      if (next_free_stack_ptr == NULL)
		{
		  next_free_stack_ptr = bond_ptr;
		  next_stack_ptr = bond_ptr;
		}

	      /* Otherwise, make last bond on stack point to this one
		 and set this one as nopw the last on the stack */
	      else
		{
		  next_free_stack_ptr->next_stack_ptr = bond_ptr;
		  next_free_stack_ptr = bond_ptr;
		}
	    }

	  /* Mark bond as added to the stack */
	  bond_ptr->checked = TRUE;
	}
      
      /* Go to the next bond-link in the linked list */
      bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
    }

  /* Loop until stack of pointers to be processed has been exhausted */
  while (next_stack_ptr != NULL)
    {
      /* Pick up the next bond from the stack, and move the stack
       pointer onto the next position */
      bond_ptr = next_stack_ptr;
      next_stack_ptr = bond_ptr->next_stack_ptr;

      /* Increment count of bonds encountered */
      nbonds++;

      /* Mark the bond's two atoms as checked */
      bond_ptr->first_atom_ptr->checked = TRUE;
      bond_ptr->second_atom_ptr->checked = TRUE;

     /* Loop through all the bonds connected to the current one
        to add to stack */

      /* Get the pointer to the first of the bonds coming off the 
	 current bond */
      bond_link_ptr = bond_ptr->first_bond_link_ptr;

      /* Process each of these bonds in turn */
      while (bond_link_ptr != NULL)
	{
	  /* Get the pointer to this bond */
	  other_bond_ptr = bond_link_ptr->bond_ptr;

	  /* If this bond hasn't already been processed, then move
	     its atoms (if they need moving) and add the bond to the
	     stack */
	  if (other_bond_ptr->checked == FALSE &&
	      other_bond_ptr->elastic == FALSE)
	    {
	      next_free_stack_ptr->next_stack_ptr = other_bond_ptr;
	      next_free_stack_ptr = other_bond_ptr;

	      /* If the stack has run out, restart it at the bond
		 just entered */
	      if (next_stack_ptr == NULL)
		next_stack_ptr = other_bond_ptr;

	      /* Mark bond as added to the stack */
	      other_bond_ptr->checked = TRUE;
	    }

	  /* Go to the next bond-link in the linked list */
	  bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
	}
    }
}
/***********************************************************************

remove_unconnected atoms  -  Find, and mark for deletion, all atoms in
                             the current residue that cannot be reached
			     from the given test-bond

***********************************************************************/

void remove_unconnected_atoms(struct bond *bond_ptr,
			      struct residue *residue_ptr,int *marked)
{
  int iatom, natoms;

  struct coordinate *atom_ptr;

  /* Initialise all atoms */
  initialise_atoms();

  /* Mark all atoms downstream of each end of the given bonds */
  mark_downstream_atoms(bond_ptr,FIRST);
  mark_downstream_atoms(bond_ptr,SECOND);

  /* Get pointer to the first of this residue's atoms */
  atom_ptr = residue_ptr->first_atom_ptr;

  /* Initialise atom count */
  iatom = 0;
  natoms = residue_ptr->natoms;

  /* Loop over all the residue's atoms checking which atoms in the
     current residue haven't been marked */
  while (iatom < natoms && atom_ptr != NULL)
    {
      /* Check which atoms in the current residue haven't been marked */
      if (atom_ptr->residue_ptr == residue_ptr &&
	  atom_ptr->checked == FALSE)
	{
	  /* Have an unconnected atom, so mark it for deletion */
	  atom_ptr->deleted = TRUE;
	  *marked = TRUE;
	}

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next;
    }
}
/***********************************************************************

delete_unwanted_atoms  -  Delete any atoms marked for deletion

***********************************************************************/

void delete_unwanted_atoms(void)
{
  int iatom, natoms;

  struct coordinate *atom_ptr, *last_atom_ptr;
  struct residue *residue_ptr;

  /* Get pointer to the first residue */
  residue_ptr = first_residue_ptr;
  last_atom_ptr = NULL;

  /* Loop over all the residues */
  while (residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
	      
      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, deleting any unwanted
         ones */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* If atom is marked for deletion, then make pointer of previous
	     atom by-pass it */
	  if (atom_ptr->deleted == TRUE)
	    {
	      /* If this is the first atom, then set pointer to
		 next one to bypass it */
	      if (atom_ptr == first_atom_ptr)
		first_atom_ptr = atom_ptr->next;

	      /* Otherwise, set last pointer around current atom */
	      else
		last_atom_ptr->next = atom_ptr->next;

	      /* If this is the first atom of the current residue,
		 then alter the residue's pointer */
	      if (atom_ptr == residue_ptr->first_atom_ptr)
		residue_ptr->first_atom_ptr = atom_ptr->next;

	      /* Reduce count of this residue's atoms */
	      residue_ptr->natoms--;
	    }

	  /* Otherwise, store pointer to this as the last undeleted
	     atom encountered */
	  else
	    last_atom_ptr = atom_ptr;
      
	  /* Get pointer to next atom in linked-list */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }
}
/***********************************************************************

delete_unwanted_bonds  -  Delete any bonds involving deleted atoms

***********************************************************************/

void delete_unwanted_bonds(void)
{
  struct bond *bond_ptr, *last_bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;

  /* Initialise */
  last_bond_ptr = NULL;

  /* Get pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop over all the bonds */
  while (bond_ptr != NULL)
    {
      /* Get the pointers to this bond's two atoms */
      atom1_ptr = bond_ptr->first_atom_ptr;
      atom2_ptr = bond_ptr->second_atom_ptr;

      /* If either atom has been deleted, then want to delete this bond */
      if (atom1_ptr->deleted == TRUE || atom2_ptr->deleted == TRUE ||
	   bond_ptr->bond_type == DELETED)
	{
	  /* Mark the bond as deleted */
	  bond_ptr->bond_type = DELETED;

	  /* Remove the bond from the linked list */
	 
	  /* If this is the first bond, then set pointer to
	     next one to bypass it */
	  if (bond_ptr == first_bond_ptr)
	    first_bond_ptr = bond_ptr->next_bond_ptr;

	  /* Otherwise, set last pointer around current bond */
	  else
	    last_bond_ptr->next_bond_ptr = bond_ptr->next_bond_ptr;
	}

      /* Otherwise, store pointer to this as the last undeleted
	 bond encountered */
      else
	last_bond_ptr = bond_ptr;

      /* Get pointer to the next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/***********************************************************************

delete_unwanted_bond_links  -  Delete any unwanted bond links

***********************************************************************/

void delete_unwanted_bond_links(void)
{
  struct bond *bond_ptr;
  struct object *object_ptr;
  struct object_bond *object_bond_ptr, *last_object_bond_ptr;

  /* Initialise */
  last_object_bond_ptr = NULL;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Get pointer to the first of this object's bonds */
      object_bond_ptr = object_ptr->first_object_bond_ptr;
      
      /* Loop through all this object's bonds, one by one */
      while (object_bond_ptr != NULL)
	{
	  /* Get the current bond */
	  bond_ptr = object_bond_ptr->bond_ptr;

	  /* If this bond has been deleted, then need to delete
	     this link also */
	  if (bond_ptr->bond_type == DELETED)
	    {
	      /* If this is the first link, then set pointer to
		 next one to bypass it */
	      if (object_bond_ptr == object_ptr->first_object_bond_ptr)
		object_ptr->first_object_bond_ptr
		  = object_bond_ptr->next_object_bond_ptr;

	      /* Otherwise, set last pointer around current bond */
	      else
		last_object_bond_ptr->next_object_bond_ptr
		  = object_bond_ptr->next_object_bond_ptr;
	    }

	  /* Otherwise, store pointer to this as the last undeleted
	     link encountered */
	  else
	    last_object_bond_ptr = object_bond_ptr;

	  /* Get pointer to the next link for this object */
	  object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
	}

      /* Get pointer to the next object */
      object_ptr = object_ptr->next_object_ptr;
    }
}
/***********************************************************************

delete_unwanted_bond_connections  -  Delete any unwanted bond-to-bond
                                     connections used in chasing downstream
				     atoms

***********************************************************************/

void delete_unwanted_bond_connections(void)
{
  struct bond *bond_ptr, *other_bond_ptr;
  struct bond_link *bond_link_ptr, *last_bond_link_ptr;

  /* Initialise */
  last_bond_link_ptr = NULL;

  /* Initialise pointer to the first stored bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all bonds */
  while (bond_ptr != NULL)
    {
      /* Get pointer to the first of the bonds springing off this one */
      bond_link_ptr = bond_ptr->first_bond_link_ptr;

      /* Loop through all this bond's adjacent bonds, one by one */
      while (bond_link_ptr != NULL)
	{
	  /* Get the connected bond */
	  other_bond_ptr = bond_link_ptr->bond_ptr;

	  /* If this bond has been deleted, or has been made into
	     an elastic bond, then need to delete this link */
	  if (other_bond_ptr->bond_type == DELETED ||
	      other_bond_ptr->elastic == TRUE)
	    {
	      /* If this is the first link, then set pointer to
		 next one to bypass it */
	      if (bond_link_ptr == bond_ptr->first_bond_link_ptr)
		bond_ptr->first_bond_link_ptr
		  = bond_link_ptr->next_bond_link_ptr;

	      /* Otherwise, set last pointer around current bond */
	      else
		last_bond_link_ptr->next_bond_link_ptr
		  = bond_link_ptr->next_bond_link_ptr;
	    }

	  /* Otherwise, store pointer to this as the last undeleted
	     link encountered */
	  else
	    last_bond_link_ptr = bond_link_ptr;

	  /* Get pointer to the next link for this bond */
	  bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
	}

      /* Get pointer to the next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/***********************************************************************

delete_atoms_and_bonds  -  Delete any atoms marked for deletion and any
                           bonds and bond-links that involve them

***********************************************************************/

void delete_atoms_and_bonds(void)
{
  /* Delete any unwanted atoms, previously marked for deletion */
  delete_unwanted_atoms();

  /* Delete any unwanted bonds */
  delete_unwanted_bonds();

  /* Delete any unwanted links */
  delete_unwanted_bond_links();

  /* Delete any unwanted bond-to-bond connections */
  delete_unwanted_bond_connections();
}
/***********************************************************************

check_connections  -  Check that all atoms within each residue are
                      connected

***********************************************************************/

void check_connections(void)
{
  int iresid, nresid;
  int iobject;
  int marked;

  struct bond *test_bond_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  marked = FALSE;

  /* Initialise pointer to the first stored object */
  printf("\nChecking for unconnected atoms ...\n");
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects to find a representative bond */
  while (object_ptr != NULL)
    {
      /* Skip hydrophobic groups, waters and simple H-bonded groups */
      if (object_ptr->object_type != HYDROPHOBIC &&
	  object_ptr->object_type != SIMPLE_HGROUP &&
	  object_ptr->object_type != WATER)
	{
	  /* Set the pointer to the first of this object's residues */
	  residue_ptr = object_ptr->first_residue_ptr;
	  nresid = object_ptr->nresidues;
	  iresid = 0;

	  /* Loop over all this object's residues */
	  while (iresid < nresid && residue_ptr != NULL)
	    {
	      /* Get a test bond for this residue */
	      test_bond_ptr = get_test_bond(residue_ptr);

	      /* If have a valid test bond then check that all atoms in
		 the current residue are reachable from it */
	      if (test_bond_ptr != NULL)
		{
		  /* Mark any unwanted atoms for deletion */
		  remove_unconnected_atoms(test_bond_ptr,residue_ptr,
					   &marked);
		}

	      /* Get pointer to the next residue */
	      residue_ptr = residue_ptr->next_residue_ptr;
	      iresid++;
	    }
	}

      /* Get pointer for next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

  /* If any atoms have been marked for deletion, then remove them and
     their associated bonds, etc */
  if (marked == TRUE)
    {
      /* Delete any unwanted atoms */
      delete_atoms_and_bonds();
    }
}
/***********************************************************************

check_for_attachment  -  Check whether the given residue is an attachment
                         to a previous residue in the given object

***********************************************************************/

int check_for_attachment(struct object *object_ptr,
			 struct residue *current_residue_ptr)
{
  int iresid, nresid;
  int connected;

  struct bond *bond_ptr;
  struct residue *residue_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise flag */
  connected = FALSE;

  /* Get number of residues making up this object */
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Get pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;

  /* Loop through all the residues until hit the current one */
  while (iresid < nresid && residue_ptr != NULL &&
	 residue_ptr != current_residue_ptr)
    {
      /* Get pointer to first bond in the list */
      bond_ptr = first_bond_ptr;

      /* Loop through all the bonds in the linked list */
      while (bond_ptr != NULL && connected == FALSE)
	{
	  /* If this is a covalent bond, check whether
	     it links the two residues */
/* v.3.1.1--> */
/*	  if (bond_ptr->bond_type == COVALENT) */
	  if (bond_ptr->bond_type == COVALENT &&
	      bond_ptr->elastic == FALSE)
/* <--v.3.1.1 */
	    {
	      /* Get this bond's two residues */
	      residue1_ptr = bond_ptr->first_atom_ptr->residue_ptr;
	      residue2_ptr = bond_ptr->second_atom_ptr->residue_ptr;

	      /* If it spans the two residues in question,
		 then have a link between them */
	      if ((residue1_ptr == residue_ptr &&
		   residue2_ptr == current_residue_ptr) ||
		  (residue2_ptr == residue_ptr &&
		   residue1_ptr == current_residue_ptr))
		connected = TRUE;
	    }

	  /* Get pointer to next atom in linked-list */
	  bond_ptr = bond_ptr->next_bond_ptr;
	}

      /* Get pointer to the next residue in the list */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Return flag indicating whether the given residue is an attachment */
  return(connected);
}
/* v.3.1.1--> */
/***********************************************************************

check_for_cyclics  -  Check whether any objects are cyclic (ie residue
                      m is covalent bound to residue n where m and n are
                      not adjacent)

***********************************************************************/

void check_for_cyclics(void)
{
  int iobject;
  int cyclic;

  struct coordinate *atom1_ptr, *atom2_ptr;
  struct bond *bond_ptr;
  struct object *object_ptr, *object1_ptr, *object2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Get pointer to the first stored bond */
      bond_ptr = first_bond_ptr;

      /* Loop through all the bonds, identifying any belonging
	 to the current object */
      while (bond_ptr != NULL)
	{
	  /* Get the two atoms on either end of this bond */
	  atom1_ptr = bond_ptr->first_atom_ptr;
	  atom2_ptr = bond_ptr->second_atom_ptr;

	  /* Get the objects containing these atoms via their respective
	     residues */
	  object1_ptr = atom1_ptr->residue_ptr->object_ptr;
	  object2_ptr = atom2_ptr->residue_ptr->object_ptr;
	  
	  /* If both atoms belong to the current object, then check
	     which residues the atoms belong to */
	  if (object1_ptr == object_ptr && object2_ptr == object_ptr)
	    {
	      /* Get which residues the two atoms belong to */
	      residue1_ptr = atom1_ptr->residue_ptr;
	      residue2_ptr = atom2_ptr->residue_ptr;

	      /* If the atoms belong to two different residues, then
		 check for cyclic peptides */
	      if (bond_ptr->bond_type == COVALENT &&
		  bond_ptr->elastic == FALSE &&
		  residue1_ptr != residue2_ptr)
		{
		  cyclic = FALSE;

		  /* If residues are not adjacent, then looks like a
		     cyclic peptide */
		  if (residue1_ptr->next_residue_ptr != residue2_ptr &&
		      residue2_ptr->next_residue_ptr != residue1_ptr)
		    {
		      cyclic = TRUE;

		      /* If they are the same residue (ie one is an
			 attachment to the other) then not cyclic */
		      if (!strncmp(residue1_ptr->res_num,
				  residue2_ptr->res_num,5) &&
			  residue1_ptr->chain == residue2_ptr->chain)
			cyclic = FALSE;
		    }

		  /* If have what appears to be a cyclic peptide, then
		     need to make the bond elastic */
		  if (cyclic == TRUE)
		    {
		      /* Print warning message and make the bond elastic */
		      printf("*** Warning. Seems to be a cyclic peptide\n");
		      printf("***          Bond made elastic: %s %s %s %c",
			     atom1_ptr->atom_type,residue1_ptr->res_name,
			     residue1_ptr->res_num,residue1_ptr->chain);
		      printf(" ->  %s %s %s %c\n",
			     atom2_ptr->atom_type,residue2_ptr->res_name,
			     residue2_ptr->res_num,residue2_ptr->chain);
		      bond_ptr->elastic = TRUE;
		      Nwarnings++;
		    }
		}
	    }

	  /* Get pointer to point to next bond in the list */
	  bond_ptr = bond_ptr->next_bond_ptr;
	}

      /* Get pointer to point to next object in the list and increment
	 object count */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }
}
/***********************************************************************

split_objects  -  Check for breaks between adjacent residues in objects
                  having more than one residue, and split into two, if
                  necessary

***********************************************************************/

void split_objects(int *nobjects,int *nligands)
{
  int iobject, iresid, nresid;
  int connected;

  struct bond *bond_ptr;
  struct object *last_object_ptr, *object_ptr, *new_object_ptr;
  struct residue *residue_ptr, *next_residue_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise variables */
  *nligands = 0;

  /* Loop through all the objects to get the last one in case need
     to create any new ones in the routine */
  last_object_ptr = NULL;
  object_ptr = first_object_ptr;
  while (object_ptr != NULL)
    {
      /* Save this object's pointer as the last encountered */
      last_object_ptr = object_ptr;

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;
  
  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Get number of residues making up this object */
      nresid = object_ptr->nresidues;
      iresid = 0;

      /* If this is a ligand, increment count */
      if (object_ptr->object_type == LIGAND)
	(*nligands)++;

      /* Get pointer to the first of this object's residues */
      residue_ptr = object_ptr->first_residue_ptr;

      /* Loop through all pairs of adjacent residues */
      while (iresid < nresid && residue_ptr != NULL)
	{
	  /* Update residue's pointer to its parent object */
	  residue_ptr->object_ptr = object_ptr;

	  /* Get pointer to the next residue */
	  next_residue_ptr = residue_ptr->next_residue_ptr;

	  /* If next residue exists, then check for covalent bonds
	     common to both residues */
	  if (next_residue_ptr != NULL && iresid < nresid - 1 &&
	      object_ptr->object_type == LIGAND)
	    {
	      /* Update residue's pointer to its parent object */
	      next_residue_ptr->object_ptr = object_ptr;

	      /* Initialise flag */
	      connected = FALSE;

	      /* Get pointer to first bond in the list */
	      bond_ptr = first_bond_ptr;

	      /* Loop through all the bonds in the linked list */
	      while (bond_ptr != NULL && connected == FALSE)
		{
		  /* If this is a covalent bond, check whether
		     it links the two residues */
/* v.3.1.1--> */
/*		  if (bond_ptr->bond_type == COVALENT)  */
		  if (bond_ptr->bond_type == COVALENT &&
		      bond_ptr->elastic == FALSE)
/* <--v.3.1.1 */
		    {
		      /* Get this bond's two residues */
		      residue1_ptr
			= bond_ptr->first_atom_ptr->residue_ptr;
		      residue2_ptr
			= bond_ptr->second_atom_ptr->residue_ptr;

		      /* If it spans the two residues in question,
			 then have a link between them */
		      if ((residue1_ptr == residue_ptr &&
			   residue2_ptr == next_residue_ptr) ||
			  (residue2_ptr == residue_ptr &&
			   residue1_ptr == next_residue_ptr))
			connected = TRUE;
		    }

		  /* Get pointer to next atom in linked-list */
		  bond_ptr = bond_ptr->next_bond_ptr;
		}

	      /* If have no link between the two residues, then may
		 be an attachment */
	      if (connected == FALSE)
		{
		  /* Check whether the residue is an attachment
		     to a previous residue in this object */
		  connected = check_for_attachment(object_ptr,
						   next_residue_ptr);
		}

	      /* If either residue is a heme, then split anyway */
	      if (!strncmp(residue_ptr->res_name,"HEM",3) ||
		  !strncmp(next_residue_ptr->res_name,"HEM",3))
		connected = FALSE;

/* v.4.0--> */
	      /* If plotting interactions across an interface, then want
		 all residues to be separate objects */
	      if (Interface_Plot == TRUE)
		connected = FALSE;
/* <--v.4.0 */

	      /* If still no link between the two residues, then may have
		 to split object into two */
	      if (connected == FALSE)
		{
/* v.4.0--> */
		  if (Interface_Plot == FALSE)
		    {
/* <--v.4.0 */
		      printf("\n");
		      printf("*** Splitting ligand into two between");
		      printf("  %s %s %c and %s %s %c\n",
			     residue_ptr->res_name,
			     residue_ptr->res_num,
			     residue_ptr->chain,
			     next_residue_ptr->res_name,
			     next_residue_ptr->res_num,
			     next_residue_ptr->chain);
/* v.3.1.2--> */
		      printf("\n");
/* <--v.3.1.2 */
		      Nwarnings++;
/* v.4.0--> */
		    }
/* <--v.4.0 */

		  /* Allocate memory for structure to hold new object's
		     details */
		  new_object_ptr = (struct object*)
		    malloc(sizeof(struct object));
		  if (new_object_ptr == NULL)
		    {
		      printf("*** ERROR. Unable to allocate memory for");
		      printf(" struct object\n");
		      printf("***        Program ligplot terminated with");
		      printf(" error.\n");
		      exit (1);
		    }

		  /* Fill in initial details for the new object */
		  new_object_ptr->object_type = object_ptr->object_type;
		  new_object_ptr->nresidues = nresid - iresid - 1;
		  new_object_ptr->first_residue_ptr = next_residue_ptr;
		  new_object_ptr->first_object_bond_ptr = NULL;
		  new_object_ptr->next_object_ptr = NULL;

		  /* Initialise remaining details (which will be
		     filled in later) */
/* v.3.2--> */
		  new_object_ptr->nbonds = 0.0;
/* <--v.3.2 */
		  new_object_ptr->nrot_bonds = 0.0;
		  new_object_ptr->nmain_chain = 0.0;
		  new_object_ptr->n_ca = 0;
		  new_object_ptr->weight = 0.0;
		  new_object_ptr->minx = 0.0;
		  new_object_ptr->miny = 0.0;
		  new_object_ptr->maxx = 0.0;
		  new_object_ptr->maxy = 0.0;
/* v.4.0--> */
		  new_object_ptr->interface = 0;
/* <--v.4.0 */
		  new_object_ptr->max_atom_size = 0.0;
		  new_object_ptr->place_x = 0.0;
		  new_object_ptr->place_y = 0.0;
		  new_object_ptr->internal_energy = 0.0;
		  new_object_ptr->total_energy = 0.0;

		  /* Update number of residues now belonging to
		     the object that has been split */
		  object_ptr->nresidues = iresid + 1;

		  /* Add link from last object to the new one just
		     created */
		  last_object_ptr->next_object_ptr = new_object_ptr;
		  last_object_ptr = new_object_ptr;

		  /* Increment count of objects */
		  (*nobjects)++;

		  /* Finish with the current object */
		  iresid = nresid + 1;
		}
	    }

	  /* Get pointer to the next residue in the list */
	  residue_ptr = residue_ptr->next_residue_ptr;
	  iresid++;
	}

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }
}
/***********************************************************************

delete_object  -  Delete the given object, marking all its atoms for
                  deletion

***********************************************************************/

void delete_object(struct object *object_ptr,
/* v.3.2--> */
/*		   struct object *previous_object_ptr) */
		   struct object *previous_object_ptr,int show_deletions)
/* <--v.3.2 */
{
  int iatom, iresid, natoms, nresid;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Show message that object to be deleted */
/* v.3.2--> */
  if (show_deletions == TRUE)
    {
/* <--v.3.2 */
      printf("+++ Residue %s %s %c not connected to any other",
	     object_ptr->first_residue_ptr->res_name,
	     object_ptr->first_residue_ptr->res_num,
	     object_ptr->first_residue_ptr->chain);
      printf("  ... deleting it!\n");
      Nwarnings++;
/* v.3.2--> */
    }
/* <--v.3.2 */

  /* Mark all the object's atoms for deletion */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop through all the object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
		  
      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, deleting any unwanted
	 ones */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* Mark the atom for deletion */
	  atom_ptr->deleted = TRUE;

	  /* Get pointer to next atom in linked-list */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

/* v.3.2--> */
      /* Mark the residue itself as deleted */
      residue_ptr->deleted = TRUE;
/* <--v.3.2 */

      /* Get pointer to the next residue in the list */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Delete the current object by making the object pointers
     bypass it */
  if (previous_object_ptr == NULL)
    first_object_ptr = object_ptr->next_object_ptr;
	      
  /* Otherwise, make pointers bypass the current link */
  else
    previous_object_ptr->next_object_ptr = object_ptr->next_object_ptr;
}
/***********************************************************************

assign_bonds_to_objects  -  Determine which bonds relate to each object

***********************************************************************/

/* <--v.3.2 */
/* void assign_bonds_to_objects(int *ndeleted) */
void assign_bonds_to_objects(int *ndeleted,int show_deletions)
/* v.3.2--> */
{
  int iobject, nbonds;
  int ncovalent, nhbonds, ncontacts;
  int cyclic;

  struct coordinate *atom1_ptr, *atom2_ptr;
  struct bond *bond_ptr;
  struct object *object_ptr, *object1_ptr, *object2_ptr, *next_object_ptr,
  *previous_object_ptr;
  struct object_bond *object_bond_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  previous_object_ptr = NULL;
  iobject = 0;
  *ndeleted = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Initialise counts of bonds between this object and all others */
      nbonds = 0;
      ncovalent = 0;
      nhbonds = 0;
      ncontacts = 0;

      /* Save pointer to the next object */
      next_object_ptr = object_ptr->next_object_ptr;

      /* Get pointer to the first stored bond */
      bond_ptr = first_bond_ptr;

      /* Loop through all the bonds, identifying any belonging
	 to the current object */
      while (bond_ptr != NULL)
	{
	  /* Get the two atoms on either end of this bond */
	  atom1_ptr = bond_ptr->first_atom_ptr;
	  atom2_ptr = bond_ptr->second_atom_ptr;

	  /* Get the objects containing these atoms via their respective
	     residues */
	  object1_ptr = atom1_ptr->residue_ptr->object_ptr;
	  object2_ptr = atom2_ptr->residue_ptr->object_ptr;
	  
	  /* If both atoms belong to the current object, then assign
	     current bond, if covalent, to current object */
	  if (object1_ptr == object_ptr && object2_ptr == object_ptr)
	    {
	      /* Get which residues the two atoms belong to */
	      residue1_ptr = atom1_ptr->residue_ptr;
	      residue2_ptr = atom2_ptr->residue_ptr;

	      /* If the atoms belong to two different residues, then
		 check for cyclic peptides */
	      if (bond_ptr->bond_type == COVALENT &&
		  residue1_ptr != residue2_ptr)
		{
		  cyclic = FALSE;

		  /* If residues are not adjacent, then looks like a
		     cyclic peptide */
		  if (residue1_ptr->next_residue_ptr != residue2_ptr &&
		      residue2_ptr->next_residue_ptr != residue1_ptr)
		    {
		      cyclic = TRUE;

		      /* If they are the same residue (ie one is an
			 attachment to the other) then not cyclic */
		      if (!strncmp(residue1_ptr->res_num,
				  residue2_ptr->res_num,5) &&
			  residue1_ptr->chain == residue2_ptr->chain)
			cyclic = FALSE;
		    }

		  /* If have what appears to be a cyclic peptide, then
		     need to make the bond elastic */
		  if (cyclic == TRUE)
		    {
		      /* Print warning message and make the bond elastic */
		      printf("*** Warning. Seems to be a cyclic peptide\n");
		      printf("***          Bond made elastic: %s %s %s %c",
			     atom1_ptr->atom_type,residue1_ptr->res_name,
			     residue1_ptr->res_num,residue1_ptr->chain);
		      printf(" ->  %s %s %s %c\n",
			     atom2_ptr->atom_type,residue2_ptr->res_name,
			     residue2_ptr->res_num,residue2_ptr->chain);
		      bond_ptr->elastic = TRUE;
		      Nwarnings++;
		    }
		}

	      /* If bond is a covalent bond then perform assignment */
	      if (bond_ptr->bond_type == COVALENT &&
		  bond_ptr->elastic == FALSE)
		{
		  /* Increment count of covalent bonds belonging to
		     the current object */
		  object_ptr->nbonds++;

		  /* Allocate memory for structure to hold this
		     object-bond link */
		  object_bond_ptr
		    = (struct object_bond*)malloc(sizeof(struct object_bond));
		  if (object_bond_ptr == NULL)
		    {
		      printf("*** ERROR. Unable to allocate memory for");
		      printf(" struct object_bond\n");
		      printf("***        Program ligplot terminated with");
		      printf(" error.\n");
		      exit (1);
		    }

		  /* Store the link and appropriate data */
		  object_bond_ptr->bond_ptr = bond_ptr;
		  object_bond_ptr->next_object_bond_ptr = NULL;
		  
		  /* If this is not the first link for current bond, make
		     new link point to previous link */
		  if (object_ptr->first_object_bond_ptr != NULL)
		    object_bond_ptr->next_object_bond_ptr
		      = object_ptr->first_object_bond_ptr;
	      
		  /* Point from current bond to current link */
		  object_ptr->first_object_bond_ptr = object_bond_ptr;
		}

	      /* If this is a hydrogen bond within the object, mark it
		 as internal */
	      else if (bond_ptr->bond_type == HBOND)
		bond_ptr->bond_type = INTERNAL;

	      /* Otherwise, if it is an internal non-bonded contact
		 interaction, then delete it altogether */
	      else if (bond_ptr->bond_type != COVALENT)
		bond_ptr->bond_type = DELETED;
	    }

	  /* Otherwise, if only one of the ends of the bond belongs to
	     this object, need to make the bond an elastic one */
	  else if (object1_ptr == object_ptr || object2_ptr == object_ptr)
	    {
	      bond_ptr->elastic = TRUE;

	      /* Increment counts of different bond types between this
		 object and any others */
	      nbonds++;
	      if (bond_ptr->bond_type == HBOND)
/* v.3.1.1--> */
/*		nhbonds++; */
		{
		  nhbonds++;

		  /* If both objects are ligand objects, then make this an
		     internal hydrogen bond */
		  if (object1_ptr->object_type == LIGAND &&
		      object2_ptr->object_type == LIGAND)
		    bond_ptr->bond_type = INTERNAL;
		}
/* <--v.3.1.1 */
	      else if (bond_ptr->bond_type == CONTACT)
		ncontacts++;
	      else if (bond_ptr->bond_type == COVALENT)
/* v.3.1--> */
/*		ncovalent++; */
		{
		  /* If this covalent bond is between the current object
		     and the ligand, then increment count of covalent
		     bonds */
		  if ((object1_ptr == object_ptr &&
		       object2_ptr->object_type == LIGAND) ||
		      (object2_ptr == object_ptr &&
		       object1_ptr->object_type == LIGAND))
		    ncovalent++;
		}
/* <--v.3.1 */
	    }

	  /* Get pointer to point to next bond in the list */
	  bond_ptr = bond_ptr->next_bond_ptr;
	}

      /* If this is not the ligand object, then determine what type it
	 is */
/* v.4.0--> */
/*      if (object_ptr->object_type != LIGAND) */
      if (object_ptr->object_type != LIGAND || Interface_Plot == TRUE)
/* <--v.4.0 */
	{
	  /* If there are no bonds whatsoever between this object and
	     any others, then might as well delete it */
	  if (nbonds == 0)
	    {
	      /* Delete the current object */
/* v.3.2--> */
/*	      delete_object(object_ptr,previous_object_ptr); */
/* v.4.0--> */
	      if (Print_as_is == FALSE)
		{
/* <--v.4.0 */
		  delete_object(object_ptr,previous_object_ptr,show_deletions);
/* <--v.3.2 */
		  (*ndeleted)++;
/* v.4.0--> */
		}
/* <--v.4.0 */

	      /* Set current pointer to the previous object */
	      object_ptr = previous_object_ptr;
	    }

	  /* Otherwise, if no hydrogen bonds or covalent bonds, then
	     mark this object as being a hydrophobic group */
	  else if (ncontacts > 0 && nhbonds == 0 && ncovalent == 0)
/* v.4.0--> */
	    {
	      /* Make object hydrophobic, or set interface flag */
	      if (Interface_Plot == FALSE)
/* <--v.4.0 */
		object_ptr->object_type = HYDROPHOBIC;
/* v.4.0--> */
	      else
		object_ptr->interface = 9;
	    }
/* <--v.4.0 */
	}

      /* Save the current object pointer */
      previous_object_ptr = object_ptr;

      /* Get pointer to point to next object in the list and increment
	 object count */
      object_ptr = next_object_ptr;
      iobject++;
    }

/* v.4.0--> */
  /* For interface plots, loop through the object again setting the
     appropriate ones as hydrophobic objects */
  if (Interface_Plot == TRUE)
    {
      /* Reinitialise to first object */
      object_ptr = first_object_ptr;

      /* Loop through all objects */
      while (object_ptr != NULL)
	{
	  /* If this should be a hydrophobic object, then make it so */
	  if (object_ptr->interface == 9)
	    {
	      /* Mark the appropriate interface */
	      if (object_ptr->object_type == LIGAND)
		object_ptr->interface = 1;
	      else
		object_ptr->interface = 2;

	      object_ptr->object_type = HYDROPHOBIC;
	    }

	  /* Otherwise set interface type according to the object type */
	  else
	    {
	      /* Mark the appropriate interface */
	      if (object_ptr->object_type == LIGAND)
		object_ptr->interface = 1;
	      else if (object_ptr->object_type != WATER)
		object_ptr->interface = 2;
	    }

	  /* Get pointer to point to next object in the list */
	  object_ptr = object_ptr->next_object_ptr;
	}
    }
/* <--v.4.0 */
}
/***********************************************************************

show_deleted_residues_message  -  Show explanatory note about residues
                                  deleted

***********************************************************************/

void show_deleted_residues_message(int ndeleted)
{
  /* Print the message */
  printf("\n");
  if (ndeleted == 1)
    {
      printf("NOTE. The above residue has been deleted as it is not ");
      printf("connected to \n");
      printf("any otherin the diagram. It may have originally been ");
      printf("read in because it \n");
      printf("appears in the CONECT records of the PDB file.\n");
    }
  else
    {
      printf("NOTE. The above residues have been deleted as they ");
      printf("are not connected to \n");
      printf("any others in the diagram. They may have originally ");
      printf("been read in because\n");
      printf("they appear in the CONECT records of the PDB file.\n");
    }
}
/***********************************************************************

mark_reachable_atoms  -  Mark all atoms that are reachable from the
                         ligand residue(s) via covalent bonds, H-bonds
                         or hydrophobic contacts

***********************************************************************/

void mark_reachable_atoms(void)
{
  int iatom, iresid, natoms, nresid;
  int checked;
/* v.3.1--> */
  int wanted;
/* <--v.3.1 */

  struct coordinate *atom_ptr, *atom1_ptr, *atom2_ptr;
  struct bond *bond_ptr;
/* v.3.1--> */
/*  struct object *object_ptr; */
  struct object *object_ptr, *object1_ptr, *object2_ptr;
/* <--v.3.1 */
  struct residue *residue_ptr;

  /* Inialise all atoms */
  initialise_atoms();

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects to mark the ligand atoms as checked */
  while (object_ptr != NULL)
    {
      /* If the object is a ligand, mark all its atoms */
      if (object_ptr->object_type == LIGAND)
	{
	  /* Get pointer to the first of this object's residues */
	  residue_ptr = object_ptr->first_residue_ptr;
	  nresid = object_ptr->nresidues;
	  iresid = 0;

	  /* Loop over all this object's residues */
	  while (iresid < nresid && residue_ptr != NULL)
	    {
	      /* Get pointer to this residue's first atom */
	      atom_ptr = residue_ptr->first_atom_ptr;
	      
	      /* Initialise atom count */
	      iatom = 0;
	      natoms = residue_ptr->natoms;

	      /* Loop over all the residue's atoms, checking whether
		 they are involved in interactions */
	      while (iatom < natoms && atom_ptr != NULL)
		{
		  /* Mark this atom as wanted */
		  atom_ptr->checked = TRUE;

		  /* Get pointer to the next atom */
		  atom_ptr = atom_ptr->next;
		  iatom++;
		}

	      /* Get pointer to the next residue */
	      residue_ptr = residue_ptr->next_residue_ptr;
	      iresid++;
	    }
	}

      /* Get pointer for next object */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* Initialise flag */
  checked = TRUE;

  /* Loop until no new atoms checked as being connected */
  while (checked == TRUE)
    {
      /* reset flag */
      checked = FALSE;

      /* Get pointer to the first stored bond */
      bond_ptr = first_bond_ptr;

      /* Loop through all the bonds, marking any unchecked atoms 
	 connected to checked atoms */
      while (bond_ptr != NULL)
	{
	  /* Get the two atoms on either end of this bond */
	  atom1_ptr = bond_ptr->first_atom_ptr;
	  atom2_ptr = bond_ptr->second_atom_ptr;

/* v.3.1--> */
	  /* Get the two objects these atoms belong to */
	  object1_ptr = atom1_ptr->residue_ptr->object_ptr;
	  object2_ptr = atom2_ptr->residue_ptr->object_ptr;

	  /* If the link is to a HYDROPHOBIC group, then ignore this
	     bond */
	  wanted = TRUE;
	  if ((object1_ptr->object_type == HYDROPHOBIC &&
	       object2_ptr->object_type != LIGAND) ||
	      (object2_ptr->object_type == HYDROPHOBIC &&
	       object1_ptr->object_type != LIGAND))
	    wanted = FALSE;
/* <--v.3.1 */

	  /* Process only if this is a valid bond */
/* v.3.1--> */
/*	  if (bond_ptr->bond_type == COVALENT ||
	      bond_ptr->bond_type == HBOND || 
	      bond_ptr->bond_type == CONTACT) */
	  if (wanted == TRUE && (bond_ptr->bond_type == COVALENT ||
	      bond_ptr->bond_type == HBOND || 
	      bond_ptr->bond_type == CONTACT))
/* <--v.3.1 */
	    {
	      /* If one is checked and other not, then mark the one that
		 isn't as reachable */
	      if (atom1_ptr->checked == TRUE &&
		  atom2_ptr->checked == FALSE)
		{
		  atom2_ptr->checked = TRUE;
		  checked = TRUE;
		}
	      else if (atom2_ptr->checked == TRUE &&
		       atom1_ptr->checked == FALSE)
		{
		  atom1_ptr->checked = TRUE;
		  checked = TRUE;
		}
	    }

	  /* Get pointer to point to next bond in the list */
	  bond_ptr = bond_ptr->next_bond_ptr;
	}
    }
}
/***********************************************************************

delete_unreachable_residues  -  Delete any unreachable residues

***********************************************************************/

/* v.3.2--> */
/* void delete_unreachable_residues(int *ndeleted) */
void delete_unreachable_residues(int *ndeleted,int show_deletions)
/* <--v.3.2 */
{
  int iatom, iobject, iresid, natoms, nresid;
  int got_atom;

  struct coordinate *atom_ptr;
  struct object *object_ptr;
  struct object *next_object_ptr, *previous_object_ptr;
  struct residue *residue_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  previous_object_ptr = NULL;
  iobject = 0;

  /* Loop through all objects to check which are reachable */
  while (object_ptr != NULL)
    {
      /* Save pointer to the next object */
      next_object_ptr = object_ptr->next_object_ptr;

      /* Only want to look at non-ligand objects */
      if (object_ptr->object_type != LIGAND)
	{
	  /* Set flag indicating whether any reachable atoms located */
	  got_atom = FALSE;

	  /* Get pointer to the first of this object's residues */
	  residue_ptr = object_ptr->first_residue_ptr;
	  nresid = object_ptr->nresidues;
	  iresid = 0;

	  /* Loop over all this object's residues */
	  while (iresid < nresid && residue_ptr != NULL &&
		 got_atom == FALSE)
	    {
	      /* Get pointer to this residue's first atom */
	      atom_ptr = residue_ptr->first_atom_ptr;
	      
	      /* Initialise atom count */
	      iatom = 0;
	      natoms = residue_ptr->natoms;

	      /* Loop over all the residue's atoms, checking whether
		 they are involved in interactions */
	      while (iatom < natoms && atom_ptr != NULL &&
		     got_atom == FALSE)
		{
		  /* If atom is reachable, then set flag */
		  if (atom_ptr->checked == TRUE)
		    got_atom = TRUE;

		  /* Get pointer to the next atom */
		  atom_ptr = atom_ptr->next;
		  iatom++;
		}

	      /* Get pointer to the next residue */
	      residue_ptr = residue_ptr->next_residue_ptr;
	      iresid++;
	    }

	  /* If object is unreachable, then need to delete it and all its
	     residues */
	  if (got_atom == FALSE)
	    {
	      /* Delete the current object */
/* v.3.2--> */
/*	      delete_object(object_ptr,previous_object_ptr); */
	      delete_object(object_ptr,previous_object_ptr,show_deletions);
/* <--v.3.2 */
	      (*ndeleted)++;

	      /* Set current pointer to the previous object */
	      object_ptr = previous_object_ptr;
	    }
	}

      /* Save the current object pointer */
      previous_object_ptr = object_ptr;

      /* Get pointer to point to next object in the list and increment
	 object count */
      object_ptr = next_object_ptr;
      iobject++;
    }
}
/***********************************************************************

get_atom_sizes  -  Get all the atom-sizes, according to the type of
                   object the atoms belong to

***********************************************************************/

void get_atom_sizes(void)
{
  int object_type;

  float atom_size;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;
  struct object *object_ptr;

  /* Initialise all the maximum object- and residue- atom sizes */
  
  /* Get pointer to the first object */
  object_ptr = first_object_ptr;

  /* Loop through all the stored objects to initialise */
  while (object_ptr != NULL)
    {
      /* Initialise */
      object_ptr->max_atom_size = 0.0;

      /* Get pointer to the next object in the linked list */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* Get pointer to the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop through all the stored residues to initialise */
  while (residue_ptr != NULL)
    {
      /* Initialise */
      residue_ptr->max_atom_size = 0.0;

      /* Get pointer to the next residue in the linked list */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Get pointer to the first atom */
  atom_ptr = first_atom_ptr;

  /* Loop over all the atoms */
  while (atom_ptr != NULL)
    {
      /* Get the object and residue to which the atom belongs */
      residue_ptr = atom_ptr->residue_ptr;
      object_ptr = residue_ptr->object_ptr;

      /* Get the object type */
      object_type = object_ptr->object_type;

      /* Initialise size */
      atom_size = 0.0;

      /* Water atom */
      if (object_type == WATER)
/* v.4.0--> */
/*	if (Include->Nonligand_Atoms == TRUE) */
	if (Include->Water_Atoms == TRUE)
/* <--v.4.0 */
	  atom_size = Size_Val->Waters;
	else
	  atom_size = Size_Val->Nonligand_Bonds;

      /* Hydrophobic contact residue */
      else if (object_type == HYDROPHOBIC)
	atom_size = SPOKE_EXTENT * Size_Val->Hydrophobics;

      /* Simplified ligand representation */
      else if (residue_ptr->residue_type == SIMPLE_LIGAND &&
	       !strncmp(atom_ptr->atom_type," CA ",4))
	atom_size = Size_Val->Simple_Residues;

      /* Ligand atom */
      else if (object_type == LIGAND)
	if (Include->Ligand_Atoms == TRUE)
	  atom_size = Size_Val->Ligand_Atoms;
	else
	  atom_size = Size_Val->Ligand_Bonds;

      /* Non-ligand atom */
      else
	if (Include->Nonligand_Atoms == TRUE)
	  atom_size = Size_Val->Nonligand_Atoms;
	else
	  atom_size = Size_Val->Nonligand_Bonds;

      /* Store the atom's size */
      atom_ptr->atom_size = atom_size;

      /* Update maximum atom size for this atom's residue and object */
      if (atom_size > object_ptr->max_atom_size)
	object_ptr->max_atom_size = atom_size;
      if (atom_size > residue_ptr->max_atom_size)
	residue_ptr->max_atom_size = atom_size;

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next;
    }
}
/***********************************************************************
      
atom_linkage  -  Determine which atoms are covalently bonded to each other

***********************************************************************/

void atom_linkage(void)
{
  int loop;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr, *other_atom_ptr;
  struct atom_link *atom_link_ptr;

  /* Initialise pointer to start of linked list of bonds */
  bond_ptr = first_bond_ptr;
    
  /* Loop through all the bonds in the linked list, excluding any elastic
     bonds */
  while (bond_ptr != NULL)
    {
      /* Loop over this bond's two atoms */
      for (loop = 0; loop < 2; loop++)
	{
	  if (loop == 0)
	    {
	      atom_ptr = bond_ptr->first_atom_ptr;
	      other_atom_ptr = bond_ptr->second_atom_ptr;
	    }
	  else
	    {
	      atom_ptr = bond_ptr->second_atom_ptr;
	      other_atom_ptr = bond_ptr->first_atom_ptr;
	    }

	  /* Add the other atom to the current atom's list of covalent
	     partners */

	  /* Allocate memory for structure to hold this atom link */
	  atom_link_ptr
	    = (struct atom_link*)malloc(sizeof(struct atom_link));
	  if (atom_link_ptr == NULL)
	    {
	      printf("*** ERROR. Unable to allocate memory for");
	      printf(" struct atom_link\n");
	      printf("***        Program ligplot terminated with");
	      printf(" error.\n");
	      exit (1);
	    }

	  /* Store the link, pointing at the other atom */
	  atom_link_ptr->atom_ptr = other_atom_ptr;
	  atom_link_ptr->next_atom_link_ptr = NULL;

	  /* Increment count of this atom's partners */
	  atom_ptr->natom_links++;
		  
	  /* If this is not the first link for current atom, make
	     new link point to previous link */
	  if (atom_ptr->first_atom_link_ptr != NULL)
	    atom_link_ptr->next_atom_link_ptr
		      = atom_ptr->first_atom_link_ptr;

	  /* Point from current bond to current link */
	  atom_ptr->first_atom_link_ptr = atom_link_ptr;
	}

      /* Get pointer to next bond in linked-list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/***********************************************************************

ring_check  -  Checks whether a given bond is in a ring - a path traced
               from one end of the bond eventually returns to the bond
	       via the other end

***********************************************************************/

int ring_check(struct bond *test_bond_ptr)
{
  int got_ring;

  struct bond *bond_ptr, *next_free_stack_ptr, *other_bond_ptr,
  *next_stack_ptr;
  struct bond_link *bond_link_ptr;

  /* Initialise variables */
  got_ring = FALSE;
  next_stack_ptr = NULL;
  next_free_stack_ptr = NULL;

  /* Initialise all the bond flags and stack pointers */
  initialise_bonds();

  /* Get the pointer to the first of the bonds coming off the 
     test_bond */
  bond_link_ptr = test_bond_ptr->first_bond_link_ptr;

  /* Initialise stack pointers for ring-search by putting the bonds
     coming off one end of the test-bond onto the stack of bonds to
     be searched */
  while (bond_link_ptr != NULL)
    {
      /* Add this bond to the stack only if it comes off the first
	 atom of the test-bond */
      if (bond_link_ptr->bond_end == FIRST)
	{
	  /* Get the bond this link is pointing to and add to stack */
	  bond_ptr = bond_link_ptr->bond_ptr;

	  /* If this is the first to be added to the stack, set both
	     stack-pointers to point to it */
	  if (next_free_stack_ptr == NULL)
	    {
	      next_free_stack_ptr = bond_ptr;
	      next_stack_ptr = bond_ptr;
	    }

	  /* Otherwise, make last bond on stack point to this one
	     and set this one as nopw the last on the stack */
	  else
	    {
	      next_free_stack_ptr->next_stack_ptr = bond_ptr;
	      next_free_stack_ptr = bond_ptr;
	    }

	  /* Mark this bond as checked */
	  bond_ptr->checked = TRUE;
	}
      
      /* Go to the next bond-link in the linked list */
      bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
    }

  /* Loop until stack of pointers to be processed has been exhausted
     or have determined that we have a ring */
  while (next_stack_ptr != NULL && got_ring == FALSE)
    {
      /* Pick up the next bond from the stack, and move the stack
       pointer onto the next position */
      bond_ptr = next_stack_ptr;
      next_stack_ptr = bond_ptr->next_stack_ptr;

      /* Loop through all the bonds connected to the current one */

      /* Get the pointer to the first of the bonds coming off the 
	 current bond */
      bond_link_ptr = bond_ptr->first_bond_link_ptr;

      /* Process each of these bonds in turn */
      while (bond_link_ptr != NULL)
	{
	  /* Get the pointer to this bond */
	  other_bond_ptr = bond_link_ptr->bond_ptr;
	  
	  /* If this is the test-bond, check whether it is connected
	     to the current bond via the atom opposite to the one we
	     started with, implying we have gone round in a ring */
	  if (other_bond_ptr == test_bond_ptr)
	    {
	      /* Check if the current bond has the atom we're after */
	      if (bond_ptr->first_atom_ptr
		  == test_bond_ptr->second_atom_ptr ||
		  bond_ptr->second_atom_ptr
		  == test_bond_ptr->second_atom_ptr)
		{
		  /* Have gone round, so must have a ring */
		  got_ring = TRUE;
		}
	    }

	  /* Otherwise, if this bond hasn't already been processed
	     or added to the stack, add it now */
	  else if (other_bond_ptr->checked == FALSE)
	    {
	      next_free_stack_ptr->next_stack_ptr = other_bond_ptr;
	      next_free_stack_ptr = other_bond_ptr;

	      /* If the stack has run out, restart it at the bond
		 just entered */
	      if (next_stack_ptr == NULL)
		next_stack_ptr = other_bond_ptr;

	      /* Mark this bond as checked */
	      other_bond_ptr->checked = TRUE;
	    }

	  /* Go to the next bond-link in the linked list */
	  bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
	}
    }
    
  /* Return the flag showing whether a ring has been found or not */
  return(got_ring);
}
/***********************************************************************

identify_rotatable_bonds  -  Determine which bonds are rotatable in
                             each object

***********************************************************************/

void identify_rotatable_bonds(void)
{
  int in_ring, iobject;

  struct bond *bond_ptr;
  struct object *object_ptr;
  struct object_bond *object_bond_ptr;

  /* Initialise pointer to the first stored object */
  printf("\nIdentifying rotatable bonds ...\n");
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Get pointer to the first of this object's bonds */
      object_bond_ptr = object_ptr->first_object_bond_ptr;

      /* Loop through all this object's bonds, one by one */
      while (object_bond_ptr != NULL)
	{
	  /* Get the current bond */
	  bond_ptr = object_bond_ptr->bond_ptr;

	  /* To be a rotatable bond, has to have other bonds linked to
	     it at either end */
	  if (bond_ptr->nfirst_atom_links > 0 &&
	      bond_ptr->nsecond_atom_links > 0)
	    {
	      /* Check whether the current bond is in some sort of ring
		 - ie can start at one end and trace through the bonds to
		 the other */
	      in_ring = ring_check(bond_ptr);

	      /* If bond is a ring-bond, then set marker */
	      if (in_ring == TRUE)
		bond_ptr->ring_bond = TRUE;

	      /* Otherwise, must be a rotatable bond */
	      else
		{
		  bond_ptr->rotatable_bond = TRUE;

		  /* Increment count of rotatable bonds in the current
		     object */
		  object_ptr->nrot_bonds++;
		}
	    }

	    /* If the current bond has other bonds linked only to one
	       of its ends, then mark it as an "end" bond */
	  else if (bond_ptr->nfirst_atom_links > 0 ||
		   bond_ptr->nsecond_atom_links > 0)
	    bond_ptr->end_bond = TRUE;
	  
	  /* Get pointer to this object's next bond entry */
	  object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
	}

      /* If object has bonds, but none are rotatable, make the first
	 one a rotatable one */
      if (object_ptr->nbonds > 0 && object_ptr->nrot_bonds == 0)
	{
	  bond_ptr = object_ptr->first_object_bond_ptr->bond_ptr;
	  bond_ptr->rotatable_bond = TRUE;
	  object_ptr->nrot_bonds++;
	}

      /* Get pointer for next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }
}
/***********************************************************************

write_frm_file  -  Read through the PDB file, writing out all the
                   residues (complete atom records) that will be used
                   in the final LIGPLOT picture

***********************************************************************/

void write_frm_file(char pdb_name[FILENAME_LEN])
{
  char line[LINELEN + 1];
  char res_number[6], residue[4], last_residue[4], last_resnum[6];
/* v.3.2--> */
/*  char atom_counter, atom_type[5], chain, last_chain; */
  char atom_type[5], chain, last_chain;
/* <--v.3.2 */

  int keep_reading;
  int cpos, ipos, want_residue;
/* v.3.2--> */
  int atom_counter;
  int got_residue;
/* <--v.3.2 */

  struct residue *residue_ptr;
    
  FILE *fa;

  /* Open PDB file */
  if ((fa = fopen(pdb_name,"r")) ==  NULL)
    {
      printf("\n*** Unable to open your PDB file [%s]\n",pdb_name);
      exit(1);
    }
    
  /* Initialise variables */
  atom_counter = 0;
  keep_reading = TRUE;
  last_residue[0] = '\0';
  last_resnum[0] = '\0';
  last_chain = '\0';
  want_residue = FALSE;

  /* Read in the data from the PDB file */
  printf("\nWriting out the ligplot.frm file ...\n");

  /* Loop through the PDB file, picking up all the atoms that
     are required */
  while (fgets(line,LINELEN,fa)!= NULL)
    {
      /* Check for the end of an NMR structure model */
      if (!strncmp("ENDMDL",line,6))
	keep_reading = FALSE;

      /* If this is an ATOM or HETATM record, process it */
      else if (keep_reading == TRUE &&
	       (!strncmp(line,"ATOM",4) || !strncmp(line,"HETATM",6)))
	{
	  /* Get the residue name, number, chain-ID and atom type */
	  strncpy(residue,line+17,3);
	  residue[3] = '\0';
	  strncpy(res_number,line+22,5);
	  res_number[5] = '\0';
	  chain = line[21];
	  strncpy(atom_type,line+12,4);
	  atom_type[4] = '\0';

	  /* Check whether atom belongs to a new residue */
/* v.3.2--> */
/*	  if (chain != last_chain || strncmp(res_number,last_resnum,5)) */
	  if (chain != last_chain || strncmp(res_number,last_resnum,5) ||
	      strncmp(residue,last_residue,3))
/* <--v.3.2 */
	    {
	      /* Residue is a new one, need to check whether it
		 is one of those used for the LIGPLOT and hence is to
		 be written out in full */
/* v.3.2--> */
	      got_residue = FALSE;
/* <--v.3.2 */
	      want_residue = FALSE;

	      /* Initialise search pointer to first residue in the
		 linked list */
	      residue_ptr = first_residue_ptr;

	      /* Loop through all the stored residues to check for a
		 match */
/* v.3.2--> */
/*	      while (residue_ptr != NULL && want_residue == FALSE) */
	      while (residue_ptr != NULL && got_residue == FALSE)
/* <--v.3.2 */
		{
		  /* Check whether this is the correct residue */
		  if (!strncmp(res_number,residue_ptr->res_num,5) &&
		      !strncmp(residue,residue_ptr->res_name,3) &&
		      chain == residue_ptr->chain)
/* v.3.2--> */
/*		    want_residue = TRUE; */
		    got_residue = TRUE;

		  /* If have the residue, check that it hasn't been
		     deleted */
		  if (got_residue == TRUE && residue_ptr->deleted == FALSE)
		    want_residue = TRUE;

		  /* If residue is wanted and the ligplot.res file is
		     required, then write out the residue details */
		  if (want_residue == TRUE && Write_Res_File == TRUE)
/* v.3.2.1--> */
/*		    fprintf(ligplot_res,"%s%s%c\n", */
		    {
		      fprintf(ligplot_res,"%s%s%c",
/* <--v.3.2.1 */
			      residue_ptr->res_name,
			      residue_ptr->res_num,
			      residue_ptr->chain);
/* v.3.2.1--> */
		      /* If this is a ligand residue, then write out
			 marker */
		      if (residue_ptr->inligand == TRUE)
			fprintf(ligplot_res," *\n");
		      else
			fprintf(ligplot_res,"  \n");
		    }
/* <--v.3.2.1 */
/* <--v.3.2 */

		  /* Get pointer to the next residue in the linked list */
		  residue_ptr = residue_ptr->next_residue_ptr;
		}
/* v.3.2--> */
	      /* Save the current residue number and chain */
	      strncpy(last_residue,residue,3);
	      last_residue[3] = '\0';
	      strncpy(last_resnum,res_number,5);
	      last_resnum[5] = '\0';
	      last_chain = chain;
/* <--v.3.2 */
	    }

	  /* If this residue is wanted, then write the atom-record
	     out to the .frm file */
	  if (want_residue == TRUE)
	    {
	      /* Truncate the record at last non-blank character */
	      ipos = 0;
	      for (cpos = 0; cpos < 101 && line[cpos] != '\n'; cpos++)
		{
		  if (line[cpos] != ' ')
		    {
		      ipos = cpos + 1;
		    }
		}
	      line[ipos] = '\0';
	      
	      /* Write the record out */
	      fprintf(ligplot_frm,"%s\n",line);
	      atom_counter++;
	    }
	}
    }

    /* Print counts of records written out */
    printf("   Number of atoms written out       = %7d\n",atom_counter);
}
/* v.4.0--> */
/***********************************************************************

write_lig_prot_bonds  -  Write out any covalent bonds between the ligand
                         and the protein to the .res file

***********************************************************************/

void write_lig_prot_bonds(void)
{
  int nbonds;

  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct object *object1_ptr, *object2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;    

  /* Initialise variables */
  nbonds = 0;

  /* Initialise bond pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds */
  while (bond_ptr != NULL)
    {
      /* Process only if this is a covalent, elastic bond */
      if (bond_ptr->bond_type == COVALENT && bond_ptr->elastic == TRUE)
	{
	  /* Get pointers to the two atoms at either end of the
	     bond, and their respective residues and object */
	  atom1_ptr = bond_ptr->first_atom_ptr;
	  atom2_ptr = bond_ptr->second_atom_ptr;
	  residue1_ptr = atom1_ptr->residue_ptr;
	  residue2_ptr = atom2_ptr->residue_ptr;
	  object1_ptr = residue1_ptr->object_ptr;
	  object2_ptr = residue2_ptr->object_ptr;

	  /* Write out the bond to the .res file only if it is across
	     two different objects, of which one is a ligand */
	  if (object1_ptr != object2_ptr &&
	      (object1_ptr->object_type == LIGAND ||
	       object2_ptr->object_type == LIGAND))
	    {
	      /* Write out the bond */
	      fprintf(ligplot_res,"#Covalent links: ");
	      if (object1_ptr->object_type == LIGAND)
		fprintf(ligplot_res,"%s %s %s %c -> %s %s %s %c\n",
			atom1_ptr->atom_type,
			residue1_ptr->res_name,
			residue1_ptr->res_num,
			residue1_ptr->chain,
			atom2_ptr->atom_type,
			residue2_ptr->res_name,
			residue2_ptr->res_num,
			residue2_ptr->chain);
	      else
		fprintf(ligplot_res,"%s %s %s %c -> %s %s %s %c\n",
			atom2_ptr->atom_type,
			residue2_ptr->res_name,
			residue2_ptr->res_num,
			residue2_ptr->chain,
			atom1_ptr->atom_type,
			residue1_ptr->res_name,
			residue1_ptr->res_num,
			residue1_ptr->chain);

	      /* Increment count of bonds written out */
	      nbonds++;
	    }
	}

      /* Get pointer for next object */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Print counts of records written out */
  printf("   Number of covalent ligand-protein bonds = %7d\n",nbonds);
}
/* <--v.4.0 */
/***********************************************************************

pare_down_hydrophobics  -  Mark any unwanted atoms in hydrophobic groups
                           for deletion

***********************************************************************/

void pare_down_hydrophobics(void)
{
  int iatom, iobject, iresid, natoms, nresid;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr;
  struct object *object_ptr;
/* v.3.2--> */
  struct object *object1_ptr, *object2_ptr;
/* <--v.3.2 */
  struct residue *residue_ptr;

  /* Initialise all atoms */
  printf("\nMarking unwanted atoms ...\n");
  initialise_atoms();

  /* Initialise bond pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds */
  while (bond_ptr != NULL)
    {
      /* Process only if this is a hydrophobic contact */
      if (bond_ptr->bond_type == CONTACT)
	{
	  /* Mark both atoms as being involved in the bond */
	  bond_ptr->first_atom_ptr->checked = TRUE;
	  bond_ptr->second_atom_ptr->checked = TRUE;
	}

/* v.3.2--> */
      /* Get the two objects these atoms belong to */
      object1_ptr = bond_ptr->first_atom_ptr->residue_ptr->object_ptr;
      object2_ptr = bond_ptr->second_atom_ptr->residue_ptr->object_ptr;

      /* If this is a bond between a hydrophobic group and anything
	 other than a ligand residue, then delete (except for case of
	 interface plots) */
      if (object1_ptr != object2_ptr)
	{
/* v.4.0--> */
	  if (Interface_Plot == FALSE)
	    {
/* <--v.4.0 */
	      if ((object1_ptr->object_type == HYDROPHOBIC &&
		   object2_ptr->first_residue_ptr->inligand == FALSE) ||
		  (object2_ptr->object_type == HYDROPHOBIC &&
		   object1_ptr->first_residue_ptr->inligand == FALSE))
		{
		  /* Mark this bond for deletion */
		  bond_ptr->bond_type = DELETED;
		}
/* v.4.0--> */
	    }
/* <--v.4.0 */
	}
/* <--v.3.2 */
      /* Get pointer for next object */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Process only if this is a hydrophobic group */
      if (object_ptr->object_type == HYDROPHOBIC)
	{
	  /* Get pointer to the first of this object's residues */
	  residue_ptr = object_ptr->first_residue_ptr;
	  nresid = object_ptr->nresidues;
	  iresid = 0;

	  /* Loop over all this object's residues */
	  while (iresid < nresid && residue_ptr != NULL)
	    {
	      /* Get pointer to this residue's first atom */
	      atom_ptr = residue_ptr->first_atom_ptr;
	      
	      /* Initialise atom count */
	      iatom = 0;
	      natoms = residue_ptr->natoms;

	      /* Loop over all the residue's atoms, checking whether
		 they are involved in interactions */
	      while (iatom < natoms && atom_ptr != NULL)
		{
		  /* If atom not involved in any hydrophobic contacts
		     then can mark it for deletion */
		  if (atom_ptr->checked == FALSE)
		    atom_ptr->deleted = TRUE;

		  /* Get pointer to the next atom */
		  atom_ptr = atom_ptr->next;
		  iatom++;
		}

	      /* Get pointer to the next residue */
	      residue_ptr = residue_ptr->next_residue_ptr;
	      iresid++;
	    }
	}

      /* Get pointer for next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }
}
/***********************************************************************

centre_of_mass  -  Find the centre of mass of the supplied coordinates

***********************************************************************/

void centre_of_mass(float coord_store[MAXATOMS][3],int natoms,
		    float *centre_x,float *centre_y,float *centre_z)
{
  int iatom;
  float mass_x = 0.0, mass_y = 0.0, mass_z = 0.0, centre_mass_x,
  centre_mass_y, centre_mass_z;

  /* Loop through all the coordinates supplied */ 
  for (iatom = 0; iatom < natoms; iatom++)
    {
      mass_x = mass_x + coord_store[iatom][0];
      mass_y = mass_y + coord_store[iatom][1];
      mass_z = mass_z + coord_store[iatom][2];
    }

  /* Calculate centre of mass */
  centre_mass_x = mass_x / natoms;
  centre_mass_y = mass_y / natoms;
  centre_mass_z = mass_z / natoms;
  *centre_x = centre_mass_x;
  *centre_y = centre_mass_y;
  *centre_z = centre_mass_z;
}
/***********************************************************************

adjust_stored_coords  -  Adjust the stored coordinates to origin

***********************************************************************/

void adjust_stored_coords(int natoms,float coord_store[MAXATOMS][3],
			  float centre_x,float centre_y,float centre_z)
{
  int iatom;

  for (iatom = 0; iatom < natoms; iatom++)
    {
      coord_store[iatom][0] = coord_store[iatom][0] - centre_x;
      coord_store[iatom][1] = coord_store[iatom][1] - centre_y;
      coord_store[iatom][2] = coord_store[iatom][2] - centre_z;
    }
}
/***********************************************************************
  
calc_eigen_values  -  Subroutine to compute eigenvalues & eigenvectors
                      of a real symmetric matrix, from IBM SSP manual
		      (see p165).

		      Ian Tickle April 1992
 
		      (modified by David Moss February 1993)

Eigenvalues & vectors of real symmetric matrix stored in triangular form.
 
Arguments:
 
mv = 0 to compute eigenvalues only.
mv = 1 to compute eigenvalues & eigenvectors.
 
n - supplied dimension of n*n matrix.
 
a - supplied array of size n*(n+1)/2 containing n*n matrix in the order:
 
 	   1      2      3    ...
     1    a[0]
     2    a[1]   a[2]
     3    a[3]   a[4]   a[5]   ...
     ...

     NOTE a is used as working space and is overwritten on return.
     Eigenvalues are written into diagonal elements of a
     i.e.  a[0]  a[2]  a[5]  for a 3*3 matrix.
 
r - Resultant matrix of eigenvectors stored columnwise in the same
     order as eigenvalues.

***********************************************************************/

void calc_eigen_values(int mv, int n, float* a, float* r)
{
  int i, il, ilq, ilr, im, imq, imr, ind, iq, j, l, ll, lm, lq, m, mm, mq;
  float anorm, anrmx, cosx, cosx2, sincs, sinx, sinx2, thr, x, y;

  if (mv)
    { for (i=1; i<n*n; i++) r[i]=0.;
      for (i=0; i<n*n; i+=n+1) r[i]=1.;
    }

  /* Initial and final norms (anorm & anrmx). */
  anorm=0.;
  iq=0;
  for (i=0; i<n; i++) for (j=0; j<=i; j++)
    { if (j!=i) anorm+=a[iq]*a[iq];
      ++iq;
    }

  if (anorm>0.)
    { anorm=sqrt(2.0*anorm);
      anrmx=1.0e-6*anorm/n;

      /* Compute threshold and initialise flag. */
      thr=anorm;
      do
	{ thr/=n;
      do
	{ ind=0;
	  l=0;
	  do
	    { lq=l*(l+1)/2;
	      ll=l+lq;
	      m=l+1;
	      ilq=n*l;
	      do

		/* Compute sin & cos. */
		{ mq=m*(m+1)/2;
		  lm=l+mq;
		  if (fabs(a[lm])>=thr)
		    { ind=1;
		      mm=m+mq;
		      x=0.5*(a[ll]-a[mm]);
		      y=-a[lm]/sqrt(a[lm]*a[lm]+x*x);
		      if (x<0.0) y=-y;
		      sinx=y/sqrt(2.0*(1.0+(sqrt(1.0-y*y))));
		      sinx2=sinx*sinx;
		      cosx=sqrt(1.0-sinx2);
		      cosx2=cosx*cosx;
		      sincs=sinx*cosx;

		      /* Rotate l & m columns. */
		      imq=n*m;
		      for (i=0; i<n; i++)
			{ iq=i*(i+1)/2;
			  if (i!=l && i!=m)
			    { if (i<m) im=i+mq;
			    else im=m+iq;
			      if (i<l) il=i+lq;
			      else il=l+iq;
			      x=a[il]*cosx-a[im]*sinx;
			      a[im]=a[il]*sinx+a[im]*cosx;
			      a[il]=x;
			    }
			  if (mv)
			    { ilr=ilq+i;
			      imr=imq+i;
			      x=r[ilr]*cosx-r[imr]*sinx;
			      r[imr]=r[ilr]*sinx+r[imr]*cosx;
			      r[ilr]=x;
			    }
			}

		      x=2.*a[lm]*sincs;
		      y=a[ll]*cosx2+a[mm]*sinx2-x;
		      x=a[ll]*sinx2+a[mm]*cosx2+x;
		      a[lm]=(a[ll]-a[mm])*sincs+a[lm]*(cosx2-sinx2);
		      a[ll]=y;
		      a[mm]=x;
		    }

		  /* Tests for completion.
		     Test for m = last column. */
		} while (++m!=n);

	      /* Test for l =penultimate column. */
	    } while (++l!=n-1);
	} while (ind);
	  
	  /* Compare threshold with final norm. */
	} while (thr>anrmx);
    }
}

/***********************************************************************

eigen_sort  -  Sort eigenvalues & eigenvectors (if mv=1) in order of
               descending eigenvalue.
	       Arguments are as for eigen, which must have been called.

***********************************************************************/

void eigen_sort(int mv, int n, float* a, float* r)
{
  int i, im, j, k, km, l, m;
  float am;

  k=0;
  for (i=0; i<n-1; i++)
    { im=i;
      km=k;
      am=a[k];
      l=0;
      for (j=0; j<n; j++)
	{ if (j>i && a[l]>am)
	    { im=j;
	      km=l;
	      am=a[l];
	    }
	  l+=j+2;
	}
      if (im!=i)
	{ a[km]=a[k];
	  a[k]=am;
	  if (mv)
	    { l=n*i;
	      m=n*im;
	      for (j=0; j<n; j++)
		{ am=r[l];
		  r[l++]=r[m];
		  r[m++]=am;
		}
	    }
	}
      k+=i+2;
    }
}
/***********************************************************************

vecprd  -  Calculates the vector product (vx,vy,vz) of two vectors 
           (px,py,pz) and (qx,qy,qz)

***********************************************************************/

void vecprd(float px, float py, float pz, float qx, float qy, float qz,
            float *vx, float *vy, float *vz)
{
  *vx = py * qz - pz * qy;
  *vy = pz * qx - px * qz;
  *vz = px * qy - py * qx;
}
/***********************************************************************

dot3  -  Returns the scalar product of the two vectors (x1,y1,z1) 
         and (x2,y2,z2)

***********************************************************************/

float dot3(float x1, float y1, float z1, float x2, float y2, float z2)
{
  float dot_prod;

  dot_prod = x1 * x2 + y1 * y2 + z1 * z2;
  return(dot_prod);
}
/***********************************************************************

orien3  -  Returns the orientation of the polygon with the consecutive
           vertices (x1,y1,z1), (x2,y2,z2), and (x3,y3,z3) as viewed
           from (ex,ey,ez).
                      -1  :  clockwise orientation
                      +1  :  anti-clockwise orientation (desired order)
                       0  :  degenerate (line or point)

 [Adapted from Angell I O & Griffith G (1987). High-resolution Computer
  Graphics using FORTRAN 77. Macmillan Education Ltd, Basingstoke, UK.]

***********************************************************************/

int orien3(float x1, float y1, float z1, float x2, float y2, float z2,
           float x3, float y3, float z3, float ex, float ey, float ez)
{
    float a, b, c, dx1, dy1, dz1, dx2, dy2, dz2, fe;
    int orientation;

    /* Initialise variables */
    dx1 = x2 - x1;
    dy1 = y2 - y1;
    dz1 = z2 - z1;
    dx2 = x3 - x2;
    dy2 = y3 - y2;
    dz2 = z3 - z2;

    /* Calculate triple scalar product */
    vecprd(dx1,dy1,dz1,dx2,dy2,dz2,&a,&b,&c);
    fe = dot3(a,b,c,ex-x1,ey-y1,ez-z1);
    if (fe < 0.0)
      orientation = -1;
    else
      orientation = 1;
    return(orientation);
}
/***********************************************************************

principal_components  -  Calculate the transformation to align the
                          supplied atom coords with their principal axes
                          along the x-, y-, and z-directions

***********************************************************************/

void principal_components(int natoms,float transformation[3][3],
			  float coord_store[MAXATOMS][3])
{
  int iatom, iorder;
  int nzeroy, nzeroz;
  int array_size, i, ipos, j, nsize;
  float amatrix[6], matrix[3][3], rmatrix[9], x, y, z;
  float eigen_value[3];
  float xxsum, xysum, xzsum, yysum, yzsum, zzsum;

  /* Initialise counts */
  nzeroy = 0;
  nzeroz = 0;
  xxsum = 0.0;
  xysum = 0.0;
  xzsum = 0.0;
  yysum = 0.0;
  yzsum = 0.0;
  zzsum = 0.0;

/* v.3.2--> */
  /* Initialise transformation matrix */
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      if (i == j)
	transformation[i][j] = 1.0;
      else
	transformation[i][j] = 0.0;

  /* If only have a single atom, then return */
  if (natoms < 2)
    return;
/* <--v.3.2 */

  /* Loop through all data points */
  for (iatom = 0; iatom < natoms; iatom++)
    {
      /* Retrieve the coordinates for this point */
      x = coord_store[iatom][0];
      y = coord_store[iatom][1];
      z = coord_store[iatom][2];

      /* Count instances where the y- and z-coords are zero */
      if (fabs(y) < 0.000001)
	nzeroy++;
      if (fabs(z) < 0.000001)
	nzeroz++;

      /* Compute the terms of the scalar products matrix */
      xxsum = xxsum+x*x;
      xysum = xysum+x*y;
      xzsum = xzsum+x*z;
      yysum = yysum+y*y;
      yzsum = yzsum+y*z;
      zzsum = zzsum+z*z;
    }
    
  /* Form the scalar products matrix */
  matrix[0][0] = xxsum;
  matrix[0][1] = xysum;
  matrix[0][2] = xzsum;
  matrix[1][0] = xysum;
  matrix[1][1] = yysum;
  matrix[1][2] = yzsum;
  matrix[2][0] = xzsum;
  matrix[2][1] = yzsum;
  matrix[2][2] = zzsum;

  /* If there are only two points, and both their z-coords are zero, or
     both their y- and z-coords are zero, adjust the order of the scalar
     products matrix whose eigenvalues are required */
  nsize = 3;
  if ((natoms == 2) && (nzeroz == natoms))
    {
      if (nzeroy == natoms)
	nsize = 1;
      else
	nsize = 2;
    }

  /* Calculate the eigen vectors and eigen values providing order of the
     matrix is not 1 */
  if (nsize > 1)
    {
      /* Convert matrix into form required by the eigen routine */
      ipos = 0;
      for (i = 0; i < nsize; i++)
	for (j = 0; j <= i; j++)
	  {
	    amatrix[ipos] = matrix[i][j];
	    ipos++;
	  }

      /* Calculate the eigen values */
      array_size = nsize * (nsize + 1) / 2;
      calc_eigen_values(1,nsize,amatrix,rmatrix);

      /* Sort the eigen vectors according to the eigen values */
      eigen_sort(1,nsize,amatrix,rmatrix);

      /* Retrieve the eigen values */
      i = j = 0;
      for (ipos = 0; ipos < array_size; ipos++)
	{
	  /* If this is a diagonal term, then is one of the eigen values */
	  if (i == j)
	    eigen_value[i] = amatrix[ipos];

	  /* Increment array indices */
	  j++;
	  if (j > i)
	    {
	      i++;
	      j = 0;
	    }
	}

      /* Retrieve the eigen vectors */
      i = j = 0;
      for (ipos = 0; ipos < nsize * nsize; ipos++)
	{
	  /* Store the transformation matrix */
	  transformation[j][i] = rmatrix[ipos];

	  /* Increment array indices */
	  j++;
	  if (j > nsize - 1)
	    {
	      i++;
	      j = 0;
	    }
	}

      /* If the matrix was of order 2, fix up row 3 and column 3 of the
	 output transformation matrix */
      if (nsize == 2)
	{
	  for (i = 0; i < 3; i++)
	    {
	      for (j = 0; j < 3; j++)
		{
		  if ((i == 2) && (j == 2))
		    transformation[i][j] = 1.0;
		  else if ((i == 2) || (j == 2))
		    transformation[i][j] = 0.0;
		}
	    }
	}

      /* If the matrix was of order 2, fix up row 3 and column 3 of the
	 output transformation matrix */
      if (nsize == 2)
	{
	  for (i = 0; i < 3; i++)
	    {
	      for (j = 0; j < 3; j++)
		{
		  if ((i == 3) && (j == 3))
		    transformation[i][j] = 1.0;
		  else if ((i == 3) || (j == 3))
		    transformation[i][j] = 0.0;
		}
	    }
	}
    }

    /* If the matrix was of order 1, make the transormation matrix the
       identity matrix */
    else
      {
	for (i = 0; i < 3; i++)
	  {
	    for (j = 1; j < 4; j++)
	      {
		if (i == j)
		  transformation[i][j] = 1.0;
		else
		  transformation[i][j] = 0.0;
	      }
	  }
      }
    
  /* Check whether the sorted eigenvectors make a right-handed set */
  iorder = orien3(transformation[0][0],transformation[0][1],
		  transformation[0][2],0.0,0.0,0.0,
		  transformation[1][0],transformation[1][1],
		  transformation[1][2],transformation[2][0],
		  transformation[2][1],transformation[2][2]);
  
  /* If the sorted vectors don't make a right-handed set, then reverse the
     signs on the third vector */
  if (iorder == 1)
    {
      transformation[0][2] = - transformation[0][2];
      transformation[1][2] = - transformation[1][2];
      transformation[2][2] = - transformation[2][2];
    }
}
/***********************************************************************

apply_rotation  -  Apply the given rotation to the supplied coordinates

***********************************************************************/

void apply_rotation(float x,float y,float z,float *xnew,
		    float *ynew,float *znew,float mat[3][3])
{
  /* Apply the tranformation */
  *xnew = mat[0][0] * x + mat[1][0] * y + mat[2][0] * z;
  *ynew = mat[0][1] * x + mat[1][1] * y + mat[2][1] * z;
  *znew = mat[0][2] * x + mat[1][2] * y + mat[2][2] * z;
}
/***********************************************************************

calc_fit_coordinates  -  Align principal axes of original coords with
                         the x-, y- and z-axes. Repeat for each object

***********************************************************************/

void calc_fit_coordinates(int nligands)
{
  int iobject, iatom, iresid, natoms, nresid, nstored;
  int wanted;

  float coord_store[MAXATOMS][3], mass_x, mass_y, mass_z;
  float new_x, new_y, new_z, x, y, z;
  float transformation[3][3];

  struct coordinate *atom_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  printf("\nCalculating fits ...\n");
  nstored = 0;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms and store their original coords */
  while (atom_ptr != NULL)
    {
      wanted = FALSE;

      /* If there is no ligand, then use all atoms */
      if (nligands == 0)
	wanted = TRUE;

      /* Otherwise, only want this atom if it belongs to the ligand */
      else
	{
	  object_ptr = atom_ptr->residue_ptr->object_ptr;
	  if (object_ptr->object_type == LIGAND)
	    wanted = TRUE;
	}

      /* Store this atom's coords */
      if (wanted == TRUE && nstored < MAXATOMS)
	{
	  coord_store[nstored][0] = atom_ptr->original_x;
	  coord_store[nstored][1] = atom_ptr->original_y;
	  coord_store[nstored][2] = atom_ptr->original_z;
	  nstored++;
	}

      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next;
    }	       

  /* If have some atoms to fit onto, then do the fit */
  if (nstored > 0)
    {
      /* Calculate the centre of mass of these atoms */
      centre_of_mass(coord_store,nstored,&mass_x,&mass_y,&mass_z);

      /* Adjust coordinates to place centre of mass at the origin */
      adjust_stored_coords(nstored,coord_store,mass_x,mass_y,mass_z);

      /* Calculate transformation to orient the principal axes of the
	 original structure along the x-, y- and z-axes */
      principal_components(nstored,transformation,coord_store);

      /* Get pointer to first atom again */
      atom_ptr = first_atom_ptr;

      /* Loop through all the atoms to compute their fitted coords */
      while (atom_ptr != NULL)
	{
	  /* Get the atom's original coords */
	  x = atom_ptr->original_x - mass_x;
	  y = atom_ptr->original_y - mass_y;
	  z = atom_ptr->original_z - mass_z;

	  /* Apply rotation matrix */
	  apply_rotation(x,y,z,&new_x,&new_y,&new_z,transformation);
	  atom_ptr->fit_x = new_x;
	  atom_ptr->fit_y = new_y;
	  atom_ptr->fit_z = new_z;

	  /* Get pointer to next atom in linked-list */
	  atom_ptr = atom_ptr->next;
	}
    }

  /* Get the translation that puts each object's fitted coords at the
     origin */

  /* Get pointer to first object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all object's residues to calculate its centre of mass */
  while (object_ptr != NULL)
    {
      /* Initialise count of stored coordinates */
      nstored = 0;

      /* Set the pointer to the first of this object's residues */
      residue_ptr = object_ptr->first_residue_ptr;
      nresid = object_ptr->nresidues;
      iresid = 0;

      /* Loop over all this object's residues */
      while (iresid < nresid && residue_ptr != NULL)
	{
	  /* Get pointer to first residue's first atom */
	  atom_ptr = residue_ptr->first_atom_ptr;
	  iatom = 0;
	  natoms = residue_ptr->natoms;

	  /* Loop over all the residue's atoms */
	  while (iatom < natoms && atom_ptr != NULL)
	    {
	      /* Store this atom's fitted coords */
	      if (nstored < MAXATOMS)
		{
		  coord_store[nstored][0] = atom_ptr->fit_x;
		  coord_store[nstored][1] = atom_ptr->fit_y;
		  coord_store[nstored][2] = 0.0;
		  nstored++;
		}

	      /* Get pointer to the next atom */
	      atom_ptr = atom_ptr->next;
	      iatom++;
	    }

	  /* Get pointer to the next residue */
	  residue_ptr = residue_ptr->next_residue_ptr;
	  iresid++;
	}

      /* Calculate the centre of mass of these atoms */
      centre_of_mass(coord_store,nstored,&mass_x,&mass_y,&mass_z);

      /* Store this centre of mass */
      object_ptr->place_x = mass_x;
      object_ptr->place_y = mass_y;

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }
}
/***********************************************************************

get_other_ring_bonds  -  Get all the other bonds in the current ring

***********************************************************************/

void get_other_ring_bonds(struct bond *test_bond_ptr)
{
  struct bond *bond_ptr, *next_free_stack_ptr, *other_bond_ptr,
  *next_stack_ptr, *next_link_ptr;
  struct bond_link *bond_link_ptr;

  /* Initialise variables */
  next_link_ptr = NULL;
  next_stack_ptr = NULL;
  next_free_stack_ptr = NULL;

  /* Initialise all the bond flags and stack pointers */
  initialise_bonds();

  /* Add the test-bond onto the start of the stack */
  bond_ptr = test_bond_ptr;
  next_stack_ptr = bond_ptr;
  next_free_stack_ptr = bond_ptr;
  bond_ptr->checked = TRUE;

  /* Initialise linked-list that will eventually contain all the
     bonds in the current ring */
  next_link_ptr = bond_ptr;

  /* Loop until stack of pointers to be processed has been exhausted
     or have determined that we have a ring */
  while (next_stack_ptr != NULL)
    {
      /* Pick up the next bond from the stack, and move the stack
       pointer onto the next position */
      bond_ptr = next_stack_ptr;
      next_stack_ptr = bond_ptr->next_stack_ptr;

      /* Loop through all the bonds connected to the current one */

      /* Get the pointer to the first of the bonds coming off the 
	 current bond */
      bond_link_ptr = bond_ptr->first_bond_link_ptr;

      /* Process each of these bonds in turn */
      while (bond_link_ptr != NULL)
	{
	  /* Get the pointer to this bond */
	  other_bond_ptr = bond_link_ptr->bond_ptr;

	  /* Process only if this is a ring-bond or end-bond which
	     hasn't already been processed */
	  if (other_bond_ptr->checked == FALSE &&
	      (other_bond_ptr->ring_bond == TRUE ||
	       other_bond_ptr->end_bond == TRUE))
	    {
	      /* If this is a ring-bond , then add it to the stack */
	      if (other_bond_ptr->ring_bond == TRUE)
		{
		  /* Add to end of stack of rings to be processed */
		  next_free_stack_ptr->next_stack_ptr = other_bond_ptr;
		  next_free_stack_ptr = other_bond_ptr;

		  /* If the stack has run out, restart it at the bond
		     just entered */
		  if (next_stack_ptr == NULL)
		    next_stack_ptr = other_bond_ptr;
		}

	      /* Add to the linked-list of bonds belonging to the
		 current ring (and end-bonds sprouting off it) */
	      next_link_ptr->next_link_ptr = other_bond_ptr;
	      next_link_ptr = other_bond_ptr;

	      /* Mark this bond as checked */
	      other_bond_ptr->checked = TRUE;
	    }

	  /* Go to the next bond-link in the linked list */
	  bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
	}
    }
}
/***********************************************************************

get_ring_atom_coords  -  Extract all the coordinates of the current
                         ring's atoms, so that can squash to make planar

***********************************************************************/

void get_ring_atom_coords(struct bond *ring_bond_ptr,
			  float coord_store[MAXATOMS][3],
			  int *nring_atoms,int *nring_bonds)
{
  int loop, natoms, nbonds;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr;

  /* Initialise variables */
  natoms = 0;
  nbonds = 0;

  /* Initialise all the atom flags */
  initialise_atoms();

  /* Get pointer to the first bond in the current ring */
  bond_ptr = ring_bond_ptr;

  /* Loop through all the bonds in the current ring */
  while (bond_ptr != NULL)
    {
      /* Only want coordinates if this is a ring-bond (as opposed
	 to an associated end-bond) */
      if (bond_ptr->ring_bond == TRUE)
	{
	  /* Increment count of bonds in this ring */
	  nbonds++;

	  /* Loop over this bond's two atoms */
	  for (loop = 0; loop < 2; loop++)
	    {
	      if (loop == 0)
		atom_ptr = bond_ptr->first_atom_ptr;
	      else
		atom_ptr = bond_ptr->second_atom_ptr;

	      /* If this atom's coordinates have not already been stored,
		 then add to list for flattening */
	      if (atom_ptr->checked == FALSE)
		{
		  /* Store the atom's coordinates */

		  /* Check that maximum storage not exceeded */
		  if (natoms > MAXSTACK)
		    {
		      printf("*** Maximum number of atoms");
		      printf(" exceeded in flatten_ring");
		      printf("  %d\n",MAXSTACK);
		      exit(1);
		    }
		  coord_store[natoms][0] = atom_ptr->x;
		  coord_store[natoms][1] = atom_ptr->y;
		  coord_store[natoms][2] = atom_ptr->z;
		  natoms++;

		  /* Mark atom as already included */
		  atom_ptr->checked = TRUE;
		}
	    }
	}

      /* Mark this bond as already flattened */
	  bond_ptr->flattened = TRUE;

      /* Get the next bond in the current ring */
      bond_ptr = bond_ptr->next_link_ptr;
    }

  /* Return the number of atoms stored */
  *nring_atoms = natoms;
  *nring_bonds = nbonds;
}
/***********************************************************************

move_object  -  Move the given object by the supplied translation

***********************************************************************/

void move_object(struct object *object_ptr,float move_x,
		 float move_y,float move_z)
{
  int iatom, iresid, natoms, nresid;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, performing the
	 transformation */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* Apply the move */
	  atom_ptr->x = atom_ptr->x - move_x;
	  atom_ptr->y = atom_ptr->y - move_y;
	  atom_ptr->z = atom_ptr->z - move_z;

	  /* Get pointer to the next atom */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

      /* Adjust the residue boundaries */
      residue_ptr->minx = residue_ptr->minx - move_x;
      residue_ptr->maxx = residue_ptr->maxx - move_x;
      residue_ptr->miny = residue_ptr->miny - move_y;
      residue_ptr->maxy = residue_ptr->maxy - move_y;

/* v.4.0--> */
      /* Adjust the residue's mean position */
      residue_ptr->flattened_mean_x = residue_ptr->flattened_mean_x - move_x;
      residue_ptr->flattened_mean_y = residue_ptr->flattened_mean_y - move_y;
/* <--v.4.0 */

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Adjust the object boundaries */
  object_ptr->minx = object_ptr->minx - move_x;
  object_ptr->maxx = object_ptr->maxx - move_x;
  object_ptr->miny = object_ptr->miny - move_y;
  object_ptr->maxy = object_ptr->maxy - move_y;
}
/***********************************************************************

place_group_at_origin  -  Place all the atoms of the current hydrophobic
                          group at the origin

***********************************************************************/

void place_group_at_origin(struct object *object_ptr)
{
  int iatom, iresid, natoms, nresid;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Initialise object boundaries */
  object_ptr->minx = 0.0;
  object_ptr->maxx = 0.0;
  object_ptr->miny = 0.0;
  object_ptr->maxy = 0.0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, zeroing the coordinates */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* Zero the coordinates */
	  atom_ptr->x = 0.0;
	  atom_ptr->y = 0.0;
	  atom_ptr->z = 0.0;

	  /* Get pointer to the next atom */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

      /* Initialise residue boundaries */
      residue_ptr->minx = 0.0;
      residue_ptr->maxx = 0.0;
      residue_ptr->miny = 0.0;
      residue_ptr->maxy = 0.0;

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

rotate_object  -  Apply the given rotation to the given object

***********************************************************************/

void rotate_object(struct object *object_ptr,float v[3][3])
{
  int iatom, icoord, iresid, natoms, nresid;
  float x, y, z;
  float new_coords[3];

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, performing the
	 transformation */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* Extract this atom's coordinates */
	  x = atom_ptr->x;
	  y = atom_ptr->y;
	  z = atom_ptr->z;
	
	  /* Apply the transformation */
	  for (icoord = 0; icoord < 3; icoord++)
	    {
	      /* Calculate transformed coordinates */
	      new_coords[icoord]
		=  ((x * v[0][icoord])
		    + (y * v[1][icoord])
		    + (z * v[2][icoord]));
	    }

	  /* Save the new coordinates */
	  atom_ptr->x = new_coords[0];
	  atom_ptr->y = new_coords[1];
	  atom_ptr->z = new_coords[2];

	  /* Get pointer to the next atom */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

squash_z  -  Set the z-coord of the ring-atoms (and attached end-bond
             atoms) to zero to completely flatten current ring

***********************************************************************/

void squash_z(struct bond *ring_bond_ptr)
{
  struct bond *bond_ptr;

  /* Get pointer to the first bond in the current ring */
  bond_ptr = ring_bond_ptr;

  /* Loop through all the bonds in the current ring */
  while (bond_ptr != NULL)
    {
      /* Set z-coordinate of both this bond's atoms to zero */
      bond_ptr->first_atom_ptr->z = 0.0;
      bond_ptr->second_atom_ptr->z = 0.0;

      /* Get the next bond in the current ring */
      bond_ptr = bond_ptr->next_link_ptr;
    }
}
/***********************************************************************

flatten_ring  -  Flatten the current ring, given the starting bond

***********************************************************************/

void flatten_ring(struct object *object_ptr,struct bond *ring_bond_ptr)
{
  int natoms, nbonds;

  float coord_store[MAXATOMS][3], mass_x, mass_y, mass_z;
  float transformation[3][3];

  struct residue *residue_ptr;

  /* Extract all the coordinates of the current ring's atoms, so that
     can squash to make whole ring perfectly planar */
  get_ring_atom_coords(ring_bond_ptr,coord_store,&natoms,&nbonds);

  /* Have the coords of all the atoms making up the current ring,
     so calculate the centre of mass of these atoms */
  centre_of_mass(coord_store,natoms,&mass_x,&mass_y,&mass_z);

  /* Transform the object to place the ring's centre of mass at the
     origin */
  move_object(object_ptr,mass_x,mass_y,mass_z);

  /* Apply the same transformation to the stored ring-atom coordinates
     prior to calculating the best-fit plane */
  adjust_stored_coords(natoms,coord_store,mass_x,mass_y,mass_z);

  /* Calculate transformation to orient the plane of the ring
     approximately in the x-y plane */
  principal_components(natoms,transformation,coord_store);

  /* Apply this transformation to the entire object */
  rotate_object(object_ptr,transformation);

  /* Zero the z-coordinates of the atoms in the ring and any end bonds
     springing off it */
  squash_z(ring_bond_ptr);

  /* Display message to show which ring has been flattened */
  residue_ptr = ring_bond_ptr->first_atom_ptr->residue_ptr;
  printf("   Flattened %2d-bond ring of:  %s %s %c\n",
	 nbonds,residue_ptr->res_name,residue_ptr->res_num,
	 residue_ptr->chain);
}
/***********************************************************************

flatten_all_rings  -  Flatten any rings in the current object

***********************************************************************/

int flatten_all_rings(struct object *object_ptr)
{
  int n_flattened;

  struct bond *bond_ptr;
  struct object_bond *object_bond_ptr;

  /* Initialise variables */
  n_flattened = 0;

  /* Get pointer to the first of this object's bonds */
  object_bond_ptr = object_ptr->first_object_bond_ptr;

  /* Loop through all object's bonds */
  while (object_bond_ptr != NULL)
    {
      /* Get the bond pointer */
      bond_ptr = object_bond_ptr->bond_ptr;

      /* If this is a ring-bond which isn't part of a ring already
	 flattened, then process */
      if (bond_ptr->ring_bond == TRUE && bond_ptr->flattened == FALSE)
	{
	  /* Get all the other bonds involved in this ring */
	  get_other_ring_bonds(bond_ptr);

	  /* Flatten the ring */
	  flatten_ring(object_ptr,bond_ptr);
	  n_flattened++;
	}

      /* Get pointer to point to next bond for this object */
      object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
    }

  /* Return the number of rings flattened */
  return(n_flattened);
}
/***********************************************************************

get_angle  -  Function to calculate the angle between line from origin
              to supplied point (x,y), and the x-axis. The value
              returned is between 0 and 2 pi radians.

***********************************************************************/

float get_angle(float x, float y)
{
  float ang;

  /* if x is very close to zero, then return pi/2 or 3pi/2, depending
     on y-value */
  if (fabs(x) < 0.0000001)
    {
      if (y < 0.0)
        ang = 3.0 * PI / 2.0;
      else
        ang =       PI / 2.0;
    }

  /* If y is very close to zero, then return 0 or pi, depending on
     x-value */
  else if (fabs(y) < 0.0000001)
    {
      if (x < 0.0)
        ang = PI;
      else
        ang = 0.0;
    }

  /* Otherwise, calculate angle from arctan y/x */
  else
    {
      ang = atan((double) (y / x));
      if (ang > 0.0)
        {
          if (y < 0.0) 
            ang = ang + PI;
        }
      else
        {
          if (y < 0.0)
            ang = 2.0 * PI + ang;
          else
            ang =       PI + ang;
        }
    }
  return ang;
}
/***********************************************************************

calculate_matrix  -  Calculate the 4 x 3 transformation matrix from the
                     supplied translation and angles to rotate about
                     the z-x-z axes

***********************************************************************/

void calculate_matrix(float transl[3], float cos_theta[2], float sin_theta[2],
            float mat[4][3])
{
  float xcys, xsyc;

  /* Calculate each of the terms in the transformation matrix */
  xcys = (transl[0] * cos_theta[0] + transl[1] * sin_theta[0]);
  xsyc = (transl[0] * sin_theta[0] - transl[1] * cos_theta[0]);
  mat[0][0] =   cos_theta[0] * cos_theta[2]
              - sin_theta[0] * cos_theta[1] * sin_theta[2];
  mat[1][0] =   sin_theta[0] * cos_theta[2]
              + cos_theta[0] * cos_theta[1] * sin_theta[2];
  mat[2][0] =   sin_theta[1] * sin_theta[2];
  mat[3][0] = - xcys * cos_theta[2]
              + xsyc * cos_theta[1] * sin_theta[2]
              - transl[2] * sin_theta[1] * sin_theta[2];
  mat[0][1] = -(cos_theta[0] * sin_theta[2]
              + sin_theta[0] * cos_theta[1] * cos_theta[2]);
  mat[1][1] = -(sin_theta[0] * sin_theta[2]
              - cos_theta[0] * cos_theta[1] * cos_theta[2]);
  mat[2][1] =   sin_theta[1] * cos_theta[2];
  mat[3][1] =   xcys * sin_theta[2]
              + xsyc * cos_theta[1] * cos_theta[2]
              - transl[2] * sin_theta[1] * cos_theta[2];
  mat[0][2] =   sin_theta[0] * sin_theta[1];
  mat[1][2] = - cos_theta[0] * sin_theta[1];
  mat[2][2] =   cos_theta[1];
  mat[3][2] = - xsyc * sin_theta[1] - transl[2] * cos_theta[1];
}              
/***********************************************************************

line_transformation  -  Calculate the full transformation that will
                        put the given 2 atoms along the x-axis, with the
		        first atom at the origin and the other atom
		        along the +ve (or -ve) x-axis

***********************************************************************/

void line_transformation(float coord_store[3][3], float matrix[4][3],
			 float splay_angle)
{
  int   icoord, itheta;
  float  calc, cos_theta[3],
        nx, ny, nz, sin_theta[3], theta[3],
        transl[3], xyz[3][3];

  /* Initialise variables */
  for (itheta = 0; itheta < 3; itheta++)
    {
      theta[itheta] = 0.0;
      sin_theta[itheta] = 0.0;
      cos_theta[itheta] = 1.0;
    }

  /* Translation to origin is just the reverse of the coordinates of the
     first atom */
  transl[0] = coord_store[0][0];
  transl[1] = coord_store[0][1];
  transl[2] = coord_store[0][2];

  /* Translate the 2 atoms to the origin */
  for (icoord = 0; icoord < 3; icoord++)
    {
      xyz[0][icoord] = coord_store[0][icoord] - transl[icoord];
      xyz[1][icoord] = coord_store[1][icoord] - transl[icoord];
    }

  /* Calculate theta1 - the angle to rotate about the z-axis */
  nx = xyz[1][0];
  ny = xyz[1][1];
  nz = xyz[1][2];
  theta[0] = get_angle(ny,-nx);

  /* Calculate theta2 - the angle to rotate about the x-axis to put the
     line along the +ve y-axis */
  calc = nx * nx + ny * ny;
  theta[1] = - get_angle(nz,sqrt(calc));
  theta[1] = theta[1] + PI / 2.0;

  /* Third rotation is thro' 90 degrees about z-axis to put from +ve
     y-axis to +ve x-axis */
  theta[2] = PI / 2.0;

  /* If want the line pointing at some splay-angle round from the +ve
     x-axis, then add that to the final rotation */
  theta[2] = theta[2] - splay_angle;

  /* Precalculate sine and cosine values for theta1 and theta2 */
  cos_theta[0] = cos(theta[0]);
  cos_theta[1] = cos(theta[1]);
  sin_theta[0] = sin(theta[0]);
  sin_theta[1] = sin(theta[1]);
  if (fabs(cos_theta[0]) < 0.000001)
    cos_theta[0] = 0.0;
  if (fabs(cos_theta[1]) < 0.000001)
    cos_theta[1] = 0.0;
  if (fabs(sin_theta[0]) < 0.000001)
    sin_theta[0] = 0.0;
  if (fabs(sin_theta[1]) < 0.000001)
    sin_theta[1] = 0.0;

  /* Precalculate sine and cosine values for theta3 */
  cos_theta[2] = cos(theta[2]);
  sin_theta[2] = sin(theta[2]);
  if (fabs(cos_theta[2]) < 0.000001)
    cos_theta[2] = 0.0;
  if (fabs(sin_theta[2]) < 0.000001)
    sin_theta[2] = 0.0;

  /* Calculate the transformation matrix */
  calculate_matrix(transl,cos_theta,sin_theta,matrix);
}
/***********************************************************************

apply_transformation  -  Apply the tranformation to the supplied
                         coordinates

***********************************************************************/

void apply_transformation(float x,float y,float z,float *xnew,
                          float *ynew,float *znew,float mat[4][3])
{
  /* Apply the tranformation */
  *xnew = mat[0][0] * x + mat[1][0] * y + mat[2][0] * z + mat[3][0];
  *ynew = mat[0][1] * x + mat[1][1] * y + mat[2][1] * z + mat[3][1];
  *znew = mat[0][2] * x + mat[1][2] * y + mat[2][2] * z + mat[3][2];
}
/* v.3.1--> */
/***********************************************************************

apply_xy_transformation  -  Apply the tranformation to the supplied
                            x and y coordinates

***********************************************************************/

void apply_xy_transformation(float x,float y,float *xnew,float *ynew,
                          float mat[4][3])
{
  /* Apply the tranformation */
  *xnew = mat[0][0] * x + mat[1][0] * y + mat[3][0];
  *ynew = mat[0][1] * x + mat[1][1] * y + mat[3][1];
}
/* <--v.3.1 */
/***********************************************************************

transform_object  -  Apply the given transformation to the given object

***********************************************************************/

void transform_object(struct object *object_ptr,float matrix[4][3])
{
  int iatom, iresid, natoms, nresid;
  float x, y, z;
  float new_x, new_y, new_z;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, performing the
	 transformation */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* Extract this atom's coordinates */
	  x = atom_ptr->x;
	  y = atom_ptr->y;
	  z = atom_ptr->z;
	
	  /* Apply the transformation */
	  apply_transformation(x,y,z,&new_x,&new_y,&new_z,matrix);

	  /* Save the new coordinates */
	  atom_ptr->x = new_x;
	  atom_ptr->y = new_y;
	  atom_ptr->z = new_z;

	  /* Get pointer to the next atom */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

stretch_bond  -  Stretch the current bond, moving all atoms off one of
                 its ends by the required stretch distance

***********************************************************************/

void stretch_bond(struct bond *stretch_bond_ptr,float move_dist)
{
  int loop, nbonds, end;

  struct bond *bond_ptr, *next_free_stack_ptr, *other_bond_ptr,
  *next_stack_ptr;
  struct bond_link *bond_link_ptr;
  struct coordinate *atom_ptr;

  /* Initialise variables */
  next_stack_ptr = NULL;
  next_free_stack_ptr = NULL;

  /* Initialise all the bond flags and stack pointers */
  initialise_bonds();

  /* Initialise all the atom flags */
  initialise_atoms();

  /* If first atom has fewer bonds downstream of it, then want to
     process the bonds coming off it */
  if (stretch_bond_ptr->nbonds_from_first_atom
      < stretch_bond_ptr->nbonds_from_second_atom)
    {
      end = FIRST;
      move_dist = -1 * move_dist;
    }

  /* Otherwise, use the other end of the bond */
  else
    end = SECOND;

  /* Set pointer to the first of the bonds sprouting off the stretch
     bond */
  bond_link_ptr = stretch_bond_ptr->first_bond_link_ptr;

  /* Set the stretch-bond as already checked */
  stretch_bond_ptr->checked = TRUE;

  /* Initialise count of bonds encountered */
  nbonds = 0;

  /* Initialise stack pointers for bond-search by putting the bonds
     coming off the selected end of the stretch-bond onto the stack
     of bonds to be searched */
  while (bond_link_ptr != NULL)
    {
      /* Add this bond to the stack only if it comes off the relevant
	 atom of the test-bond */
      if (bond_link_ptr->bond_end == end)
	{
	  /* Get the bond this link is pointing to and add to stack */
	  bond_ptr = bond_link_ptr->bond_ptr;

	  /* Add to stack provided it is not an elastic bond */
	  if (bond_ptr->elastic == FALSE)
	    {
	      /* If this is the first to be added to the stack, set both
		 stack-pointers to point to it */
	      if (next_free_stack_ptr == NULL)
		{
		  next_free_stack_ptr = bond_ptr;
		  next_stack_ptr = bond_ptr;
		}

	      /* Otherwise, make last bond on stack point to this one
		 and set this one as nopw the last on the stack */
	      else
		{
		  next_free_stack_ptr->next_stack_ptr = bond_ptr;
		  next_free_stack_ptr = bond_ptr;
		}
	    }

	  /* Mark bond as added to the stack */
	  bond_ptr->checked = TRUE;
	}
      
      /* Go to the next bond-link in the linked list */
      bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
    }

  /* Loop until stack of pointers to be processed has been exhausted */
  while (next_stack_ptr != NULL)
    {
      /* Pick up the next bond from the stack, and move the stack
       pointer onto the next position */
      bond_ptr = next_stack_ptr;
      next_stack_ptr = bond_ptr->next_stack_ptr;

      /* Increment count of bonds encountered */
      nbonds++;

      /* Loop over the bond's two atoms */
      for (loop = 0; loop < 2; loop++)
	{
	  /* Get pointer to appropriate atom */
	  if (loop == 0)
	    atom_ptr = bond_ptr->first_atom_ptr;
	  else
	    atom_ptr = bond_ptr->second_atom_ptr;

	  /* If first atom has not been moved, then move it */
	  if (atom_ptr->checked == FALSE)
	    {
	      atom_ptr->x = atom_ptr->x + move_dist;
	      atom_ptr->checked = TRUE;
	    }
	}

     /* Loop through all the bonds connected to the current one
        to add to stack */

      /* Get the pointer to the first of the bonds coming off the 
	 current bond */
      bond_link_ptr = bond_ptr->first_bond_link_ptr;

      /* Process each of these bonds in turn */
      while (bond_link_ptr != NULL)
	{
	  /* Get the pointer to this bond */
	  other_bond_ptr = bond_link_ptr->bond_ptr;

	  /* If this bond hasn't already been processed, then move
	     its atoms (if they need moving) and add the bond to the
	     stack */
	  if (other_bond_ptr->checked == FALSE &&
	      other_bond_ptr->elastic == FALSE)
	    {
	      next_free_stack_ptr->next_stack_ptr = other_bond_ptr;
	      next_free_stack_ptr = other_bond_ptr;

	      /* If the stack has run out, restart it at the bond
		 just entered */
	      if (next_stack_ptr == NULL)
		next_stack_ptr = other_bond_ptr;

	      /* Mark bond as added to the stack */
	      other_bond_ptr->checked = TRUE;
	    }

	  /* Go to the next bond-link in the linked list */
	  bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
	}
    }
    
  /* Update the counter of bonds downstream on stretched bond */
  if (end == FIRST)
    stretch_bond_ptr->nbonds_from_first_atom = nbonds;
  else
    stretch_bond_ptr->nbonds_from_second_atom = nbonds;
}
/***********************************************************************

correct_bond_lengths  -  Correct bond-lengths in case any have become
                         deformed during the flattening of the ring groups

***********************************************************************/

void correct_bond_lengths(struct object *object_ptr)
{
  float bond_length, coord_store[3][3], diff, matrix[4][3], move_dist;
  float x1, x2, y1, y2, z1, z2;

  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct object_bond *object_bond_ptr;

  /* Get pointer to the first of this object's bonds */
  object_bond_ptr = object_ptr->first_object_bond_ptr;

  /* Loop through all object's bonds */
  while (object_bond_ptr != NULL)
    {
      /* Get the bond pointer */
      bond_ptr = object_bond_ptr->bond_ptr;

      /* If this is a rotatable bond, or an end bond, then process */
      if (bond_ptr->rotatable_bond == TRUE || bond_ptr->end_bond == TRUE)
	{
	  /* Get the pointers to the two atoms at either end of the bond */
	  atom1_ptr = bond_ptr->first_atom_ptr;
	  atom2_ptr = bond_ptr->second_atom_ptr;

	  /* Get the atom coordinates */
	  x1 = atom1_ptr->x;
	  y1 = atom1_ptr->y;
	  z1 = atom1_ptr->z;
	  x2 = atom2_ptr->x;
	  y2 = atom2_ptr->y;
	  z2 = atom2_ptr->z;

	  /* Calculate length of bond */
	  bond_length = sqrt((x2 - x1) * (x2 - x1)
			     + (y2 - y1) * (y2 - y1)
			     + (z2 - z1) * (z2 - z1));

	  /* If this bond is significantly different from its original
	     value, need to correct it */
	  diff = fabs(bond_length - bond_ptr->bond_length);
	  if (diff > 0.01)
	    {
	      /* Get the distance by which the bond needs to be stretched
		 or shrunk */
	      move_dist = bond_ptr->bond_length - bond_length;

	      /* Get the coordinates of the two atoms defining this bond */
	      coord_store[0][0] = x1;
	      coord_store[0][1] = y1;
	      coord_store[0][2] = z1;
	      coord_store[1][0] = x2;
	      coord_store[1][1] = y2;
	      coord_store[1][2] = z2;
	      
	      /* Calculate the transformation that will put this bond
		 along the x-axis, with atom 1 at the origin and atom 2
		 pointing in the +ve x-direction */
	      line_transformation(coord_store,matrix,0.0);

	      /* Apply the transformation to all atoms in the current
	         object */
	      transform_object(object_ptr,matrix);

	      /* If this is an end-bond, then only need to move a single
		 atom */
	      if (bond_ptr->end_bond == TRUE)
		{
		  /* If the bond's first atom is at the free end */
		  if (bond_ptr->nfirst_atom_links == 0)
		    {
		      /* Move first atom (which should currently be
		         at the origin) back in the -ve x-direction */
		      atom1_ptr->x = atom1_ptr->x - move_dist;
		    }

		  /* Otherwise, move the other atom */
		  else
		    {
		      /* Move second atom (which should currently be
		         on the +ve x-axis) in the +ve x-direction */
		      atom2_ptr->x = atom2_ptr->x + move_dist;
		    }
		}
	      /* Otherwise, if this is a rotatable bond, need to
		 move all atoms reachable from one of its ends
		 (the one with fewest connections) */
	      else
		stretch_bond(bond_ptr,move_dist);
	    }
	}
      
      /* Get pointer to point to next bond for this object */
      object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
    }
}
/***********************************************************************

plane_transformation - Calculate the full transformation that will
                       put the given 3 atoms in the x-y plane, with
		       the middle atom at the origin and the other
		       two atoms placed either side of the positive
		       x-axis s.t. the angle 1-2-3 is bisected by
		       that axis

***********************************************************************/

void plane_transformation(float coord[3][3], float matrix[4][3],
			  int in_neg_direction)
{
  int   icoord, itheta;
  float bisect[2], calc, cos_theta[3], mag1, mag2, normal[3],
        newx, newy, newz, nx, ny, nz, sin_theta[3], theta[3],
        transl[3], xyz[3][3];

  /* Initialise variables */
  for (itheta = 0; itheta < 3; itheta++)
    {
      theta[itheta] = 0.0;
      sin_theta[itheta] = 0.0;
      cos_theta[itheta] = 1.0;
    }

  /* Translation to origin is just the reverse of the coordinates of the
     middle atom */
  transl[0] = coord[1][0];
  transl[1] = coord[1][1];
  transl[2] = coord[1][2];

  /* Translate the 3 atoms to the origin */
  for (icoord = 0; icoord < 3; icoord++)
    {
      xyz[0][icoord] = coord[0][icoord] - transl[icoord];
      xyz[1][icoord] = coord[1][icoord] - transl[icoord];
      xyz[2][icoord] = coord[2][icoord] - transl[icoord];
    }

  /* Calculate normal to plane holding the 3 atoms */
  normal[0] = xyz[2][1] * xyz[0][2] - xyz[2][2] * xyz[0][1];
  normal[1] = xyz[2][2] * xyz[0][0] - xyz[2][0] * xyz[0][2];
  normal[2] = xyz[2][0] * xyz[0][1] - xyz[2][1] * xyz[0][0];
  nx = normal[0];
  ny = normal[1];
  nz = normal[2];

  /* Calculate theta1 - the angle to rotate about the z-axis */
  theta[0] = get_angle(ny,-nx);

  /* Calculate theta2 - the angle to rotate about the x-axis */
  calc = nx * nx + ny * ny;
  theta[1] = - get_angle(nz,sqrt(calc));

  /* Precalculate sine and cosine values for theta1 and theta2 */
  cos_theta[0] = cos(theta[0]);
  cos_theta[1] = cos(theta[1]);
  sin_theta[0] = sin(theta[0]);
  sin_theta[1] = sin(theta[1]);
  if (fabs(cos_theta[0]) < 0.000001)
    cos_theta[0] = 0.0;
  if (fabs(cos_theta[1]) < 0.000001)
    cos_theta[1] = 0.0;
  if (fabs(sin_theta[0]) < 0.000001)
    sin_theta[0] = 0.0;
  if (fabs(sin_theta[1]) < 0.000001)
    sin_theta[1] = 0.0;

  /* Calculate the transformation matrix */
  calculate_matrix(transl,cos_theta,sin_theta,matrix);

  /* Transform coordinates of atoms 1 and 3 */
  apply_transformation(coord[0][0],coord[0][1],coord[0][2],
                       &newx,&newy,&newz,matrix);
  xyz[0][0] = newx;
  xyz[0][1] = newy;
  xyz[0][2] = newz;
  apply_transformation(coord[2][0],coord[2][1],coord[2][2],
                       &newx,&newy,&newz,matrix);
  xyz[2][0] = newx;
  xyz[2][1] = newy;
  xyz[2][2] = newz;

  /* Calculate theta3 - the get_angle to rotate about the z-axis s.t.
     the get_angle between the fragment's atoms 1-2-3 is bisected by
     the +ve x-axis */

  /* Find vector in direction of the bisector */
  mag1 = xyz[0][0] * xyz[0][0] + xyz[0][1] * xyz[0][1];
  mag2 = xyz[2][0] * xyz[2][0] + xyz[2][1] * xyz[2][1];
  if (mag1 > 0.000001 && mag2 > 0.000001)
    {
      bisect[0] = xyz[0][0] / sqrt(mag1) + xyz[2][0] / sqrt(mag2);
      bisect[1] = xyz[0][1] / sqrt(mag1) + xyz[2][1] / sqrt(mag2);
    }
  else
    printf("+WARNING: Bond-length zero!!!\n");

  /* Calculate get_angle the bisector makes with the x-axis */
  theta[2] = get_angle(bisect[0],bisect[1]);

  /* If want the bisector pointing along the -ve, rather than the +ve,
     x-axis, then turn through a further 180 degrees */
  if (in_neg_direction == TRUE)
    theta[2] = theta[2] - PI;

  /* Precalculate sine and cosine values for theta3 */
  cos_theta[2] = cos(theta[2]);
  sin_theta[2] = sin(theta[2]);
  if (fabs(cos_theta[2]) < 0.000001)
    cos_theta[2] = 0.0;
  if (fabs(sin_theta[2]) < 0.000001)
    sin_theta[2] = 0.0;

  /* Calculate the transformation matrix */
  calculate_matrix(transl,cos_theta,sin_theta,matrix);
}
/***********************************************************************

transform_downstream_atoms  -  Apply given transformation to all the
                               atoms downstream of the given bond

***********************************************************************/

void transform_downstream_atoms(struct bond *reference_bond_ptr,int end,
				float matrix[4][3])
{
  int loop, nbonds;

  float newx, newy, newz;

  struct bond *bond_ptr, *next_free_stack_ptr, *other_bond_ptr,
  *next_stack_ptr;
  struct bond_link *bond_link_ptr;
  struct coordinate *atom_ptr;

  /* Initialise variables */
  next_stack_ptr = NULL;
  next_free_stack_ptr = NULL;

  /* Initialise all the bond flags and stack pointers */
  initialise_bonds();

  /* Initialise all the atom flags */
  initialise_atoms();

  /* Mark both this bond's atoms as already checked */
  reference_bond_ptr->first_atom_ptr->checked = TRUE;
  reference_bond_ptr->second_atom_ptr->checked = TRUE;

  /* Set pointer to the first of the bonds sprouting off the reference
     bond */
  bond_link_ptr = reference_bond_ptr->first_bond_link_ptr;

  /* Set the reference-bond as already checked */
  reference_bond_ptr->checked = TRUE;

  /* Initialise count of bonds encountered */
  nbonds = 0;

  /* Initialise stack pointers for bond-search by putting the bonds
     coming off the selected end of the reference-bond onto the stack
     of bonds to be searched */
  while (bond_link_ptr != NULL)
    {
      /* Add this bond to the stack only if it comes off the relevant
	 end */
      if (bond_link_ptr->bond_end == end)
	{
	  /* Get the bond this link is pointing to and add to stack */
	  bond_ptr = bond_link_ptr->bond_ptr;

	  /* Add to stack provided it is not an elastic bond */
	  if (bond_ptr->elastic == FALSE)
	    {
	      /* If this is the first to be added to the stack, set both
		 stack-pointers to point to it */
	      if (next_free_stack_ptr == NULL)
		{
		  next_free_stack_ptr = bond_ptr;
		  next_stack_ptr = bond_ptr;
		}

	      /* Otherwise, make last bond on stack point to this one
		 and set this one as nopw the last on the stack */
	      else
		{
		  next_free_stack_ptr->next_stack_ptr = bond_ptr;
		  next_free_stack_ptr = bond_ptr;
		}
	    }

	  /* Mark bond as added to the stack */
	  bond_ptr->checked = TRUE;
	}
      
      /* Go to the next bond-link in the linked list */
      bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
    }

  /* Loop until stack of pointers to be processed has been exhausted */
  while (next_stack_ptr != NULL)
    {
      /* Pick up the next bond from the stack, and move the stack
       pointer onto the next position */
      bond_ptr = next_stack_ptr;
      next_stack_ptr = bond_ptr->next_stack_ptr;

      /* Increment count of bonds encountered */
      nbonds++;

      /* Loop over the bond's two atoms */
      for (loop = 0; loop < 2; loop++)
	{
	  /* Get pointer to appropriate atom */
	  if (loop == 0)
	    atom_ptr = bond_ptr->first_atom_ptr;
	  else
	    atom_ptr = bond_ptr->second_atom_ptr;

	  /* If first atom has not been transformed, then apply the
	     matrix transformation */
	  if (atom_ptr->checked == FALSE)
	    {
	      /* Apply the transformation to the current atom */
	      apply_transformation(atom_ptr->x,atom_ptr->y,atom_ptr->z,
				   &newx,&newy,&newz,matrix);
	      atom_ptr->x = newx;
	      atom_ptr->y = newy;
	      atom_ptr->z = newz;
	      atom_ptr->checked = TRUE;
	    }
	}

     /* Loop through all the bonds connected to the current one
        to add to stack */

      /* Get the pointer to the first of the bonds coming off the 
	 current bond */
      bond_link_ptr = bond_ptr->first_bond_link_ptr;

      /* Process each of these bonds in turn */
      while (bond_link_ptr != NULL)
	{
	  /* Get the pointer to this bond */
	  other_bond_ptr = bond_link_ptr->bond_ptr;

	  /* If this bond hasn't already been processed, then move
	     its atoms (if they need moving) and add the bond to the
	     stack */
	  if (other_bond_ptr->checked == FALSE &&
	      other_bond_ptr->elastic == FALSE)
	    {
	      next_free_stack_ptr->next_stack_ptr = other_bond_ptr;
	      next_free_stack_ptr = other_bond_ptr;

	      /* If the stack has run out, restart it at the bond
		 just entered */
	      if (next_stack_ptr == NULL)
		next_stack_ptr = other_bond_ptr;

	      /* Mark bond as added to the stack */
	      other_bond_ptr->checked = TRUE;
	    }

	  /* Go to the next bond-link in the linked list */
	  bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
	}
    }
}
/***********************************************************************

offshoot_splay  -  Check bonds springing off each ring to set them at
                   appropriate splay angles away from the ring

***********************************************************************/

void offshoot_splay(struct object *object_ptr,struct bond *ring_bond_ptr,
		    struct bond *other_ring_bond_ptr,
		    struct bond *connect_bond_ptr[CONNECTIONS],
		    int nstore,int end)
{
  int istore;

  float coord_store[3][3], matrix[4][3], splay_angle;
  float newx, newy, newz;
  float splay_start, splay_diff;

  struct coordinate *atom_ptr, *atom1_ptr, *atom2_ptr, *common_atom_ptr;
  struct bond *bond_ptr;

/*    int counter, istore, nzeroz,
    stored_coords, ibond;
    int other_atom[CONNECTIONS + 1], other_bond[CONNECTIONS + 1],
    other_sign[CONNECTIONS + 1];
    struct link_bond *other_bond_ptr; */

  /* Get pointers to the two atoms of the ring-bond of interest */
  if (end == FIRST)
    {
      atom1_ptr = ring_bond_ptr->first_atom_ptr;
      atom2_ptr = ring_bond_ptr->second_atom_ptr;
    }
  else
    {
      atom1_ptr = ring_bond_ptr->second_atom_ptr;
      atom2_ptr = ring_bond_ptr->first_atom_ptr;
    }

  /* Store the pointer to the atom that is common to both ring bonds
     and hence is the atom about which all the transformations will be
     performed */
  common_atom_ptr = atom1_ptr;

  /* Store the atom coordinates, with the common atom between the
     two ring bonds as the second (middle) atom */
  coord_store[0][0] = atom2_ptr->x;
  coord_store[0][1] = atom2_ptr->y;
  coord_store[0][2] = atom2_ptr->z;
  coord_store[1][0] = atom1_ptr->x;
  coord_store[1][1] = atom1_ptr->y;
  coord_store[1][2] = atom1_ptr->z;

  /* Get the pointer to whichever of the other ring-bond's atoms is
     the third one of interest */
  if (other_ring_bond_ptr->first_atom_ptr == common_atom_ptr)
    atom_ptr = other_ring_bond_ptr->second_atom_ptr;
  else
    atom_ptr = other_ring_bond_ptr->first_atom_ptr;

  /* Store the coordinates of this third atom */
  coord_store[2][0] = atom_ptr->x;
  coord_store[2][1] = atom_ptr->y;
  coord_store[2][2] = atom_ptr->z;

  /* Calculate transformation that will place atoms 1 and 3 in the
     x-y plane, with their bisector pointing in the -ve x-direction */
  plane_transformation(coord_store,matrix,TRUE);

  /* Apply the transformation to all atoms in the object */
  transform_object(object_ptr,matrix);

  /* Calculate the splay-angle for the non-ring bonds */
  if (nstore == 1)
    {
      splay_start = 0.0;
      splay_diff = 0.0;
    }
  else if (nstore > 1)
    {
      splay_start = - PI / (nstore + 1);
      splay_diff = - 2.0 * splay_start;
    }

  /* Initialise values for splaying out the non-ring bonds */
  coord_store[0][0] = 0.0;
  coord_store[0][1] = 0.0;
  coord_store[0][2] = 0.0;

  /* Loop over all the stored non-ring bonds that are to be splayed */
  for (istore = 0; istore < nstore; istore++)
    {
      /* Get pointer to this bond */
      bond_ptr = connect_bond_ptr[istore];

      /* Get pointers to the bond's two atoms */
      atom1_ptr = bond_ptr->first_atom_ptr;
      atom2_ptr = bond_ptr->second_atom_ptr;

      /* Determine which is at the other end of the bond */
      if (bond_ptr->first_atom_ptr == common_atom_ptr)
	{
	  atom_ptr = bond_ptr->second_atom_ptr;
	  end = SECOND;
	}
      else
	{
	  atom_ptr = bond_ptr->first_atom_ptr;
	  end = FIRST;
	}

      /* Store ths atom's coords */
      coord_store[1][0] = atom_ptr->x;
      coord_store[1][1] = atom_ptr->y;
      coord_store[1][2] = atom_ptr->z;

      /* Calculate the transformation that will place the current
	 bond along the +ve x-axis */
      splay_angle = splay_start + istore * splay_diff;

      line_transformation(coord_store,matrix,splay_angle);

      /* Apply the transformation to the current atom */
      apply_transformation(atom_ptr->x,atom_ptr->y,atom_ptr->z,
			   &newx,&newy,&newz,matrix);
      atom_ptr->x = newx;
      atom_ptr->y = newy;
      atom_ptr->z = newz;

      /* If the current bond is not an end-bond, then need to apply
	 the transformation to all the atoms downstream of it */
      if (bond_ptr->end_bond == FALSE)
	{
	  /* Transform all atoms downstream of the current bond */
	  transform_downstream_atoms(bond_ptr,end,matrix);
	}
    }

}
/***********************************************************************

flatten_ring_offshoots  -  Check bonds springing off each ring in the
                           current object to set them at appropriate
                           splay angles away from the ring

***********************************************************************/

void flatten_ring_offshoots(struct object *object_ptr)
{
  int end, loop;
  int non_ring_bonds, nring_bonds, nstore;

  struct bond *bond_ptr, *other_bond_ptr, *other_ring_bond_ptr,
  *connect_bond_ptr[CONNECTIONS];
  struct bond_link *bond_link_ptr;
  struct object_bond *object_bond_ptr;

  /* Get pointer to the first of this object's bonds */
  object_bond_ptr = object_ptr->first_object_bond_ptr;

  /* Loop through all object's bonds */
  while (object_bond_ptr != NULL)
    {
      /* Get the bond pointer */
      bond_ptr = object_bond_ptr->bond_ptr;

      /* If this is a ring bond, then process */
      if (bond_ptr->ring_bond == TRUE)
	{
	  /* Loop over the bond's two atoms */
	  for (loop = 0; loop < 2; loop++)
	    {
	      /* Get pointer to appropriate atom */
	      if (loop == 0)
		end = FIRST;
	      else
		end = SECOND;

	      /* See how many bonds, and of what types, are attached to
		 this atom */

	      /* Initalise counts and pointers */
	      nstore = 0;
	      non_ring_bonds = 0;
	      nring_bonds = 0;
	      other_ring_bond_ptr = NULL;
	      
	      /* Get the pointer to the first of the bonds coming off the 
		 current bond */
	      bond_link_ptr = bond_ptr->first_bond_link_ptr;
	      
	      /* Check the bond-types of all the bonds coming off the first
		 atom*/
	      while (bond_link_ptr != NULL)
		{
		  /* If it comes off the first atom the current bond, then
		     check its type */
		  if (bond_link_ptr->bond_end == end)
		    {
		      /* Get pointer to the other bond and store */
		      other_bond_ptr = bond_link_ptr->bond_ptr;
		      
		      /* Increment counts of ring- and non-ring bonds, and
			 store pointers to these bonds */
		      if (other_bond_ptr->ring_bond == TRUE)
			{
			  other_ring_bond_ptr = other_bond_ptr;
			  nring_bonds++;
			}
		      else
			{
			  non_ring_bonds++;
			  connect_bond_ptr[nstore] = other_bond_ptr;
			  nstore++;
			}
		    }
	      
		  /* Go to the next bond-link in the linked list */
		  bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
		}

	      /* If have both ring-bonds and non-ring bonds at this
		 junction, need to splay the non-ring bonds */
	      if (nring_bonds == 1 && non_ring_bonds > 0)
		{
		  /* Perform the splay */
		  offshoot_splay(object_ptr,bond_ptr,other_ring_bond_ptr,
				 connect_bond_ptr,nstore,end);
		}
	      else if (nring_bonds > 1 && non_ring_bonds > 0)
		{
		  printf("*** WARNING. %d ring bonds and %d ",
			 nring_bonds + 1,non_ring_bonds);
		  printf(" non-ring bonds. May cause problems!\n");
		  Nwarnings++;
		}
	    }
	}

      /* Get pointer to this object's next bond entry */
      object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
    }		
}
/***********************************************************************

flatten_bonds  -  Flatten out the bonds on the given side of the current
                  bond so that they lie in the x-y plane

***********************************************************************/

void flatten_bonds(struct bond *current_bond_ptr,int end)
{
  int have_ring, have_nonring, ibond, nstored, other_end;

  float coord_store[3][3], stored_coords[CONNECTIONS][3], matrix[4][3];
  float max_z, splay_angle;
  float new_x, new_y, new_z, x, y, z;

  struct bond *bond_ptr, *stored_bond_ptr[CONNECTIONS];
  struct bond_link *bond_link_ptr;
  struct coordinate *atom_ptr, *common_atom_ptr;

  /* Initialise variables */
  have_ring = FALSE;
  have_nonring = FALSE;
  nstored = 0;
  max_z = 0.0;

  /* Get pointer to the pivot atom (currently at the origin) */
  if (end == FIRST)
    common_atom_ptr = current_bond_ptr->first_atom_ptr;
  else
    common_atom_ptr = current_bond_ptr->second_atom_ptr;

  /* Get the pointer to the first of the bonds coming off the 
     test_bond */
  bond_link_ptr = current_bond_ptr->first_bond_link_ptr;

  /* Loop over all the bonds attached to the current one at this end */
  while (bond_link_ptr != NULL)
    {
      /* Check this bond comes off the appropriate end */
      if (bond_link_ptr->bond_end == end)
	{
	  /* Get pointer to this bond */
	  bond_ptr = bond_link_ptr->bond_ptr;

	  /* Check if this is a ring bond. If it is, then don't
	     want to apply the splaying */
	  if (bond_ptr->ring_bond == TRUE)
	    have_ring = TRUE;
	  else
	    have_nonring = TRUE;

	  /* Get the pointer to the atom at the other end of the bond */
	  if (bond_ptr->first_atom_ptr == common_atom_ptr)
	    atom_ptr = bond_ptr->second_atom_ptr;
	  else
	    atom_ptr = bond_ptr->first_atom_ptr;

	  /* Store ths atom's coords */
	  if (nstored < CONNECTIONS)
	    {
	      stored_bond_ptr[nstored] = bond_ptr;
	      stored_coords[nstored][0] = atom_ptr->x;
	      stored_coords[nstored][1] = atom_ptr->y;
	      stored_coords[nstored][2] = atom_ptr->z;
	      nstored++;
	    }

	  /* Get maximum z-deviation */
	  z = atom_ptr->z;
	  if (fabs(z) > max_z)
	    max_z = fabs(z);
	}

      /* Go to the next bond-link in the linked list */
      bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
    }

  /* If the bonds require flattening, then do it */
  if (max_z > 0.001)
    {
      /* If have both a ring, and other bonds springing off this end,
	 then flatten as for one bond */
      if (have_ring == TRUE && have_nonring == TRUE)
	nstored = 1;

      /* If have just a ring, then flatten as for one bond */
      if (have_ring == TRUE)
	nstored = 1;
    
      /* If there is just the one bond coming off this one, swivel
         it, and all bonds beyond it, into the x-y plane */
      if (nstored == 1)
	{
	  /* Create two dummy atoms, one at the origin, and the other
	     having y- and z-coords opposite in sign to the first atom
	     stored */
	  coord_store[0][0] = stored_coords[0][0];
	  coord_store[0][1] = stored_coords[0][1];
	  coord_store[0][2] = stored_coords[0][2];
	  coord_store[1][0] = 0.0;
	  coord_store[1][1] = 0.0;
	  coord_store[1][2] = 0.0;
	  coord_store[2][0] = stored_coords[0][0];
	  coord_store[2][1] = - stored_coords[0][1];
	  coord_store[2][2] = - stored_coords[0][2];

	  /* Calculate transformation that will this place bond (and its
	     mirror image in the x-y plane, with their bisector pointing
	     in the +ve x-direction */
	  plane_transformation(coord_store,matrix,FALSE);
	      
	  /* Apply the transformation to all atoms springing off
	     the current bond */
	  transform_downstream_atoms(current_bond_ptr,end,matrix);
	}

  /* If there are two bonds coming off this one, then calculate the
     transformation that will swivel the plane defined by the two bonds
     into the x-y plane */
      else if (nstored == 2)
	{
	  /* Move the first of the atoms to position 1, and replace
	     position two with a dummy atom at the origin */
	  coord_store[0][0] = stored_coords[0][0];
	  coord_store[0][1] = stored_coords[0][1];
	  coord_store[0][2] = stored_coords[0][2];
	  coord_store[1][0] = 0.0;
	  coord_store[1][1] = 0.0;
	  coord_store[1][2] = 0.0;
	  coord_store[2][0] = stored_coords[1][0];
	  coord_store[2][1] = stored_coords[1][1];
	  coord_store[2][2] = stored_coords[1][2];

	  /* Calculate transformation that will this place atoms 1 and 3
	     in the x-y plane, with their bisector pointing in the +ve
	     x-direction */
	  plane_transformation(coord_store,matrix,FALSE);
	      
	  /* Apply the transformation to all atoms springing off
	     the current bond */
	  transform_downstream_atoms(current_bond_ptr,end,matrix);
	}

      /* If there are more than two bonds coming off this one, then need
	 to treat each one separately, splaying them out in the x-y
	 plane at regularly-spaced angular intervals */
      else if (nstored > 2)
	{
	  /* Loop through each of the bonds in turn */
	  for (ibond = 0; ibond < nstored; ibond++)
	    {
	      /* Get the coordinates of this bond */
	      coord_store[0][0] = 0.0;
	      coord_store[0][1] = 0.0;
	      coord_store[0][2] = 0.0;
	      coord_store[1][0] = stored_coords[ibond][0];
	      coord_store[1][1] = stored_coords[ibond][1];
	      coord_store[1][2] = stored_coords[ibond][2];

	      /* Calculate the splay-angle round from x-axis that want
		 this bond to be placed at in the x-y plane */
	      splay_angle = PI - (ibond + 1) * (2.0 * PI / (nstored + 1));
	      
	      /* Calculate the transformation that will place the current
		 bond at the correct splay angle */
	      line_transformation(coord_store,matrix,splay_angle);

	      /* Get the bond to apply the transformation to */
	      bond_ptr = stored_bond_ptr[ibond];

	      /* Get the pointer to the atom at the other end of the bond */
	      if (bond_ptr->first_atom_ptr == common_atom_ptr)
		{
		  atom_ptr = bond_ptr->second_atom_ptr;
		  other_end = SECOND;
		}
	      else
		{
		  atom_ptr = bond_ptr->first_atom_ptr;
		  other_end = FIRST;
		}

	      /* Extract this atom's coordinates */
	      x = atom_ptr->x;
	      y = atom_ptr->y;
	      z = atom_ptr->z;
	
	      /* Apply the transformation */
	      apply_transformation(x,y,z,&new_x,&new_y,&new_z,matrix);

	      /* Save the new coordinates */
	      atom_ptr->x = new_x;
	      atom_ptr->y = new_y;
	      atom_ptr->z = new_z;

	      /* Apply the transformation to all atoms springing off
		 the current bond */
	      transform_downstream_atoms(bond_ptr,other_end,matrix);
	    }
	}
    }
}
/***********************************************************************

flip_object  -  Flip the current object about the x- or y-axis

***********************************************************************/

void flip_object(struct object *object_ptr,int xflip,int yflip)
{
  int iatom, iresid, natoms, nresid;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, performing the
	 transformation */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* Apply the transformation */
	  atom_ptr->x = xflip * atom_ptr->x;
	  atom_ptr->y = yflip * atom_ptr->y;

	  /* Mark atom as being done */
	  atom_ptr->checked = TRUE;

	  /* Get pointer to the next atom */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

flatten_object  -  Transform each rotatable bond in the current object
                   to the origin and flatten bonds either side of it

***********************************************************************/

void flatten_object(struct object *object_ptr)
{
  float coord_store[3][3], matrix[4][3];

  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct object_bond *object_bond_ptr;


  /* Get pointer to the first of this object's bonds */
  object_bond_ptr = object_ptr->first_object_bond_ptr;

  /* Loop through all object's bonds */
  while (object_bond_ptr != NULL)
    {
      /* Get the bond pointer */
      bond_ptr = object_bond_ptr->bond_ptr;

      /* If this is a rotatable bond perform unfolding of the object
	 around this bond */
      if (bond_ptr->rotatable_bond == TRUE)
	{
	  /* Get pointers to the bond's two atoms */
	  atom1_ptr = bond_ptr->first_atom_ptr;
	  atom2_ptr = bond_ptr->second_atom_ptr;

	  /* Store the two atoms' coordinates */
	  coord_store[0][0] = atom1_ptr->x;
	  coord_store[0][1] = atom1_ptr->y;
	  coord_store[0][2] = atom1_ptr->z;
	  coord_store[1][0] = atom2_ptr->x;
	  coord_store[1][1] = atom2_ptr->y;
	  coord_store[1][2] = atom2_ptr->z;

	  /* Calculate the transformation that will put this bond along
	     the x-axis, with atom 1 at the origin and atom 2 pointing
	     in the -ve x-direction */
	  line_transformation(coord_store,matrix,-PI);

	  /* Apply the transformation to all atoms in the current
	     object */
	  transform_object(object_ptr,matrix);

	  /* Perform flattening of bonds attached to atom 1 */
	  flatten_bonds(bond_ptr,FIRST);
	  
	  /* Shift the entire structure onto second of this bond's
	     two atoms and flip right over about the y-axis */
	  move_object(object_ptr,atom2_ptr->x,atom2_ptr->y,
				  atom2_ptr->z);
	  flip_object(object_ptr,-1,1);
	  
	  /* Perform flattening of bonds attached to atom 2 */
	  flatten_bonds(bond_ptr,SECOND);
	}
      
      /* Get pointer to this object's next bond entry */
      object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
    }
}
/***********************************************************************

unroll_structure  -  Unroll the structure of corrent object until
                     z-coords have been flattened completely

***********************************************************************/

void unroll_structure(struct object *object_ptr)
{
  int flat, flat_loop, loop;
  int natoms;
  float centre_x, centre_y, centre_z, max_off, x, y, z;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr;
  struct object_bond *object_bond_ptr;


  /* Initialise variables */
  flat = FALSE;

  /* Loop until object is reasonably flat */
  for (flat_loop = 0; flat_loop < MAX_FLAT_LOOP && flat == FALSE;
       flat_loop++)
    {
      /* Initialise variables for this loop */
      centre_x = 0.0;
      centre_y = 0.0;
      centre_z = 0.0;
      natoms = 0;
      flat = TRUE;
      max_off = 0.0;

      /* Initialise all the atom flags */
      initialise_atoms();

      /* Get pointer to the first of this object's bonds */
      object_bond_ptr = object_ptr->first_object_bond_ptr;

      /* Loop through all object's bonds */
      while (object_bond_ptr != NULL)
	{
	  /* Get the bond pointer */
	  bond_ptr = object_bond_ptr->bond_ptr;

	  /* Loop over the bond's two atoms */
	  for (loop = 0; loop < 2; loop++)
	    {
	      /* Get pointer to appropriate atom */
	      if (loop == 0)
		atom_ptr = bond_ptr->first_atom_ptr;
	      else
		atom_ptr = bond_ptr->second_atom_ptr;

	      /* Process only if this atom hasn't already been checked */
	      if (atom_ptr->checked == FALSE)
		{
		  /* Get the coordinates of this atom */
		  x = atom_ptr->x;
		  y = atom_ptr->y;
		  z = atom_ptr->z;

		  /* Accumulate centre-of-mass data */
		  centre_x = centre_x + x;
		  centre_y = centre_y + y;
		  centre_z = centre_z + z;
		  natoms++;

		  /* If the z-coordinate is close to zero, set it
		     to exactly zero */
		  if (fabs(z) < 0.001)
		    atom_ptr->z = 0.0;

		  /* Otherwise, object is not flat, so check if this is
		     the highest deviant from the x-y plane so far */
		  else
		    {
		      flat = FALSE;
		      if (fabs(z) > max_off)
			max_off = fabs(z);
		    }

		  /* Set this atom's flag as checked */
		  atom_ptr->checked = TRUE;
		}
	    }
      
	  /* Get pointer to this object's next bond entry */
	  object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
	}

      /* If the object needs flattening, then do it */
      if (flat == FALSE)
	{
	  /* If object contains rotatable bonds, then can unroll it */
	  if (object_ptr->nrot_bonds > 0)
	    flatten_object(object_ptr);

	  /* Otherwise, just move object to origin */
	  else
	    if (natoms > 0)
	      {
		centre_x = centre_x / natoms;
		centre_y = centre_y / natoms;
		centre_z = centre_z / natoms;
		move_object(object_ptr,centre_x,centre_y,centre_z);
	      }
	}
    }
}
/***********************************************************************

adjust_to_origin  -  Translate all the atoms by the given amount

***********************************************************************/

void adjust_to_origin(float centre_x,float centre_y,float centre_z)
{
  struct coordinate *atom_ptr;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms to initialise the flag indicating whether
     coords have already been stored */
  while (atom_ptr != NULL)
    {
      /* Apply the transformation */
      atom_ptr->x = atom_ptr->x - centre_x;
      atom_ptr->y = atom_ptr->y - centre_y;
      atom_ptr->z = atom_ptr->z - centre_z;
      
      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next;
    }	       
}
/***********************************************************************

align_object  -  Align object's principal axes with the x-y plane

***********************************************************************/

void align_object(struct object *object_ptr)
{
/* v.4.0--> */
/*  int iobject, iatom, iresid, natoms, nresid, nstored;
  int wanted; */
  int iatom, iresid, natoms, nresid, nstored;
/* <--v.4.0 */

  float coord_store[MAXATOMS][3], mass_x, mass_y, mass_z;
/* v.4.0--> */
/*  float new_x, new_y, new_z, x, y, z;
  float x, y, z; */
/* <--v.4.0 */
  float transformation[3][3];

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  nstored = 0;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to first residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* Store this atom's fitted coords */
	  if (nstored < MAXATOMS)
	    {
	      coord_store[nstored][0] = atom_ptr->x;
	      coord_store[nstored][1] = atom_ptr->y;
	      coord_store[nstored][2] = atom_ptr->z;
	      nstored++;
	    }

	  /* Get pointer to the next atom */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Calculate the centre of mass of these atoms */
  centre_of_mass(coord_store,nstored,&mass_x,&mass_y,&mass_z);

  /* Place the centre of mass at the origin */
  adjust_to_origin(mass_x,mass_y,mass_z);

  /* Apply the same translation to the coordinates stored in
     coord_store */
  adjust_stored_coords(nstored,coord_store,mass_x,mass_y,mass_z);

  /* Calculate transformation to orient the principal axis along
     the x-direction */
  principal_components(nstored,transformation,coord_store);

  /* Apply this transformation to the entire object */
  rotate_object(object_ptr,transformation);
}
/***********************************************************************

force_flat  -  Zero all the z-coordinates (whether or not the structure
               has been properly flattened) so that the flipping
	       operations aren't thrown out by any non-planarity in
	       the structure

***********************************************************************/

void force_flat(struct object *object_ptr)
{
  int iatom, iresid, natoms, nresid;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, performing the
	 zeroing the z-coords */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* Set the z-coord to zero */
	  atom_ptr->z = 0.0;
      
	  /* Get pointer to the next atom */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

create_simple_hgroup  -  Adjust all the non-ligand residues for the
                         schematic plot, marking any unwanted atoms
			 for deletion

***********************************************************************/

void create_simple_hgroup(struct object *object_ptr)
{
  int iatom, iresid, natoms, nresid;
  float hgroup_size;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise all atoms */
  initialise_atoms();

  /* Calculate size for the simple H-group */
  hgroup_size = 5.0 * Text_Size_Val->Nonligand_Residue_Names * CHAR_ASPECT;

  /* Initialise bond pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds */
  while (bond_ptr != NULL)
    {
      /* Process only if this is an H-bond */
      if (bond_ptr->bond_type == HBOND)
	{
	  /* Mark both atoms as being involved in the bond */
	  bond_ptr->first_atom_ptr->checked = TRUE;
	  bond_ptr->second_atom_ptr->checked = TRUE;
	}

      /* Get pointer for next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Set object type to simplified residue */
  object_ptr->object_type = SIMPLE_HGROUP;

  /* Adjust object- and residue- maximum atom sizes */
  object_ptr->max_atom_size = hgroup_size;
  residue_ptr->max_atom_size = hgroup_size;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
      
      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, marking any not involved
	 in H-bonds for deletion */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* If atom not involved in any interactions with other
	     residues, then can mark it for deletion */
	  if (atom_ptr->checked == FALSE)
	    atom_ptr->deleted = TRUE;

	  /* Get pointer to the next atom */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

remove_mainchains  -  Remove any non-ligand mainchain atoms if mainchain
                      not involved in any interactions

***********************************************************************/

void remove_mainchains(struct object *object_ptr)
{
  int iatom, iresid, natoms, nresid;
  int ninteract;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise all atoms */
  initialise_atoms();

  /* Initialise bond pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds */
  while (bond_ptr != NULL)
    {
      /* Process only if this is an H-bond or a contact bond */
      if (bond_ptr->bond_type == HBOND || bond_ptr->bond_type == CONTACT)
	{
	  /* Mark both atoms as being involved in the bond */
	  bond_ptr->first_atom_ptr->checked = TRUE;
	  bond_ptr->second_atom_ptr->checked = TRUE;
	}

      /* Get pointer for next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Initialise count of mainchain atoms involved in interactions */
      ninteract = 0;

      /* Loop over all the residue's atoms counting the number of
         mainchain atoms involved in interactions */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* If this is a mainchain atom and it is involved in
	     an interaction, then increment count */
	  if (atom_ptr->side_chain == FALSE && atom_ptr->checked == TRUE)
	    ninteract++;

	  /* Get pointer to the next atom */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

      /* If no interactions with any of the mainchain atoms, then
	 delete them all */
      if (ninteract == 0)
	{
	  /* Re-initialise to this residue's first atom */
	  atom_ptr = residue_ptr->first_atom_ptr;
	  iatom = 0;

	  /* Loop over all the residue's atoms, marking mainchain
	     atoms for deletion */
	  while (iatom < natoms && atom_ptr != NULL)
	    {
	      /* If this is a mainchain atom, mark it for deletion */
	      if (atom_ptr->side_chain == FALSE)
		atom_ptr->deleted = TRUE;

	      /* Get pointer to the next atom */
	      atom_ptr = atom_ptr->next;
	      iatom++;
	    }
	}

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

mark_sidechain_atoms  -  Mark all the current residue's sidechain atoms
                         for deletion

***********************************************************************/

void mark_sidechain_atoms(struct residue *residue_ptr)
{
  int iatom, natoms;

  struct coordinate *atom_ptr;

  /* Get pointer to this residue's first atom */
  atom_ptr = residue_ptr->first_atom_ptr;

  /* Initialise atom count */
  iatom = 0;
  natoms = residue_ptr->natoms;

  /* Loop over all the residue's atoms, to mark all the sidechain atoms  */
  while (iatom < natoms && atom_ptr != NULL)
    {
      /* Initialise the side-chain flag */
      if (atom_ptr->side_chain == TRUE)
	atom_ptr->deleted = TRUE;

      /* Get pointer to the next atom in this residue */
      atom_ptr = atom_ptr->next;
      iatom++;
    }
}
/* v.3.1--> */
/***********************************************************************

transform_downlist_atoms  -  Apply the given transformation to all the
                             atoms in the list from the given CA atom
                             onwards by finding the first downlist bond

***********************************************************************/

void transform_downlist_atoms(struct object *object_ptr,
			      struct coordinate *current_atom_ptr,
			      float matrix[4][3])
{
  int end, iatom, iresid, natoms, nresid;
  int got_bond, transform, wanted;
/* v.3.1.1--> */
/*  float x, y, z;
  float new_x, new_y, new_z; */
/* <--v.3.1.1 */

  struct bond *bond_ptr, *use_bond_ptr;
  struct coordinate *atom_ptr;
  struct residue *current_residue_ptr, *residue_ptr;

  /* Initialise variables */
  got_bond = FALSE;
  transform = FALSE;
  use_bond_ptr = NULL;
  current_residue_ptr = current_atom_ptr->residue_ptr;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL && got_bond == FALSE)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, performing the
	 transformation where required */
      while (iatom < natoms && atom_ptr != NULL && got_bond == FALSE)
	{
	  /* Check whether this is the starting atom */
	  if (atom_ptr == current_atom_ptr)
	    transform = TRUE;

	  /* Check if this is a sidechain atom of the starting
	     residue */
	  wanted = TRUE;
	  if (residue_ptr == current_residue_ptr &&
	      atom_ptr->side_chain == TRUE)
	    wanted = FALSE;

	  /* If have a potential starting bond, then try to find it */
	  if (transform == TRUE && wanted == TRUE)
	    {
	      /* Initialise to first bond */
	      bond_ptr = first_bond_ptr;

	      /* Loop through all the bonds to find a covalent bond
		 between the two atoms */
	      while (bond_ptr != NULL && got_bond == FALSE)
		{
		  /* Check if this is a covalent bond between the 
		     two atoms */
		  if (bond_ptr->bond_type == COVALENT)
		    {
		      /* Atoms in bond, with the far end the one
			 we're after */
		      if (bond_ptr->first_atom_ptr == current_atom_ptr &&
			  bond_ptr->second_atom_ptr == atom_ptr)
			{
			  got_bond = TRUE;
			  use_bond_ptr = bond_ptr;
			  end = SECOND;
			}

		      /* Atoms in bond, with the far end the one
			 we're after */
 		      else if (bond_ptr->first_atom_ptr == atom_ptr &&
			       bond_ptr->second_atom_ptr == current_atom_ptr)
			{
			  got_bond = TRUE;
			  use_bond_ptr = bond_ptr;
			  end = FIRST;
			}
		    }

		  /* Go to the next bond in the linked list */
		  bond_ptr = bond_ptr->next_bond_ptr;
		}
	    }

	  /* Get pointer to the next atom */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* If have a bond that can be used for the transformation, use it */
  if (use_bond_ptr != NULL)
    transform_downstream_atoms(use_bond_ptr,end,matrix);
}
/* <--v.3.1 */
/***********************************************************************

stretch_mainchain  -  Stretch the mainchain of the current residue out
                      into a straight line for the simplified ligand
		      representation

***********************************************************************/

void stretch_mainchain(struct object *object_ptr,
		       struct coordinate *n_ptr,struct coordinate *ca_ptr,
		       struct coordinate *c_ptr,struct coordinate *o_ptr,
		       struct coordinate *previous_ca_ptr,
		       struct coordinate *last_n_ptr,
		       struct coordinate *last_ca_ptr,
		       struct coordinate *last_c_ptr,
		       struct coordinate *last_o_ptr)
{
  int ncovalent, negside;
  int lose_C, nhbonds_O;
  int iatom, natoms;

  float coord_store[3][3], matrix[4][3];
  float x, x1, x2, y, y1, y2;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr;

  /* Initialise variables */
  lose_C = FALSE;

  /* Adjust the atom-sizes of the main-chain atoms */
  ca_ptr->atom_size = Size_Val->Simple_Residues;
  n_ptr->atom_size = Size_Val->Ligand_Bonds;
  c_ptr->atom_size = Size_Val->Ligand_Bonds;
  ca_ptr->plot_atom = FALSE;
  n_ptr->plot_atom = FALSE;
  c_ptr->plot_atom = FALSE;

  /* Update residue and object maximum sizes */
  if (ca_ptr->atom_size > ca_ptr->residue_ptr->object_ptr->max_atom_size)
    ca_ptr->residue_ptr->object_ptr->max_atom_size = ca_ptr->atom_size;
  if (ca_ptr->atom_size > ca_ptr->residue_ptr->max_atom_size)
    ca_ptr->residue_ptr->max_atom_size = ca_ptr->atom_size;

  /* Increment count of object's CA atoms */
  object_ptr->n_ca++;

  /* If have the previous CA or last N atoms, then can stretch the mainchain
     by pivoting at the last CA and getting the previous CA (or last N) in
     line with the current CA */
  if (last_ca_ptr != NULL && (previous_ca_ptr != NULL || last_n_ptr != NULL))
    {
      /* Store the coords of the last CA */
      coord_store[0][0] = last_ca_ptr->x;
      coord_store[0][1] = last_ca_ptr->y;
      coord_store[0][2] = last_ca_ptr->z;

      /* Store the coords of the previous CA, or last N, depending on
         which are available */
      if (previous_ca_ptr != NULL)
	{
	  coord_store[1][0] = previous_ca_ptr->x;
	  coord_store[1][1] = previous_ca_ptr->y;
	  coord_store[1][2] = previous_ca_ptr->z;
	}
      else
	{
	  coord_store[1][0] = last_n_ptr->x;
	  coord_store[1][1] = last_n_ptr->y;
	  coord_store[1][2] = last_n_ptr->z;
	}

      /* Calculate the transformation that will place the last CA at the
	 origin, with the previous CA (or last N) in the negative
	 x-direction */
      line_transformation(coord_store,matrix,-PI);

      /* Apply the transformation to all atoms in the current object */
      transform_object(object_ptr,matrix);

      /* Store the coords of the current CA atom */
      coord_store[0][0] = 0.0;
      coord_store[0][1] = 0.0;
      coord_store[0][2] = 0.0;
      coord_store[1][0] = ca_ptr->x;
      coord_store[1][1] = ca_ptr->y;
      coord_store[1][2] = ca_ptr->z;

      /* Calculate transformation to rotate about the last CA
	 to put the current CA atom along the +ve x-axis */
      line_transformation(coord_store,matrix,0.0);

      /* Apply the transformation to all the atoms in the object
         from the last CA onwards */
      transform_downlist_atoms(object_ptr,last_ca_ptr,matrix);
    }

  /* If have a previous CA and C atoms, then stretch the N and O residues
     equidistantly between the two */
  if (last_ca_ptr != NULL)
    {
      /* Get the coordinates of the 2 CA atoms */
      x1 = last_ca_ptr->x;
      y1 = last_ca_ptr->y;
      x2 = ca_ptr->x;
      y2 = ca_ptr->y;

      /* Get point one third of the way along */
      x = x1 + (x2 - x1) / 3.0;
      y = y1 + (y2 - y1) / 3.0;

      /* Place the previous carbon there */
      if (last_c_ptr != NULL)
	{
	  last_c_ptr->x = x;
	  last_c_ptr->y = y;

	  /* If there was an associated carbonyl oxygen, then place it 
	     on one side */
	  if (last_o_ptr != NULL)
	    {
	      last_o_ptr->x = x;
	      last_o_ptr->y = y + 1.23;
	    }
	}

      /* Get point two thirds of the way along */
      x = x1 + 2.0 * (x2 - x1) / 3.0;
      y = y1 + 2.0 * (y2 - y1) / 3.0;

      /* Place the current residue's nitrogen there */
      n_ptr->x = x;
      n_ptr->y = y;

      /* If residue is a proline, then need to tidy it up a little */
      if (!strncmp(n_ptr->residue_ptr->res_name,"PRO",3))
	{
	  /* Get the Pro's first atom */
	  atom_ptr = n_ptr->residue_ptr->first_atom_ptr;
	  iatom = 0;
	  natoms = n_ptr->residue_ptr->natoms;
	  negside = 0;

	  /* Loop over the residue's atoms to adjust the sidechain
	     coords */
	  while (iatom < natoms && atom_ptr != NULL)
	    {
	      /* If this is a sidechain atom, then fix its coords */
	      if (!strncmp(atom_ptr->atom_type," CB ",4))
		{
		  /* Check which side of the x-axis this atom lies on */
		  if (negside == 0)
		    if (atom_ptr->y < 1)
		      negside = -1;
		    else
		      negside = 1;
		  atom_ptr->x = ca_ptr->x + 0.400;
		  atom_ptr->y = negside * 1.500;
		  atom_ptr->z = 0.0;
		}

	      if (!strncmp(atom_ptr->atom_type," CG ",4))
		{
		  /* Check which side of the x-axis this atom lies on */
		  if (negside == 0)
		    if (atom_ptr->y < 1)
		      negside = -1;
		    else
		      negside = 1;
		  atom_ptr->x = (n_ptr->x + ca_ptr->x) / 2.0;
		  atom_ptr->y = negside * 2.300;
		  atom_ptr->z = 0.0;
		}

	      if (!strncmp(atom_ptr->atom_type," CD ",4))
		{
		  /* Check which side of the x-axis this atom lies on */
		  if (negside == 0)
		    if (atom_ptr->y < 1)
		      negside = -1;
		    else
		      negside = 1;
		  atom_ptr->x = n_ptr->x - 0.400;
		  atom_ptr->y = negside * 1.500;
		  atom_ptr->z = 0.0;
		}

	      /* Get pointer to the next atom */
	      atom_ptr = atom_ptr->next;
	      iatom++;
	    }
	}
    }

  /* Place the carbonyl oxygen at the carbonyl carbon position and 
     alter the bonds, if required */
  if (lose_C == TRUE && o_ptr != NULL)
    {
      /* Replace the carbonyl carbon by the oxygen and delete the
	 oxygen */
      strncpy(c_ptr->atom_type," O  ",4);
      strncpy(c_ptr->print_name," O  ",4);
      c_ptr->atom_type[4] = '\0';
      c_ptr->print_name[4] = '\0';
      o_ptr->deleted = TRUE;
    }

  /* Blank out the name of the CA atom */
  ca_ptr->print_name[0] = '\0';
  ncovalent = 0;

  /* Loop through all bonds to transfer any H-bonds from the oxygen
     to its new position on the carbonyl carbon (if required), or delete
     the oxygen if not involved in any H-bonds, and to transfer any
     non-bonded contacts from the carbon to the CA */
  nhbonds_O = 0;

  /* Get pointer to first bond in the linked list */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds */
  while (bond_ptr != NULL)
    {
      /* If this is an H-bond, check whether the O is involved */
      if (bond_ptr->bond_type == HBOND)
	{
	  /* Check for H-bonds to the current O atom */
	  if (bond_ptr->first_atom_ptr == o_ptr ||
	      bond_ptr->second_atom_ptr == o_ptr)
	    {
	      /* Increment count of H-bonds to the O */
	      nhbonds_O++;

	  /* If either atom of the bond is the oxygen which
	     has just been transferred, then transfer pointer, too */
	      if (lose_C == TRUE)
		{
		  if (bond_ptr->first_atom_ptr == o_ptr)
		    bond_ptr->first_atom_ptr = c_ptr;
		  if (bond_ptr->second_atom_ptr == o_ptr)
		    bond_ptr->second_atom_ptr = c_ptr;
		}
	    }
	}

      /* Transfer any contacts with the C to the CA */
      else if (bond_ptr->bond_type == CONTACT && lose_C == TRUE)
	{
	  /* If either atom of the bond is the carbonyl carbon
	     then transfer to the CA */
	  if (bond_ptr->first_atom_ptr == c_ptr)
	    bond_ptr->first_atom_ptr = ca_ptr;
	  else if (bond_ptr->second_atom_ptr == c_ptr)
	    bond_ptr->second_atom_ptr = ca_ptr;
	}

      /* Count the number of covalent bonds to the C */
      else if (bond_ptr->bond_type == COVALENT)
	{
	  /* If either atom of the bond is the carbon which
	     has just been overwritten, then increment count of
	     its covalent bonds */
	  if (bond_ptr->first_atom_ptr == c_ptr ||
	      bond_ptr->second_atom_ptr == c_ptr)
	    ncovalent++;
	}

      /* Get pointer for next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* If oxygen has no hydrogen bonds, then can be deleted */
/* v.3.2--> */
/*  if (lose_C == FALSE && nhbonds_O == 0) */
  if (lose_C == FALSE && nhbonds_O == 0 && o_ptr != NULL)
/* <--v.3.2 */
    o_ptr->deleted = TRUE;
}
/***********************************************************************

create_simple_ligand  -  Adjust all the ligand residues for the
                         schematic plot, marking any unwanted atoms
			 for deletion

***********************************************************************/

void create_simple_ligand(struct object *object_ptr)
{
  int iatom, iresid, natoms, nresid;
  int sidechain_interaction;

  struct bond *bond_ptr;
  struct coordinate *c_ptr, *ca_ptr, *n_ptr, *o_ptr;
  struct coordinate *last_ca_ptr, *last_n_ptr, *last_c_ptr, *last_o_ptr;
  struct coordinate *previous_ca_ptr;
  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise all atoms */
  initialise_atoms();

  /* Initialise bond pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds */
  while (bond_ptr != NULL)
    {
      /* Process only if this is an H-bond or hydrophobic contact */
      if (bond_ptr->bond_type == HBOND || bond_ptr->bond_type == CONTACT)
	{
	  /* Mark both atoms as being involved in the bond */
	  bond_ptr->first_atom_ptr->checked = TRUE;
	  bond_ptr->second_atom_ptr->checked = TRUE;
	}

      /* Get pointer for next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Loop over the object's residues to mark any unwanted sidechain
     atoms for deletion */

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Initialise flag */
      sidechain_interaction = FALSE;

      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, to locate all the mainchain
	 and sidechain atoms  */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* Initialise sidechain marker */
	  atom_ptr->side_chain = FALSE;

	  /* Check if atom is one of the mainchain atoms */
	  if (!strncmp(atom_ptr->atom_type," N  ",4))
	    n_ptr = atom_ptr;
	  else if (!strncmp(atom_ptr->atom_type," CA ",4))
	    ca_ptr = atom_ptr;
	  else if (!strncmp(atom_ptr->atom_type," C  ",4))
	    c_ptr = atom_ptr;
	  else if (!strncmp(atom_ptr->atom_type," O  ",4))
	    o_ptr = atom_ptr;

	  /* Otherwise, if it is a sidechain atom, then see whether
	     it is involved in an interaction */
	  else
	    {
	      atom_ptr->side_chain = TRUE;
	      if (atom_ptr->checked == TRUE)
		sidechain_interaction = TRUE;
	    }

	  /* Get pointer to the next atom in this residue */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

      /* If have all the mainchain atoms, then can stretch this
	 residue out into a straight line */
      if (n_ptr != NULL && ca_ptr != NULL && c_ptr != NULL)
	{
	  /* If none of the sidechain atoms are involved in interactions,
	     then mark the whole sidechain for deletion */
	  if (sidechain_interaction == FALSE)
	    mark_sidechain_atoms(residue_ptr);
	}

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Loop to stretch out the mainchain */

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Initialise pointer to previous residue's CA */
  previous_ca_ptr = NULL;
  last_ca_ptr = NULL;
  last_n_ptr = NULL;
  last_c_ptr = NULL;
  last_o_ptr = NULL;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Initialise flags */
      n_ptr = NULL;
      ca_ptr = NULL;
      c_ptr = NULL;
      o_ptr = NULL;

      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, to locate all the mainchain
	 and sidechain atoms  */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* Check if atom is one of the mainchain atoms */
	  if (!strncmp(atom_ptr->atom_type," N  ",4))
	    n_ptr = atom_ptr;
	  else if (!strncmp(atom_ptr->atom_type," CA ",4))
	    ca_ptr = atom_ptr;
	  else if (!strncmp(atom_ptr->atom_type," C  ",4))
	    c_ptr = atom_ptr;
	  else if (!strncmp(atom_ptr->atom_type," O  ",4))
	    o_ptr = atom_ptr;

	  /* Get pointer to the next atom in this residue */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

      /* If have all the mainchain atoms, then can stretch this
	 residue out into a straight line */
      if (n_ptr != NULL && ca_ptr != NULL && c_ptr != NULL)
	{
	  /* Stretch the mainchain out */
	  stretch_mainchain(object_ptr,n_ptr,ca_ptr,c_ptr,o_ptr,
			    previous_ca_ptr,last_n_ptr,last_ca_ptr,
			    last_c_ptr,last_o_ptr);

	  /* Mark the residue as being of the simplified type */
	  residue_ptr->residue_type = SIMPLE_LIGAND;

	  /* Save the pointers (note that O's position is saved as the
	     current carbonyl C) */
	  previous_ca_ptr = last_ca_ptr;
	  last_ca_ptr = ca_ptr;
	  last_n_ptr = n_ptr;
	  last_c_ptr = c_ptr;
	  last_o_ptr = o_ptr;
	}

      /* Otherwise, no stretching to be done */
      else
	{
	  last_n_ptr = NULL;
	  last_c_ptr = NULL;
	  last_o_ptr = NULL;
	}

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

Andrew Martin's matfit and qikfit routines for least-squares fitting of
one set of coordinates onto another (both sets centred around the origin)

***********************************************************************/
/*>static void qikfit(REAL umat[3][3], REAL rm[3][3], BOOL column)
   ---------------------------------------------------------------
   Input:   REAL  umat[3][3]     The U matrix
            BOOL  column         TRUE: Create a column-wise matrix
                                 (other way round from normal).
   Output:  REAL  rm[3][3]       The output rotation matrix
  
   Does the actual fitting for matfit().
   04.02.91 Original
   01.06.92 ANSIed & doc'd
   11.03.94 column changed to BOOL
*/
static void qikfit(REAL  umat[3][3],
                   REAL  rm[3][3],
                   BOOL  column)
{
   
   REAL  rot[3][3],
         turmat[3][3],
         c[3][3],
         coup[3],
         dir[3],
         step[3],
         v[3],
         rtsum,rtsump,
         rsum,
         stp,stcoup,
         ud,tr,ta,cs,sn,ac,
         delta,deltap,
         gfac,
         cle,clep;
   int   i,j,k,l,m,
         jmax,
         ncyc,
         nsteep,
         nrem;

   /* Rotate repeatedly to reduce couple about initial direction
      to zero.
      Clear the rotation matrix */
   for(l=0;l<3;l++)
   {
      for(m=0;m<3;m++)
         rot[l][m] = 0.0;
      rot[l][l] = 1.0;
   }

   /* Copy vmat[][] (sp) into umat[][] (dp) */
   jmax = 30;
   rtsum = umat[0][0] + umat[1][1] + umat[2][2];
   delta = 0.0;

   for(ncyc=0;ncyc<jmax;ncyc++)
   {
      /* Modified CG. For first and every NSTEEP cycles, set previous
         step as zero and do an SD step */
      nsteep = 3;
      nrem = ncyc-nsteep*(int)(ncyc/nsteep);

      if(!nrem)
      {
         for(i=0;i<3;i++) step[i]=0.0;
         clep = 1.0;
      }
      
      /* Couple */
      coup[0] = umat[1][2]-umat[2][1];
      coup[1] = umat[2][0]-umat[0][2];
      coup[2] = umat[0][1]-umat[1][0];
      cle = sqrt(coup[0]*coup[0] + coup[1]*coup[1] + coup[2]*coup[2]);

      /* Gradient vector is now -coup */
      gfac = (cle/clep)*(cle/clep);

      /* Value of rtsum from previous step */
      rtsump = rtsum;
      deltap = delta;
      clep   = cle;
      if(cle < SMALL) break;

      /* Step vector conjugate to  previous */
      stp = 0.0;
      for(i=0;i<3;i++)
      {
         step[i]=coup[i]+step[i]*gfac;
         stp   += (step[i] * step[i]);
      }
      stp = 1.0/sqrt(stp);
         
      /* Normalised step */
      for(i=0;i<3;i++) dir[i] = stp*step[i];

      /* Couple resolved along step direction */
      stcoup = coup[0]*dir[0] + coup[1]*dir[1] + coup[2]*dir[2];

      /* Component of UMAT along direction */
      ud = 0.0;
      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
            ud += umat[l][m]*dir[l]*dir[m];


      tr = umat[0][0]+umat[1][1]+umat[2][2]-ud;
      ta = sqrt(tr*tr + stcoup*stcoup);
      cs=tr/ta;
      sn=stcoup/ta;
         
      /* If cs<0 then posiiton is unstable, so don't stop */
      if((cs>0.0) && (ABS(sn)<SMALSN)) break;
            
      /* Turn matrix for correcting rotation:

         Symmetric part */
      ac = 1.0-cs;
      for(l=0;l<3;l++)
      {
         v[l] = ac*dir[l];
         for(m=0;m<3;m++)
            turmat[l][m] = v[l]*dir[m];
         turmat[l][l] += cs;
         v[l]=dir[l]*sn;
      }

      /* Asymmetric part */
      turmat[0][1] -= v[2];
      turmat[1][2] -= v[0];
      turmat[2][0] -= v[1];
      turmat[1][0] += v[2];
      turmat[2][1] += v[0];
      turmat[0][2] += v[1];

      /* Update total rotation matrix */
      for(l=0;l<3;l++)
      {
         for(m=0;m<3;m++)
         {
            c[l][m] = 0.0;
            for(k=0;k<3;k++)
               c[l][m] += turmat[l][k]*rot[k][m];
         }
      }

      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
            rot[l][m] = c[l][m];

      /* Update umat tensor */
      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
         {
            c[l][m] = 0.0;
            for(k=0;k<3;k++)
               c[l][m] += turmat[l][k]*umat[k][m];
         }

      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
            umat[l][m] = c[l][m];

      rtsum = umat[0][0] + umat[1][1] + umat[2][2];
      delta = rtsum - rtsump;

      /* If no improvement in this cycle then stop */
      if(ABS(delta)<SMALL) break;

      /* Next cycle */
   }

   rsum = rtsum;

   /* Copy rotation matrix for output */
   if(column)
   {
      for(i=0;i<3;i++)
         for(j=0;j<3;j++)
            rm[j][i] = rot[i][j];
   }
   else
   {
      for(i=0;i<3;i++)
         for(j=0;j<3;j++)
            rm[i][j] = rot[i][j];
   }
}
/**********************************************************************/
/**********************************************************************/
/*>BOOL matfit(COOR *x1, COOR *x2, REAL rm[3][3], int n,
               REAL *wt1, BOOL column)
   -----------------------------------------------------
   Input:   COOR  *x1         First array of coordinates
            COOR  *x2         Second array of coordinates
            int   n           Number of coordinates
            REAL  *wt1        Weight array or NULL
            int   column      TRUE: Output a column-wise matrix (as used
                                 by FRODO)
                              FALSE: Output a standard row-wise matrix.
   Output:  REAL  rm[3][3]    Returned rotation matrix
   Returns: rmsd              RMS distance between the two sets of
                              coordinates after fitting

   Fit coordinate array x1 to x2 both centred around the origin and of 
   length n. Optionally weighted with the wt1 array if iw flag is set.
   If flag is set, error messages will be issued internally. 
   If column is set the matrix will be returned column-wise rather
   than row-wise.

   04.02.91 Original
   01.06.92 ANSIed & doc'd
   17.06.93 various changes for release (including parameters)
   11.03.94 column changed to BOOL
   04.03.97 Calc of rmsd added (RAL)
*/
float matfit(COOR    *x1,        /* First coord array    */
	     COOR    *x2,        /* Second coord array   */
	     REAL    rm[3][3],   /* Rotation matrix      */
	     int     n,          /* Number of points     */
	     REAL    *wt1,       /* Weight array         */
	     BOOL    column)     /* Column-wise output   */
{
   int  i, ipoint, j;
   float dist, rmsd, x, y, z;

   REAL umat[3][3];

   if(wt1)
   {
      for(i=0;i<3;i++)
      {
         for(j=0;j<3;j++) umat[i][j] = 0.0;

         for(j=0;j<n;j++)
         {
            switch(i)
            {
               case 0:
                  umat[i][0] += wt1[j] * x1[j].x * x2[j].x;
                  umat[i][1] += wt1[j] * x1[j].x * x2[j].y;
                  umat[i][2] += wt1[j] * x1[j].x * x2[j].z;
                  break;
               case 1:
                  umat[i][0] += wt1[j] * x1[j].y * x2[j].x;
                  umat[i][1] += wt1[j] * x1[j].y * x2[j].y;
                  umat[i][2] += wt1[j] * x1[j].y * x2[j].z;
                  break;
               case 2:
                  umat[i][0] += wt1[j] * x1[j].z * x2[j].x;
                  umat[i][1] += wt1[j] * x1[j].z * x2[j].y;
                  umat[i][2] += wt1[j] * x1[j].z * x2[j].z;
                  break;
            }
         }
      }
   } 
   else
   {
      for(i=0;i<3;i++)
      {
         for(j=0;j<3;j++) umat[i][j] = 0.0;

         for(j=0;j<n;j++)
         {
            switch(i)
            {
               case 0:
                  umat[i][0] += x1[j].x * x2[j].x;
                  umat[i][1] += x1[j].x * x2[j].y;
                  umat[i][2] += x1[j].x * x2[j].z;
                  break;
               case 1:
                  umat[i][0] += x1[j].y * x2[j].x;
                  umat[i][1] += x1[j].y * x2[j].y;
                  umat[i][2] += x1[j].y * x2[j].z;
                  break;
               case 2:
                  umat[i][0] += x1[j].z * x2[j].x;
                  umat[i][1] += x1[j].z * x2[j].y;
                  umat[i][2] += x1[j].z * x2[j].z;
                  break;
            }
         }
      }
   }
   qikfit(umat,rm,column);

   /* Calculate rmsd between the two fitted sets of coordinates */
   rmsd = 0.0;

   /* Loop over all the pairs of coordinates */
   for (ipoint = 0; ipoint < n; ipoint++)
     {
       /* Calculate the transformed coordinates of the current point */
       if (column)
	 {
	   x = rm[0][0] * x1[ipoint].x + rm[1][0] * x1[ipoint].y
	     + rm[2][0] * x1[ipoint].z;
	   y = rm[0][1] * x1[ipoint].x + rm[1][1] * x1[ipoint].y
	     + rm[2][1] * x1[ipoint].z;
	   z = rm[0][2] * x1[ipoint].x + rm[1][2] * x1[ipoint].y
	     + rm[2][2] * x1[ipoint].z;
	 }
       else
	 {
	   x = rm[0][0] * x1[ipoint].x + rm[0][1] * x1[ipoint].y
	     + rm[0][2] * x1[ipoint].z;
	   y = rm[1][0] * x1[ipoint].x + rm[1][1] * x1[ipoint].y
	     + rm[1][2] * x1[ipoint].z;
	   z = rm[2][0] * x1[ipoint].x + rm[2][1] * x1[ipoint].y
	     + rm[2][2] * x1[ipoint].z;
	 }

       /* Calculate squared distance between these two points */
       dist = (x - x2[ipoint].x) * (x - x2[ipoint].x)
	 + (y - x2[ipoint].y) * (y - x2[ipoint].y)
	 + (z - x2[ipoint].z) * (z - x2[ipoint].z);
       rmsd = rmsd + dist;
     }

   /* Calculate the rmsd */
   if (n > 0)
     {
       rmsd = rmsd / n;
       if (rmsd > 0.0)
	 rmsd = sqrt(rmsd);
       else
	 rmsd = 0.0;
     }

   /* Return the calculated rmsd */
   return(rmsd);
}
/***********************************************************************

get_object_distance_match  -  Calculate the rms deviation between the
                              atom-atom distances within the current
			      object and the equivalent distances in the
			      original coordinates system

***********************************************************************/

float get_object_distance_match(struct object *object_ptr)
{
  int iatom, iresid, jatom, jresid, natoms1, nresid1, nresid2, natoms2;
  int nterms;
/* v.3.1--> */
  int atom_interval;
/* <--v.3.1 */

  float distance, dist2, orig_distance, rmsd, total_rmsd;
  float x1, x2, y1, y2;
  float orig_x1, orig_x2, orig_y1, orig_y2, orig_z1, orig_z2;

  struct coordinate *atom1_ptr, *atom2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise variables */
  total_rmsd = 0.0;

  /* Set the pointer to the object's first residue */
  residue1_ptr = object_ptr->first_residue_ptr;
  nresid1 = object_ptr->nresidues;
  iresid = 0;
  nterms = 0;

/* v.3.1--> */
  /* If have only one residue, use every atom, for two residues use every
     other atom, and so on */
  atom_interval = nresid1;
/* <--v.3.1 */

  /* Loop over all this object's residues */
  while (iresid < nresid1 && residue1_ptr != NULL)
    {
      /* Set second pointer to the object's first residue */
      residue2_ptr = object_ptr->first_residue_ptr;
      nresid2 = object_ptr->nresidues;
      jresid = 0;

      /* Loop over all the object's residues for the second atoms */
      while (jresid < nresid2 && residue2_ptr != NULL)
	{
	  /* Get pointer to residue's first atom */
	  atom1_ptr = residue1_ptr->first_atom_ptr;

	  /* Initialise atom count */
	  iatom = 0;
	  natoms1 = residue1_ptr->natoms;

	  /* Loop over all the residue's atoms */
	  while (iatom < natoms1 && atom1_ptr != NULL)
	    {
/* v.3.1--> */
	      /* If this is a wanted atom, then carry on */
	      if (iatom % atom_interval == 0)
		{
/* <--v.3.1 */
		  /* Get the atom's coordinates */
		  x1 = atom1_ptr->x;
		  y1 = atom1_ptr->y;
		  orig_x1 = atom1_ptr->original_x;
		  orig_y1 = atom1_ptr->original_y;
		  orig_z1 = atom1_ptr->original_z;

		  /* Get pointer to first compare atom */
		  atom2_ptr = residue2_ptr->first_atom_ptr;

		  /* Initialise atom count */
		  jatom = 0;
		  natoms2 = residue2_ptr->natoms;

		  /* Loop over all the atoms making up the compare atom */
		  while (jatom < natoms2 && atom2_ptr != NULL)
		    {
/* v.3.1--> */
		      /* If this is a wanted atom, then carry on */
		      if (jatom % atom_interval == 0)
			{
/* <--v.3.1 */
			  /* Skip if dealing with the same atom */
			  if (atom1_ptr != atom2_ptr)
			    {
			      /* Get atom's coordinates */
			      x2 = atom2_ptr->x;
			      y2 = atom2_ptr->y;
			      orig_x2 = atom2_ptr->original_x;
			      orig_y2 = atom2_ptr->original_y;
			      orig_z2 = atom2_ptr->original_z;

			      /* Calculate the distance between the
				 two atoms */
			      dist2 = (x1 - x2) * (x1 - x2)
				+ (y1 - y2) * (y1 - y2);
			      distance = sqrt((double) dist2);

			      /* Calculate the original distance between
				 these two atoms in the structure */
			      dist2 = (orig_x1 - orig_x2) * (orig_x1 - orig_x2)
				+ (orig_y1 - orig_y2) * (orig_y1 - orig_y2)
				  + (orig_z1 - orig_z2) * (orig_z1 - orig_z2);
			      orig_distance = sqrt((double) dist2);

			      /* Calculate the rmsd between the two
				 distances */
			      rmsd = distance - orig_distance;
			      rmsd = rmsd * rmsd;

			      /* Add to total rmsd value */
			      total_rmsd = total_rmsd + rmsd;
			      nterms++;
			    }
			}

		      /* Get pointer to other residue's next atom */
		      atom2_ptr = atom2_ptr->next;
		      jatom++;
		    }
		}

	      /* Get pointer to the next atom */
	      atom1_ptr = atom1_ptr->next;
	      iatom++;
	    }

          /* Get pointer to the next residue */
          residue2_ptr = residue2_ptr->next_residue_ptr;
          jresid++;
	}

      /* Get pointer to the next residue */
      residue1_ptr = residue1_ptr->next_residue_ptr;
      iresid++;
    }

  /* Return the calculated energy */
  if (nterms > 0)
    total_rmsd = total_rmsd / nterms;
  if (total_rmsd > 0.0)
    total_rmsd = sqrt(total_rmsd);
  else
    total_rmsd = 0.0;
  return(total_rmsd);
}
/***********************************************************************

fit_object  -  Fit the given object as best as possible to the 
               original version (to give approximate orientation and
	       relative disposition of the component atoms)

***********************************************************************/

float fit_object(struct object *object_ptr,int no_move)
{
  int i, iatom, iresid, istored, j;
  int natoms, nresid, nstored;
  int best_flipped;
/* v.4.0--> */
  int internal_fit, wanted;
/* <--v.4.0 */

  float transform[3][3];
  float mean_x, mean_y, mean_fit_x, mean_fit_y;
  float flipped_rmsd, rmsd;
  float distance_match, flipped_distance_match;
/* v.4.0--> */
  float residue_mean_x, residue_mean_y;
/* <--v.4.0 */

  BOOL column;

  REAL  flip_matrix[3][3], matrix[3][3];

  COOR coords1[MAXATOMS], coords1_flip[MAXATOMS], coords2[MAXATOMS];

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
/* v.4.0--> */
  if (no_move == TRUE)
    internal_fit = TRUE;
  else
    internal_fit = FALSE;
/* <--v.4.0 */
  best_flipped = FALSE;
  mean_x = 0.0;
  mean_y = 0.0;
  mean_fit_x = 0.0;
  mean_fit_y = 0.0;
  nstored = 0;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
/* v.4.0--> */
      /* Initialise residue mean positions */
      residue_mean_x = 0.0;
      residue_mean_y = 0.0;
/* <--v.4.0 */

      /* Get pointer to first residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms */
      while (iatom < natoms && atom_ptr != NULL)
	{
/* v.4.0--> */
	  /* Determine whether this atom's coords are to be included
	     in the fitting */
	  if (object_ptr->have_anchors == FALSE || internal_fit == TRUE)
	    wanted = TRUE;

	  /* If have any anchor points, then only use if this atom is
	     one of them */
	  else if (atom_ptr->have_anchor == TRUE)
	    wanted = TRUE;
	  else
	    wanted = FALSE;
/* <--v.4.0 */

	  /* Store the current coords and the fit coords */
/* v.4.0--> */
/*	  if (nstored < MAXATOMS) */
	  if (nstored < MAXATOMS && wanted == TRUE)
/* <--v.4.0 */
	    {
/* v.4.0--> */
	      /* If this is an anchor-atom, use its coordinates instead */
	      if (atom_ptr->have_anchor == TRUE)
		{
		  atom_ptr->fit_x = atom_ptr->anchor_pstn_x;
		  atom_ptr->fit_y = atom_ptr->anchor_pstn_y;
		}
/* <--v.4.0 */

	      /* Store the two sets of coordinates to be
		 fitted */
	      coords1[nstored].x = atom_ptr->x;
	      coords1[nstored].y = atom_ptr->y;
	      coords1[nstored].z = 0.0;
	      coords2[nstored].x = atom_ptr->fit_x;
	      coords2[nstored].y = atom_ptr->fit_y;
	      coords2[nstored].z = 0.0;

	      /* Increment count of stored coords */
	      nstored++;

	      /* Accumulate mean values */
	      mean_x = mean_x + atom_ptr->x;
	      mean_y = mean_y + atom_ptr->y;
	      mean_fit_x = mean_fit_x + atom_ptr->fit_x;
	      mean_fit_y = mean_fit_y + atom_ptr->fit_y;
	    }

/* v.4.0--> */
	  /* Accumulate residue means */
	  residue_mean_x = residue_mean_x + atom_ptr->x;
	  residue_mean_y = residue_mean_y + atom_ptr->y;
/* <--v.4.0 */

	  /* Get pointer to the next atom */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

/* v.4.0--> */
      /* If this is an anchor residue, calculate residue mean and
         use this as a fitting point */
      if (residue_ptr->have_anchor == TRUE && natoms > 0 &&
	  internal_fit == FALSE)
	{
	  /* Get means */
	  residue_mean_x = residue_mean_x / natoms;
	  residue_mean_y = residue_mean_y / natoms;

	  /* Store residue mean position and anchor position */
	  coords1[nstored].x = residue_mean_x;
	  coords1[nstored].y = residue_mean_y;
	  coords1[nstored].z = 0.0;
	  coords2[nstored].x = residue_ptr->anchor_pstn_x;
	  coords2[nstored].y = residue_ptr->anchor_pstn_y;
	  coords2[nstored].z = 0.0;

	  /* Increment count of stored coords */
	  nstored++;

	  /* Accumulate mean values */
	  mean_x = mean_x + residue_mean_x;
	  mean_y = mean_y + residue_mean_y;
	  mean_fit_x = mean_fit_x + residue_ptr->anchor_pstn_x;
	  mean_fit_y = mean_fit_y + residue_ptr->anchor_pstn_y;
	}
/* <--v.4.0 */

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Check that have some atoms belonging to this object */
  if (nstored > 0)
    {
      /* Calculate the mean values of the two sets of coords */
      mean_x = mean_x / nstored;
      mean_y = mean_y / nstored;
      mean_fit_x = mean_fit_x / nstored;
      mean_fit_y = mean_fit_y / nstored;

      /* Loop through all the stored coordinates to centre them at
	 the origin */
      for (istored = 0; istored < nstored; istored++)
	{
	  /* Adjust both sets of coords */
	  coords1[istored].x = coords1[istored].x - mean_x;
	  coords1[istored].y = coords1[istored].y - mean_y;
	  coords2[istored].x = coords2[istored].x - mean_fit_x;
	  coords2[istored].y = coords2[istored].y - mean_fit_y;

	  /* Store the flipped version of the object also,
	     to see if this gives a better fit */
	  coords1_flip[istored].x = - coords1[istored].x;
	  coords1_flip[istored].y = coords1[istored].y;
	  coords1_flip[istored].z = 0.0;
	}

      /* Calculate transformation that gives the best fit between
	 the flattened object and the centred original structure */
      if (nstored > 1)
	{
	  /* Set flag to return column-wise rotation matrix */
	  column = TRUE;

	  /* Try fitting the flipped version of the object first */
	  flipped_rmsd = matfit(coords1_flip,coords2,flip_matrix,
				nstored,NULL,column);

	  /* Calculate agreement in terms of all atom-atom distances
	     within the object */
	  flipped_distance_match = get_object_distance_match(object_ptr);
	  flipped_rmsd = flipped_rmsd + flipped_distance_match;

	  /* Now fit object straight */
	  rmsd = matfit(coords1,coords2,matrix,nstored,NULL,column);

	  /* Repeat calculation of atom-atom distance matches within
	     object */
	  distance_match = get_object_distance_match(object_ptr);
	  rmsd = rmsd + distance_match;

	  /* Determine which gives the best fit */
	  if (flipped_rmsd < rmsd)
	    {
	      best_flipped = TRUE;
	      rmsd = flipped_rmsd;
	    }

	  /* Get required rotation matrix */
	  for (i = 0; i < 3; i++)
	    for (j = 0; j < 3; j++)
	      if (best_flipped == TRUE)
		transform[i][j] = flip_matrix[i][j];
	      else
		transform[i][j] = matrix[i][j];
	}

      /* If only a single atom, then make matrix the unit matrix */
      else
	{
	  for (i = 0; i < 3; i++)
	    for (j = 0; j < 3; j++)
	      if (i == j)
/* v.3.2--> */
/*		matrix[i][j] = 1.0; */
		transform[i][j] = 1.0;
/* <--v.3.2 */
	      else
/* v.3.2--> */
/*		matrix[i][j] = 0.0; */
		transform[i][j] = 0.0;
/* <--v.3.2 */
	}

      /* Move object to its desired position, if required */
      if (no_move == FALSE)
	{
	  /* Translate object to centre of mass */
	  move_object(object_ptr,mean_x,mean_y,0.0);

	  /* If object is best flipped, then do so */
	  if (best_flipped == TRUE)
	    flip_object(object_ptr,-1,1);

	  /* Apply the rotation matrix to coords */
	  rotate_object(object_ptr,transform);

	  /* Move object back to where the fit coords came from */
	  move_object(object_ptr,- mean_fit_x,- mean_fit_y,0.0);
	}
    }

  /* Return the calculated rms deviation */
  return(rmsd);
}
/***********************************************************************

calc_internal_energy  -  Calculate the object's total energy from any
                         close contacts between its atoms

***********************************************************************/

float calc_internal_energy(struct object *object_ptr,int *clash)
{
  int conatoms, iatom, in_range, iresid, jatom, jresid, katom, natoms,
  nresid, other_natoms, wanted;

  float atom_size1, atom_size2, distance, dist2, interact_dist;
  float clash_energy, energy, fit_energy, rmsd, total_energy;
  float x1, x2, y1, y2;

  struct atom_link *atom_link_ptr;
  struct coordinate *atom_ptr, *other_atom_ptr;
  struct residue *residue_ptr, *other_residue_ptr;

  /* Initialise variables */
  *clash = FALSE;
  clash_energy = 0.0;
  fit_energy = 0.0;
  total_energy = 0.0;

  /* Calculate the maximum and minimum coords of atoms in each
     residue to speed up computation */
  update_residue_boundaries(object_ptr);

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get a second pointer to this residue for calculating
	 inter-residue contacts */
      other_residue_ptr = residue_ptr;
      jresid = iresid;

      /* Loop over all the possible other residues */
      while (jresid < nresid && other_residue_ptr != NULL)
	{
	  /* Get the maximum atom-size of each residue */
	  atom_size1 = residue_ptr->max_atom_size;
	  atom_size2 = other_residue_ptr->max_atom_size;
	  interact_dist = atom_size1 + atom_size2 + INTERNAL_INTERACT_DIST;

	  /* Check whether the residues are within range of one another */
	  in_range = TRUE;
	  if ((other_residue_ptr->minx >
	       (residue_ptr->maxx + interact_dist)) ||
	      (residue_ptr->minx >
	       (other_residue_ptr->maxx + interact_dist)) ||
	      (other_residue_ptr->miny >
	       (residue_ptr->maxy + interact_dist)) ||
	      (residue_ptr->miny >
	       (other_residue_ptr->maxy + interact_dist)))
	    in_range = FALSE;

	  /* If residues are within range, then check all distances
	     between them */
	  if (in_range == TRUE)
	    {
	      /* Get pointer to first residue's first atom */
	      atom_ptr = residue_ptr->first_atom_ptr;

	      /* Initialise atom count */
	      iatom = 0;
	      natoms = residue_ptr->natoms;

	      /* Loop over all the first residue's atoms */
	      while (iatom < natoms && atom_ptr != NULL)
		{
		  /* Get atom's coordinates */
		  x1 = atom_ptr->x;
		  y1 = atom_ptr->y;

		  /* Get the atom's size */
		  atom_size1 = atom_ptr->atom_size;

		  /* Get pointer to other residue's first atom */
		  other_atom_ptr = other_residue_ptr->first_atom_ptr;
		  
		  /* Initialise atom count */
		  jatom = 0;
		  other_natoms = other_residue_ptr->natoms;
		  
		  /* Loop over all the other residue's atoms */
		  while (jatom < other_natoms && other_atom_ptr != NULL)
		    {
		      /* Initialise flag */
		      wanted = TRUE;

		      /* If both are the same atom, then don't want them */
		      if (other_atom_ptr == atom_ptr)
			wanted = FALSE;
		      
		      /* Otherwise, check that they are not covalently
		         bonded */
		      else
			{
			  /* Get pointer to the first atom covalently
			     bonded to current atom of interest */
			  atom_link_ptr = atom_ptr->first_atom_link_ptr;
			  katom = 0;
			  conatoms = atom_ptr->natom_links;

			  /* Loop over all the connected atoms */
			  while (katom < conatoms && atom_link_ptr != NULL
				 && wanted == TRUE)
			    {
			      /* Check whether the connected atom is
				 the same as the other atom currently
				 under consideration */
			      if (atom_link_ptr->atom_ptr ==
				  other_atom_ptr)
				wanted = FALSE;
	
			      /* Get pointer to next connected atom */
			      atom_link_ptr
				= atom_link_ptr->next_atom_link_ptr;
			      katom++;
			    }
			}

		      /* If atoms are not bonded, calculate the distance
			 between them */
			if (wanted == TRUE)
			{
			  /* Get other atom's coordinates */
			  x2 = other_atom_ptr->x;
			  y2 = other_atom_ptr->y;

			  /* Get the other atom's size */
			  atom_size2 = other_atom_ptr->atom_size;

			  /* Calculate the distance between the
			     two atoms */
			  dist2 = (x1 - x2) * (x1 - x2)
			    + (y1 - y2) * (y1 - y2);
			  distance = sqrt((double) dist2);

			  /* Adjust distance to allow for the atom sizes */
			  if (distance > (atom_size1 + atom_size2))
			    distance = distance - (atom_size1 + atom_size2);
			  else
			    distance = 0.0;
			  dist2 = distance * distance;

			  /* Add energy of this interaction to the
			     overall energy sum */
			  if (dist2 < 0.001)
			    dist2 = 0.001;
			  if (dist2 < INTERNAL_DIST2)
			    energy = Atom_Atom_Clash / dist2;
			  else
			    energy = 0.0;
			  clash_energy = clash_energy + energy;
			}

		      /* Get pointer to other residue's next atom */
		      other_atom_ptr = other_atom_ptr->next;
		      jatom++;
		    }

		  /* Get pointer to the next atom */
		  atom_ptr = atom_ptr->next;
		  iatom++;
		}
	    }

          /* Get pointer to the next residue */
          other_residue_ptr = other_residue_ptr->next_residue_ptr;
          jresid++;
	}

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* See how well the current object can be fitted onto the original */
  rmsd = fit_object(object_ptr,TRUE);

  /* Calculate an equivalent energy according to the rms deviation of the
     best fit */
  fit_energy = rmsd;
  
  /* Get total object energy from its components */
  total_energy = clash_energy + fit_energy;
  if (fabs(clash_energy) < 0.001)
    *clash = FALSE;
  else
    *clash = TRUE;

  /* Scale the energy up by its weight, as read in from the run parameters */
  total_energy = Internal_Energy_Weight * total_energy;

  /* Return the calculated energy */
  return(total_energy);
}
/***********************************************************************

flip_about_bond  -  Flip part of the current object about the given bond,
                    transforming all the appropriate atoms' coordinates
                    from +y to -y

***********************************************************************/

void flip_about_bond(struct bond *flip_bond_ptr,int end)
{
  int loop, nbonds;

  struct bond *bond_ptr, *next_free_stack_ptr, *other_bond_ptr,
  *next_stack_ptr;
  struct bond_link *bond_link_ptr;
  struct coordinate *atom_ptr;

  /* Initialise variables */
  next_stack_ptr = NULL;
  next_free_stack_ptr = NULL;

  /* Initialise all the bond flags and stack pointers */
  initialise_bonds();

  /* Initialise all the atom flags */
  initialise_atoms();

  /* If haven't been given which end to use for the flip, then use
     whichever has fewer downstream of it */
  if (end != FIRST && end != SECOND)
    {
      /* If first atom has fewer bonds downstream of it, then want to
	 process the bonds coming off it */
      if (flip_bond_ptr->nbonds_from_first_atom
	  < flip_bond_ptr->nbonds_from_second_atom)
	end = FIRST;
  
      /* Otherwise, use the other end of the bond */
      else
	end = SECOND;
    }

  /* Set pointer to the first of the bonds sprouting off the flip
     bond */
  bond_link_ptr = flip_bond_ptr->first_bond_link_ptr;

  /* Set the flip-bond as already checked */
  flip_bond_ptr->checked = TRUE;

  /* Initialise count of bonds encountered */
  nbonds = 0;

  /* Initialise stack pointers for bond-search by putting the bonds
     coming off the selected end of the flip-bond onto the stack
     of bonds to be searched */
  while (bond_link_ptr != NULL)
    {
      /* Add this bond to the stack only if it comes off the relevant
	 atom of the test-bond */
      if (bond_link_ptr->bond_end == end)
	{
	  /* Get the bond this link is pointing to and add to stack */
	  bond_ptr = bond_link_ptr->bond_ptr;

	  /* Add to stack provided it is not an elastic bond */
	  if (bond_ptr->elastic == FALSE)
	    {
	      /* If this is the first to be added to the stack, set both
		 stack-pointers to point to it */
	      if (next_free_stack_ptr == NULL)
		{
		  next_free_stack_ptr = bond_ptr;
		  next_stack_ptr = bond_ptr;
		}

	      /* Otherwise, make last bond on stack point to this one
		 and set this one as nopw the last on the stack */
	      else
		{
		  next_free_stack_ptr->next_stack_ptr = bond_ptr;
		  next_free_stack_ptr = bond_ptr;
		}
	    }

	  /* Mark bond as added to the stack */
	  bond_ptr->checked = TRUE;
	}
      
      /* Go to the next bond-link in the linked list */
      bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
    }

  /* Loop until stack of pointers to be processed has been exhausted */
  while (next_stack_ptr != NULL)
    {
      /* Pick up the next bond from the stack, and move the stack
       pointer onto the next position */
      bond_ptr = next_stack_ptr;
      next_stack_ptr = bond_ptr->next_stack_ptr;

      /* Increment count of bonds encountered */
      nbonds++;

      /* Loop over the bond's two atoms */
      for (loop = 0; loop < 2; loop++)
	{
	  /* Get pointer to appropriate atom */
	  if (loop == 0)
	    atom_ptr = bond_ptr->first_atom_ptr;
	  else
	    atom_ptr = bond_ptr->second_atom_ptr;

	  /* If first atom has not been flipped, then flip it */
	  if (atom_ptr->checked == FALSE)
	    {
	      atom_ptr->y = - atom_ptr->y;
	      atom_ptr->checked = TRUE;
	    }
	}

     /* Loop through all the bonds connected to the current one
        to add to stack */

      /* Get the pointer to the first of the bonds coming off the 
	 current bond */
      bond_link_ptr = bond_ptr->first_bond_link_ptr;

      /* Process each of these bonds in turn */
      while (bond_link_ptr != NULL)
	{
	  /* Get the pointer to this bond */
	  other_bond_ptr = bond_link_ptr->bond_ptr;

	  /* If this bond hasn't already been processed, then move
	     its atoms (if they need moving) and add the bond to the
	     stack */
	  if (other_bond_ptr->checked == FALSE &&
	      other_bond_ptr->elastic == FALSE)
	    {
	      next_free_stack_ptr->next_stack_ptr = other_bond_ptr;
	      next_free_stack_ptr = other_bond_ptr;

	      /* If the stack has run out, restart it at the bond
		 just entered */
	      if (next_stack_ptr == NULL)
		next_stack_ptr = other_bond_ptr;

	      /* Mark bond as added to the stack */
	      other_bond_ptr->checked = TRUE;
	    }

	  /* Go to the next bond-link in the linked list */
	  bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
	}
    }
    
  /* Update the counter of bonds downstream on fliped bond */
  if (end == FIRST)
    flip_bond_ptr->nbonds_from_first_atom = nbonds;
  else
    flip_bond_ptr->nbonds_from_second_atom = nbonds;
}
/***********************************************************************

flip_back  -  Flip all the flagged atoms back the way they were

***********************************************************************/

void flip_back(void)
{
  struct coordinate *atom_ptr;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms to initialise the flag indicating whether
     coords have already been stored */
  while (atom_ptr != NULL)
    {
      /* If the atom is flagged as flipped, then flip it back */
      if (atom_ptr->checked == TRUE)
	atom_ptr->y = - atom_ptr->y;

      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next;
    }	       
}
/***********************************************************************

get_random_number  -  Routine for generating a random number from the
                      supplied seed. Returns a uniform random deviate
		      between 0.0 and 1.0. On the first call, the
		      seed must be set to any negative value to
		      initialise or reinitialise the sequence.
		      (Algorithm is that described for routine ran2
		      in Numerical Recipes).

***********************************************************************/

float get_random_number(int *iseed,int random_start)
{
  static int ir[97], iy, j;

  static float random_number, seed;

  static int first_pass = TRUE;
  static int m = 714025;
  static int ia = 1366;
  static int ic = 150889;
  static float rm = 1.4005112E-6;

  FILE *file_ptr;

  /* Initialise the random number */
  random_number = 0.0;

  /* If this is the first pass, or sequence to be re-initialised, then
     pick up seed from file */
  if (*iseed < 0 || first_pass == TRUE)
    {
      /* Initialise the seed */
      seed = 0.0;

      /* If a completely random start is required, then read in the
	 number written out to the file ranseed.dat last time this
	 routine was run */
      if (random_start == TRUE)
	{
	  /* Open the random seed file and read in seed */
	  if ((file_ptr = fopen("ranseed.dat","r")) !=  NULL)
	    {
	      /* If file exists, then read in the seed */
	      fscanf(file_ptr,"%f",&seed);

	      /* Close the file */
	      fclose(file_ptr);
	    }

	  /* Store integer version of the seed */
	  *iseed = seed;
	}

      /* Otherwise, use a standard starting point */
      else
	*iseed = 54813;

      /* Check that the seed is non-zero */
      if (*iseed == 0) 
	*iseed = 1;

      /* Generate random number from seed */
      *iseed = (ic - *iseed) % m;
      for (j = 0; j < 97; j++)
	{
	  *iseed = (ia * (*iseed) + ic) % m;
	  ir[j] = *iseed;
	}
      *iseed = (ia * (*iseed) + ic) % m;
      iy = *iseed;
    }

  /* Check for error */
  j = (97 * iy) / m;
  if (j >= 97 || j < 0)
    {
      printf("*** Random number error\n");
      exit(1);
    }

  /* Get the random number */
  iy = ir[j];
  random_number = iy * rm;
  *iseed = (ia * (*iseed) + ic) % m;
  ir[j] = *iseed;

  /* If this is the first call, use the generated random number to
     create a new seed and write to disk */
  if (first_pass == TRUE)
    {
      first_pass = FALSE;
      seed = iy * rm * 100000;

      /* Open the random seed file */
      if ((file_ptr = fopen("ranseed.dat","w")) !=  NULL)
	{
	  /* If no file error, then write out the seed */
	  fprintf(file_ptr,"%f\n",seed);

	  /* Close the file */
	  fclose(file_ptr);
	}
    }

  /* Return the random number generated */
  return(random_number);
}
/***********************************************************************

save_restore_coords  -  Save the given object's atom positions

***********************************************************************/

void save_restore_coords(struct object *object_ptr,int action)
{
  int iatom, iresid, natoms, nresid;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, performing the save */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* Save the x- and y-coords */
	  if (action == SAVE)
	    {
	      atom_ptr->save_x = atom_ptr->x;
	      atom_ptr->save_y = atom_ptr->y;
	    }

	  /* Restore the coords */
	  else
	    {
	      atom_ptr->x = atom_ptr->save_x;
	      atom_ptr->y = atom_ptr->save_y;
	    }

	  /* Get pointer to the next atom */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

save_restore_all_objects  -  Loop through all objects to save/restore
                             their coordinates

***********************************************************************/

void save_restore_all_objects(int restore)
{
  int iobject;

  struct object *object_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Save/restore this object's coords */
      save_restore_coords(object_ptr,restore);

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }
}
/***********************************************************************

get_matrix  -  Calculate the rotation matrix for a rotation about the
               z-axi through the given angle

***********************************************************************/

void get_matrix(float angle,float matrix[4][3])
{
  int i, j;

  float cos_theta, sin_theta;

  /* Initialise the matrix */
  for (i = 0; i < 4; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  if (i == j && i < 3)
	    matrix[i][j] = 1.0;
	  else
	    matrix[i][j] = 0.0;
	}
    }

  /* Calculate sines and cosines of the angle */
  cos_theta = cos(angle);
  sin_theta = sin(angle);

  /* Form the appropriate rotation matrix */
  matrix[0][0] = cos_theta;
  matrix[0][1] = sin_theta;
  matrix[1][0] = - sin_theta;
  matrix[1][1] = cos_theta;
}
/* v.3.2--> */
/***********************************************************************

calc_rotation_matrix  -  Calculate the rotation matrix for a given
                         rotation about the z-axis

***********************************************************************/

void calc_rotation_matrix(float angle,float matrix[3][3])
{
  int i, j;

  float cos_theta, sin_theta;

  /* Initialise the matrix */
  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  if (i == j)
	    matrix[i][j] = 1.0;
	  else
	    matrix[i][j] = 0.0;
	}
    }

  /* Calculate sines and cosines of the angle */
  cos_theta = cos(- angle / RADDEG);
  sin_theta = sin(- angle / RADDEG);

  /* Form the appropriate rotation matrix */
  matrix[0][0] = cos_theta;
  matrix[0][1] = sin_theta;
  matrix[1][0] = - sin_theta;
  matrix[1][1] = cos_theta;
}
/* <--v.3.2 */
/***********************************************************************

untangle_object  -  Untangle atom and bond overlaps by flipping structure
                    about its rotatable bonds

***********************************************************************/

void untangle_object(struct object *object_ptr,int *rseed)
{
  int end, iseed;
  int cycle, nflips;
  int clash, last_clash, swivel;

  float coord_store[3][3], matrix[4][3];
  float energy, last_energy;
  float rotate_angle, x;

  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr, *swap_ptr;
  struct object_bond *object_bond_ptr;

  /* Get the current random number generator seed */
  iseed = *rseed;

  /* Initialise */
  nflips = 1;

  /* Perform the energy minimization process */
  for (cycle = 0; cycle < MAXCYCLE && nflips != 0; cycle++)
    {
      /* Initialise count of flips made to structure */
      nflips = 0;
      if (cycle > 0 && cycle % 10 == 0)
	{
	  if (object_ptr->object_type == LIGAND)
	    printf("   Untangling ligand,           cycle %4d\n",cycle);
	  else
	    printf("   Untangling residue %s %s %c, cycle %4d\n",
		   object_ptr->first_residue_ptr->res_name,
		   object_ptr->first_residue_ptr->res_num,
		   object_ptr->first_residue_ptr->chain,cycle);
	}

      /* Calculate current energy of the object from atom clashes */
      energy = calc_internal_energy(object_ptr,&clash);
      last_clash = clash;
      last_energy = energy;

      /* Get pointer to the first of this object's bonds */
      object_bond_ptr = object_ptr->first_object_bond_ptr;
      
      /* Loop through all object's bonds */
      while (object_bond_ptr != NULL)
	{
	  /* Get the bond pointer */
	  bond_ptr = object_bond_ptr->bond_ptr;

	  /* If this is a rotatable bond see if object will look
	     better if flipped about this bond */
	  if (bond_ptr->rotatable_bond == TRUE)
	    {
	      /* Initialise swivel flag */
	      swivel = FALSE;

	      /* Get pointers to the bond's two atoms */
	      atom1_ptr = bond_ptr->first_atom_ptr;
	      atom2_ptr = bond_ptr->second_atom_ptr;

	      /* If have been flipping a long time, and this is a
		 side-chain bond, may want to try sivelling the
		 sidechain about this bond */
	      if (cycle + 1 >= START_SWIVELS &&
		  (atom1_ptr->side_chain == TRUE ||
		  atom2_ptr->side_chain == TRUE))
		{
		  /* Determine whether to flip or swivel */
		  x = get_random_number(&iseed,Random_Start);
		  if (x <= 0.5)
		    swivel = TRUE;
		}

	      /* If swivelling, then store all the object's coordinates
	         ane determine about which end of the bond to perform the
		 swivel */
	      if (swivel == TRUE)
		{
		  /* Store all the object's atoms */
		  save_restore_coords(object_ptr,SAVE);

		  /* Determine which end to use for the swivel */
		  x = get_random_number(&iseed,Random_Start);
		  if (x <= 0.5)
		    end = FIRST;
		  else
		    {
		      /* To swivel about the second end, swap the
			 bond's two atoms around */
		      end = SECOND;
		      swap_ptr = atom1_ptr;
		      atom1_ptr = atom2_ptr;
		      atom2_ptr = swap_ptr;
		    }
		}

	      /* Store the two atoms' coordinates */
	      coord_store[0][0] = atom1_ptr->x;
	      coord_store[0][1] = atom1_ptr->y;
	      coord_store[0][2] = atom1_ptr->z;
	      coord_store[1][0] = atom2_ptr->x;
	      coord_store[1][1] = atom2_ptr->y;
	      coord_store[1][2] = atom2_ptr->z;

	      /* Calculate the transformation that will put this bond
		 along the x-axis, with atom 1 at the origin and atom 2
		 pointing in the -ve x-direction */
	      line_transformation(coord_store,matrix,0.0);

	      /* Apply the transformation to all atoms in the current
		 object */
	      transform_object(object_ptr,matrix);

	      /* If swivelling about the bonds, then do so */
	      if (swivel == TRUE)
		{
		  /* Get swivel angle */
		  x = get_random_number(&iseed,Random_Start);
		  if (x <= 0.5)
		    rotate_angle = SWIVEL_ANGLE / RADDEG;
		  else
		    rotate_angle = - SWIVEL_ANGLE / RADDEG;

		  /* Form the transformation matrix to perform the
		     swivel */
		  get_matrix(rotate_angle,matrix);

		  /* Apply the transformation to all atoms springing off
		     the current bond */
		  transform_downstream_atoms(bond_ptr,end,matrix);
		}

	      /* Otherwise, flip one side of the object around this
		 rotatable bond (currently lying along the x-axis) by
		 changing the signs of all the y-coords */
	      else
		flip_about_bond(bond_ptr,EITHER);

	      /* Calculate energy of the flipped configuration */
	      energy = calc_internal_energy(object_ptr,&clash);

	      /* If the energy is greater than that of the previous
		 configuration, then restore to previous conformation */
	      if (energy > (last_energy - 0.000001))
		{
		  /* If performed a swivel, then restore all the
		     coordinates as they were before */
		  if (swivel == TRUE)
		    save_restore_coords(object_ptr,RESTORE);

		  /* Otherwise, just flip the y-values back */
		  else
		    flip_back();

		  /* Retrieve the previous energy */
		  energy = last_energy;
		  clash = last_clash;
		}

	      /* Otherwise, increment count of flips made */
	      else
		nflips++;

	      /* Save the current energy */
	      last_energy = energy;
	      object_ptr->internal_energy = energy;
	    }

	  /* Get pointer to this object's next bond entry */
	  object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
	}

      /* Keep going until have started svivelling */
      if (clash == TRUE && nflips == 0)
	nflips = 1;
    }

  /* Return the current random number seed */
  *rseed = iseed;
}
/***********************************************************************

flatten_all_objects  -  Flatten all objects into 2D by squashing rings
                        and unravelling along all rotatable bonds

***********************************************************************/

void flatten_all_objects(int *iseed)
{
  int iobject, n_flattened;
  int rseed;

  struct object *object_ptr;

  /* Initialise pointer to the first stored object */
  printf("\nFlattening and untangling objects ...\n");
  object_ptr = first_object_ptr;
  iobject = 0;
  rseed = *iseed;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* If plotting schematic diagram, prepare the ligand or H-group */
      if (Include->Simple_Nonligand_Residues == TRUE)
	{
	  /* If this is an H-group, remove all unwanted atoms */
	  if (object_ptr->object_type == HGROUP)
	    create_simple_hgroup(object_ptr);

	  /* Delete any unwanted atoms */
	  delete_atoms_and_bonds();
	}

      /* Otherwise, if not plotting non-ligand mainchain atoms, then
	 get rid of mainchains not involved in interactions */
      else if (Include->Mainchain_Atoms == FALSE &&
	       object_ptr->object_type == HGROUP)
	{
	  /* Remove unwanted mainchain atoms */
	  remove_mainchains(object_ptr);

	  /* Delete any unwanted atoms */
	  delete_atoms_and_bonds();
	}

      /* Update each residue's maximum and minimum coords */
      update_residue_boundaries(object_ptr);

      /* If object is just a hydrophobic group, (or a H-bonded group
	 in the simplified representation) then place at origin */
      if (object_ptr->object_type == HYDROPHOBIC ||
	  object_ptr->object_type == SIMPLE_HGROUP ||
	  object_ptr->object_type == WATER)
	place_group_at_origin(object_ptr);

      /* Otherwise, perform the full flattening of rings, unrolling
	 and untangling procedure */
      else
	{
	  /* Flatten any rings in the structure so that their non-planarity
	     doesn't mess up the process of unrolling the structure */
	  n_flattened = flatten_all_rings(object_ptr);

	  /* If any rings have been flattened, then may need to correct
	     the lengths of bonds attached to them as these may have
	     become deformed */
	  if (n_flattened > 0)
	    correct_bond_lengths(object_ptr);

	  /* Fix any bonds springing off the rings */
	  flatten_ring_offshoots(object_ptr);

	  /* Unroll the structure by flattening out the bonds either
	     side of each rotatable bond */
	  unroll_structure(object_ptr);

/* v.3.2--> */
	  /* Prior to forcing flat, align object as best as possible
	     with x-y plane */
	  align_object(object_ptr);
/* <--v.3.2 */

	  /* Force the object into x-y plane, whether or not it's quite
	     there yet, so that flipping processes don't mess it up
	     completely */
	  force_flat(object_ptr);

	  /* If plotting schematic diagram, prepare this object if
	     it is a ligand */
	  if (Include->Simple_Ligand_Residues == TRUE &&
	      object_ptr->object_type == LIGAND)
	    {
	      /* Stretch the ligand residues out along the mainchain
		 and delete any unwanted sidechain atoms */
	      create_simple_ligand(object_ptr);

	      /* Delete any unwanted atoms */
	      delete_atoms_and_bonds();
	    }

	  /* Flip about rotatable bonds to minimize atom and
	     bond overlaps as far as is possible for this object */
	  untangle_object(object_ptr,&rseed);
	}

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

  /* Retrieve the current random-number seed */
  *iseed = rseed;
}
/***********************************************************************
 
score_objects - Score all objects according to their importance

***********************************************************************/

void score_objects(void)
{
  int iresid, natoms, nresid;
  int weight;

  struct bond *bond_ptr;
  struct object *object_ptr, *object1_ptr, *object2_ptr;
  struct residue *residue_ptr, *residue1_ptr, *residue2_ptr;

  /* Get pointer to the first of the objects */
  object_ptr = first_object_ptr;

  /* Loop through all the objects */
  while (object_ptr != NULL)
    {
      /* Initialise weight */
      weight = 0;

      /* Determine weight according to object type */

      /* Highest weight goes to ligand itself */
      if (object_ptr->object_type == LIGAND)
	weight = 100;

      /* Single-atom groups weight the lowest */
      else if (object_ptr->object_type == HYDROPHOBIC ||
	       object_ptr->object_type == SIMPLE_HGROUP ||
	       object_ptr->object_type == WATER)
	weight = 1;

      /* Otherwise, for ordinary H-groups, weight according to their
	 number of atoms */
      else
	{
	  /* Set the pointer to the first of this object's residues */
	  residue_ptr = object_ptr->first_residue_ptr;
	  nresid = object_ptr->nresidues;
	  iresid = 0;

	  /* Loop over all this object's residues */
	  while (iresid < nresid && residue_ptr != NULL)
	    {
	      /* Get the number of atoms in this residue and add to weight */
	      natoms = residue_ptr->natoms;
	      weight = weight + natoms;

	      /* Get pointer to the next residue */
	      residue_ptr = residue_ptr->next_residue_ptr;
	      iresid++;
	    }
	}

      /* Store object's weight */
      object_ptr->weight = weight;

      /* Get pointer to the next object */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* Get pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through the stored bond-list to see if already have bond */
  while (bond_ptr != NULL)
    {
      /* Initialise the weight */
      weight = 0;

      /* Get the two objects at either end of the bond */
      residue1_ptr = bond_ptr->first_atom_ptr->residue_ptr;
      residue2_ptr = bond_ptr->second_atom_ptr->residue_ptr;
      object1_ptr = residue1_ptr->object_ptr;
      object2_ptr = residue2_ptr->object_ptr;

      /* If this is a covalent bond, then only want it if it is between
	 ligand and non-ligand residues */
      if (bond_ptr->bond_type == COVALENT)
	{
	  /* If not the same object, then check if one is the ligand */
	  if (object1_ptr != object2_ptr)
	    {
	      if (residue1_ptr->inligand == TRUE ||
		  residue2_ptr->inligand == TRUE)
		weight = 20;
	    }
	}

      /* Check for H-bond */
      else if (bond_ptr->bond_type == HBOND)
	{
	  /* Assign H-bond weight */
	  weight = 15;
	}

      /* Check for non-bonded contact */
      else if (bond_ptr->bond_type == CONTACT)
	{
	  /* Assign weight for non-bonded contact */
	  weight = 5;
	}

      /* Add the weight for the current bond to both objects involved */
      object1_ptr->weight = object1_ptr->weight + weight;
      if (object1_ptr != object2_ptr)
	object2_ptr->weight = object2_ptr->weight + weight;

      /* Go to the next bond in the linked list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/***********************************************************************
 
sort_objects - Sort all objects according to their importance using their
               computed weights

***********************************************************************/
void sort_objects(void)
{
  int flipped;

  float weight, other_weight;

  struct object *object_ptr, *last_object_ptr, *other_object_ptr;

  /* Initialise variables */
  flipped = TRUE;

  /* Loop until all the values have been sorted */
  while (flipped == TRUE)
    {
      /* Initialise flag */
      flipped = FALSE;

      /* Get pointer to the first of the objects */
      object_ptr = first_object_ptr;
      last_object_ptr = NULL;

      /* Loop through all the objects */
      while (object_ptr != NULL)
	{
	  /* Get the object's weight */
	  weight = object_ptr->weight;

	  /* Get pointer to the next object */
	  other_object_ptr = object_ptr->next_object_ptr;

	  /* If have the next object, then proceed */
	  if (other_object_ptr != NULL)
	    {
	      /* Get the other object's weight */
	      other_weight = other_object_ptr->weight;

	      /* If other object's weight is greater, then swap the
		 order of the objects */
	      if (other_weight > weight)
		{
		  /* Set flag to indicate that swap occurred */
		  flipped = TRUE;

		  /* If this is the first object, get the first pointer
		     to point to the second object */
		  if (last_object_ptr == NULL)
		      first_object_ptr = other_object_ptr;

		  /* Otherwise, set the last object's pointer to the
		     other object */
		  else
		    last_object_ptr->next_object_ptr = other_object_ptr;

		  /* Swap the pointers */
		  object_ptr->next_object_ptr
		    = other_object_ptr->next_object_ptr;
		  other_object_ptr->next_object_ptr = object_ptr;

		  /* Set the other object as the last object */
		  last_object_ptr = other_object_ptr;
		}

	      /* Otherwise, set the current object as the last object */
	      else
		last_object_ptr = object_ptr;
	    }
	  else
	    last_object_ptr = object_ptr;

	  /* Get pointer to the next object in the stack */
	  object_ptr = last_object_ptr->next_object_ptr;
	}
    }
}
/* v.4.0--> */
/***********************************************************************

calc_residue_means  -  Calculate mean positions of each residue in the
                       plot

***********************************************************************/

void calc_residue_means(void)
{
  int iatom, iresid, natoms, ncoords, nresid;

  float mean_x, mean_y, mean_z;

  struct coordinate *atom_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Set the pointer to the first of this object's residues */
      residue_ptr = object_ptr->first_residue_ptr;
      nresid = object_ptr->nresidues;
      iresid = 0;

      /* Loop over all this object's residues */
      while (iresid < nresid && residue_ptr != NULL)
	{
	  /* Initialise the residue's average position */
	  mean_x = 0.0;
	  mean_y = 0.0;
	  mean_z = 0.0;
	  ncoords = 0;

	  /* Get pointer to first residue's first atom */
	  atom_ptr = residue_ptr->first_atom_ptr;
	  iatom = 0;
	  natoms = residue_ptr->natoms;

	  /* Loop over all this residue's atoms */
	  while (iatom < natoms && atom_ptr != NULL)
	    {
	      /* Update the mean coordinates */
	      mean_x = mean_x + atom_ptr->original_x;
	      mean_y = mean_y + atom_ptr->original_y;
	      mean_z = mean_z + atom_ptr->original_z;
	      ncoords++;

	      /* Get pointer to the next atom */
	      atom_ptr = atom_ptr->next;
	      iatom++;
	    }

	  /* Calculate this residue's mean coordinates */
	  if (ncoords > 0)
	    {
	      /* Get mean coords */
	      mean_x = mean_x / ncoords;
	      mean_y = mean_y / ncoords;
	      mean_z = mean_z / ncoords;

	      /* Store the object's original and flattened mean coords */
	      residue_ptr->original_mean_x = mean_x;
	      residue_ptr->original_mean_y = mean_y;
	      residue_ptr->original_mean_z = mean_z;
	    }

	  /* Get pointer to the next residue */
	  residue_ptr = residue_ptr->next_residue_ptr;
	  iresid++;
	}

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }
}
/***********************************************************************

read_rcm_file  -  Read in the .rcm file, if present and match up the
                  relevant objects

***********************************************************************/

void read_rcm_file(char pdb_name[FILENAME_LEN])
{
  int iatom, iresid, natoms, nresid;

  char atom_type[5];
  char input_line[LINELEN + 1], rcm_name[FILENAME_LEN];
  char chain, res_name[4], res_num[6];

  int done, have_file;
  int line, natom_anchors, nres_anchors;

  float x, y;

  struct coordinate *atom_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  have_file = TRUE;
  line = 0;
  natom_anchors = 0;
  nres_anchors = 0;

  /* Get the name of the .rcm file */
  get_filename(pdb_name,rcm_name,".rcm");

  /* Open .rcm file to see if it is present */
  if ((file_ptr = fopen(rcm_name,"r")) ==  NULL)
    have_file = FALSE;

  /* If file is present, read through all the records */
  if (have_file == TRUE)
    {
      /* Read in the data from the PDB file */
      printf("\nReading in the residue restraints file [%s] ...\n",
	     rcm_name);

      /* Loop through the PDB file, picking up all the atoms that
	 are required */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL)
	{
	  /* If not a header record, get the details */
	  if (!strncmp(input_line,"ATOM  ",6) ||
	      !strncmp(input_line,"HETATM",6) ||
	      !strncmp(input_line,"RESDUE",6))
	    {
	      /* Get the details from this line */
	      strncpy(res_name,input_line+17,3);
	      res_name[3] = '\0';
	      strncpy(res_num,input_line+22,5);
	      res_num[5] = '\0';
	      chain = input_line[21];
	      if (chain == '-')
		chain = ' ';

	      /* Get the x- and y-coordinates */
	      sscanf(input_line+30," %f %f ",&x,&y);

	      /* Initialise pointer to the first stored object */
	      object_ptr = first_object_ptr;
	      done = FALSE;

	      /* Loop through all objects to find the residue corresponding
	         to the anchor-position just read in */
	      while (object_ptr != NULL && done == FALSE)
		{
		  /* Get pointer to this object's first residue */
		  residue_ptr = object_ptr->first_residue_ptr;
		  nresid = object_ptr->nresidues;
		  iresid = 0;

		  /* Loop over all this object's residues */
		  while (iresid < nresid && residue_ptr != NULL &&
			 done == FALSE)
		    {
		      /* See if this is the relevant residue */
		      if (!strncmp(residue_ptr->res_name,res_name,3) &&
			  !strncmp(residue_ptr->res_num,res_num,5) &&
			  residue_ptr->chain == chain)
			{
			  /* If the anchor position refers to the CofM
			     of the reaisue, then store it */
			  if (!strncmp(input_line+12,"CofM",4))
			    {
			      residue_ptr->have_anchor = TRUE;
			      residue_ptr->anchor_pstn_x = x;
			      residue_ptr->anchor_pstn_y = y;
			      object_ptr->have_anchors = TRUE;
			      done = TRUE;
			      nres_anchors++;
			    }

			  /* Otherwise, take this to be the anchor for one
			     of this residue's atom-positions */
			  else
			    {
			      /* Get the atom name */
			      strncpy(atom_type,input_line+12,4);
			      atom_type[4] = '\0';

			      /* Get pointer to first residue's first atom */
			      atom_ptr = residue_ptr->first_atom_ptr;
			      iatom = 0;
			      natoms = residue_ptr->natoms;

			      /* Loop over all this residue's atoms */
			      while (iatom < natoms && atom_ptr != NULL &&
				     done == FALSE)
				{
				  /* If this is the right atom, then
				     store its anchor position */
				  if (!strncmp(atom_ptr->atom_type,
					       atom_type,4))
				    {
				      /* Increment count of this residue's
					 atom-anchor positions */
				      if (atom_ptr->have_anchor == FALSE)
					residue_ptr->nanchored_atoms++;

				      /* Store anchor position */
				      atom_ptr->have_anchor = TRUE;
				      atom_ptr->anchor_pstn_x = x;
				      atom_ptr->anchor_pstn_y = y;
				      object_ptr->have_anchors = TRUE;
				      done = TRUE;
				      natom_anchors++;
				    }

				  /* Get pointer to the next atom */
				  atom_ptr = atom_ptr->next;
				  iatom++;
				}

			      /* If appropriate atom not found, then
				 display error message */
			      if (done == FALSE)
				{
				  /* Print message */
				  printf("*** WARNING. Atom [%s] for ",
					 atom_type);
				  printf("residue [%s %s %c] in .rcm file",
					 res_name,res_num,chain);
				  printf(" not found.\n");
				  printf("***          Restraint ignored.\n");
				  Nwarnings++;
				  done = TRUE;
				}
			    }
			}

		      /* Get pointer to the next residue */
		      residue_ptr = residue_ptr->next_residue_ptr;
		      iresid++;
		    }

		  /* Get pointer to next object */
		  object_ptr = object_ptr->next_object_ptr;
		}

	      /* If residue not found, then print warning message */
	      if (done == FALSE)
		{
		  /* Print message */
		  printf("*** WARNING. Residue [%s %s %c] in .rcm file",
			 res_name,res_num,chain);
		  printf(" not found. Restraint ignored.\n");
		  Nwarnings++;
		}
	    }

	  /* Increment line-count */
	  line++;
	}
    }

  /* If have at least two anchor positions, then use anchoring */
  if ((nres_anchors + natom_anchors) > 1)
    {
      /* Set flag and show number of anchor-points read in */
      Have_Anchors = TRUE;
      printf("   Number of anchoring residues identified:    %d\n",
	     nres_anchors);
      printf("   Number of atom-anchor positions identified: %d\n",
	     natom_anchors);
    }
}
/***********************************************************************

rotate_centres_of_mass  -  Apply the given rotation to the residue centres
                           of mass

***********************************************************************/

void rotate_centres_of_mass(int natoms, float coord_store[MAXATOMS][3],
			    float v[3][3])
{
  int iatom, icoord;
  float x, y, z;
  float new_coords[3];

  /* Loop over all the centres of mass */
  for (iatom = 0; iatom < natoms; iatom++)
    {
      /* Get the coordinates of this centre of mass */
      x = coord_store[iatom][0];
      y = coord_store[iatom][1];
      z = coord_store[iatom][2];

      /* Apply the transformation */
      for (icoord = 0; icoord < 3; icoord++)
	{
	  /* Calculate transformed coordinates */
	      new_coords[icoord]
		=  ((x * v[0][icoord])
		    + (y * v[1][icoord])
		    + (z * v[2][icoord]));
	}

      /* Save the new coordinates */
      coord_store[iatom][0] = new_coords[0];
      coord_store[iatom][1] = new_coords[1];
      coord_store[iatom][2] = new_coords[2];
    }
}
/***********************************************************************

get_object_means  -  Calculate mean positions of eaach object's original
                     and flattened residues

***********************************************************************/

void get_object_means(struct object *store_object_ptr[MAXATOMS],
		      float coord_store[MAXATOMS][3],
		      float flattened_mean[MAXATOMS][2],int *nstored)
{
  int iatom, iresid, natoms, ncoords, nobj_coords, nresid;

  float flattened_mean_x, flattened_mean_y, mean_x, mean_y, mean_z;
  float flattened_object_mean_x, flattened_object_mean_y;
  float object_mean_x, object_mean_y, object_mean_z;

  struct coordinate *atom_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Initialise the object's average position */
  object_mean_x = 0.0;
  object_mean_y = 0.0;
  object_mean_z = 0.0;
  flattened_object_mean_x = 0.0;
  flattened_object_mean_y = 0.0;
  nobj_coords = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Set the pointer to the first of this object's residues */
      residue_ptr = object_ptr->first_residue_ptr;
      nresid = object_ptr->nresidues;
      iresid = 0;

      /* Loop over all this object's residues */
      while (iresid < nresid && residue_ptr != NULL)
	{
	  /* Initialise the residue's average position */
	  mean_x = 0.0;
	  mean_y = 0.0;
	  mean_z = 0.0;
	  flattened_mean_x = 0.0;
	  flattened_mean_y = 0.0;
	  ncoords = 0;

	  /* Get pointer to first residue's first atom */
	  atom_ptr = residue_ptr->first_atom_ptr;
	  iatom = 0;
	  natoms = residue_ptr->natoms;

	  /* Loop over all this residue's atoms */
	  while (iatom < natoms && atom_ptr != NULL)
	    {
	      /* Update the residue's mean coordinates */
	      mean_x = mean_x + atom_ptr->original_x;
	      mean_y = mean_y + atom_ptr->original_y;
	      mean_z = mean_z + atom_ptr->original_z;
	      ncoords++;

	      /* Update mean coords of the current x-y coords of the
		 atom in the residue's flattened state */
	      flattened_mean_x = flattened_mean_x + atom_ptr->x;
	      flattened_mean_y = flattened_mean_y + atom_ptr->y;

	      /* Repeat for current object */
	      object_mean_x = object_mean_x + atom_ptr->original_x;
	      object_mean_y = object_mean_y + atom_ptr->original_y;
	      object_mean_z = object_mean_z + atom_ptr->original_z;
	      flattened_object_mean_x = flattened_object_mean_x + atom_ptr->x;
	      flattened_object_mean_y = flattened_object_mean_y + atom_ptr->y;
	      nobj_coords++;

	      /* Get pointer to the next atom */
	      atom_ptr = atom_ptr->next;
	      iatom++;
	    }

	  /* Get this residue's mean coordinates */
	  if (ncoords > 0)
	    {
	      /* Get mean coords */
	      mean_x = mean_x / ncoords;
	      mean_y = mean_y / ncoords;
	      mean_z = mean_z / ncoords;

	      /* Calculate the flattened-out residue's mean position, and
		 store the coords */
	      flattened_mean_x = flattened_mean_x / ncoords;
	      flattened_mean_y = flattened_mean_y / ncoords;

	      /* Store the residue's original and flattened mean coords */
	      residue_ptr->flattened_mean_x = flattened_mean_x;
	      residue_ptr->flattened_mean_y = flattened_mean_y;
	    }

	  /* Get pointer to the next residue */
	  residue_ptr = residue_ptr->next_residue_ptr;
	  iresid++;
	}

      /* Get this object's mean coordinates */
      if (nobj_coords > 0)
	{
	  /* Get mean coords */
	  object_mean_x = object_mean_x / nobj_coords;
	  object_mean_y = object_mean_y / nobj_coords;
	  object_mean_z = object_mean_z / nobj_coords;

	  /* Store the current object pointer */
	  store_object_ptr[*nstored] = object_ptr;

	  /* Store the mean coords */
	  coord_store[*nstored][0] = object_mean_x;
	  coord_store[*nstored][1] = object_mean_y;
	  coord_store[*nstored][2] = object_mean_z;

	  /* Calculate the flattened-out object's mean position, and
	     store the coords */
	  flattened_object_mean_x = flattened_object_mean_x / nobj_coords;
	  flattened_object_mean_y = flattened_object_mean_y / nobj_coords;
	  flattened_mean[*nstored][0] = flattened_object_mean_x;
	  flattened_mean[*nstored][1] = flattened_object_mean_y;

	  /* Increment count of stored sets of coords */
	  (*nstored)++;
	}

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }
}
/***********************************************************************

place_interacting_residues  -  Lay out the current object's residues on
                               the appropriate side of the interface

***********************************************************************/

void place_interacting_residues(float coord_store[MAXATOMS][3],
				struct object *store_object_ptr[MAXATOMS],
				int nstored)
{
  int iatom, iresid, istored, natoms, nresid;

  float dx, dy, flattened_mean_x, flattened_mean_y, spread_factor, x_pos;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;
  struct object *object_ptr;

  /* Calculate the spread-factor */
  spread_factor = 2.0;

  /* Loop through all the stored positions to lay out the objects
     according to the principal component analysis */
  for (istored = 0; istored < nstored; istored++)
    {
      /* Get the corresponding object */
      object_ptr = store_object_ptr[istored];

      /* Get the mean position of the flattened object's coords */
      flattened_mean_x = object_ptr->first_residue_ptr->flattened_mean_x;
      flattened_mean_y = object_ptr->first_residue_ptr->flattened_mean_y;

      /* Calculate the adjustment to be made to each of this objects'
	 flattened coords */

      /* The x-coordinate is just the transformed z-coord, with the
	 residues from one domain/chain on the left side of it, and those
	 from the other on the right */
      x_pos = spread_factor * fabs(coord_store[istored][2]);
      if (object_ptr->object_type == WATER)
	x_pos = 0.0;
      else if (object_ptr->interface == 1)
	x_pos = - x_pos;
      dx = x_pos - flattened_mean_x;

      /* The y-coordinate is just the transformed x-coord */
      dy = spread_factor * coord_store[istored][0] - flattened_mean_y;

      /* Set the pointer to the first of this object's residues */
      residue_ptr = object_ptr->first_residue_ptr;
      nresid = object_ptr->nresidues;
      iresid = 0;

      /* Loop over all this object's residues */
      while (iresid < nresid && residue_ptr != NULL)
	{
	  /* Get pointer to first residue's first atom */
	  atom_ptr = residue_ptr->first_atom_ptr;
	  iatom = 0;
	  natoms = residue_ptr->natoms;

	  /* Loop over all this residue's atoms */
	  while (iatom < natoms && atom_ptr != NULL)
	    {
	      /* Transform the atom's coordinates */
	      atom_ptr->x = atom_ptr->x + dx;
	      atom_ptr->y = atom_ptr->y + dy;
	      atom_ptr->z = 0.0;

	      /* Get pointer to the next atom */
	      atom_ptr = atom_ptr->next;
	      iatom++;
	    }

	  /* Shift the flattened means by the same distance that
	     the individual atomic coords have been shifted */
	  residue_ptr->flattened_mean_x = residue_ptr->flattened_mean_x + dx;
	  residue_ptr->flattened_mean_y = residue_ptr->flattened_mean_y + dy;

	  /* Get pointer to the next residue */
	  residue_ptr = residue_ptr->next_residue_ptr;
	  iresid++;
	}
    }
}
/***********************************************************************

lay_interface_objects_out  -  Lay out the interacting objects on either
                              side of an interface

***********************************************************************/

void lay_interface_objects_out(void)
{
  int nstored;

  float coord_store[MAXATOMS][3], flattened_mean[MAXATOMS][2];
  float mass_x, mass_y, mass_z;
  float transformation[3][3];

  struct object *store_object_ptr[MAXATOMS];

  /* Initialise variables */
  nstored = 0;

  /* Get the mean coordinates of the flattened objects */
  get_object_means(store_object_ptr,coord_store,flattened_mean,
		   &nstored);

  /* Calculate the centre of mass of the stored coordinates */
  centre_of_mass(coord_store,nstored,&mass_x,&mass_y,&mass_z);

  /* Translate to mean coordinate to the origin */
  adjust_stored_coords(nstored,coord_store,mass_x,mass_y,mass_z);

  /* Calculate transformation to orient the principal axis along
     the x-direction */
  principal_components(nstored,transformation,coord_store);

  /* Apply this transformation to the centres of mass */
  rotate_centres_of_mass(nstored,coord_store,transformation);

  /* Place all the residues of this object on the appropriate
     side of the interface */
  place_interacting_residues(coord_store,store_object_ptr,nstored);
}
/* <--v.4.0 */
/***********************************************************************

lay_objects_out  -  Lay the objects out on the page in a starting
                    arrangement that approximately reflects the
		    relative interactions

***********************************************************************/

void lay_objects_out(void)
{
  int iatom, iresid, natoms, nresid;

/* v.3.2--> */
  float spread_factor;
/* <--v.3.2 */
  float rmsd;
/* v.4.0--> */
  float x, y;
/* <--v.4.0 */

  struct coordinate *atom_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
/* v.3.2--> */
      /* Calculate the spread-factor according to the object type */
      spread_factor = SPREAD_FACTOR;
      if (object_ptr->object_type == HYDROPHOBIC)
	spread_factor = 10.0 * spread_factor;
/* <--v.3.2 */

      /* If this if a hydrophobic group or simplified H-group, use
         only a single atom to represent it */
      if (object_ptr->object_type == HYDROPHOBIC ||
	  object_ptr->object_type == SIMPLE_HGROUP ||
	  object_ptr->object_type == WATER)
	{
	  /* Set object's position using its placement position */
/* v.3.2--> */
/*	  object_ptr->minx = SPREAD_FACTOR * object_ptr->place_x;
	  object_ptr->maxx = SPREAD_FACTOR * object_ptr->place_x;
	  object_ptr->miny = SPREAD_FACTOR * object_ptr->place_y;
	  object_ptr->maxy = SPREAD_FACTOR * object_ptr->place_y; */
/* v.4.0--> */
/*	  object_ptr->minx = spread_factor * object_ptr->place_x;
	  object_ptr->maxx = spread_factor * object_ptr->place_x;
	  object_ptr->miny = spread_factor * object_ptr->place_y;
	  object_ptr->maxy = spread_factor * object_ptr->place_y; */

	  /* Use the residue's centre-of-mass coords in the original
	     structure (times the spread factor) as starting coords
	     for object */ 
	  x = spread_factor * object_ptr->place_x;
	  y = spread_factor * object_ptr->place_y;

	  /* If have an anchor-position for this object, use that instead */
	  if (object_ptr->have_anchors == TRUE)
	    {
	      /* Extract the anchor position */
	      x = object_ptr->first_residue_ptr->anchor_pstn_x;
	      y = object_ptr->first_residue_ptr->anchor_pstn_y;
	    }

	  /* Place the object */
	  object_ptr->minx = x;
	  object_ptr->maxx = x;
	  object_ptr->miny = y;
	  object_ptr->maxy = y;
/* <--v.4.0 */
/* <--v.3.2 */
	}
      else
	{
	  /* Perform the fit for the current object */
	  rmsd = fit_object(object_ptr,FALSE);
	}

      /* Set the pointer to the first of this object's residues */
      residue_ptr = object_ptr->first_residue_ptr;
      nresid = object_ptr->nresidues;
      iresid = 0;

      /* Loop over all this object's residues */
      while (iresid < nresid && residue_ptr != NULL)
	{
	  /* Get pointer to first residue's first atom */
	  atom_ptr = residue_ptr->first_atom_ptr;
	  iatom = 0;
	  natoms = residue_ptr->natoms;

	  /* Loop over all this residue's atoms */
	  while (iatom < natoms && atom_ptr != NULL)
	    {
	      /* If this if a hydrophobic group or simplified H-group,
		 set all atom coords to the placement position */
	      if (object_ptr->object_type == HYDROPHOBIC ||
		  object_ptr->object_type == SIMPLE_HGROUP ||
		  object_ptr->object_type == WATER)
		{
		  /* Update the atom's coords */
		  atom_ptr->x = object_ptr->minx;
		  atom_ptr->y = object_ptr->miny;
		}

	      /* Otherwise, move the object to its placement position */
	      else
		{
		  /* Move object to its placement position */
		  atom_ptr->x = atom_ptr->x
		    + SPREAD_FACTOR * object_ptr->place_x;
		  atom_ptr->y = atom_ptr->y
		    + SPREAD_FACTOR * object_ptr->place_y;
		}

	      /* Force z-coord to be zero */
	      atom_ptr->z = 0.0;

	      /* Get pointer to the next atom */
	      atom_ptr = atom_ptr->next;
	      iatom++;
	    }

	  /* Get pointer to the next residue */
	  residue_ptr = residue_ptr->next_residue_ptr;
	  iresid++;
	}

      /* Update residue- and object-boundaries */
      update_residue_boundaries(object_ptr);

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }
}
/***********************************************************************

get_atom_clash_energy  -  Calculate the energy for atom-clashes between
                          the given objects

***********************************************************************/

float get_atom_clash_energy(struct object *object1_ptr,
			    struct object *object2_ptr)
{
  int iatom, iresid, jatom, jresid, natoms1, nresid1, nresid2, natoms2;
  int object_type1, object_type2;
  int in_range;

  float atom_size1, atom_size2, distance, dist2, interact_dist;
  float energy, total_energy;
  float x1, x2, y1, y2;

  struct coordinate *atom1_ptr, *atom2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise variables */
  total_energy = 0.0;

  /* Get the two object types */
  object_type1 = object1_ptr->object_type;
  object_type2 = object2_ptr->object_type;

  /* Set the pointer to the first object's first residue */
  residue1_ptr = object1_ptr->first_residue_ptr;
  nresid1 = object1_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid1 && residue1_ptr != NULL)
    {
      /* Set the pointer to the second object's first residue */
      residue2_ptr = object2_ptr->first_residue_ptr;
      nresid2 = object2_ptr->nresidues;
      jresid = 0;

      /* Loop over other object's residues */
      while (jresid < nresid2 && residue2_ptr != NULL)
	{
	  /* Get the maximum atom-size of each residue */
	  atom_size1 = residue1_ptr->max_atom_size;
	  atom_size2 = residue2_ptr->max_atom_size;
	  interact_dist = atom_size1 + atom_size2 + ATOM_INTERACT_DIST;

	  /* Check whether the residues are within range of one another */
	  in_range = TRUE;
	  if ((residue2_ptr->minx > (residue1_ptr->maxx + interact_dist)) ||
	      (residue1_ptr->minx > (residue2_ptr->maxx + interact_dist)) ||
	      (residue2_ptr->miny > (residue1_ptr->maxy + interact_dist)) ||
	      (residue1_ptr->miny > (residue2_ptr->maxy + interact_dist)))
	    in_range = FALSE;

/* v.3.2--> */
	  /* If either residue doesn't have any atoms, then skip */
	  if (residue1_ptr->natoms == 0 || residue2_ptr->natoms == 0)
	    in_range = FALSE;
/* <--v.3.2 */

	  /* If residues are within range, then check all distances
	     between them */
	  if (in_range == TRUE)
	    {
	      /* Get pointer to first residue's first atom */
	      atom1_ptr = residue1_ptr->first_atom_ptr;

	      /* Initialise atom count */
	      iatom = 0;
	      natoms1 = residue1_ptr->natoms;

	      /* If this object is a hydrophobic group or a simplified
		 residue, then need only consider a single position */
	      if (object_type1 == HYDROPHOBIC ||
		  object_type1 == SIMPLE_HGROUP)
		{
		  x1 = object1_ptr->minx;
		  y1 = object1_ptr->miny;
		  atom_size1 = object1_ptr->max_atom_size;
		  natoms1 = 1;
		}

	      /* Loop over all the first residue's atoms */
	      while (iatom < natoms1 && atom1_ptr != NULL)
		{
		  /* Get atom's coordinates and size */
		  if (object_type1 != HYDROPHOBIC &&
		      object_type1 != SIMPLE_HGROUP)
		    {
		      /* Get atom's coordinates */
		      x1 = atom1_ptr->x;
		      y1 = atom1_ptr->y;

		      /* Get the atom's size */
		      atom_size1 = atom1_ptr->atom_size;
		    }

		  /* Get pointer to other residue's first atom */
		  atom2_ptr = residue2_ptr->first_atom_ptr;
		  
		  /* Initialise atom count */
		  jatom = 0;
		  natoms2 = residue2_ptr->natoms;
		  
		  /* If this object is a hydrophobic group or a simplified
		     residue, then need only consider a single position */
		  if (object_type2 == HYDROPHOBIC ||
		      object_type2 == SIMPLE_HGROUP)
		    {
		      x2 = object2_ptr->minx;
		      y2 = object2_ptr->miny;
		      atom_size2 = object2_ptr->max_atom_size;
		      natoms2 = 1;
		    }

		  /* Loop over all the other residue's atoms */
		  while (jatom < natoms2 && atom2_ptr != NULL)
		    {
		      /* Get atom's coordinates and size */
		      if (object_type2 != HYDROPHOBIC &&
			  object_type2 != SIMPLE_HGROUP)
			{
			  /* Get other atom's coordinates */
			  x2 = atom2_ptr->x;
			  y2 = atom2_ptr->y;

			  /* Get the other atom's size */
			  atom_size2 = atom2_ptr->atom_size;
			}

		      /* Calculate the distance between the
			 two atoms */
		      dist2 = (x1 - x2) * (x1 - x2)
			+ (y1 - y2) * (y1 - y2);
		      distance = sqrt((double) dist2);

		      /* Adjust distance to allow for the atom sizes */
		      if (distance > (atom_size1 + atom_size2))
			distance = distance - (atom_size1 + atom_size2);
		      else
			distance = 0.0;
		      dist2 = distance * distance;

		      /* Add energy of this interaction to the
			 overall energy sum */
		      if (dist2 < 0.001)
			dist2 = 0.001;
		      if (dist2 < ATOM_DIST2)
			energy = Atom_Atom_Clash / dist2;
		      else
			energy = 0.0;
		      total_energy = total_energy + energy;

		      /* Get pointer to other residue's next atom */
		      atom2_ptr = atom2_ptr->next;
		      jatom++;
		    }

		  /* Get pointer to the next atom */
		  atom1_ptr = atom1_ptr->next;
		  iatom++;
		}
	    }

          /* Get pointer to the next residue */
          residue2_ptr = residue2_ptr->next_residue_ptr;
          jresid++;
	}

      /* Get pointer to the next residue */
      residue1_ptr = residue1_ptr->next_residue_ptr;
      iresid++;
    }

  /* Return the calculated energy */
  return(total_energy);
}
/***********************************************************************

get_distance_energy  -  Calculate the energy given the bond-atom distance
                        and the atom-size

***********************************************************************/

float get_distance_energy(float distance, float atom_size)
{
  float energy;
  
  /* Adjust distance to allow for the atom sizes */
  if (distance > atom_size)
    distance = distance - atom_size;
  else
    distance = 0.0;

  /* If the distance is within the interaction distance, then
     compute bond-atom clash energy */
  if (distance < 0.03)
    distance = 0.03;
  if (distance < BOND_INTERACT_DIST)
    energy = Bond_Atom_Clash / (distance * distance);
  else
    energy = 0.0;

  /* Return the calculated energy */
  return(energy);
}
/***********************************************************************

get_bond_clash_energy  -  Calculate the energy for atom-bond clashes
                          between the given bond and all objects/atoms

***********************************************************************/

float get_bond_clash_energy(struct bond *current_bond_ptr,
			    float matrix[4][3],int external_bond,
			    struct object *current_object_ptr)
{
  int iatom, iobject, iresid, natoms, nresid;
  int internal_bond, object_type;
  int in_range, want_atom, wanted;

  float atom_size, distance, interact_dist, length;
  float energy, total_energy;
  float max_x, max_y, min_x, min_y, x, x1, x2, y, y1, y2;
/* v.3.1--> */
/*  float conv_x, conv_x1, conv_x2, conv_y, conv_y1, conv_y2, conv_z; */
  float conv_x, conv_x1, conv_x2, conv_y, conv_y1, conv_y2;
/* <--v.3.1 */

  struct coordinate *atom_ptr, *atom1_ptr, *atom2_ptr;
  struct object *object_ptr, *object1_ptr, *object2_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  total_energy = 0.0;

  /* Get the pointers to the two atoms at either end of the bond */
  atom1_ptr = current_bond_ptr->first_atom_ptr;
  atom2_ptr = current_bond_ptr->second_atom_ptr;

  /* Get the corresponding objects */
  object1_ptr = atom1_ptr->residue_ptr->object_ptr;
  object2_ptr = atom2_ptr->residue_ptr->object_ptr;

  /* Check whether this is an internal bond */
  internal_bond = FALSE;
  if (object1_ptr == object2_ptr)
    internal_bond = TRUE;

  /* Get the coordinates of the two atoms */
  x1 = atom1_ptr->x;
  y1 = atom1_ptr->y;
  x2 = atom2_ptr->x;
  y2 = atom2_ptr->y;

  /* Determine max and min extents of the bond */
  min_x = x1;
  max_x = x1;
  min_y = y1;
  max_y = y1;
  if (x2 < min_x)
    min_x = x2;
  if (x2 > max_x)
    max_x = x2;
  if (y2 < min_y)
    min_y = y2;
  if (y2 > max_y)
    max_y = y2;

  /* Apply transformation to the atom coords */
/* v.3.1--> */
/*  apply_transformation(x1,y1,0.0,&conv_x1,&conv_y1,&conv_z,matrix);
  apply_transformation(x2,y2,0.0,&conv_x2,&conv_y2,&conv_z,matrix); */
  apply_xy_transformation(x1,y1,&conv_x1,&conv_y1,matrix);
  apply_xy_transformation(x2,y2,&conv_x2,&conv_y2,matrix);
/* <--v.3.1 */

  /* Calculate the bond length */
  length = conv_x2;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Set the pointer to the object's first residue */
      residue_ptr = object_ptr->first_residue_ptr;
      nresid = object_ptr->nresidues;
      iresid = 0;
      wanted = TRUE;

      /* If this if a hydrophobic group or simplified H-group, then
         is represented by a single atom */
      object_type = object_ptr->object_type;
      if (object_type == HYDROPHOBIC || object_type == SIMPLE_HGROUP)
	{
	  /* Get atom size and coordinates */
	  x = object_ptr->minx;
	  y = object_ptr->miny;
	  nresid = 1;
	  natoms = 1;

	  /* If bond is attached to this object at either end, then
	     don't want to consider it */
	  if (object_ptr == object1_ptr || object_ptr == object2_ptr)
	    wanted = FALSE;
	}

      /* If the current bond is an internal one to a single object,
	 then don't need to consider its clashes with other atoms in
	 this object */
      if (internal_bond == TRUE && object_ptr == object1_ptr)
	wanted = FALSE;

      /* If the bond we're looking at is an external bond, then are
	 only interested in its clashes with the object of interest */
      if (external_bond == TRUE && object_ptr != current_object_ptr)
	wanted = FALSE;

      /* Loop over all the object's residues */
      while (iresid < nresid && residue_ptr != NULL && wanted == TRUE)
	{
	  /* Get the maximum atom-size of this residue */
	  atom_size = residue_ptr->max_atom_size;
	  interact_dist = atom_size + BOND_INTERACT_DIST;

	  /* Check whether the residues is within range of the bond */
	  in_range = TRUE;
	  if ((min_x > (residue_ptr->maxx + interact_dist)) ||
	      (max_x < (residue_ptr->minx - interact_dist)) ||
	      (min_y > (residue_ptr->maxy + interact_dist)) ||
	      (max_y < (residue_ptr->miny - interact_dist)))
	    in_range = FALSE;

	  /* If residue within range of the bond, then check all its
	     atoms against the bond */
	  if (in_range == TRUE)
	    {
	      /* Get pointer to residue's first atom */
	      atom_ptr = residue_ptr->first_atom_ptr;

	      /* Initialise atom count */
	      iatom = 0;
	      natoms = residue_ptr->natoms;

	      /* If object is a hydrophobic group or a simplified
		 residue, then only have a single position */
	      if (object_type == HYDROPHOBIC ||
		  object_type == SIMPLE_HGROUP)
		natoms = 1;

	      /* Loop over all the residue's atoms */
	      while (iatom < natoms && atom_ptr != NULL)
		{
		  /* Get atom's coordinates and size */
		  if (object_type != HYDROPHOBIC &&
		      object_type != SIMPLE_HGROUP)
		    {
		      /* Get atom's coordinates */
		      x = atom_ptr->x;
		      y = atom_ptr->y;

		      /* Get the atom's size */
		      atom_size = atom_ptr->atom_size;
		    }

		  /* If current atom is at either end of the bond, then
		     don't want it */
		  want_atom = TRUE;
		  if (atom_ptr == atom1_ptr || atom_ptr == atom2_ptr)
		    want_atom = FALSE;

		  /* If atom is wanted, then calculate bond-clash energy */
		  if (want_atom == TRUE)
		    {
		      /* Initialise energy */
		      energy = 0.0;

		      /* Apply transformation to atom's coordinates */
/* v.3.1--> */
/*		      apply_transformation(x,y,0.0,&conv_x,&conv_y,
					   &conv_z,matrix); */
		      apply_xy_transformation(x,y,&conv_x,&conv_y,
					      matrix);
/* <--v.3.1 */
		      x = conv_x;
		      y = conv_y;

		      /* If atom's x-value lies within the length of
			 the bond, take its y-value as the distance */
		      if (x >= 0.0 && x <= length)
			distance = fabs(y);

		      /* If atom is to the left of the bond, then get
			 its distance to the origin (where one end of
			 the bond currently lies) */
		      else if (x < 0)
			distance = sqrt(x * x + y * y);

		      /* Otherwise, calculate distance to the other
			 end of the bond */
		      else
			distance = sqrt((x - length) * (x - length)
					+ y * y);

		      /* Calculate the energy based on this distance 
			 and the atom size */
		      energy = get_distance_energy(distance,atom_size);

		      /* Add to total */
		      total_energy = total_energy + energy;
		    }

		  /* Get pointer to the next atom */
		  atom_ptr = atom_ptr->next;
		  iatom++;
		}
	    }
	  
          /* Get pointer to the next residue */
          residue_ptr = residue_ptr->next_residue_ptr;
          iresid++;
	}

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

  /* Return the calculated energy */
  return(total_energy);
}
/***********************************************************************

get_bond_overlap_energy  -  Calculate the energy for bond overlaps
                            between the given bond and all other bonds

***********************************************************************/

float get_bond_overlap_energy(struct bond *current_bond_ptr,
			      float matrix[4][3],int external_bond,
			      struct object *current_object_ptr)
{
  int bond_type, current_bond_type;
  int in_range, wanted;
/* v.3.1--> */
  int point_object1, point_object2;
/* <--v.3.1 */

  float length;
  float overlap_energy, total_energy;
/* v.3.1.1--> */
  /*  float max_x, max_y, min_x, min_y, x1, x2, y1, y2; */
  float x1, x2, y1, y2;
/* <--v.3.1.1 */
  float other_x1, other_x2, other_y1, other_y2;
/* v.3.1--> */
/*  float conv_x, conv_y, conv_z, cut_x; */
  float conv_x, conv_y, cut_x;
/* <--v.3.1 */
  float energy;

  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct coordinate *other_atom1_ptr, *other_atom2_ptr;
  struct object *object1_ptr, *object2_ptr;
  struct object *other_object1_ptr, *other_object2_ptr;

  /* Initialise variables */
  overlap_energy = Overlap_Score;
  total_energy = 0.0;

  /* Get the current bond's type */
  current_bond_type = current_bond_ptr->bond_type;

  /* Get the pointers to the two atoms at either end of the bond */
  atom1_ptr = current_bond_ptr->first_atom_ptr;
  atom2_ptr = current_bond_ptr->second_atom_ptr;

  /* Get the corresponding objects */
  object1_ptr = atom1_ptr->residue_ptr->object_ptr;
  object2_ptr = atom2_ptr->residue_ptr->object_ptr;

/* v.3.1--> */
  /* Check whether either object at the end of the bond is a point
     object */
  point_object1 = FALSE;
  if (object1_ptr->object_type == HYDROPHOBIC ||
      object1_ptr->object_type == WATER ||
      object1_ptr->object_type == SIMPLE_HGROUP)
    point_object1 = TRUE;
  point_object2 = FALSE;
  if (object2_ptr->object_type == HYDROPHOBIC ||
      object2_ptr->object_type == WATER ||
      object2_ptr->object_type == SIMPLE_HGROUP)
    point_object2 = TRUE;
/* <--v.3.1 */

  /* Get the transformed coordinates of the two atoms */
  x1 = atom1_ptr->x;
  y1 = atom1_ptr->y;
  x2 = atom2_ptr->x;
  y2 = atom2_ptr->y;

  /* Apply transformation to the atom coords */
/* v.3.1--> */
/*  apply_transformation(x1,y1,0.0,&conv_x,&conv_y,&conv_z,matrix); */
  apply_xy_transformation(x1,y1,&conv_x,&conv_y,matrix);
/* <--v.3.1 */
  x1 = conv_x;
  y1 = conv_y;
/* v.3.1--> */
/*  apply_transformation(x2,y2,0.0,&conv_x,&conv_y,&conv_z,matrix); */
  apply_xy_transformation(x2,y2,&conv_x,&conv_y,matrix);
/* <--v.3.1 */
  x2 = conv_x;
  y2 = conv_y;

  /* Calculate the bond length */
  length = x2;

  /* Get pointer to first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds, one by one */
  while (bond_ptr != NULL)
    {
      /* Initialise flag giving whether this bond is wanted */
      wanted = TRUE;

      /* Get the bond type */
      bond_type = bond_ptr->bond_type;

      /* Check that this isn't the current bond */
      if (bond_ptr == current_bond_ptr)
	wanted = FALSE;

      /* Not interested in internal H-bonds  */
      if (bond_type == INTERNAL || bond_type == DELETED)
	wanted = FALSE;

      /* Get the pointers to the two atoms at either end of the bond */
      other_atom1_ptr = bond_ptr->first_atom_ptr;
      other_atom2_ptr = bond_ptr->second_atom_ptr;
	  
      /* Check that the two bonds don't have any atoms in common */
      if (atom1_ptr == other_atom1_ptr || atom1_ptr == other_atom2_ptr ||
	  atom2_ptr == other_atom1_ptr || atom2_ptr == other_atom2_ptr)
	wanted = FALSE;

      /* Get the corresponding objects */
      other_object1_ptr = other_atom1_ptr->residue_ptr->object_ptr;
      other_object2_ptr = other_atom2_ptr->residue_ptr->object_ptr;

      /* If bonds both spring from a hydrophobic group, or
	 a simplified residue, then don't want to check for
	 an overlap */
/* v.3.1--> */
/*      if ((object1_ptr == other_object1_ptr || 
	   object1_ptr == other_object2_ptr) &&
	  (object1_ptr->object_type == HYDROPHOBIC ||
	   object1_ptr->object_type == WATER ||
	   object1_ptr->object_type == SIMPLE_HGROUP))
	wanted = FALSE;
      if ((object2_ptr == other_object1_ptr || 
	   object2_ptr == other_object2_ptr) &&
	  (object2_ptr->object_type == HYDROPHOBIC ||
	   object2_ptr->object_type == WATER ||
	   object2_ptr->object_type == SIMPLE_HGROUP))
	wanted = FALSE; */
      if (point_object1 == TRUE &&
	  (object1_ptr == other_object1_ptr || 
	   object1_ptr == other_object2_ptr))
	wanted = FALSE;
      if (point_object2 == TRUE &&
	  (object2_ptr == other_object1_ptr || 
	   object2_ptr == other_object2_ptr))
	wanted = FALSE;
/* <--v.3.1 */

/* v.3.2--> */
      /* If either bond is a contact bond which doesn't spring from a
	 hydrophobic group (eg is between ligand and an H-group), then
	 not interested in it */
      if (current_bond_type == CONTACT &&
	  object1_ptr->object_type != HYDROPHOBIC &&
	  object2_ptr->object_type != HYDROPHOBIC)
	wanted = FALSE;
      if (bond_type == CONTACT &&
	  other_object1_ptr->object_type != HYDROPHOBIC &&
	  other_object2_ptr->object_type != HYDROPHOBIC)
	wanted = FALSE;
/* <--v.3.2 */

      /* If the bond we're looking at is an external bond,
	 then are only interested in its overlaps with bonds
	 involving the object of interest */
/* v.3.1--> */
/*      if (external_bond == TRUE && */
      if (external_bond == TRUE &&
/* <--v.3.1 */
	  other_object1_ptr != current_object_ptr &&
	  other_object2_ptr != current_object_ptr)
	wanted = FALSE;

      /* If bond is still wanted, then check if within range */
      if (wanted == TRUE)
	{
	  /* Get the coordinates of the two atoms */
	  other_x1 = other_atom1_ptr->x;
	  other_y1 = other_atom1_ptr->y;
	  other_x2 = other_atom2_ptr->x;
	  other_y2 = other_atom2_ptr->y;

	  /* Apply transformation to the atom coords */
/* v.3.1--> */
/*	  apply_transformation(other_x1,other_y1,0.0,&conv_x,&conv_y,
			       &conv_z,matrix); */
	  apply_xy_transformation(other_x1,other_y1,&conv_x,&conv_y,
				  matrix);
/* <--v.3.1 */
	  other_x1 = conv_x;
	  other_y1 = conv_y;
/* v.3.1--> */
/*	  apply_transformation(other_x2,other_y2,0.0,&conv_x,&conv_y,
			       &conv_z,matrix); */
	  apply_xy_transformation(other_x2,other_y2,&conv_x,&conv_y,
				  matrix);
/* <--v.3.1 */
	  other_x2 = conv_x;
	  other_y2 = conv_y;

	  /* Determine max and min extents of the bond */
/* v.3.1--> */
/*	  min_x = other_x1;
	  max_x = other_x1;
	  min_y = other_y1;
	  max_y = other_y1;
	  if (other_x2 < min_x)
	    min_x = other_x2;
	  if (other_x2 > max_x)
	    max_x = other_x2;
	  if (other_y2 < min_y)
	    min_y = other_y2;
	  if (other_y2 > max_y)
	    max_y = other_y2; */
/* <--v.3.1 */

	  /* Check whether the bonds are within range of one another */
	  in_range = TRUE;
/* v.3.1--> */
/*	  if (min_x > x2 || max_x < 0.0 || min_y > 0.0 || max_y < 0.0)
	    in_range = FALSE; */
	  if (other_y1 < 0.0 && other_y2 < 0.0)
	    in_range = FALSE;
	  else if (other_y1 > 0.0 && other_y2 > 0.0)
	    in_range = FALSE;
	  else if (other_x1 < 0.0 && other_x2 < 0.0)
	    in_range = FALSE;
	  else if (other_x1 > x2 && other_x2 > x2)
	    in_range = FALSE;
/* <--v.3.1 */

	  /* If the two bonds are within range of one another, calculate
	     their cross-over point */
	  if (in_range == TRUE)
	    {
	      /* If second bond is not parallel to the x-axis calculate
		 where it cuts the x-axis */
	      if (fabs(other_y2 - other_y1) > 0.01)
		{	    
		  cut_x = other_x2 - other_y2 * (other_x2 - other_x1)
		      / (other_y2 - other_y1);

		  /* If the cut-point is within the length of the first
		     bond, then the bonds cross one another */
		  if (cut_x >= 0.0 && cut_x <= length)
		    {
		      if (current_bond_type == CONTACT ||
			  bond_type == CONTACT)
			energy = 0.0005;
		      else
			energy = overlap_energy;
		      total_energy = total_energy + energy;
		    }
		}
	    }
	}

      /* Get pointer to the next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Return the calculated energy */
  return(total_energy);
}
/***********************************************************************

calculate_energy  -  Calculate the energy between the current object and
                     all the others

***********************************************************************/

float calculate_energy(struct object *current_object_ptr,
		       float *atom_clash,float *bond_length,
/* v.4.0--> */
/*		       float *bond_clash,float *bond_overlap) */
		       float *bond_clash,float *bond_overlap,
		       float *boundary_energy,float *rel_pstn_energy)
/* <--v.4.0 */
{
  int bond_type;
  int external_bond, in_range, to_ligand, wanted;

  float coord_store[3][3], matrix[4][3];
  float atom_size1, atom_size2, energy, interact_dist;
  float atom_clash_energy, total_atom_clash_energy;
  float bond_length_energy, total_bond_length_energy;
  float bond_clash_energy, total_bond_clash_energy;
  float bond_overlap_energy, total_bond_overlap_energy;
  float internal_energy;
  float diff, length, x1, x2, y1, y2;

  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct object *object_ptr, *object1_ptr, *object2_ptr;

  /* Initialise variables */
  energy = 0.0;
  total_atom_clash_energy = 0.0;
  total_bond_length_energy = 0.0;
  total_bond_clash_energy = 0.0;
  total_bond_overlap_energy = 0.0;
/* v.4.0--> */
  *rel_pstn_energy = 0.0;
  *boundary_energy = 0.0;
/* <--v.4.0 */

  /* Get this object's internal energy */
  internal_energy = current_object_ptr->internal_energy;  

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects to calculate the atom-atom energies */
  while (object_ptr != NULL)
    {
      /* If this is not the same as the current object, calculate
	 the energy between them */
      if (object_ptr != current_object_ptr)
	{
	  /* Get the maximum atom-size of each object */
	  atom_size1 = current_object_ptr->max_atom_size;
	  atom_size2 = object_ptr->max_atom_size;
	  interact_dist = atom_size1 + atom_size2 + ATOM_INTERACT_DIST;

	  /* Check whether the objects are within range of one another */
	  in_range = TRUE;
	  if ((object_ptr->minx >
	       (current_object_ptr->maxx + interact_dist)) ||
	      (current_object_ptr->minx >
	       (object_ptr->maxx + interact_dist)) ||
	      (object_ptr->miny >
	       (current_object_ptr->maxy + interact_dist)) ||
	      (current_object_ptr->miny >
	       (object_ptr->maxy + interact_dist)))
	    in_range = FALSE;

	  /* If objects are within range, then check all distances
	     between their constituent residues */
	  if (in_range == TRUE)
	    {
	      /* Get the energy between the current pair of objects */
	      atom_clash_energy
		= get_atom_clash_energy(current_object_ptr,object_ptr);

	      /* Add the energy to the current object's total energy */
	      total_atom_clash_energy = total_atom_clash_energy
		+ atom_clash_energy;
	    }
	}

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* Loop through all the elastic bonds involving the current object
     to get their contributions to the energy */

  /* Get pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds, one by one */
  while (bond_ptr != NULL)
    {
      /* Initialise flag giving whether this bond is wanted */
      to_ligand = FALSE;
      wanted = TRUE;

      /* Get the bond type */
      bond_type = bond_ptr->bond_type;

      /* Only interested in external elastic bonds */
      if (bond_ptr->elastic == FALSE || bond_type == INTERNAL)
	wanted = FALSE;

      /* Get the pointers to the two atoms at either end of the bond */
      atom1_ptr = bond_ptr->first_atom_ptr;
      atom2_ptr = bond_ptr->second_atom_ptr;

      /* Get the objects these two atoms belong to */
      object1_ptr = atom1_ptr->residue_ptr->object_ptr;
      object2_ptr = atom2_ptr->residue_ptr->object_ptr;

      /* If this bond does not involve the current object, then
	 mark it as an external one */
      external_bond = FALSE;
      if (object1_ptr != current_object_ptr &&
	  object2_ptr != current_object_ptr)
	external_bond = TRUE;

      /* If this a wanted bond, get coords of its two ends */
      if (wanted == TRUE)
	{
	  /* Get the atom coordinates */
	  x1 = atom1_ptr->x;
	  y1 = atom1_ptr->y;
	  x2 = atom2_ptr->x;
	  y2 = atom2_ptr->y;

/* v.3.1--> */
	  /* If this is an external bond, check whether it is too far
	     from the current object to ever interact with it */
	  if (external_bond == TRUE)
	    {
	      /* Get maximum interaction distance */
	      interact_dist = current_object_ptr->max_atom_size
		+ BOND_INTERACT_DIST;

	      /* Check whether bond might be within range of the current
		 object */
	      if (x1 > current_object_ptr->maxx + interact_dist &&
		  x2 > current_object_ptr->maxx + interact_dist)
		wanted = FALSE;
	      else if (x1 < current_object_ptr->minx - interact_dist &&
		       x2 < current_object_ptr->minx - interact_dist)
		wanted = FALSE;
	      else if (y1 > current_object_ptr->maxy + interact_dist &&
		       y2 > current_object_ptr->maxy + interact_dist)
		wanted = FALSE;
	      else if (y1 < current_object_ptr->miny - interact_dist &&
		       y2 < current_object_ptr->miny - interact_dist)
		wanted = FALSE;
	    }
	}

      /* If bond still wanted, then process it */
      if (wanted == TRUE)
	{
/* <--v.3.1 */
	  /* Store the coordinates of the two atoms defining this bond */
	  coord_store[0][0] = x1;
	  coord_store[0][1] = y1;
	  coord_store[0][2] = 0.0;
	  coord_store[1][0] = x2;
	  coord_store[1][1] = y2;
	  coord_store[1][2] = 0.0;
	      
	  /* Calculate the transformation that will put this bond
	     along the x-axis, with atom 1 at the origin and atom 2
	     pointing in the +ve x-direction */
	  line_transformation(coord_store,matrix,0.0);

	  /* If this is any bond other than an internal one, calculate
	     its length and compare with actual length */
	  if ((bond_type == HBOND || bond_type == CONTACT ||
	      (bond_type == COVALENT && bond_ptr->elastic == TRUE)) &&
	      external_bond == FALSE)
	    {
	      /* Calculate length of bond */
	      length = sqrt((x2 - x1) * (x2 - x1)
				 + (y2 - y1) * (y2 - y1));

	      /* Get the difference from the bond's original length (times
	         the required "stretch factor" which determines how
		 closely the H-bond groups are to be packed in around
		 the ligand) */
	      if (bond_type == COVALENT)
		diff = length - bond_ptr->bond_length;
	      else
		diff = length - Bond_Stretch * bond_ptr->bond_length;
	      diff = diff * diff;

	      /* Calculate the energy due to the discrepancy in the
		 bond's length */
	      if (bond_type == COVALENT)
		bond_length_energy = 10.0 * diff;
	      else if (bond_type == HBOND)
		bond_length_energy = HB_Weight * diff;
	      else
/* v.3.2--> */
/*		bond_length_energy = Nonbond_Weight * diff; */
	      {
		bond_length_energy = Nonbond_Weight * diff;
		if (object1_ptr->object_type != HYDROPHOBIC &&
		    object2_ptr->object_type != HYDROPHOBIC)
		  bond_length_energy = bond_length_energy / 100.0;
	      }
/* <--v.3.2 */

	      /* Check if either atom belongs to the ligand */
	      if (atom1_ptr->residue_ptr->inligand == TRUE ||
		  atom2_ptr->residue_ptr->inligand == TRUE)
		to_ligand = TRUE;

	      /* If the bond involves the ligand, then double its
		 weight */
	      if (bond_type == HBOND && to_ligand == TRUE)
		bond_length_energy = 3 * bond_length_energy;

	      /* Add to total bond-length energy for this object */
	      total_bond_length_energy = total_bond_length_energy
		+ bond_length_energy;
	    }

	  /* Calculate energy from any clashes between the bond and atoms */
	  if (bond_type == CONTACT)
	      bond_clash_energy = 0.0;
	  else
	    {
	      bond_clash_energy
		= get_bond_clash_energy(bond_ptr,matrix,external_bond,
					current_object_ptr);
	    }

	  /* Calculate energy from bond overlaps */
	  bond_overlap_energy
	    = get_bond_overlap_energy(bond_ptr,matrix,external_bond,
				      current_object_ptr);

	  /* Update energy totals */
	  total_bond_clash_energy = total_bond_clash_energy
	    + bond_clash_energy;
	  total_bond_overlap_energy = total_bond_overlap_energy
	    + bond_overlap_energy;
	}

      /* Get pointer to the next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

/* v.4.0--> */
  /* Calculate object's anchor energy or its relative position energy */
  if (Have_Anchors == TRUE && current_object_ptr->have_anchors == TRUE)
    *rel_pstn_energy = get_anchor_energy(current_object_ptr);
  else if (Rel_Pstn_Energy_Weight > 0.0)
    *rel_pstn_energy = get_rel_pstn_energy(current_object_ptr);

  /* For interface plot, calculate boundary energy of object */
  if (Interface_Plot == TRUE)
    *boundary_energy = get_boundary_energy(current_object_ptr);
/* <--v.4.0 */

  /* Calculate the total energy from the sum of its contributions */
  energy = total_atom_clash_energy + total_bond_length_energy
/* v.4.0--> */
/*    + total_bond_clash_energy + total_bond_overlap_energy; */
    + total_bond_clash_energy + total_bond_overlap_energy
      + (*rel_pstn_energy);
/* <--v.4.0 */

  /* Store this object's total energy */
/* v.4.0--> */
/*  current_object_ptr->total_energy = internal_energy + energy / 2.0; */
  current_object_ptr->total_energy = internal_energy + (*boundary_energy)
    + energy / 2.0;
/* <--v.4.0 */
  energy = current_object_ptr->total_energy;

  /* Return the calculated energies */
  *atom_clash = total_atom_clash_energy;
  *bond_length = total_bond_length_energy;
  *bond_clash = total_bond_clash_energy;
  *bond_overlap = total_bond_overlap_energy;
  return(energy);
}
/***********************************************************************

make_random_move  -  Apply a random translation to the given object

***********************************************************************/

void make_random_move(struct object *object_ptr,float max_move_size,
		      int *rseed)
{
  int iseed;
  int got_dir;

  float length, length2, move_size, move_x, move_y;
  float x, y, z;

  /* Get the current random number generator seed */
  iseed = *rseed;

  /* Get the random number direction for the move */
  got_dir = FALSE;
  while (got_dir == FALSE)
    {
      /* Get trial vector for direction of move */
      x = get_random_number(&iseed,Random_Start);
      y = get_random_number(&iseed,Random_Start);

      /* Convert so that values lie between -1 and +1 */
      x = 2.0 * x - 1.0;
      y = 2.0 * y - 1.0;

      /* Check that the point representing this direction lies within
	 a unit circle (to ensure direction is genuinely random */
      length2 = x * x + y * y;
      if (length2 <= 1.0)
	got_dir = TRUE;
    }

  /* Calculate the move size */
  z = get_random_number(&iseed,Random_Start);
  move_size = z * max_move_size;

  /* Make the size of the move equal to the required size */
  length = sqrt(length2);
  move_x = x * move_size / length;
  move_y = y * move_size / length;

  /* Apply the move to the current object */
  move_object(object_ptr,move_x,move_y,0.0);

  /* Return the current random number seed */
  *rseed = iseed;
}
/***********************************************************************

make_random_rotation  -  Apply a random rotation or swing to the given
                         object

***********************************************************************/

void make_random_rotation(struct object *object_ptr,float max_angle,
                          int *rseed,int swing)
{
  int iatom, natoms;
  int i, iseed, j, nstored;
  int swing_about_origin, wanted;

  float angle, cos_theta, sin_theta, matrix[3][3], x, y;
  float random_no;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr, *atom1_ptr, *atom2_ptr;
  struct coordinate *stored_atom_ptr[MAXATOMS];
  struct object *object1_ptr, *object2_ptr;
  struct residue *residue_ptr;

  /* Get the current random number generator seed */
  iseed = *rseed;
  x = y = 0.0;
  nstored = 0;

  /* Initialise both matrices */
  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  if (i == j)
	    matrix[i][j] = 1.0;
	  else
	    matrix[i][j] = 0.0;
	  /*	  inverse_matrix[i][j] = matrix[i][j];  */
	}
    }

  /* Generate a random number */
  random_no = get_random_number(&iseed,Random_Start);

  /* Convert so that values lie between plus or minus of the 
     maximum rotation angle */
  angle = 2.0 * max_angle * random_no - max_angle;

  /* Determine whether to swing about the origin */
  swing_about_origin = FALSE;
  if (swing == TRUE)
    {
      random_no = get_random_number(&iseed,Random_Start);
      if (random_no < 0.1)
	swing_about_origin = TRUE;
    }

  /* Calculate sines and cosines of the angle */
  cos_theta = cos(angle);
  sin_theta = sin(angle);

  /* Form the appropriate rotation matrix */
  matrix[0][0] = cos_theta;
  matrix[1][0] = sin_theta;
  matrix[0][1] = - sin_theta;
  matrix[1][1] = cos_theta;

  /* Form the inverse matrix */
  /*  inverse_matrix[0][0] = cos_theta;
  inverse_matrix[0][1] = - sin_theta;
  inverse_matrix[1][0] = sin_theta;
  inverse_matrix[1][1] = cos_theta; */

  /* If performing a swing, get all the bonds from this object to
     other objects */
  if (swing == TRUE && swing_about_origin == FALSE)
    {
      /* Initialise pointer to start of linked list of bonds */
      bond_ptr = first_bond_ptr;

      /* Loop through all the bonds in the linked list */
      while (bond_ptr != NULL)
	{
	  /* Check if this is an elastic bond */
	  if (bond_ptr->elastic == TRUE)
	    {
	      /* Get the pointers to the two atoms at either end
		 of the bond, and the objects they belong to */
	      atom1_ptr = bond_ptr->first_atom_ptr;
	      atom2_ptr = bond_ptr->second_atom_ptr;
	      object1_ptr = atom1_ptr->residue_ptr->object_ptr;
	      object2_ptr = atom2_ptr->residue_ptr->object_ptr;

	      /* Check whether one atom is on the current object and
		 the other is not */
	      wanted = FALSE;
	      if ((object1_ptr == object_ptr && object2_ptr != object_ptr) ||
		  (object2_ptr == object_ptr && object1_ptr != object_ptr))
		wanted = TRUE;

	      /* For H-group, prefer to swing about an H-bond */
	      if (object_ptr->object_type == HGROUP &&
		  bond_ptr->bond_type != HBOND)
		wanted = FALSE;

	      /* If this is a possible bond, then store pointer to the
		 atom at the other end of it */
	      if (wanted == TRUE && nstored < MAXATOMS)
		{
		  if (object1_ptr == object_ptr)
		    stored_atom_ptr[nstored] = atom2_ptr;
		  else
		    stored_atom_ptr[nstored] = atom1_ptr;
		  nstored++;
		}
	    }

	  /* Get pointer to next atom in linked-list */
	  bond_ptr = bond_ptr->next_bond_ptr;
	}

      /* If have some possible bonds, then pick one of them to be the
	 pivotal one */
      if (nstored > 0)
	{
	  /* Generate a random number */
	  random_no = get_random_number(&iseed,Random_Start);

	  /* Get the corresponding atom position */
	  iatom = nstored * random_no;
	  if (iatom < 0)
	    iatom = 0;
	  else if (iatom > nstored - 1)
	    iatom = nstored - 1;

	  /* Get the corresponding atom pointer and its coordinates */
	  atom_ptr = stored_atom_ptr[iatom];
	  x = atom_ptr->x;
	  y = atom_ptr->y;
	}
    }

  /* If performing rotation about one of its atoms, transform the object
     to place that atom at the origin */
  else if (swing == FALSE)
    {
      /* Set the pointer to the first of this object's residues */
      residue_ptr = object_ptr->first_residue_ptr;

      /* Get the number of atoms and decide which is to be the pivot
	 about which to perform the swivel */
      natoms = residue_ptr->natoms;
      random_no = get_random_number(&iseed,Random_Start);
      iatom = natoms * random_no + 1.0;
      if (iatom < 1)
	natoms = 1;
      else if (iatom < natoms)
	natoms = iatom;

      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
      iatom = 0;

      /* Loop through all the atoms until get to the one about
	 which the swivel is to be performed */
      while (iatom < natoms && atom_ptr != NULL)
	{
	  /* Store the coordinates of this atom */
	  x = atom_ptr->x;
	  y = atom_ptr->y;
	
	  /* Get pointer to the next atom */
	  atom_ptr = atom_ptr->next;
	  iatom++;
	}
    }

  /* Move object to place the pivotal atom at the origin */
  move_object(object_ptr,x,y,0.0);

  /* Apply the rotation matrix to the current object */
  rotate_object(object_ptr,matrix);

  /* Move object back to where it came from */
  move_object(object_ptr,-x,-y,0.0);

  /* Return the current random number seed */
  *rseed = iseed;
}
/***********************************************************************

flip_attached_groups  -  Flip all the groups attached by H-bonds and
                         non-bonded contacts to the ligand atoms that
                         have just been flipped

***********************************************************************/

void flip_attached_groups(struct object *flipped_object_ptr,
			  float move_x,float move_y,float matrix[4][3])
{
  int i, j;
  int wanted;
  float inverse_matrix[4][3];

  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct object *object_ptr, *object1_ptr, *object2_ptr;

  /* Initialise pointer to start of linked list of bonds */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds in the linked list */
  while (bond_ptr != NULL)
    {
      /* Check if this is an elastic bond */
      if (bond_ptr->elastic == TRUE)
	{
	  /* Get the pointers to the two atoms at either end of the bond */
	  atom1_ptr = bond_ptr->first_atom_ptr;
	  atom2_ptr = bond_ptr->second_atom_ptr;

	  /* Check whether one of the atoms has been flipped and the
	     other not */
	  wanted = TRUE;
	  if ((atom1_ptr->checked == TRUE && atom2_ptr->checked == TRUE) ||
	      (atom1_ptr->checked == FALSE && atom2_ptr->checked == FALSE))
	    wanted = FALSE;

	  /* Check that one atom belongs to the flipped object and
	     the other one doesn't */
	  else
	    {
	      /* Get the objects these two atoms belong to */
	      object1_ptr = atom1_ptr->residue_ptr->object_ptr;
	      object2_ptr = atom2_ptr->residue_ptr->object_ptr;

	      /* Check that one is the flipped and the other isn't */
	      if (object1_ptr == flipped_object_ptr &&
		  object2_ptr != flipped_object_ptr)
		{
		  /* Confirm that the atom on the flipped object is the
		     checked one */
		  if (atom1_ptr->checked == TRUE)
		    {
		      /* Save the pointer to the other object */
		      object_ptr = object2_ptr;
		    }
		  else
		    wanted = FALSE;
		}		
	      else if (object2_ptr == flipped_object_ptr &&
		       object1_ptr != flipped_object_ptr)
		{
		  /* Confirm that the atom on the flipped object is the
		     checked one */
		  if (atom2_ptr->checked == TRUE)
		    {
		      /* Save the pointer to the other object */
		      object_ptr = object1_ptr;
		    }
		  else
		    wanted = FALSE;
		}		

	      /* Otherwise, don't want to flip anything */
	      else
		wanted = FALSE;
	    }

	  /* If the bond is a contact bond, then only want it it the
	     group involved is a hydrophonic group */
	  if (wanted == TRUE && bond_ptr->bond_type == CONTACT &&
	      object_ptr->object_type != HYDROPHOBIC)
	    wanted = FALSE;

	  /* If have a valid bond, perform the flip of this object */
	  if (wanted == TRUE)
	    {
	      /* Transform the current object into the reference
		 frame of the flipped object */
	      move_object(object_ptr,move_x,move_y,0.0);
	      transform_object(object_ptr,matrix);

	      /* Perform the flip */
	      flip_object(object_ptr,1,-1);
	      
	      /* Generate the inverse rotation matrix to that used
		 above and perform the inverse rotation */
	      for (i = 0; i < 4; i++)
		for (j = 0; j < 3; j++)
		  inverse_matrix[i][j] = matrix[i][j];
	      inverse_matrix[0][1] = - matrix[0][1];
	      inverse_matrix[1][0] = - matrix[1][0];
	      transform_object(object_ptr,inverse_matrix);

	      /* Move object back to where it came from */
	      move_object(object_ptr,- move_x,- move_y,0.0);
	    }
	}

      /* Get pointer to next atom in linked-list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/***********************************************************************

make_random_flip  -  Make a random flip of the object about one of its
                     rotatable bonds

***********************************************************************/

void make_random_flip(struct object *object_ptr,int *rseed,int flip_all)
{
  int iseed;
  int flip_bond, ibond;
  int both_ends, clash, done, use_end;

  float coord_store[MAXATOMS][3], matrix[4][3], move_x, move_y, x;
  float internal_energy;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr;
  struct object_bond *object_bond_ptr;

  /* Get the current random number generator seed */
  iseed = *rseed;

  /* Save object's current internal energy */
  internal_energy = object_ptr->internal_energy;

  /* Get a random number */
  x = get_random_number(&iseed,Random_Start);

  /* Convert into the sequential number of the rotatable bond to
     be flipped */
  flip_bond = object_ptr->nrot_bonds * x;

  /* Initialise bond count */
  ibond = 0;
  both_ends = FALSE;
  done = FALSE;

  /* Get pointer to the first of this object's bonds */
  object_bond_ptr = object_ptr->first_object_bond_ptr;
      
  /* Loop through all object's bonds */
  while (object_bond_ptr != NULL && done == FALSE)
    {
      /* Get the bond pointer */
      bond_ptr = object_bond_ptr->bond_ptr;

      /* Check whether this bond is a rotatable one */
      if (bond_ptr->rotatable_bond == TRUE)
	{
	  /* If this is the bond selected for the flip, then perform
	     the flip */
	  if (ibond == flip_bond)
	    {
	      /* Determine which end of this bond to flip */
	      x = get_random_number(&iseed,Random_Start);
	      if (x <= 0.5)
		use_end = FIRST;
	      else
		use_end = SECOND;

	      /* Determine whether to flip whole object */
	      if (flip_all == FALSE)
		{
		  /* Use previously-generated random number to decide
		     whether to flip the whole object */
		  if (fabs(x - 0.5) < 0.05)
		    both_ends = TRUE;
		}

	      /* Get the coordinates of the atom at the required end */
	      if (use_end == FIRST)
		{
		  move_x = bond_ptr->first_atom_ptr->x;
		  move_y = bond_ptr->first_atom_ptr->y;
		  atom_ptr = bond_ptr->second_atom_ptr;
		}
	      else
		{
		  move_x = bond_ptr->second_atom_ptr->x;
		  move_y = bond_ptr->second_atom_ptr->y;
		  atom_ptr = bond_ptr->first_atom_ptr;
		}

	      /* Move the current object to place this atom
		 at the origin */
	      move_object(object_ptr,move_x,move_y,0.0);

	      /* Store the coordinates of the bond two ends */
	      coord_store[0][0] = 0.0;
	      coord_store[0][1] = 0.0;
	      coord_store[0][2] = 0.0;
	      coord_store[1][0] = atom_ptr->x;
	      coord_store[1][1] = atom_ptr->y;
	      coord_store[1][2] = 0.0;

	      /* Calculate the transformation that will put this bond
		 along the x-axis, with atom 1 at the origin and atom 2
		 pointing in the -ve x-direction */
	      line_transformation(coord_store,matrix,0.0);

	      /* Apply the transformation to all atoms in the current
		 object */
	      transform_object(object_ptr,matrix);

	      /* Flip one side of the object around this rotatable bond
		 (currently lying along the x-axis) by changing the signs
		 of all the y-coords */
	      flip_about_bond(bond_ptr,use_end);

	      /* If flipping both ends, then flip the other end too */
	      if (both_ends == TRUE)
		flip_about_bond(bond_ptr,SECOND);

	      /* Calculate energy of the flipped configuration */
 	      internal_energy = calc_internal_energy(object_ptr,&clash);
	      object_ptr->internal_energy = internal_energy;

	      /* If are flipping all the attached groups, as well,
		 then do those */
	      if (flip_all == TRUE)
		{
		  bond_ptr->first_atom_ptr->checked = FALSE;
		  bond_ptr->second_atom_ptr->checked = FALSE;
		  flip_attached_groups(object_ptr,move_x,move_y,matrix);
		}

	      /* Generate the inverse rotation matrix to that used
		 above and perform the inverse rotation */
	      matrix[0][1] = - matrix[0][1];
	      matrix[1][0] = - matrix[1][0];
	      transform_object(object_ptr,matrix);

	      /* Move object back to where it came from */
	      move_object(object_ptr,- move_x,- move_y,0.0);
	      done = TRUE;
	    }

	  /* Increment count of rotatable bonds */
	  ibond++;
	}
      
      /* Get pointer to this object's next bond entry */
      object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
    }

  /* Return the current random number seed */
  *rseed = iseed;
}
/***********************************************************************

make_random_atom_flip  -  Make a random flip of the object about one of
                          its atoms

***********************************************************************/

void make_random_atom_flip(struct object *object_ptr,int *rseed,
			   int flip_about_x)
{
  int flip_atom, iseed;
  int iatom, natoms;

  float x, y;
  float random_no;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Get the current random number generator seed */
  iseed = *rseed;

  /* Set the pointer to the first of this object's residues and get the
     number of atoms */
  residue_ptr = object_ptr->first_residue_ptr;
  natoms = residue_ptr->natoms;

  /* Generate a random number */
  random_no = get_random_number(&iseed,Random_Start);

  /* Convert into the sequential number of the atom to be flipped */
  flip_atom = natoms * random_no + 1.0;
  if (flip_atom < 1)
    natoms = 1;
  else if (flip_atom < natoms)
    natoms = flip_atom;

  /* Initialise atom count */
  iatom = 0;

  /* Get pointer to the first of this residue's atoms */
  atom_ptr = residue_ptr->first_atom_ptr;

  /* Loop through all the atoms until get to the one about
     which the swivel is to be performed */
  while (iatom < natoms && atom_ptr != NULL)
    {
      /* Store the coordinates of this atom */
      x = atom_ptr->x;
      y = atom_ptr->y;
	
      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next;
      iatom++;
    }

  /* Move object to the origin */
  move_object(object_ptr,x,y,0.0);

  /* Apply the flip to the current object */
  if (flip_about_x == TRUE)
    flip_object(object_ptr,-1,1);
  else
    flip_object(object_ptr,1,-1);

  /* Move object back */
  move_object(object_ptr,-x,-y,0.0);

  /* Return the current random number seed */
  *rseed = iseed;
}
/***********************************************************************

get_total_energy  -  Calculate and print total energy for all objects

***********************************************************************/

float get_total_energy(int loop,int print_headings,int print_data)
{
  int i, iobject;

  float object_energy, total_energy;
  float atom_clash_energy, bond_length_energy, bond_clash_energy;
  float bond_overlap_energy, internal_energy;
/* v.4.0--> */
  float boundary_energy, rel_pstn_energy;
/* <--v.4.0 */
  float energy_sum[NTERMS], energy_term[NTERMS];

  struct object *object_ptr;

  /* If headings to be printed, then do so */
  if (print_headings == TRUE)
    {
      printf("\n");
/* v.4.0--> */
/*      printf("                      Atom     Bond     Bond     Bond\n");
      printf(" Loop     Internal    Clash   Length    Clash  overlap");
      printf("       Total\n"); */
      printf("                 Atom    Bond    Bond    Bond ");
      printf("  Inter-   Resid           \n");
      printf("Loop  Internal   Clash  Length   Clash overlap");
      printf("   face    pstns      Total\n");
/* <--v.4.0 */
    }

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;
  for (i = 0; i < NTERMS; i++)
    energy_sum[i] = 0.0;
  total_energy = 0.0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Get the object's internal energy */
      internal_energy = object_ptr->internal_energy;

      /* Calculate the energy between this object and all
	 the others */
      object_energy = calculate_energy(object_ptr,
				       &atom_clash_energy,
				       &bond_length_energy,
				       &bond_clash_energy,
/* v.4.0--> */
/*				       &bond_overlap_energy); */
				       &bond_overlap_energy,
				       &boundary_energy,
				       &rel_pstn_energy);
/* <--v.4.0 */

      /* Store the individual energy terms */
      energy_term[0] = internal_energy;
      energy_term[1] = atom_clash_energy / 2.0;
      energy_term[2] = bond_length_energy / 2.0;
      energy_term[3] = bond_clash_energy / 2.0;
      energy_term[4] = bond_overlap_energy / 2.0;
/* v.4.0--> */
      energy_term[5] = boundary_energy;
      energy_term[6] = rel_pstn_energy / 2.0;
/* <--v.4.0 */

      /* Increment energy sums */
      for (i = 0; i < NTERMS; i++)
	{
	  energy_sum[i] = energy_sum[i] + energy_term[i];
	  total_energy = total_energy + energy_term[i];
	}

      /* Get pointer for next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

  /* Print out summary of the total energy and its components */
  if (print_data == TRUE)
    {
/* v.4.0--> */
/*      printf("%5d.   ",loop); */
      printf("%4d. ",loop);
/* <--v.4.0 */
      for (i = 0; i < NTERMS; i++)
/* v.4.0--> */
/*	printf(" %8.2f",energy_sum[i]);
      printf("  %10.2f\n",total_energy); */
	printf("%8.2f",energy_sum[i]);
      printf(" %10.2f\n",total_energy);
/* <--v.4.0 */
    }

  /* Return the total energy */
  return(total_energy);
}
/***********************************************************************

energy_minimize  -  Minimize overall energy of the entire structure to
                    make numbers of atom and bond clashes a minimum

***********************************************************************/

void energy_minimize(int *iseed,int nligands)
{
  int loop, nlines, object_type;
  int accept_move, flip_about_x, flip_all, got_move, move_type, swing;
  int rseed;
  int energy_limit_reached;
  int n_accept, naccept_move, ntrial_move, total_accept;
  int naccept_swivel, ntrial_swivel;
  int naccept_swing, ntrial_swing;
  int naccept_flip, ntrial_flip;
  int nloops, nmove;
  int print_data, print_headings;

  float atom_clash_energy, bond_length_energy, bond_clash_energy;
  float bond_overlap_energy, internal_energy;
  float energy, previous_energy;
  float energy_drop, last_total_energy;
  float move_size, x;
  float max_swivel_angle;
  float max_swing_angle;
/* v.4.0--> */
  float boundary_energy;
  float rel_pstn_energy;
/* <--v.4.0 */

  struct object *move_object_ptr;

  /* Initialise variables */
  rseed = *iseed;
  total_accept = 0;
  nlines = 0;
  last_total_energy = 0.0;
  print_headings = TRUE;
  energy_limit_reached = FALSE;

  /* Update all the residue boundaries */
  update_all_boundaries();

  /* Initialise minimization variables */
  loop = 0;
  max_swivel_angle = PI;
  max_swing_angle = PI;
  move_size = Max_HBdist;
  move_type = MOVE_TYPES;
  n_accept = 0;
  nloops = Max_Loops;
  printf("\n");
  printf("\n");
  printf("Energy minimization ... %d loops",nloops);
  if (Energy_Drop_Maximum > 0.0)
    printf(", terminating if energy-drop < %f\n",Energy_Drop_Maximum);
  else
    printf("\n");

  /* Loop until energy reaches an acceptable level, or minimization
     has gone on for the required number of loops */
  while (loop < nloops && energy_limit_reached == FALSE)
    {
      /* Initialise counts */
      naccept_move = 0;
      ntrial_move = 0;
      naccept_swivel = 0;
      ntrial_swivel = 0;
      naccept_swing = 0;
      ntrial_swing = 0;
      naccept_flip = 0;
      ntrial_flip = 0;
      if (loop % 4 == 0)
	n_accept = 0;
      nmove = 0;

      /* Determine the move_type */
      move_type++;
      if (move_type >= MOVE_TYPES)
	move_type = 0;

      /* Initialise pointer to the first stored object */
      move_object_ptr = first_object_ptr;

      /* Loop through all objects, moving each one in turn */
      while (move_object_ptr != NULL)
	{
	  /* Initialise flags */
	  accept_move = FALSE;
	  flip_all = FALSE;
	  object_type = move_object_ptr->object_type;

	  /* Make sure the move type is compatible with the object
	     type */
	  got_move = FALSE;
	  if (object_type == LIGAND)
	    {
	      if (nligands == 1 && move_type == FLIP_OBJECT)
		got_move = TRUE;
	      else if (nligands > 1)
		got_move = TRUE;
	    }
	  else if (object_type == HYDROPHOBIC || object_type == WATER ||
		   object_type == SIMPLE_HGROUP)
	    {
	      if (move_type == MOVE_OBJECT)
		got_move = TRUE;
	    }
	  else if (object_type == HGROUP)
	    got_move = TRUE;

	  /* If flipping about a ligand bond, determine whether to take
	     all connected group when performing the flip */
	  if (object_type == LIGAND && move_type == FLIP_OBJECT)
	    {
	      /* Generate a random number and use to determine whether
		 to take all connected objects with the flipped part
		 of the ligand */
	      x = get_random_number(&rseed,Random_Start);
	      if (x < 0.5)
		flip_all = TRUE;
	    }

	  /* If flipping about a nonligand bond, decide whether to
	     change flip to be about an atom instead */
	  if (object_type != LIGAND && move_type == FLIP_OBJECT)
	    {
	      /* Generate a random number and determine whether to
		 flip about an atom instead */
	      x = get_random_number(&rseed,Random_Start);
	      if (x < 0.2)
		{
		  move_type = FLIP_ABOUT_ATOM;
		  flip_about_x = TRUE;
		  if (x < 0.1)
		    flip_about_x = FALSE;
		}
	    }

	  /* Switch off swings if near the end */
	  if (loop > (nloops - 1) && move_type == SWING_OBJECT)
	    got_move = FALSE;

	  /* If have a valid move for this object, then do it */
	  if (got_move == TRUE)
	    {
	      /* Increment count of valid moves on this loop */
	      nmove++;

	      /* Get the object's internal energy */
	      internal_energy = move_object_ptr->internal_energy;

	      /* If move might affect more than one object, compute
		 the total energy */
	      if (flip_all == TRUE)
		{
		  print_headings = FALSE;
		  print_data = FALSE;
		  previous_energy = get_total_energy(loop,print_headings,
						     print_data);

		  /* Save the coordinates of all objects */
		  save_restore_all_objects(SAVE);
		}

	      /* Otherwise, just calculate the energy between this
		 object and all the others */
	      else
		{
		  /* Calculate this object's energy */
		  previous_energy = calculate_energy(move_object_ptr,
						     &atom_clash_energy,
						     &bond_length_energy,
						     &bond_clash_energy,
/* v.4.0--> */
/*						     &bond_overlap_energy); */
						     &bond_overlap_energy,
						     &boundary_energy,
						     &rel_pstn_energy);
/* <--v.4.0 */

		  /* Save its coordinates */
		  save_restore_coords(move_object_ptr,SAVE);
		}		  

	      /* If moving object, make a random move to the current
		 object */
	      if (move_type == MOVE_OBJECT)
		{
		  /* Make the move */
		  make_random_move(move_object_ptr,move_size,&rseed);

		  /* Increment count of trial moves */
		  ntrial_move++;
		}

	      /* Or, make a random rotation or swing of the
		 current object */
	      else if (move_type == SWIVEL_OBJECT ||
		       move_type == SWING_OBJECT)
		{
		  /* Set flag about whether this is a swing or
		     a rotation */
		  swing = FALSE;
		  if (move_type == SWING_OBJECT)
		    swing = TRUE;

		  /* Apply a trial rotation/swing */
		  make_random_rotation(move_object_ptr,max_swivel_angle,
				       &rseed,swing);

		  /* Increment count of trial rotations */
		  if (move_type == SWIVEL_OBJECT)
		    ntrial_swivel++;
		  else
		    ntrial_swing++;
		}

	      /* Or, make a random flip about one of the object's
		 rotatable bonds */
	      else if (move_type == FLIP_OBJECT)
		{
		  /* Make the flip */
		  make_random_flip(move_object_ptr,&rseed,flip_all);

		  /* Increment count of trial flips */
		  ntrial_flip++;
		}

	      /* Or, make a random flip about one of the object's
		 rotatable bonds */
	      else if (move_type == FLIP_ABOUT_ATOM)
		{
		  /* Make the flip */
		  make_random_atom_flip(move_object_ptr,&rseed,
					flip_about_x);

		  /* Increment count of trial flips */
		  ntrial_flip++;
		}

	      /* If move affected more than one object, compute the new
		 total energy */
	      if (flip_all == TRUE)
		{
		  /* Update all residue- and object-boundaries */
		  update_all_boundaries();

		  /* Get total energy */
		  print_headings = FALSE;
		  print_data = FALSE;
		  energy = get_total_energy(loop,print_headings,print_data);
		}

	      /* Otherwise, calculate new energy between object and all
		 the others */
	      else
		{
		  /* Update residue- and object-boundaries */
		  update_residue_boundaries(move_object_ptr);

		  /* Calculate energy */
		  energy = calculate_energy(move_object_ptr,
					    &atom_clash_energy,
					    &bond_length_energy,
					    &bond_clash_energy,
/* v.4.0--> */
/*					    &bond_overlap_energy); */
					    &bond_overlap_energy,
					    &boundary_energy,
					    &rel_pstn_energy);
/* <--v.4.0 */
		}

	      /* Determine whether to accept or reject this move */
	      if (energy < previous_energy)
		accept_move = TRUE;

	      /* Update total energy, depending on whether move has
		 been accepted or not */
	      if (accept_move == TRUE)
		{
		  /* Increment count of accepted moves */
		  if (move_type == MOVE_OBJECT)
		    naccept_move++;
		  else if (move_type == SWIVEL_OBJECT)
		    naccept_swivel++;
		  else if (move_type == SWING_OBJECT)
		    naccept_swing++;
		  else if (move_type == FLIP_OBJECT ||
			   move_type == FLIP_ABOUT_ATOM)
		    naccept_flip++;
		}

	      /* Otherwise, if no move, revert to previous conformation */
	      else
		{
		  /* If several objects might have been moved, restore
		     the coordinates of all the objects */
		  if (flip_all == TRUE)
		    {
		      /* Restore the coords */
		      save_restore_all_objects(RESTORE);

		      /* Update all residue- and object-boundaries */
		      update_all_boundaries();
		    }
		  
		  /* Otherwise, restore the moved object only */
		  else
		    {
		      /* Restore the coords */
		      save_restore_coords(move_object_ptr,RESTORE);

		      /* Update residue- and object-boundaries */
		      update_residue_boundaries(move_object_ptr);
		    }

		  /* Update object's internal energy */
		  move_object_ptr->internal_energy = internal_energy;
		}
	    }

	  /* Get pointer to next object */
	  move_object_ptr = move_object_ptr->next_object_ptr;
	}

      /* Increment loop count */
      loop++;

      /* Count numbers of trial moves accepted */
      n_accept = naccept_move + naccept_swivel + naccept_swing
	+ naccept_flip;
      total_accept = total_accept + n_accept;
	  
      /* Show progress of minimization */
      if (loop % 20 == 0 || loop == nloops)
	{
	  if (nlines % 20 == 0)
	    print_headings = TRUE;
	  else
	    print_headings = FALSE;
	  print_data = TRUE;
	  energy = get_total_energy(loop,print_headings,print_data);

	  /* If this is the first loop, save the current total energy */
	  if (nlines == 0)
	    last_total_energy = energy;

	  /* Otherwise, check whether the energy drop since last time is
	     smaller than the stop-limit */
	  else
	    {
	      /* Get the drop in energy since the last loop */
	      energy_drop = last_total_energy - energy;

	      /* If this drop is lower than the user-defined limit, then
		 set flag to stop minimization process */
	      if (energy_drop > 0.0 && energy_drop < Energy_Drop_Maximum)
		energy_limit_reached = TRUE;
	    }

	  /* Save the current total energy */
	  last_total_energy = energy;
	  nlines++;
	}

      /* Adjust move size, if necessary */
      if (ntrial_move > 0)
	{
	  /* accept_percen = 100.0 * (float) naccept_move
	    / (float) ntrial_move;
	  if (accept_percen < 50.0)
	    move_size = 0.9 * move_size;  */

	    /* Check that move distance doesn't fall below minimum value */
	    if (move_size < MIN_MOVE)
	      move_size = MIN_MOVE;
	}

      /* Adjust rotation angle, if necessary */
      if (ntrial_swivel > 0)
	{
	  /* accept_percen
	    = 100.0 * (float) naccept_swivel / (float) ntrial_swivel;
	  if (accept_percen < 50.0)
	    max_swivel_angle = 0.9 * max_swivel_angle;  */
	    
	    /* Check that angle doesn't fall below minimum value */
	    if (max_swivel_angle * RADDEG < MIN_ANGLE)
	      max_swivel_angle = MIN_ANGLE / RADDEG;
	}

      /* Adjust swing angle, if necessary */
      if (ntrial_swing > 0)
	{
	  /* accept_percen
	    = 100.0 * (float) naccept_swing / (float) ntrial_swing;
	  if (accept_percen < 50.0)
	    max_swing_angle = 0.9 * max_swing_angle;  */

	    /* Check that angle doesn't fall below minimum value */
	    if (max_swing_angle * RADDEG < MIN_ANGLE)
	      max_swing_angle = MIN_ANGLE / RADDEG;
	}
    }

  /* Retrieve the current random-number seed */
  *iseed = rseed;
}
/***********************************************************************

rotate_whole_structure  -  Apply the given rotation to all the atoms

***********************************************************************/

void rotate_whole_structure(float v[3][3])
{
  int icoord;
  float x, y, z;
  float new_coords[3];

  struct coordinate *atom_ptr;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms to initialise the flag indicating whether
     coords have already been stored */
  while (atom_ptr != NULL)
    {
      /* Extract this atom's coordinates */
      x = atom_ptr->x;
      y = atom_ptr->y;
      z = atom_ptr->z;
	
      /* Apply the transformation */
      for (icoord = 0; icoord < 3; icoord++)
	{
	  /* Calculate transformed coordinates */
	  new_coords[icoord]
	    =  ((x * v[0][icoord])
		+ (y * v[1][icoord])
		+ (z * v[2][icoord]));
	}

      /* Save the new coordinates */
      atom_ptr->x = new_coords[0];
      atom_ptr->y = new_coords[1];
      atom_ptr->z = new_coords[2];

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next;
    }
}
/***********************************************************************

swap_x_and_y  -  Swap the x- and y- coordinates

***********************************************************************/

void swap_x_and_y(void)
{
  float swap;

  struct coordinate *atom_ptr;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms */
  while (atom_ptr != NULL)
    {
      /* Swap the x- and y- coordinates */
      swap = atom_ptr->x;
      atom_ptr->x = atom_ptr->y;
      atom_ptr->y = swap;

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next;
    }
}
/***********************************************************************

flip_about_x  -  Flip the whole structure about the x-axis

***********************************************************************/

void flip_about_x(void)
{
  struct coordinate *atom_ptr;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms */
  while (atom_ptr != NULL)
    {
      /* Perform the flip */
      atom_ptr->y = - atom_ptr->y;

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next;
    }
}
/***********************************************************************

flip_about_y  -  Flip the whole structure about the y-axis

***********************************************************************/

void flip_about_y(void)
{
  struct coordinate *atom_ptr;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms */
  while (atom_ptr != NULL)
    {
      /* Perform the flip */
      atom_ptr->x = - atom_ptr->x;

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next;
    }
}
/***********************************************************************

check_ligand_direction  -  Test if plot needs to be flipped to maintain
	                   residue order down the page, or left-to-right
	                   according to the orientation

***********************************************************************/

void check_ligand_direction(void)
{
  int iatom, iobject, iresid, natoms, ngot, nresid;
/* v.3.1--> */
  int last_resno, resno;

  char resno_string[5];
/* <--v.3.1 */

  float C_x, C_y, N_x, N_y;

  struct coordinate *atom_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;

  /* Initialise maximum and minimum coordinates */

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;
  ngot = 0;
/* v.3.1--> */
  last_resno = -999;
/* <--v.3.1 */

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Process only if this is the ligand */
      if (object_ptr->object_type == LIGAND)
	{
	  /* Set the pointer to the first of this object's residues */
	  residue_ptr = object_ptr->first_residue_ptr;
	  nresid = object_ptr->nresidues;
	  iresid = 0;

	  /* Loop over all this object's residues */
	  while (iresid < nresid && residue_ptr != NULL)
	    {
	      /* Get pointer to this residue's first atom */
	      atom_ptr = residue_ptr->first_atom_ptr;
	      iatom = 0;
	      natoms = residue_ptr->natoms;

	      /* Loop over all the residue's atoms, checking whether
		 they are involved in interactions */
	      while (iatom < natoms && atom_ptr != NULL)
		{
		  /* If this is the very first atom, store its
		     coordinates as the N-terminus coords */
		  if (ngot == 0)
		    {
		      N_x = atom_ptr->x;
		      N_y = atom_ptr->y;
/* v.3.1--> */
		      /* Get the residue number */
		      strncpy(resno_string,atom_ptr->residue_ptr->res_num,4);
		      resno_string[4] = '\0';
		      resno = atoi(resno_string);
/* <--v.3.1 */
		    }

		  /* Otherwise, store as the C-terminal coords */
		  else
		    {
/* v.3.1--> */
		      if (resno > last_resno)
			{
/* <--v.3.1 */
			  C_x = atom_ptr->x;
			  C_y = atom_ptr->y;
/* v.3.1--> */
			}
/* <--v.3.1 */
		    }

		  /* Increment number of atoms processed */
		  ngot++;
/* v.3.1--> */
		  if (resno > last_resno)
		    last_resno = resno;
/* <--v.3.1 */

		  /* Get pointer to the next atom */
		  atom_ptr = atom_ptr->next;
		  iatom++;
		}

	      /* Get pointer to the next residue */
	      residue_ptr = residue_ptr->next_residue_ptr;
	      iresid++;
	    }
	}

      /* Get pointer for next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

  /* Determine whether whole picture needs to be flipped */

  /* If portrait plot needs to be flipped, the flip about y */
  if (Picture->Portrait == TRUE && N_y < C_y)
    flip_about_x();

  /* If landscape plot needs to be flipped, the flip about x */
  else if (Picture->Portrait == FALSE && N_x > C_x)
    flip_about_y();

}
/***********************************************************************

ligand_orientate  -  Orientate the ligand on  the y-axis on the output
                     page

***********************************************************************/

void ligand_orientate(void)
{
  int iatom, iobject, iresid, natoms, nresid, nstored;
  int store_ca, store_mainchain, wanted;

  float transformation[3][3];
  float mass_x, mass_y, mass_z, coord_store[MAXATOMS][3];

  struct coordinate *atom_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;


  /* Initialise variables */
  nstored = 0;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Process only if this is the ligand */
      if (object_ptr->object_type == LIGAND)
	{
	  /* Determine which atoms are to be stored */
	  store_ca = FALSE;
	  store_mainchain = FALSE;

	  /* For the simplified plot, then use CA, providing have at
	     least 2 of them */
	  if (Include->Simple_Ligand_Residues && object_ptr->n_ca > 1)
	    store_ca = TRUE;

	  /* Otherwise, use just the mainchain atoms, if have a
	     sufficient number */
	  else if (object_ptr->nmain_chain > 9)
	    store_mainchain = TRUE;

	  /* Get pointer to the first of this object's residues */
	  residue_ptr = object_ptr->first_residue_ptr;
	  nresid = object_ptr->nresidues;
	  iresid = 0;

	  /* Loop over all this object's residues */
	  while (iresid < nresid && residue_ptr != NULL)
	    {
	      /* Get pointer to this residue's first atom */
	      atom_ptr = residue_ptr->first_atom_ptr;
	      
	      /* Initialise atom count */
	      iatom = 0;
	      natoms = residue_ptr->natoms;

	      /* Loop over all the residue's atoms, checking whether
		 they are involved in interactions */
	      while (iatom < natoms && atom_ptr != NULL)
		{
		  /* Initialise flag */
		  wanted = FALSE;

		  /* Determine if atom is of the type required */
		  if (store_ca == TRUE)
		    {
		      /* If this is a CA, then want it */
		      if (!strncmp(atom_ptr->atom_type," CA ",4))
			wanted = TRUE;
		    }
		  else if (store_mainchain == TRUE)
		    {
		      /* If this is a mainchain, then want it */
		      if (atom_ptr->side_chain == FALSE)
			wanted = TRUE;
		    }
		  else
		    wanted = TRUE;

		  /* If atom is wanted, then store its coordinates */
		  if (wanted == TRUE)
		    {
		      /* Store the coords */
		      coord_store[nstored][0] = atom_ptr->x;
		      coord_store[nstored][1] = atom_ptr->y;
		      coord_store[nstored][2] = atom_ptr->z;
		      nstored++;
		    }

		  /* Get pointer to the next atom */
		  atom_ptr = atom_ptr->next;
		  iatom++;
		}

	      /* Get pointer to the next residue */
	      residue_ptr = residue_ptr->next_residue_ptr;
	      iresid++;
	    }
	}

      /* Get pointer for next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

  /* If have stored any coordinates, perform a principal components
     analysis */
  if (nstored > 1)
    {
      /* If plotting the simplified representation, then all residues
	 will be in a straight line, so just use the first and last
	 positions */
      if (store_ca == TRUE)
	{
	  /* Move coords of last CA position to the second position */
	  coord_store[1][0] = coord_store[nstored - 1][0];
	  coord_store[1][1] = coord_store[nstored - 1][1];
	  coord_store[1][2] = coord_store[nstored - 1][2];

	  /* Set the number of points to 2 */
	  nstored = 2;
	}

      /* Calculate the centre of mass of these atoms */
      centre_of_mass(coord_store,nstored,&mass_x,&mass_y,&mass_z);

      /* Place the centre of mass at the origin */
      adjust_to_origin(mass_x,mass_y,mass_z);

      /* Apply the same translation to the coordinates stored in
	 coord_store */
      adjust_stored_coords(nstored,coord_store,mass_x,mass_y,mass_z);

      /* Calculate transformation to orient the principal axis along
	 the x-direction */
      principal_components(nstored,transformation,coord_store);

      /* Apply this transformation to the entire structure */
      rotate_whole_structure(transformation);

      /* If the plot is to be in Portrait, rather than Landscape,
	 orientation swap the x coords with the y coords */
      if (Picture->Portrait == TRUE)
	  swap_x_and_y();

      /* Test if plot needs to be flipped to maintain residue order
	 down the page, or left-to-right according to the orientation */
      check_ligand_direction();

      /* Update all the residue- and object-boundaries */
      update_all_boundaries();
    }
}
/***********************************************************************

write_pdb_file  -  Print out the co-ordinates in a PDB fashion

***********************************************************************/

void write_pdb_file(void)
{
/* v.4.0--> */
  int interface;
/* <--v.4.0 */

  struct coordinate *atom_ptr;

  /* Write out the title (if there is one) */
  printf("\nWriting out the PDB file...\n");
  if (Print_Title[0] != '\0')
    fprintf(ligplot_pdb,"TITLE  %s\n",Print_Title);

/* v.4.0--> */
  /* If this is a DIMPLOT, then write out identifier */
  if (Interface_Plot == TRUE)
    fprintf(ligplot_pdb,"DIMPLOT\n");

  /* Initialise interface identifier */
  interface = 0;
/* <--v.4.0 */

  /* Initialise pointer to start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms in the linked list */
  while (atom_ptr != NULL)
    {
/* v.4.0--> */
      /* If this is a DIMPLOT, then check which interface is being
	 written out */
      if (Interface_Plot == TRUE)
	{
	  if (atom_ptr->residue_ptr->object_ptr->interface != interface)
	    {
	      fprintf(ligplot_pdb,"SURFACE %d\n",
		      atom_ptr->residue_ptr->object_ptr->interface);
	      interface = atom_ptr->residue_ptr->object_ptr->interface;
	    }
	}
/* <--v.4.0 */

      /* Write out the atom record */
      fprintf(ligplot_pdb,"ATOM  %5s ",atom_ptr->atom_number);
      fprintf(ligplot_pdb,"%4s ",atom_ptr->atom_type);
      fprintf(ligplot_pdb,"%3s %c%5s   %8.3f%8.3f%8.3f%6s%6s%8.3f\n",
	      atom_ptr->residue_ptr->res_name,
	      atom_ptr->residue_ptr->chain,
	      atom_ptr->residue_ptr->res_num,
	      atom_ptr->x,
	      atom_ptr->y,
	      atom_ptr->z,
	      atom_ptr->occupancy,
	      atom_ptr->bvalue,
	      atom_ptr->accessibility);

      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next;
    }

/* v.4.0--> */
  /* Write out the maximum accessibility/B-value */
  fprintf(ligplot_pdb,"ENDLIN                                        ");
  fprintf(ligplot_pdb,"                    %8.3f\n",Maximum_accessibility);
 /* <--v.4.0 */
}
/***********************************************************************

write_bonds  -  Write out all the bonds in the structure

***********************************************************************/

void write_bonds(void)
{
  int ibond;
  char bond_type_char, ebond_char, rbond_char, single_double, source[7];
  struct bond *bond_ptr;
  struct coordinate *atom_ptr1, *atom_ptr2;

  printf("\nWriting out the bonds file...\n");

  /* Initialise pointer to start of linked list of bonds */
  bond_ptr = first_bond_ptr;
  ibond = 0;

  /* Loop through all the bonds in the linked list */
  while (bond_ptr != NULL)
    {
      /* Write out provided bond hasn't been marked for deletion */
      if (bond_ptr->bond_type != DELETED)
	{
	  bond_type_char = '.';
	  if (bond_ptr->bond_type == HBOND)
	    bond_type_char = 'H';
	  else if (bond_ptr->bond_type == COVALENT)
	    bond_type_char = 'c';
	  else if (bond_ptr->bond_type == CONTACT)
	    bond_type_char = 'n';
	  else if (bond_ptr->bond_type == INTERNAL)
	    bond_type_char = 'i';

	  /* Determine whether a rotatable bond */
	  rbond_char = ' ';
	  if (bond_ptr->rotatable_bond == TRUE)
	    rbond_char = 'r';
	  else if (bond_ptr->end_bond == TRUE)
	    rbond_char = 'e';

	  /* Determine whether an elastic bond */
	  ebond_char = ' ';
	  if (bond_ptr->elastic == TRUE)
	    ebond_char = 'E';

	  /* Get the bond's two atoms */
	  atom_ptr1 = bond_ptr->first_atom_ptr;
	  atom_ptr2 = bond_ptr->second_atom_ptr;

	  /* Determine bond-order, single/double/treble */
	  single_double = '-';
	  if (bond_ptr->bond_order == 2)
	    single_double = 'd';
	  else if (bond_ptr->bond_order == 3)
	    single_double = 't';

	  /* Determine where bond info came from */
	  if (bond_ptr->bond_source == CONECT)
	    strcpy(source,"CONECT");
	  else if (bond_ptr->bond_source == CALCULATED)
	    strcpy(source,"Calc  ");
	  else if (bond_ptr->bond_source == HBPLUS)
	    strcpy(source,"HBPLUS");
	  else
	    strcpy(source,"      ");
	  
	  /* Write the bond information out to the ligplot.bonds file */
	  fprintf(ligplot_bonds_out,"Bond %3d. Atoms%s [%s %s %s %c]",ibond,
		  atom_ptr1->atom_number,atom_ptr1->atom_type,
		  atom_ptr1->residue_ptr->res_name,
		  atom_ptr1->residue_ptr->res_num,
		  atom_ptr1->residue_ptr->chain);
	  fprintf(ligplot_bonds_out," ->%s [%s %s %s %c]",
		  atom_ptr2->atom_number,atom_ptr2->atom_type,
		  atom_ptr2->residue_ptr->res_name,
		  atom_ptr2->residue_ptr->res_num,
		  atom_ptr2->residue_ptr->chain);
	  fprintf(ligplot_bonds_out,"%c%c%c%7.4f %c %s\n",bond_type_char,
		  rbond_char,ebond_char,bond_ptr->bond_length,
		  single_double,source);

	  /* If this is a hydrogen bond, write out to the ligplot.hhb
	     file */
	  if (bond_ptr->bond_type == HBOND ||
	      bond_ptr->bond_type == INTERNAL)
	    {
	      fprintf(ligplot_hhb_out,"%s %c %s %s     ",
		      atom_ptr1->residue_ptr->res_name,
		      atom_ptr1->residue_ptr->chain,
		      atom_ptr1->residue_ptr->res_num,
		      atom_ptr1->atom_type);
	      fprintf(ligplot_hhb_out,"%s %c %s %s    %4.2f\n",
		      atom_ptr2->residue_ptr->res_name,
		      atom_ptr2->residue_ptr->chain,
		      atom_ptr2->residue_ptr->res_num,
		      atom_ptr2->atom_type,bond_ptr->bond_length);
	    }

	  /* If this is a non-bonded contact, write out to the ligplot.nnb
	     file */
	  else if (bond_ptr->bond_type == CONTACT)
	    {
	      fprintf(ligplot_nnb_out,"%s %c %s %s     ",
		      atom_ptr1->residue_ptr->res_name,
		      atom_ptr1->residue_ptr->chain,
		      atom_ptr1->residue_ptr->res_num,
		      atom_ptr1->atom_type);
	      fprintf(ligplot_nnb_out,"%s %c %s %s    %4.2f\n",
		      atom_ptr2->residue_ptr->res_name,
		      atom_ptr2->residue_ptr->chain,
		      atom_ptr2->residue_ptr->res_num,
		      atom_ptr2->atom_type,bond_ptr->bond_length);
	    }

	  /* Increment bond counter */
	  ibond++;
	}

      /* Get pointer to next atom in linked-list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/* v.4.0--> */
/***********************************************************************

write_rcm_file  -  Write out the co-ordinates of the residue centres
                   of mass to the ligplot.rcm file

***********************************************************************/

void write_rcm_file(void)
{
  int iresid, nresid;

  char chain;

  struct object *object_ptr;
  struct residue *residue_ptr;

  /* Write out the header records */
  fprintf(ligplot_rcm,"                 Res.  Res       Flattened    \n");
  fprintf(ligplot_rcm,"            Atom Name  Num      coordinates   \n");
  fprintf(ligplot_rcm,"            ---- --- -----    ----------------\n");

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* If chain is blank, write out as a hyphen */
      chain = object_ptr->first_residue_ptr->chain;
      if (chain == ' ')
	chain = '-';

      /* Set the pointer to the first of this object's residues */
      residue_ptr = object_ptr->first_residue_ptr;
      nresid = object_ptr->nresidues;
      iresid = 0;

      /* Loop over all this object's residues */
      while (iresid < nresid && residue_ptr != NULL)
	{
	  /* Write out the residue's centre-of-mass coords */
	  fprintf(ligplot_rcm,"RESDUE      CofM %s %c%s   %8.3f%8.3f\n",
		  residue_ptr->res_name,
		  chain,
		  residue_ptr->res_num,
		  residue_ptr->flattened_mean_x,
		  residue_ptr->flattened_mean_y);

	  /* Get pointer to the next residue */
	  residue_ptr = residue_ptr->next_residue_ptr;
	  iresid++;
	}

      /* Get pointer for next object */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* Close the output file */
  fclose(ligplot_rcm);
}
/* <--v.4.0 */


/**************   P O S T S C R I P T   R O U T I N E S   *************/

/***********************************************************************

psarc_  -  Write arc out to PostScript file

***********************************************************************/

void psarc_(float x1,float y1,float radius,float angle_start,float angle_end)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f %7.2f %7.2f Arc\n",x1, y1, radius,
          angle_start, angle_end);
}
/***********************************************************************

psbbox_ - Function to write out bounded box to Postscript file

**********************************************************************/

void psbbox_(float x1,float y1,float x2,float y2,float x3,float y3,
	     float x4,float y4)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f Pline4\n",
	  x1,y1,x2,y2,x3,y3,x4,y4);
}
/***********************************************************************

psccol_  -  Write circle colour out to PostScript file

***********************************************************************/

void psccol_(char colour[PS_STRING_LENGTH])
{
  fprintf(ps_file,"/Circol { %s } def\n",colour);
}
/***********************************************************************

pscirc_  -  Write a circle out to PostScript file

***********************************************************************/

void pscirc_(float x,float y,char radius[PS_STRING_LENGTH])
{
  fprintf(ps_file,"%7.2f %7.2f %s Circle\n",x,y,radius);
}
/***********************************************************************

psclos_   Write final lines to PostScript file and close

***********************************************************************/

void psclos_(void)
{

    /* Write out closing lines to PostScript file*/
    fprintf(ps_file,"\n\n");
    fprintf(ps_file,"LigplotSave restore\n");
    fprintf(ps_file,"showpage\n");
    fprintf(ps_file,"%%%%Trailer\n");
    fprintf(ps_file,"%%%%BoundingBox: %d %d %d %d\n",(int) BBOXX1 - 1,
        (int) BBOXY1 - 1, (int) BBOXX2 + 1, (int) BBOXY2 + 1);
    fprintf(ps_file,"%%%%EOF");
}
/***********************************************************************

pscolb_  -  Write background colour level (for text, etc)

***********************************************************************/

void pscolb_(char text_string[PS_STRING_LENGTH])
{
    fprintf(ps_file,"%s\n",text_string);
}
/***********************************************************************

pscomm_   Write out a comment to the PostScript file

***********************************************************************/

void pscomm_(char text_string[PS_STRING_LENGTH])
{
  fprintf(ps_file,"\n");
  fprintf(ps_file,"%% %s\n",text_string);
  fprintf(ps_file,"\n");
}
/***********************************************************************

pscrgb_  -  Write RGB circle colour out to PostScript file

***********************************************************************/

void pscrgb_(float red, float green, float blue)
{
  fprintf(ps_file,"/Circol { %7.4f %7.4f %7.4f setrgbcolor } def\n",
          red, green, blue);
}
/***********************************************************************

psctxt_   Write out centred text to PostScript file

***********************************************************************/

void psctxt_(float x,float y,char size[PS_STRING_LENGTH],
    char text_string[PS_STRING_LENGTH])
{
  fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
  fprintf(ps_file,"(%s) %s Center\n",text_string,size);
  fprintf(ps_file,"(%s) %s Print\n",text_string,size);
}
/***********************************************************************

psdash_  -  Switch dashed lines on/off

***********************************************************************/

void psdash_(int ONOFF)
{
  if (ONOFF == 0)
    fprintf(ps_file,"[]  0 setdash\n");
  else
    fprintf(ps_file,"[ %2d %2d ] 0 setdash\n",ONOFF,ONOFF);
}
/***********************************************************************

pshade_  -  Write colour/shade out to PostScript file

***********************************************************************/

void pshade_(char colour[PS_STRING_LENGTH])
{
    fprintf(ps_file,"G %s\n",colour);
}
/***********************************************************************

psline_  -   Write line out to PostScript file

***********************************************************************/

void psline_(float x1,float y1,float x2,float y2)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f %7.2f L\n",x1, y1, x2, y2);
}
/***********************************************************************

pslwid_  -  Write line-width out to PostScript file

***********************************************************************/

void pslwid_(char width[PS_STRING_LENGTH])
{
    fprintf(ps_file,"%s W\n",width);
}
/***********************************************************************

psphcl_  -  Write sphere colour out to PostScript file

***********************************************************************/

void psphcl_(char colour[PS_STRING_LENGTH])
{
    fprintf(ps_file,"/Sphcol { %s } def\n",colour);
}
/***********************************************************************

pspher_  -  Write a sphere out to PostScript file

***********************************************************************/

void pspher_(float x,float y,char radius[PS_STRING_LENGTH])
{
  fprintf(ps_file,"%7.2f %7.2f %s Sphere\n",x,y,radius);
}
/***********************************************************************

psrest_  -  Restore graphics state

***********************************************************************/

void psrest_(void)
{
    fprintf(ps_file,"grestore\n");
}
/***********************************************************************

psrot_  -  Rotate graphics page through 90 degrees

***********************************************************************/

void psrot_(float new_origin_x, float new_origin_y)
{
    fprintf(ps_file," %7.2f %7.2f moveto Rot90\n",new_origin_x,new_origin_y);
}
/***********************************************************************

pssave_  -  Save graphics state

***********************************************************************/

void pssave_(void)
{
    fprintf(ps_file,"gsave\n");
}
/***********************************************************************

pstext_   Write out text to PostScript file

***********************************************************************/

void pstext_(float x,float y,char size[PS_STRING_LENGTH],
    char text[PS_STRING_LENGTH])
{
  fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
  fprintf(ps_file,"(%s) %s Print\n",text,size);
}
/***********************************************************************

psubox_  -  Write out unbounded box to PostScript file

***********************************************************************/

void psubox_(float x1,float y1,float x2,float y2,float x3,float y3,
        float x4,float y4)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f Pl4\n",
          x1, y1, x2, y2, x3, y3, x4, y4);
}
/***********************************************************************

psucir_  -  Write out a grey/coloured circle, with no border, to
            PostScript file

***********************************************************************/

void psucir_(float x,float y,float radius)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f Ucircle\n",x,y,radius);
}
/***********************************************************************

match_colour_names  -  Match up the colour names with the Colour
                       Names Table

***********************************************************************/

void match_colour_names(void)
{
    int icolour, iletter, jcolour, last_letter, match;

    /* Check that colour name can be used, adjusting where necessary */
    for (icolour = 0; icolour < MAX_COLOURS; icolour++)
    {
        last_letter = 0;

        /* Loop through all the letters of the name */
        for (iletter = 0; iletter < COL_NAME_LEN; iletter++)
        {
            /* Replace blanks by underscores */
            if (Colour_Table_Name[icolour][iletter] == ' ' ||
                Colour_Table_Name[icolour][iletter] == '\'')
                Colour_Table_Name[icolour][iletter] = '_';
            else
                last_letter = iletter;
        }

        /* Blank out any trailing spaces */
        for (iletter = last_letter + 1; iletter < COL_NAME_LEN; iletter++)
            Colour_Table_Name[icolour][iletter] = ' ';

        /* Check for duplicate colours, discarding second instance */
        for (jcolour = 0; jcolour < icolour; jcolour++)
        {
            if (!strncmp(Colour_Table_Name[jcolour],
                Colour_Table_Name[icolour],COL_NAME_LEN))
                Colour_Table_Name[jcolour][0] = '\0';
        }
    }

    /* Check for atom-colour option on ligand and non-ligand bonds */
    if (!strncmp(object_colour[1],"ATOM           ",15))
      Split_Colour_Ligand_Bonds = TRUE;
    if (!strncmp(object_colour[2],"ATOM           ",15))
      Split_Colour_Nonligand_Bonds = TRUE;

    /* Adjust the object colours */
    for (icolour = 0; Picture->In_Colour == TRUE && icolour < OBJECT_COLOURS;
	 icolour++)
    {
        last_letter = 0;

        /* Loop through all the letters of the name */
        for (iletter = 0; iletter < COL_NAME_LEN; iletter++)
        {
            /* Replace blanks by underscores */
            if (object_colour[icolour][iletter] == ' ')
                object_colour[icolour][iletter] = '_';
            else
                last_letter = iletter;
        }

        /* Blank out any trailing spaces */
        for (iletter = last_letter + 1; iletter < COL_NAME_LEN; iletter++)
            object_colour[icolour][iletter] = ' ';

	/* If the colour of ligand or non-ligand bonds has been defined
	   as ATOM colour, then don't need to match against the colour
	   table */
        match = FALSE;
	if ((icolour == 1 && Split_Colour_Ligand_Bonds == TRUE) ||
	    (icolour == 2 && Split_Colour_Nonligand_Bonds == TRUE))
	  match = TRUE;

        /* Check that this colour has a matching entry in the table
           of colour names */
        for (jcolour = 0; jcolour < MAX_COLOURS && match == FALSE; jcolour++)
        {
            if (!strncmp(object_colour[icolour],
                Colour_Table_Name[jcolour],COL_NAME_LEN))
                match = TRUE;
        }

        /* If failed to find a match, use whichever is the first colour */
        if (match == FALSE)
        {
            printf("*** Warning. Unidentified colour [%s]. Replaced by [%s]\n",
                object_colour[icolour],Colour_Table_Name[0]);
            strcpy(object_colour[icolour],Colour_Table_Name[0]);
	    Nwarnings++;
        }
    }

    /* Adjust the text colours */
    for (icolour = 0; Picture->In_Colour == TRUE && icolour < TEXT_COLOURS; icolour++)
    {
        last_letter = 0;

        /* Loop through all the letters of the name */
        for (iletter = 0; iletter < COL_NAME_LEN; iletter++)
        {
            /* Replace blanks by underscores */
            if (text_colour[icolour][iletter] == ' ')
                text_colour[icolour][iletter] = '_';
            else
                last_letter = iletter;
        }

        /* Blank out any trailing spaces */
        for (iletter = last_letter + 1; iletter < COL_NAME_LEN; iletter++)
            text_colour[icolour][iletter] = ' ';

        /* Check that this colour has a matching entry in the table
           of colour names */
        match = FALSE;
        for (jcolour = 0; jcolour < MAX_COLOURS; jcolour++)
        {
            if (!strncmp(text_colour[icolour],
                Colour_Table_Name[jcolour],COL_NAME_LEN))
                match = TRUE;
        }

        /* If failed to find a match, use whichever is the first colour */
        if (match  == FALSE)
        {
            printf("*** Warning. Unidentified colour [%s]. Replaced by [%s]\n",
                text_colour[icolour],Colour_Table_Name[0]);
            strcpy(text_colour[icolour],Colour_Table_Name[0]);
	    Nwarnings++;
        }
    }
}
/***********************************************************************

define_object_colours  -  Write the object colours and definitions
                          to the PostScript file

***********************************************************************/

void define_object_colours(void)
{
    char colour[COL_NAME_LEN + 1];
    int icolour, jcolour;

    /* Write out the object colours */
    for (icolour = 0; icolour < OBJECT_COLOURS; icolour++)
    {
        /* Retrieve the appropriate description */
        switch (icolour)
        {
            case 0 : strcpy(colour,Colour->Background); break;
            case 1 : strcpy(colour,Colour->Ligand_Bonds); break;
            case 2 : strcpy(colour,Colour->Nonligand_Bonds); break;
            case 3 : strcpy(colour,Colour->Hydrogen_Bonds); break;
            case 4 : strcpy(colour,Colour->External_Bonds); break;
            case 5 : strcpy(colour,Colour->Hydrophobics); break;
            case 6 :
                strcpy(colour,Colour->Accessibility_Min);

                /* Save the RGB values for the minimum accessibility colour */
                for (jcolour = 0; jcolour < MAX_COLOURS; jcolour++)
                {
                    if (!strncmp(object_colour[icolour],
                        Colour_Table_Name[jcolour],COL_NAME_LEN))
                    {
                        Accessibility_Min[0] = Colour_Table[jcolour][0];
                        Accessibility_Min[1] = Colour_Table[jcolour][1];
                        Accessibility_Min[2] = Colour_Table[jcolour][2];
                    }
                }
                break;
            case 7 :
                strcpy(colour,Colour->Accessibility_Max);

                /* Save the RGB values for the maximum accessibility colour */
                for (jcolour = 0; jcolour < MAX_COLOURS; jcolour++)
                {
                    if (!strncmp(object_colour[icolour],
                        Colour_Table_Name[jcolour],COL_NAME_LEN))
                    {
                        Accessibility_Max[0] = Colour_Table[jcolour][0];
                        Accessibility_Max[1] = Colour_Table[jcolour][1];
                        Accessibility_Max[2] = Colour_Table[jcolour][2];
                    }
                }
                break;
            case 8 : strcpy(colour,Colour->Nitrogen); break;
            case 9 : strcpy(colour,Colour->Oxygen); break;
            case 10 : strcpy(colour,Colour->Carbon); break;
            case 11 : strcpy(colour,Colour->Sulphur); break;
            case 12 : strcpy(colour,Colour->Water); break;
            case 13 : strcpy(colour,Colour->Phosphorus); break;
            case 14 : strcpy(colour,Colour->Iron); break;
            case 15 : strcpy(colour,Colour->Other); break;
            case 16 : strcpy(colour,Colour->Atom_Edges); break;
            case 17 : strcpy(colour,Colour->Simple_Residues); break;
            default :
	      printf("*** Program error. Invalid Colours parameter: %d\n",
		     icolour);
	      exit(1);
        }

        /* Write out the colour and its definition to the PostScript file */
	if (!strncmp(object_colour[icolour],"ATOM           ",15))
	  fprintf(ps_file,"%% %s defined as ATOM colour\n",colour);
	else
	  fprintf(ps_file,"/%s { %s } def\n",colour,object_colour[icolour]);
    }
}
/***********************************************************************

define_text_colours  -  Write the text colours and definitions
                        to the PostScript file

***********************************************************************/

void define_text_colours(void)
{
  char colour[COL_NAME_LEN + 1];
  int icolour;

  /* Write out the text colours */
  for (icolour = 0; icolour < TEXT_COLOURS; icolour++)
    {
      /* Retrieve the appropriate description */
      switch (icolour)
        {
	case 0 : strcpy(colour,Text_Colour->Title); break;
	case 1 : strcpy(colour,Text_Colour->Key_Text); break;
	case 2 : strcpy(colour,Text_Colour->Ligand_Residue_Names); break;
	case 3 : strcpy(colour,Text_Colour->Nonligand_Residue_Names); break;
	case 4 : strcpy(colour,Text_Colour->Water_Names); break;
	case 5 : strcpy(colour,Text_Colour->Hydrophobic_Names); break;
	case 6 : strcpy(colour,Text_Colour->Ligand_Atom_Names); break;
	case 7 : strcpy(colour,Text_Colour->Nonligand_Atom_Names); break;
	case 8 : strcpy(colour,Text_Colour->Hbond_Lengths); break;
	default :
	  printf("*** Program error. Invalid Text Colours parameter: %d\n",
		 icolour);
	  exit(1);
        }

      /* Write out the colour and its definition to the PostScript file */
      fprintf(ps_file,"/%s { %s } def\n",colour,text_colour[icolour]);
    }
}
/***********************************************************************

define_object_sizes  -  Write the object size definitions to the
                        PostScript file

***********************************************************************/

void define_object_sizes(float scale_factor,int print_out)
{
  char object[COL_NAME_LEN + 1];
  int isize;
  float simple_thick;

  /* Initialise maximum object sizes */
  Max_object_size = 0.0;

  /* Write out the object sizes */
  for (isize = 0; isize < OBJECT_SIZES; isize++)
    {
      /* Use scaling factor to change from Angstroms to PostScript
	 sizes */
      object_size[isize] = scale_factor * object_size[isize];
      
      /* Retrieve the appropriate description */
      switch (isize)
        {
	case 0 : 
	  strcpy(object,Size->Ligand_Atoms);
	  Size_Val->Ligand_Atoms = object_size[isize];
	  if (Include->Ligand_Atoms == TRUE &&
	      object_size[isize] > Max_object_size)
	    Max_object_size = object_size[isize];
	  break;
	case 1 : 
	  strcpy(object,Size->Nonligand_Atoms);
	  Size_Val->Nonligand_Atoms = object_size[isize];
	  if (Include->Nonligand_Atoms == TRUE &&
	      object_size[isize] > Max_object_size)
	    Max_object_size = object_size[isize];
	  break;
	case 2 :
	  strcpy(object,Size->Waters);
	  Size_Val->Waters = object_size[isize];
	  if (Include->Waters == TRUE &&
	      object_size[isize] > Max_object_size)
	    Max_object_size = object_size[isize];
	  break;
	case 3 : 
	  strcpy(object,Size->Hydrophobics);
	  Size_Val->Hydrophobics = object_size[isize];
	  if (Include->Hydrophobics == TRUE &&
	      object_size[isize] > Max_object_size)
	    Max_object_size = SPOKE_EXTENT * object_size[isize];
	  break;
	case 4 :
	  strcpy(object,Size->Simple_Residues);
	  Size_Val->Simple_Residues = object_size[isize];
	  if (Include->Simple_Ligand_Residues == TRUE &&
	      object_size[isize] > Max_object_size)
	    Max_object_size = object_size[isize];
	  simple_thick = object_size[isize] / 10.0;
	  break;
	case 5 : 
	  strcpy(object,Size->Ligand_Bonds);
	  Size_Val->Ligand_Bonds = object_size[isize];
	  break;
	case 6 : 
	  strcpy(object,Size->Nonligand_Bonds);
	  Size_Val->Nonligand_Bonds = object_size[isize];
	  break;
	case 7 : 
	  strcpy(object,Size->Hydrogen_Bonds);
	  Size_Val->Hydrogen_Bonds = object_size[isize];
	  break;
	case 8 : 
	  strcpy(object,Size->External_Bonds);
	  Size_Val->External_Bonds = object_size[isize];
	  break;
	default :
	  printf("*** Program error. Invalid Size parameter: %d\n",
		 isize);
	  exit(1);
        }

      /* Write out the size definition and its value to the PostScript
	 file */
      if (print_out == TRUE)
	fprintf(ps_file,"/%s { %8.3f } def\n",object,object_size[isize]);
    }
  /* Write out corresponding line-thicknesses */
  if (print_out == TRUE)
    {
      fprintf(ps_file,"/%s { %8.3f } def\n","Default_linewidth",0.2);
      fprintf(ps_file,"/%s { %8.3f } def\n","Arc_linewidth    ",0.75);
      fprintf(ps_file,"/%s { %8.3f } def\n","Simple_width     ",
	      simple_thick);
    }
}
/***********************************************************************

define_text_sizes  -  Write the text size definitions to the
                      PostScript file

***********************************************************************/

void define_text_sizes(float scale_factor,int print_out)
{
  char text[COL_NAME_LEN + 1];
  int isize;

  /* Initialise maximum text height and width */
  Max_Textsize = 0.0;
  Max_Textwidth = 0.0;

  /* Write out the textt sizes */
  for (isize = 0; isize < TEXT_SIZES; isize++)
    {
      /* Use scaling factor to change from Angstroms to PostScript
	 sizes */
      text_size[isize] = scale_factor * text_size[isize];
	
      /* Retrieve the appropriate description */
      switch (isize)
        {
	case 0 :
	  strcpy(text,Text_Size->Ligand_Residue_Names);
	  Text_Size_Val->Ligand_Residue_Names = text_size[isize];
	  if (Include->Residue_Names == TRUE &&
	      text_size[isize] > Max_Textsize)
	    {
	      Max_Textsize = text_size[isize];
	      Max_Textwidth = 10 * Max_Textsize * CHAR_ASPECT;
	    }
	  break;
	case 1 :
	  strcpy(text,Text_Size->Nonligand_Residue_Names);
	  Text_Size_Val->Nonligand_Residue_Names = text_size[isize];
	  if (Include->Residue_Names == TRUE &&
	      text_size[isize] > Max_Textsize)
	    {
	      Max_Textsize = text_size[isize];
	      Max_Textwidth = 10 * Max_Textsize * CHAR_ASPECT;
	    }
	  break;
	case 2 : 
	  strcpy(text,Text_Size->Water_Names);
	  Text_Size_Val->Water_Names = text_size[isize];
	  if (Include->Residue_Names == TRUE &&
	      text_size[isize] > Max_Textsize)
	    {
	      Max_Textsize = text_size[isize];
	      Max_Textwidth = 10 * Max_Textsize * CHAR_ASPECT;
	    }
	  break;
	case 3 :
	  strcpy(text,Text_Size->Hydrophobic_Names);
	  Text_Size_Val->Hydrophobic_Names = text_size[isize];
	  break;
	case 4 : 
	  strcpy(text,Text_Size->Simple_Residue_Names); 
	  Text_Size_Val->Simple_Residue_Names = text_size[isize];
	  break;
	case 5 : 
	  strcpy(text,Text_Size->Ligand_Atom_Names);
	  Text_Size_Val->Ligand_Atom_Names = text_size[isize];
	  if (Include->Atom_Names == TRUE &&
	      text_size[isize] > Max_Textsize)
	    {
	      Max_Textsize = text_size[isize];
	      Max_Textwidth = 4 * Max_Textsize * CHAR_ASPECT;
	    }
	  break;
	case 6 : 
	  strcpy(text,Text_Size->Nonligand_Atom_Names);
	  Text_Size_Val->Nonligand_Atom_Names = text_size[isize];
	  if (Include->Residue_Names == TRUE &&
	      text_size[isize] > Max_Textsize)
	    {
	      Max_Textsize = text_size[isize];
	      Max_Textwidth = 4 * Max_Textsize * CHAR_ASPECT;
	    }
	  break;
	case 7 : 
	  strcpy(text,Text_Size->Hbond_Lengths); 
	  Text_Size_Val->Hbond_Lengths = text_size[isize];
	  break;
	default :
	  printf("*** Program error. Invalid Text Size parameter: %d\n",
		 isize);
	  exit(1);
        }

      /* Write out the size definition and its value to the PostScript
	 file */
      if (print_out == TRUE)
	fprintf(ps_file,"/%s { %8.3f } def\n",text,text_size[isize]);
    }
}
/***********************************************************************

psopen_  -  Open PostScript file and write out header records

***********************************************************************/

void psopen_(void)
{
  char  colno[3], exclamation, percent;
  float grey_shade, xx1, xx2, xy1, xy2;
  int   icolour, igrey;

  /* Initialise variables */
  percent = '%';
  exclamation = '!';

  /* Define border round picture */
  xx1 = (float)BBOXX1;
  xx2 = (float)BBOXX2;
  xy1 = (float)BBOXY1;
  xy2 = (float)BBOXY2;

  /* No border round picture */
  xx1 = -1.0;
  xx2 = 650.0;
  xy1 = -1.0;
  xy2 = 951.0;

  /* Print the plot title */
  printf("\nWriting to PostScript output file: ligplot.ps\n\n");
  if (Print_Title[0] != '\0')
      printf("Title: %s\n\n",Print_Title);

  /* Open the PostScript file */
  if ((ps_file = fopen("ligplot.ps","w")) == NULL)
    {
      printf("\n*** Unable to open PostScript file: ligplot.ps\n");
      exit(1);
    }

  /* Check and adjust colour names in colour definitions */
  match_colour_names();

  /* Write out headings to PostScript file*/
  fprintf(ps_file,"%c%cPS-Adobe-3.0\n",percent,exclamation);
  fprintf(ps_file,"%%%%Creator: Ligplot v.3.0\n");
  fprintf(ps_file,"%%%%DocumentNeededResources: font Times-Roman Symbol\n");
  fprintf(ps_file,"%%%%BoundingBox: (atend)\n");
  fprintf(ps_file,"%%%%Pages: 1\n");
  if (Print_Title[0] != '\0')
      fprintf(ps_file,"%%%%Title: %s\n",Print_Title);
  fprintf(ps_file,"%%%%EndComments\n");
  fprintf(ps_file,"%%%%BeginProlog\n");
  fprintf(ps_file,"/L { moveto lineto stroke } bind def\n");
  fprintf(ps_file,"/G { gsave } bind def\n");
  fprintf(ps_file,"/W { setlinewidth } bind def\n");
  fprintf(ps_file,"/D { setdash } bind def\n");
  fprintf(ps_file,"/Col { setrgbcolor } bind def\n");
  fprintf(ps_file,"/Zero_linewidth { 0.0 } def\n");
  fprintf(ps_file,"/Sphcol { 1 setgray } def\n");
  fprintf(ps_file,"/Sphere {");
  fprintf(ps_file,"  newpath 3 copy 0 360 arc gsave Sphcol fill 0\n");
  fprintf(ps_file,"  setgray 0.5 setlinewidth\n");
  fprintf(ps_file,"  3 copy 0.94 mul 260 350 arc stroke 3 copy 0.87");
  fprintf(ps_file," mul 275 335 arc stroke\n");
  fprintf(ps_file,"  3 copy 0.79 mul 295 315 arc stroke 3 copy 0.8");
  fprintf(ps_file," mul 115 135 arc\n");
  fprintf(ps_file,"  3 copy 0.6 mul 135 115 arcn closepath gsave 1");
  fprintf(ps_file," setgray fill grestore stroke\n");
  fprintf(ps_file,"  3 copy 0.7 mul 115 135 arc stroke 3 copy 0.6");
  fprintf(ps_file," mul 124.9 125 arc\n");
  fprintf(ps_file,"  0.8 mul 125 125.1 arc stroke grestore stroke");
  fprintf(ps_file," } bind def\n");

  /* Loop to write out the different grey-shades */
  for (igrey = 0; igrey < 11; igrey++)
  {
      grey_shade = (float)igrey / 10.0;
      sprintf(colno,"%2d",igrey);
      colno[2] = '\0';

      if (colno[0] == ' ')
          colno[0] = '0';
      fprintf(ps_file,"/Grey%s { %6.2f setgray } def\n",colno,grey_shade);
  }

  /* Loop to write out colour definitions to PostScript file */
  for (icolour = 0; icolour < MAX_COLOURS; icolour++)
  {
      if (Colour_Table_Name[icolour][0] != '\0')
      {
          fprintf(ps_file,"/%s { %8.4f %8.4f %8.4f",
              Colour_Table_Name[icolour],Colour_Table[icolour][0],
              Colour_Table[icolour][1],Colour_Table[icolour][2]);
          fprintf(ps_file," setrgbcolor } def\n");
      }
  }

  /* Write out the object- and text-colours */
  define_object_colours();
  define_text_colours();
  define_object_sizes(1.0,FALSE);
  define_text_sizes(1.0,FALSE);

  /* Print the remainder of the definitions */
  fprintf(ps_file,"/Title_size { %8.3f } def\n",TITLE_SIZE);
  fprintf(ps_file,"/Poly3 { moveto lineto lineto fill grestore } bind ");
  fprintf(ps_file,"def\n");
  fprintf(ps_file,"/Pl3 { 6 copy Poly3 moveto moveto moveto closepath ");
  fprintf(ps_file,"stroke } bind def\n");
  fprintf(ps_file,"/Pline3 { 6 copy Poly3 moveto lineto lineto closepa");
  fprintf(ps_file,"th stroke } bind def\n");
  fprintf(ps_file,"/Poly4 { moveto lineto lineto lineto fill grestore } ");
  fprintf(ps_file,"bind def\n");
  fprintf(ps_file,"/Pl4 { 8 copy Poly4 moveto moveto moveto moveto ");
  fprintf(ps_file,"closepath stroke } bind def\n");
  fprintf(ps_file,"/Pline4 { 8 copy Poly4 moveto lineto lineto lineto");
  fprintf(ps_file," closepath stroke } bind def\n");
  fprintf(ps_file,"/Circol { 1 setgray } def\n");
  fprintf(ps_file,"/Circle { gsave newpath 0 360 arc gsave Circol fill ");
  fprintf(ps_file,"grestore stroke grestore } bind def\n");
  fprintf(ps_file,"/Ucircle { gsave newpath 0 360 arc gsave Circol fill ");
  fprintf(ps_file,"grestore grestore } bind def\n");
  fprintf(ps_file,"/Arc { newpath arc stroke newpath } bind def\n");
  fprintf(ps_file,"/Print { /Times-Roman findfont exch scalefont setfont");
  fprintf(ps_file," show } bind def\n");
  fprintf(ps_file,"/Gprint { /Symbol findfont exch scalefont setfont show");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/Center {\n");
  fprintf(ps_file,"  dup /Times-Roman findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop -2 div exch -3 div rmoveto\n");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/CenterRot90 {\n");
  fprintf(ps_file,"  dup /Times-Roman findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop -2 div exch 3 div exch rmoveto\n");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/UncenterRot90 {\n");
  fprintf(ps_file,"  dup /Times-Roman findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth } bind def\n");
  fprintf(ps_file,"/Rot90 { gsave currentpoint translate 90 rotate }");
  fprintf(ps_file," bind def\n");
  fprintf(ps_file,"%%%%EndProlog\n");
  fprintf(ps_file,"%%%%BeginSetup\n");
  fprintf(ps_file,"1 setlinecap 1 setlinejoin 1 setlinewidth 0 setgray");
  fprintf(ps_file," [ ] 0 setdash newpath\n");
  fprintf(ps_file,"%%%%EndSetup\n");
  fprintf(ps_file,"%%%%Page: 1 1 \n");
  fprintf(ps_file,"/LigplotSave save def\n");
  fprintf(ps_file,"%6.2f%6.2f moveto %7.2f%7.2f lineto%7.2f%7.2f lineto\n",
          xx1,xy1,xx2,xy1,xx2,xy2);
  fprintf(ps_file,"%7.2f%7.2f lineto closepath\n",xx1,xy2);
  fprintf(ps_file,"gsave 1.0000 setgray fill grestore\n");
  fprintf(ps_file,"stroke gsave\n");

  /*Draw in the background colour*/
  if (Picture->In_Colour == TRUE)
    {
      pscomm_("Background");
      pshade_(Colour->Background);
      psubox_(BBOXX1, BBOXY1, BBOXX2, BBOXY1, BBOXX2, BBOXY2, BBOXX1, BBOXY2);
    }
}

/****************   P L O T T I N G    R O U T I N E S   ***************/

/***********************************************************************

adjust_for_landscape  -  Perform boundary adjustments for printing
                         in Landscape, as opposed to Portrait,
                         orientation

***********************************************************************/

void adjust_for_landscape(void)
{

    /* Make necessary adjustments */
    pscomm_("Landscape oritentation");
    Page_Min_x = BBOXY1;
    Page_Max_x = BBOXY2;
    Page_Min_y = BBOXX1;
    Page_Max_y = BBOXX2;
    Plot_Centre_x = (Page_Min_x + Page_Max_x) / 2.0;
    Plot_Centre_y = (Page_Min_y + Page_Max_y) / 2.0;
}
/***********************************************************************

print_title  -  Print the title for the plot

***********************************************************************/

void print_title(char string[TITLE_LEN + 1])
{
  float y;

  /* Define the margin taken up by the title area */
  y = Margin_y + TITLE_FRACTION * (BBOXY2 - BBOXY1) / 2.0;
  Margin_y = Margin_y + TITLE_FRACTION * (BBOXY2 - BBOXY1);

  /* Print the title */
  pscomm_("Plot heading");
  pscolb_(Text_Colour->Title);
  psctxt_(Plot_Centre_x,y,"Title_size",string);
}
/***********************************************************************

print_key_text  -  Print the text items for the key

***********************************************************************/

void print_key_text(float x_pos, float y_pos, float xmin, float ymin,
    float width, float height, char text_size[PS_STRING_LENGTH],
    char text_string[PS_STRING_LENGTH])
{
  float x, y;

  /* Print the text item */
  x = xmin + x_pos * width;
  y = ymin + y_pos * height;
  pstext_(x,y,text_size,text_string);
}
/***********************************************************************

plot_spokes  -  Plot the spokes to signify hydrophobic contacts

***********************************************************************/

void plot_spokes(float x, float y, float centre_x, float centre_y,
                 float scale_factor, float outer_radius, float radius,
                 float arc_angle, float direction_angle, float ang_step)
{
  int done, sign;
  float angle, angle_radians;
  float tot_step;
  float xc, x1, x2, yc, y1, y2;
  float ps_coord_x1, ps_coord_x2, ps_coord_y1, ps_coord_y2;

  /* Loop twice to plot the spokes in the clockwise and anti-clockwise
     directions */
  for (sign = -1; sign < 2; sign = sign + 2)
  {
      /* Initialise start for plotting in this direction */
      angle = ang_step * (int)(0.5 + direction_angle / ang_step);
      tot_step = 0.0;
      done = FALSE;
      while (done == FALSE)
      {
          /* Calculate angle in radians and start- and end-points of spoke */
          angle_radians = angle / RADDEG;
          xc = radius * cos(angle_radians) / scale_factor;
          yc = radius * sin(angle_radians) / scale_factor;
          x1 = x + xc;
          y1 = y + yc;
          xc = outer_radius * cos(angle_radians) / scale_factor;
          yc = outer_radius * sin(angle_radians) / scale_factor;
          x2 = x + xc;
          y2 = y + yc;

          /* Convert to PostScript coordinates */
          ps_coord_x1 = scale_factor * (x1 - centre_x) + Plot_Centre_x;
          ps_coord_x2 = scale_factor * (x2 - centre_x) + Plot_Centre_x;
          ps_coord_y1 = scale_factor * (y1 - centre_y) + Plot_Centre_y;
          ps_coord_y2 = scale_factor * (y2 - centre_y) + Plot_Centre_y;

          /* Plot this spoke */
          psline_(ps_coord_x1,ps_coord_y1,ps_coord_x2,ps_coord_y2);

          /* Increment the angle for the next spoke */
          angle = angle + sign * ang_step;
          if (angle < 0.0)
              angle = angle + 360.0;
          if (angle > 360.0)
              angle = angle - 360.0;

          /* Check whether have done all the spokes required in this
             direction */
          tot_step = tot_step + ang_step;
          if (tot_step > arc_angle / 2.0)
              done = TRUE;
      }
  }
}
/***********************************************************************

plot_contact_group  -  Plot a symbol for hydrophobic contacts - the
                       direction is defined as the start and end of
                       the virtual bond joining the contact group to
                       the atom it makes contact with

***********************************************************************/

void plot_contact_group(float x, float y, float dir_x, float dir_y,
                        float centre_x, float centre_y, float scale_factor,
                        float radius, float arc_angle,
/* v.4.0--> */
/*                        int spokes_only) */
                        int spokes_only,int interface)
/* <--v.4.0 */
{
  float angle_end, angle_start, ang_step, outer_radius;
  float direction_angle;
  float x1, y1;
  float ps_coord_x, ps_coord_y, ps_radius;

  /* Initialise values */
  ang_step = SPOKE_ANGLE;
  if (spokes_only == TRUE)
  {
      radius = radius * (1.75 / 1.4);
      ang_step = SPOKE_ANGLE * 3.0 / 2.0;
  }
  outer_radius = SPOKE_EXTENT * radius;

  /* Initialise PostScript output */
  pssave_();
/* v.4.0--> */
/*  pscolb_(Colour->Hydrophobics); */
  if (interface != 1)
    pscolb_(Colour->Hydrophobics);
  else
    pscolb_(Colour->Ligand_Bonds);
/* <--v.4.0 */

  /* Calculate the angle defining the direction of this group */
  x1 = dir_x - x;
  y1 = dir_y - y;
  if (fabs(x1) > 0.0001)
  {
      direction_angle = RADDEG * atan(y1 / x1);
      if (direction_angle >= 0.0)
      {
          if (x1 < 0.0 && y1 <= 0.0)
            direction_angle = direction_angle + 180.0;
      }
      else
      {
          direction_angle = direction_angle + 360.0;
          if (y1 > 0.0)
            direction_angle = direction_angle - 180.0;
      }
  }
  else if (fabs(y1) < 0.0001)
      direction_angle = 0.0;
  else if (y1 > 0.0)
      direction_angle = 90.0;
  else
      direction_angle = 270.0;

  /* Calculate start- and end-angles of the arc */
  angle_start = direction_angle - arc_angle / 2.0;
  if (angle_start < 0.0)
      angle_start = angle_start + 360.0;
  angle_end = direction_angle + arc_angle / 2.0;
  if (angle_end > 360.0)
      angle_end = angle_end - 360.0;

  /* Plot the arc */
  ps_coord_x = scale_factor * (x - centre_x) + Plot_Centre_x;
  ps_coord_y = scale_factor * (y - centre_y) + Plot_Centre_y;
  ps_radius = radius;
  if (spokes_only == FALSE)
  {
      pslwid_("Arc_linewidth");
      psarc_(ps_coord_x,ps_coord_y,ps_radius,angle_start,angle_end);
  }
  else
      pslwid_("Default_linewidth");

  /* Plot the spokes for the interacting residue */
  plot_spokes(x,y,centre_x,centre_y,scale_factor,outer_radius,radius,
              arc_angle,direction_angle,ang_step);

  /* Reset graphics state */
  psrest_();
}
/***********************************************************************

accessibility_circle  -  Plot the coloured/shaded circle representing
                         the current atom's accessibility

***********************************************************************/

void accessibility_circle(float ps_coord_x, float ps_coord_y,
    float min_circle_size, float max_circle_size, float scale_factor,
    float accessibility)
{
    int icol;
    float angle, circle_size, col[3], ratio, shade;

    /* Calculate the accessibility ratio */
    ratio = accessibility / Maximum_accessibility;
    if (ratio > 1.0)
        ratio = 1.0;
    if (ratio < 0.0)
        ratio = 0.0;
    if (xsite_file == TRUE)
    {
        if (ratio < 0.1)
            ratio = 0.0;
        else if (ratio > 0.25)
            ratio = 1.0;
        else
            ratio = (ratio - 0.1) / 0.15;
    }

    /* Compute the corresponding shade for the circle representing the
       accessibility */
    angle = PI_BY_2 * sqrt(ratio);
    shade = sin(angle);
    for (icol = 0; icol < 3; icol++)
    {
         col[icol] = Accessibility_Min[icol]
             + shade * (Accessibility_Max[icol] - Accessibility_Min[icol]);
         if (xsite_file == TRUE)
             col[icol] = sin(angle);
    }

    /* Draw in the circle of the appropriate size and shade */
    pscrgb_(col[0],col[1],col[2]);
    circle_size = min_circle_size + (max_circle_size - min_circle_size)
        * (1.0 - ratio);
    psucir_(ps_coord_x,ps_coord_y,circle_size * scale_factor);
}
/***********************************************************************

print_symbol_key  -  Print the key to the symbols used in the
                     figure
***********************************************************************/

void print_symbol_key(void)
{
  float margin;
  float yscale, width, text_size, x, x1, x2, xmid, xmin, y, y1, y2, ymin;
  float key_size, key_x, key_y;
  float ligand_atom_size, nonligand_atom_size;
  float atom_x1_1, atom_x1_2, atom_x1_3;
  float atom_x2_1, atom_x2_2, atom_x2_3;
  float atom_x3, atom_x4, atom_x5;
  float atom_y11, atom_y12, atom_y13, atom_y21, atom_y22, atom_y23,
        atom_y3;
  float bond1, bond2, bond3, ps_scale;
  float text_x1_1, text_x1_2, text_x1_3;
  float text_x2, text_x31, text_x32, text_x33;
  float text_y11, text_y12, text_y13, text_y21, text_y22, text_y23,
        text_y3, text_drop;

  /* Define positioning of elements */
  key_x = 0.03;
  key_y = 0.85;
  key_size = 0.17;
  text_size = 0.09;
  text_drop = text_size / 4.0;

  /* Define sizes, ensuring that not too large for key-area */
  ps_scale = 20.0;
  ligand_atom_size = ps_scale *  Size_Val->Ligand_Atoms;
  if (ligand_atom_size > 5.0)
    ligand_atom_size = 5.0;
  nonligand_atom_size = ps_scale *  Size_Val->Nonligand_Atoms;
  if (nonligand_atom_size > 5.0)
    nonligand_atom_size = 5.0;
  bond1 = ps_scale * Size_Val->Ligand_Bonds;
  bond2 = ps_scale * Size_Val->Nonligand_Bonds;
  bond3 = ps_scale * Size_Val->Hydrogen_Bonds;

  /* Bond positioning */
  atom_x1_1 = 0.06;
  atom_x2_1 = 0.12;
  atom_x1_2 = 0.307;
  atom_x2_2 = 0.367;
  atom_x1_3 = 0.59;
  atom_x2_3 = 0.65;
  text_x1_1 = 0.15;
  text_x1_2 = 0.397;
  text_x1_3 = 0.68;
  atom_y11 = 0.64;
  atom_y12 = 0.48;
  atom_y13 = 0.32;
  text_y11 = atom_y11 - text_drop;
  text_y12 = atom_y12 - text_drop;
  text_y13 = atom_y13 - text_drop;

  /* Hydrophobic contacts */
  atom_x3 = 0.46;
  text_x2 = 0.50;
  atom_y21 = 0.62;
  atom_y22 = 0.52;
  atom_y23 = 0.32;
  text_y21 = atom_y21 - text_drop;
  text_y22 = atom_y22 - text_drop;
  text_y23 = atom_y23 - text_drop;

  /* Accessibility shading */
  atom_x4 = 0.332;
  atom_x5 = 0.46;
  text_x31 = 0.052;
  text_x32 = 0.375;
  text_x33 = 0.50;
  atom_y3 = 0.09;
  text_y3 = atom_y3 - text_drop;

  /* Determine size of margin according to what is to be plotted */
  margin = 0.0;
  if (Include->Accessibilities == TRUE)
      margin = MARGIN_ACC * (BBOXY2 - BBOXY1);

  /* If no hydrophobic interactions, then show bonds on a single line */
  if (Include->Hydrophobics == FALSE)
  {
      margin = margin - MARGIN_NO_HPH * (BBOXY2 - BBOXY1);
      atom_y11 = atom_y13;
      atom_y12 = atom_y13;
      text_y11 = text_y13;
      text_y12 = text_y13;
      key_y = 0.53;
  }
  else
  {
      atom_x1_2 = atom_x1_1;
      atom_x1_3 = atom_x1_1;
      atom_x2_2 = atom_x2_1;
      atom_x2_3 = atom_x2_1;
      text_x1_2 = text_x1_1;
      text_x1_3 = text_x1_1;
  }

  /* Define border round the key area */
  x1 = (float)BBOXX1;
  x2 = (float)BBOXX2;
  width = x2 - x1;
  y1 = Margin_y;
  y2 = y1 + margin + KEY_HEIGHT * (BBOXY2 - BBOXY1);
  yscale = KEY_SCALE * (BBOXY2 - BBOXY1);
  Margin_y = y2;

  /* Determine minimum points for text */
  xmin = x1;
  if (Picture->Portrait == FALSE)
      xmin = x1 + (Page_Max_x - Page_Min_x - width) / 2.0;
  ymin = y1;
  if (Include->Accessibilities == TRUE)
      ymin = y1 + MARGIN_ACC * (BBOXY2 - BBOXY1);

  /* Write comment to PostScript file */
  pscomm_("KEY TO SYMBOLS");
  pssave_();

  /* Define the various sizes for text and objects */
  fprintf(ps_file,"/%s { %8.3f } def\n","Default_linewidth",0.2);
  fprintf(ps_file,"/%s { %8.3f } def\n","Arc_linewidth",0.75);
  fprintf(ps_file,"/%s { %8.3f } def\n","Key_size",key_size * yscale);
  fprintf(ps_file,"/%s { %8.3f } def\n","Text_size",text_size * yscale);
  fprintf(ps_file,"/%s { %8.3f } def\n",Size->Ligand_Atoms,ligand_atom_size);
  fprintf(ps_file,"/%s { %8.3f } def\n",Size->Nonligand_Atoms,
	  nonligand_atom_size);
  fprintf(ps_file,"/%s { %8.3f } def\n",Size->Ligand_Bonds,bond1);
  fprintf(ps_file,"/%s { %8.3f } def\n",Size->Nonligand_Bonds,bond2);
  fprintf(ps_file,"/%s { %8.3f } def\n",Size->Hydrogen_Bonds,bond3);
  fprintf(ps_file,"/%s { %8.3f } def\n",Text_Size->Hbond_Lengths,
      text_size * yscale / 2.0);
  fprintf(ps_file,"/%s { %8.3f } def\n",Text_Size->Hydrophobic_Names,
      text_size * yscale * 3.0 / 4.0);

  /* Print key heading */
  pscolb_(Text_Colour->Key_Text);
  pslwid_("Default_linewidth");
  print_key_text(key_x,key_y,xmin,ymin,width,yscale,"Key_size","Key");

  /* Print ligand bond */
  pscomm_("Ligand bond");
  pssave_();
  x1 = xmin + atom_x1_1 * width;
  x2 = xmin + atom_x2_1 * width;
  xmid = (x1 + x2) / 2.0;
  y = ymin + atom_y11 * yscale;
  pslwid_(Size->Ligand_Bonds);

  /* If bond-split required, then print in two halves */
  if (Split_Colour_Ligand_Bonds == TRUE)
    {
      pscolb_(Colour->Nitrogen);
      psline_(x1,y,xmid,y);
      pscolb_(Colour->Carbon);
      psline_(x2,y,xmid,y);
    }
  else
    {
      pscolb_(Colour->Ligand_Bonds);
      psline_(x1,y,x2,y);
    }
  psrest_();

  /* Print the two atoms */
  pssave_();
  if (Include->Ligand_Atoms == TRUE)
    {
      pscolb_(Colour->Atom_Edges);
      psphcl_(Colour->Nitrogen);
      pspher_(x1,y,Size->Ligand_Atoms);
      psphcl_(Colour->Carbon);
      pspher_(x2,y,Size->Ligand_Atoms);
    }
  psrest_();

  /* Print the text */
  if (Include->Simple_Ligand_Residues == TRUE)
      print_key_text(text_x1_1,text_y11,xmin,ymin,width,yscale,
          "Text_size","Ligand m/c");
  else
/* v.4.0--> */
    {
      if (Interface_Plot == FALSE)
/* <--v.4.0 */
	print_key_text(text_x1_1,text_y11,xmin,ymin,width,yscale,
		       "Text_size","Ligand bond");
/* v.4.0--> */
      else
	print_key_text(text_x1_1,text_y11,xmin,ymin,width,yscale,
		       "Text_size","Residues of first surface");
    }
/* <--v.4.0 */

  /* Print non-ligand bond */
  pscomm_("Nonligand bond");
  pssave_();
  x1 = xmin + atom_x1_2 * width;
  x2 = xmin + atom_x2_2 * width;
  xmid = (x1 + x2) / 2.0;
  y = ymin + atom_y12 * yscale;
  pslwid_(Size->Nonligand_Bonds);

  /* If bond-split required, then print in two halves */
  if (Split_Colour_Nonligand_Bonds == TRUE)
    {
      pscolb_(Colour->Nitrogen);
      psline_(x1,y,xmid,y);
      pscolb_(Colour->Carbon);
      psline_(x2,y,xmid,y);
    }
  else
    {
      pscolb_(Colour->Nonligand_Bonds);
      psline_(x1,y,x2,y);
    }
  psrest_();

  /* Print the two atoms */
  pssave_();
  if (Include->Nonligand_Atoms == TRUE)
    {
      pscolb_(Colour->Atom_Edges);
      psphcl_(Colour->Nitrogen);
      pspher_(x1,y,Size->Nonligand_Atoms);
      psphcl_(Colour->Carbon);
      pspher_(x2,y,Size->Nonligand_Atoms);
    }
  psrest_();

  /* Print the text */
  if (Include->Simple_Ligand_Residues == TRUE)
      print_key_text(text_x1_2,text_y12,xmin,ymin,width,yscale,
		     "Text_size","Ligand sidechain");
  else
/* v.4.0--> */
    {
      if (Interface_Plot == FALSE)
/* <--v.4.0 */
	print_key_text(text_x1_2,text_y12,xmin,ymin,width,yscale,
		       "Text_size","Non-ligand bond");
/* v.4.0--> */
      else
	print_key_text(text_x1_2,text_y12,xmin,ymin,width,yscale,
		       "Text_size","Residues of second surface");
    }
/* <--v.4.0 */

  /* Print H-bond */
  if (Include->Hbonds == TRUE)
    {
      pscomm_("Hydrogen bond");
      pssave_();
      x1 = xmin + atom_x1_3 * width;
      x2 = xmin + atom_x2_3 * width;
      y = ymin + atom_y13 * yscale;
      pscolb_(Colour->Hydrogen_Bonds);
      pslwid_(Size->Hydrogen_Bonds);
      psdash_(2);
      psline_(x1,y,x1 + 0.3 * (x2 - x1),y);
      psline_(x2,y,x1 + 0.6 * (x2 - x1),y);
      psdash_(0);
      psrest_();

      /* Print the length */
      x = (x1 + x2) / 2.0;
      pssave_();
      pscolb_(Text_Colour->Hbond_Lengths);
      psctxt_(x,y,Text_Size->Hbond_Lengths,"3.0");
      psrest_();

      /* Print the two atoms */

      /* Print ligand nitrogen atom */
      pssave_();
      if (Include->Ligand_Atoms == TRUE)
	{
	  pscolb_(Colour->Atom_Edges);
	  psphcl_(Colour->Nitrogen);
	  pspher_(x1,y,Size->Ligand_Atoms);
	}

      /* If not plotting atoms, then show tip of the bond */
      else
	{
	  if (Split_Colour_Ligand_Bonds == TRUE)
	    psccol_(Colour->Nitrogen);
	  else
	    psccol_(Colour->Ligand_Bonds);
	  psucir_(x1,y,bond1 / 2.0);
	}
      psrest_();

      /* Print non-ligand oxygen atom */
      pssave_();
      if (Include->Nonligand_Atoms == TRUE)
	{
	  pscolb_(Colour->Atom_Edges);
	  psphcl_(Colour->Oxygen);
	  pspher_(x2,y,Size->Nonligand_Atoms);
	}

      /* If not plotting atoms, then show tip of the bond */
      else
	{
	  if (Split_Colour_Nonligand_Bonds == TRUE)
	    psccol_(Colour->Oxygen);
	  else
	    psccol_(Colour->Nonligand_Bonds);
	  psucir_(x2,y,bond2 / 2.0);
	}
      psrest_();

      /* Print explanatory text */
      if (Include->Hbond_Lengths == TRUE)
	print_key_text(text_x1_3,text_y13,xmin,ymin,width,yscale,
		       "Text_size","Hydrogen bond and its length");
      else
	print_key_text(text_x1_3,text_y13,xmin,ymin,width,yscale,
		       "Text_size","Hydrogen bond");
    }

  /* If hydrophobic contacts have been shown, explain their symbols */
  if (Include->Hydrophobics == TRUE)
  {
      /* Print residues involved in hydrophobic contacts */
      pscomm_("Hydrophobic interaction");
      pssave_();
      x = xmin + atom_x3 * width;
      y = ymin + atom_y21 * yscale;
      pscolb_(Text_Colour->Hydrophobic_Names);
      psctxt_(x,y,Text_Size->Hydrophobic_Names,"His 53");
      psrest_();

      /* Print the carbon atom from the ligand */
      x1 = xmin + atom_x3 * width;
      y1 = ymin + atom_y23 * yscale;
      pssave_();
      if (Include->Ligand_Atoms == TRUE)
	{
	  pscolb_(Colour->Atom_Edges);
	  psphcl_(Colour->Carbon);
	  pspher_(x1,y1,Size->Ligand_Atoms);
	}

      /* If not plotting atoms, then show tip of the bond */
      else
	{
	  if (Split_Colour_Ligand_Bonds == TRUE)
	    psccol_(Colour->Carbon);
	  else
	    psccol_(Colour->Ligand_Bonds);
	  psucir_(x1,y1,bond1 / 2.0);
	}
      psrest_();

      /* Plot the spokes radiating from atom */
      plot_contact_group(x1 - Plot_Centre_x,y1 - Plot_Centre_y,
                         x - Plot_Centre_x,y - Plot_Centre_y,
/* v.4.0--> */
/*                         0.0,0.0,1.0,ligand_atom_size,120.0,TRUE); */
                         0.0,0.0,1.0,ligand_atom_size,120.0,TRUE,0);
/* <--v.4.0 */

      /* Plot the spokes rediating from hydrophobic group */
      plot_contact_group(x - Plot_Centre_x,y - Plot_Centre_y,
                         x1 - Plot_Centre_x,y1 - Plot_Centre_y,
/* v.4.0--> */
/*                         0.0,0.0,1.0,2.2 * ligand_atom_size,120.0,FALSE); */
                         0.0,0.0,1.0,2.2 * ligand_atom_size,120.0,FALSE,0);
/* <--v.4.0 */

      /* Explanatory text: Line 1 */
/* v.4.0--> */
      if (Interface_Plot == FALSE)
/* <--v.4.0 */
	print_key_text(text_x2,text_y21,xmin,ymin,width,yscale,"Text_size",
		       "Non-ligand residues involved in hydrophobic");
/* v.4.0--> */
      else
	print_key_text(text_x2,text_y21,xmin,ymin,width,yscale,"Text_size",
		       "Residues involved in hydrophobic");
/* <--v.4.0 */

      /* Explanatory text: Line 2 */
      print_key_text(text_x2,text_y22,xmin,ymin,width,yscale,
          "Text_size","contact(s)");

      /* Explanatory text: Line 3 */
      print_key_text(text_x2,text_y23,xmin,ymin,width,yscale,
          "Text_size",
          "Corresponding atoms involved in hydrophobic contact(s)");

  }

  /* If accessibility-shading on the plot, explain the shading */
  if (Include->Accessibilities == TRUE)
  {
      /* Show shadings/colours for buried and accessible atoms */
      pscomm_("Accessibility shades/colours");
      pssave_();

      /* Explanatory text: Line 1 */
      print_key_text(text_x31,text_y3,xmin,ymin,width,yscale,
          "Text_size","Solvent accessibility shading:");

      /* Shading for completely buried atoms */
      x = xmin + atom_x4 * width;
      y = ymin + atom_y3 * yscale;
      accessibility_circle(x,y,2.1 * ligand_atom_size,
			   3.0 * ligand_atom_size,1.0,0.0);
      psccol_(Colour->Background);
      psucir_(x,y,2.1 * ligand_atom_size);

      /* Show the atom */
      pssave_();
      if (Include->Ligand_Atoms == TRUE)
	{
	  pscolb_(Colour->Atom_Edges);
	  psphcl_(Colour->Nitrogen);
	  pspher_(x,y,Size->Ligand_Atoms);
	}

      /* If not plotting atoms, then show tip of the bond */
      else
	{
	  if (Split_Colour_Ligand_Bonds == TRUE)
	    psccol_(Colour->Nitrogen);
	  else
	    psccol_(Colour->Ligand_Bonds);
	  psucir_(x,y,bond1 / 2.0);
	}
      psrest_();

      /* Explanatory text: Line 2 - Buried shade/colour */
      print_key_text(text_x32,text_y3,xmin,ymin,width,yscale,
          "Text_size","Buried");

      /* Shading for completely accessible atoms */
      x = xmin + atom_x5 * width;
      y = ymin + atom_y3 * yscale;
      accessibility_circle(x,y,2.1 * ligand_atom_size,
			   3.0 * ligand_atom_size,1.0,
			   0.9 * Maximum_accessibility);
      psccol_(Colour->Background);
      psucir_(x,y,2.1 * ligand_atom_size);

      /* Show the atom */
      pssave_();
      if (Include->Ligand_Atoms == TRUE)
	{
	  pscolb_(Colour->Atom_Edges);
	  psphcl_(Colour->Nitrogen);
	  pspher_(x,y,Size->Ligand_Atoms);
	}

      /* If not plotting atoms, then show tip of the bond */
      else
	{
	  if (Split_Colour_Ligand_Bonds == TRUE)
	    psccol_(Colour->Nitrogen);
	  else
	    psccol_(Colour->Ligand_Bonds);
	  psucir_(x,y,bond1 / 2.0);
	}
      psrest_();
      psrest_();

      /* Explanatory text: Line 2 - Accessible shade/colour */
      print_key_text(text_x33,text_y3,xmin,ymin,width,yscale,
          "Text_size","Highly accessible");
  }

  /* Reset graphics defaults */
  psrest_();
}
/***********************************************************************

extract_minmax  -  Determine the max and min x and y coordinates of the
                   structure

***********************************************************************/

void extract_minmax(float *min_x,float *max_x,float *min_y,float *max_y)
{
  int first;

  float atom_size, minx, maxx, miny, maxy, x, y;

  struct coordinate *atom_ptr;

  /* Initialise variables */
  first = TRUE;
  minx = miny = 0.0;
  maxx = maxy = 1.0;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms to initialise the flag indicating whether
     coords have already been stored */
  while (atom_ptr != NULL)
    {
      /* Get this atom's relative size */
      atom_size = atom_ptr->atom_size;

      /* Get atom's x- and y-coordinates */
      x = atom_ptr->x;
      y = atom_ptr->y;

      /* If this is the first atom, then use it to initialise the
	 minimum and maximum values */
      if (first == TRUE)
	{
	  minx = x - atom_size;
	  maxx = x + atom_size;
	  miny = y - atom_size;
	  maxy = y + atom_size;
	  first = FALSE;
	}

      /* Otherwise, update maximum and minimum values */
      else
	{
	  if ((x - atom_size) < minx)
	    minx = x - atom_size;
	  if ((x + atom_size) > maxx)
	    maxx = x + atom_size;
	  if ((y - atom_size) < miny)
	    miny = y - atom_size;
	  if ((y + atom_size) > maxy)
	    maxy = y + atom_size;
	}

      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next;
    }

  /* Return the computed values */
  *min_x = minx;
  *max_x = maxx;
  *min_y = miny;
  *max_y = maxy;
}
/***********************************************************************

scale  -  Scale coordinates to the PostScript page

***********************************************************************/

void scale(float *scale_factor,float min_x,float max_x,float min_y,
    float max_y)
{
  float diff_x, diff_y, picture_height, picture_width;

  /* Calculate size of picture in real coordinates */
  diff_x = max_x - min_x;
  diff_y = max_y - min_y;

  /* If residue-names being printed, extend the size by the label-width
     and label-height */
  Max_object_width = Max_Textwidth;
  Max_object_height = Max_Textsize;
  if (Max_object_size > Max_object_width)
    Max_object_width = Max_object_size;
  if (Max_object_size > Max_object_height)
    Max_object_height = Max_object_size;

  /* Compute extra margin required */
  diff_x = diff_x + 2.0 * Max_object_width;
  diff_y = diff_y + 2.0 * Max_object_height;

  /* Determine extent of PostScript page that is available for the
     picture area */
  picture_height = Page_Max_y - Margin_y;
  picture_width = Page_Max_x - Page_Min_x;
  Plot_Centre_x = Page_Min_x + picture_width / 2.0;
  Plot_Centre_y = Margin_y + picture_height / 2.0;

  /* Determine whether picture bounded in the x- or y-direction */
  if ((diff_y / diff_x) < (picture_height / picture_width))
    {
      *scale_factor = BORDER_MARGIN * picture_width / diff_x;
    }
  else
    {
      *scale_factor = BORDER_MARGIN * picture_height / diff_y;
    }
}
/***********************************************************************

centre_point  -  Find the centre x- and y-coordinate on the
                 PostScript page

***********************************************************************/

void centre_point(float *centre_x,float *centre_y,float scale_factor,
             float min_x,float max_x,float min_y,float max_y)
{
  float ps_maxx, ps_minx, ps_maxy, ps_miny;

  *centre_x = (max_x + min_x) / 2.0;
  *centre_y = (max_y + min_y) / 2.0;

  /* Calculate PostScript limits of current picture */
  ps_minx = scale_factor * (min_x - *centre_x - Max_object_width)
    + Plot_Centre_x;
  ps_maxx = scale_factor * (max_x - *centre_x + Max_object_width)
    + Plot_Centre_x;
  ps_miny = scale_factor * (min_y - *centre_y - Max_object_height)
    + Plot_Centre_y;
  ps_maxy = scale_factor * (max_y - *centre_y + Max_object_height)
    + Plot_Centre_y;
  if (ps_minx < Page_Min_x || Include->Key == TRUE)
    ps_minx = Page_Min_x;
  if (ps_maxx > Page_Max_x || Include->Key == TRUE)
    ps_maxx = Page_Max_x;
  if (ps_miny < Page_Min_y || Include->Key == TRUE || Print_Title[0] != '\0')
    ps_miny = Page_Min_y;
  if (ps_maxy > Page_Max_y)
    ps_maxy = Page_Max_y;

  /* Use full limits for landscape orientation */
  if (Picture->Portrait == FALSE)
    {
      ps_minx = BBOXX1;
      ps_maxx = BBOXX2;
      ps_miny = BBOXY1;
      ps_maxy = BBOXY2;
    }

  /* Write out "convert" command to PostScript file as a comment (used
     when using "convert" from within a script file to generate a .gif
     file from the ligplot.ps PostScript file) */
  fprintf(ps_file,"%%convert -clip %dx%d+%d-%d\n",(int)(ps_maxx - ps_minx),
	  (int)(ps_maxy - ps_miny),(int) ps_minx,(int) ps_miny);
}
/***********************************************************************
 
sort_accessibilities - Sort routine using bubble sort method. Sorts
                       the atomic accessibilities in ascending or
		       descending order

***********************************************************************/

void sort_accessibilities(struct coordinate **first_stack_ptr,int order)
{
  int flipped;

  float accessibility, other_accessibility;

  struct coordinate *atom_ptr, *last_atom_ptr, *other_atom_ptr;

  /* Initialise variables */
  flipped = TRUE;

  /* Loop until all the values have been sorted */
  while (flipped == TRUE)
    {
      /* Initialise flag */
      flipped = FALSE;

      /* Get pointer to the first of the atoms */
      atom_ptr = *first_stack_ptr;
      last_atom_ptr = NULL;

      /* Loop through all the atoms */
      while (atom_ptr != NULL)
	{
	  /* Get the atom's accessibility */
	  accessibility = order * atom_ptr->accessibility;

	  /* Get pointer to the next atom */
	  other_atom_ptr = atom_ptr->next_stack_ptr;

	  /* If have the next atom, then proceed */
	  if (other_atom_ptr != NULL)
	    {
	      /* Get the other atom's accessibility */
	      other_accessibility = order * other_atom_ptr->accessibility;

	      /* If other atom's accessibility is smaller,
		 then swap the order of the atoms */
	      if (other_accessibility < accessibility)
		{
		  /* Set flag to indicate that swap occurred */
		  flipped = TRUE;

		  /* If this is the first atom, get the first pointer
		     to point to the second atom */
		  if (last_atom_ptr == NULL)
		      *first_stack_ptr = other_atom_ptr;

		  /* Otherwise, set the last atom's pointer to the
		     other atom */
		  else
		    last_atom_ptr->next_stack_ptr = other_atom_ptr;

		  /* Swap the pointers */
		  atom_ptr->next_stack_ptr = other_atom_ptr->next_stack_ptr;
		  other_atom_ptr->next_stack_ptr = atom_ptr;

		  /* Set the other atom as the last atom */
		  last_atom_ptr = other_atom_ptr;
		}

	      /* Otherwise, set the current atom as the last atom */
	      else
		last_atom_ptr = atom_ptr;
	    }
	  else
	    last_atom_ptr = atom_ptr;

	  /* Get pointer to the next atom in the stack */
	  atom_ptr = last_atom_ptr->next_stack_ptr;
	}
    }
}
/***********************************************************************

shade_accessibilities  -  Shade in the accessibility values for each
                          atom in the ligand

***********************************************************************/

void shade_accessibilities(float centre_x, float centre_y,
			   float scale_factor)
{
  int iatom, natoms, order;
/* v.4.0--> */
  int plot_accessibility;
/* <--v.4.0 */

  float ps_coord_x, ps_coord_y;
  float min_circle_size, max_circle_size, circle_size;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr, *bond_atom_ptr[2], *first_stack_ptr,
  *last_stack_ptr;
  struct object *object_ptr;
  struct object_bond *object_bond_ptr;

  /* Initialise plot of atom accessibilities */
  printf("   Shading accessibilities ...\n");
  min_circle_size = 2.8 * Size_Val->Ligand_Atoms / scale_factor;
  max_circle_size = 4.0 * Size_Val->Ligand_Atoms / scale_factor;
  pscomm_("Atomic accessibilities");

  /* Adjust for maximum accessibility read in */
  if (Maximum_accessibility == 0.0)
    Maximum_accessibility = 1.0;

  /* For an XSITE file, set the maximum accessibility to 90.0% */
  if (xsite_file == TRUE)
    Maximum_accessibility = 100.0;

  /* Initialise all the atom flags */
  initialise_atoms();

  /* Initialise first and last stack pointers */
  first_stack_ptr = NULL;
  last_stack_ptr = NULL;
  natoms = 0;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
/* v.4.0--> */
      /* Determine whether accessibility to be plotted */
      plot_accessibility = TRUE;
      if (object_ptr->object_type == HYDROPHOBIC)
	plot_accessibility = FALSE;
      else if (object_ptr->object_type != LIGAND)
	{
	  if (object_ptr->object_type == WATER)
	    plot_accessibility = FALSE;
	  else if (Include->Ligand_Accessibilities_Only == TRUE)
	    plot_accessibility = FALSE;
	}
/* <--v.4.0 */
      /* If object is not a hydrophobic group, then plot all its atoms */
/* v.4.0--> */
/*      if (object_ptr->object_type != HYDROPHOBIC) */
      if (plot_accessibility == TRUE)
/* <--v.4.0 */
	{
	  /* Get pointer to the first of this object's bonds */
	  object_bond_ptr = object_ptr->first_object_bond_ptr;
      
	  /* Loop through all object's bonds */
	  while (object_bond_ptr != NULL)
	    {
	      /* Get the bond pointer */
	      bond_ptr = object_bond_ptr->bond_ptr;

	      /* Get pointers to the bond's two atoms */
	      bond_atom_ptr[0] = bond_ptr->first_atom_ptr;
	      bond_atom_ptr[1] = bond_ptr->second_atom_ptr;

	      /* Loop to plot the bond's two atoms */
	      for (iatom = 0; iatom < 2; iatom++)
		{
		  /* Check that this atom hasn't already been plotted */
		  if (bond_atom_ptr[iatom]->checked == FALSE)
		    {
		      /* Add atom to the stack */

		      /* If this is the first atom on the stack, set
			 first stack pointer */
		      if (first_stack_ptr == NULL)
			first_stack_ptr = bond_atom_ptr[iatom];

		      /* Otherwise, get the last atom to point to
			 the current atom */
		      else
			last_stack_ptr->next_stack_ptr = bond_atom_ptr[iatom];

		      /* Save the current atom's pointer */
		      last_stack_ptr = bond_atom_ptr[iatom];
		      natoms++;

		      /* Update flag to indicate atom has already been
			 plotted */
		      bond_atom_ptr[iatom]->checked = TRUE;
		    }
		}
      
	      /* Get pointer to this object's next bond entry */
	      object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
	    }
	}
      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* Process if there are any atoms to be done */
  if (natoms > 0)
    {
      /* Set the sort-order to be ascending order */
      order = 1;

      /* If this is an XSITE file, then plot the shades in the reverse
	 order so that the darker ones come out on top */
      if (xsite_file == TRUE)
	order = -1;

      /* Sort all the stored accessibility values in ascending order */
      sort_accessibilities(&first_stack_ptr,order);

      /* Plot all the accessibility values as shaded circles */
      
      /* Get pointer to the first of the atoms to be plotted */
      atom_ptr = first_stack_ptr;

      /* Loop through all the atoms in order of increasing accessibility */
      while (atom_ptr != NULL)
	{
	  /* Calculate PostScript coords of the atom's centre */
	  ps_coord_x = scale_factor * (atom_ptr->x - centre_x)
	    + Plot_Centre_x;
	  ps_coord_y = scale_factor * (atom_ptr->y - centre_y)
	    + Plot_Centre_y;

	  /* Set the shade/colour of the circle according
	     to this atom's accessibility */
	  accessibility_circle(ps_coord_x,ps_coord_y,min_circle_size,
			       max_circle_size,scale_factor,
			       atom_ptr->accessibility);

	  /* Get pointer to the next atom */
	  atom_ptr = atom_ptr->next_stack_ptr;
	}
    }

  /* Set the shade to the background colour */
  pscomm_("Blank out centres");
  psccol_(Colour->Background);

  /* Get pointer to the first of the atoms again */
  atom_ptr = first_stack_ptr;

  /* Loop through all the atoms a second time to blank out the
       central region */
  while (atom_ptr != NULL)
    {
      /* Calculate PostScript coords of the atom's centre */
      ps_coord_x = scale_factor * (atom_ptr->x - centre_x)
	+ Plot_Centre_x;
      ps_coord_y = scale_factor * (atom_ptr->y - centre_y)
	+ Plot_Centre_y;

      /* Draw in the circle */
      circle_size = min_circle_size;
      psucir_(ps_coord_x,ps_coord_y,circle_size * scale_factor);

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next_stack_ptr;
    }
}
/***********************************************************************

get_atom_colour  -  Get the colour of the current atom

***********************************************************************/

void get_atom_colour(char atom_type[5],char res_name[4],
		     char colour[COL_NAME_LEN + 1])
{
  /* Determine the atom-type and its colour */

  /* Nitrogen */
  if (atom_type[1] == 'N' && atom_type[0] != 'Z' && atom_type[0] != 'M')
    strcpy(colour,Colour->Nitrogen);

  /* Carbon */
  else if (atom_type[1] == 'C')
    strcpy(colour,Colour->Carbon);

  /* Water */
  else if (!strncmp(res_name,"HOH",3))
    strcpy(colour,Colour->Water);

  /* Oxygen */
  else if (atom_type[1] == 'O' && atom_type[0] != 'C' && atom_type[0] != 'M')
    strcpy(colour,Colour->Oxygen);

  /* Sulphur */
  else if (atom_type[1] == 'S')
    strcpy(colour,Colour->Sulphur);

  /* Phosphorus */
  else if (atom_type[1] == 'P')
    strcpy(colour,Colour->Phosphorus);

  /* Iron */
  else if (!strncmp(atom_type,"FE",2))
    strcpy(colour,Colour->Iron);

  /* Default colour */
  else
    strcpy(colour,Colour->Other);
  }
/***********************************************************************

double_bond_split  -  Calculate the coordinate shifts required for
                      printing a double- or triple_bond

***********************************************************************/

void double_bond_split(float x1,float y1,float x2,float y2,
		       float line_width,float scale_factor,
		       int single_double,float *ps_dx,float *ps_dy)
{
  float dx, dy, gap, len, x, y;

  /*Initialize */
  dx = 0.0;
  dy = 0.0;

  /* Check that bond is either a double- or triple-bond */
  if (single_double == 2 || single_double == 3)
    {
      /* Get the lengths of the bond in the x- and y-directions */
      x = x2 - x1;
      y = y2 - y1;

      /* Get the size of the gap between the two parallel lines defining
	 the current double-bond */
      gap = 0.75 * line_width / scale_factor;
      if (single_double == 't')
	gap = 1.5 * line_width / scale_factor;

      /* Calculate special cases where bond is close to horizontal or
	 close to vertical */
      if (fabs(x) < 0.001)
	{
	  dx = gap;
	  dy = 0.0;
	}
      else if (fabs(y) < 0.001)
	{
	  dx = 0.0;
	  dy = gap;
	}

      /* Otherwise calculate the shift-sizes in the x- and y-directions */
      else
	{
	  len = x * x + y * y;
	  len = sqrt((double) len);
	  dx = gap * y / len;
	  dy = gap * x / len; 
	}
    }

  /* Compute the shift-lengths as scaled to the PostScript coordinates */
  *ps_dx = scale_factor * dx;
  *ps_dy = scale_factor * dy;
}
/***********************************************************************

plot_bond_lines  -  Print all the bond lines in each object in turn
  
***********************************************************************/

void plot_bond_lines(float centre_x,float centre_y,float scale_factor)
{
  char atom_type[5];
  char atom1_colour[COL_NAME_LEN + 1], atom2_colour[COL_NAME_LEN + 1];

  int iobject, width_change;
  int object_type;

  float ps_dx, ps_dy, line_width;
  float ps_coord_ax, ps_coord_ay, ps_coord_bx, ps_coord_by;
  float ps_coord_midx, ps_coord_midy;
  
  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct object *object_ptr;
  struct object_bond *object_bond_ptr;

  /* Initialise print */
  atom_type[0] = '\0';
  pscomm_("Bond lines");
  pssave_();
  line_width = Size_Val->Ligand_Bonds;
  ps_dx = 0.0;
  ps_dy = 0.0;

  /* Initialise pointer to the first stored object */
  printf("   Plotting covalent bonds ...\n");
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Get this object's type */
      object_type = object_ptr->object_type;

      /* Process only if not a hydrophobic group, and not a simplified
         H-group */
      if (object_type != HYDROPHOBIC && object_type != SIMPLE_HGROUP)
	{
	  /* Set default line colour and thickness according to the
	     object type */
	  if (object_type == LIGAND)
	    {
	      if (Split_Colour_Ligand_Bonds == FALSE)
		pscolb_(Colour->Ligand_Bonds);
	      pslwid_(Size->Ligand_Bonds);
	    }

	  /* Default colours and lines thicknesses for non-ligand bonds */
	  else
	    {
	      if (Split_Colour_Nonligand_Bonds == FALSE)
		pscolb_(Colour->Nonligand_Bonds);
	      pslwid_(Size->Nonligand_Bonds);
	    }

	  /* Get pointer to the first of this object's bonds */
	  object_bond_ptr = object_ptr->first_object_bond_ptr;

	  /* Loop through all this object's bonds, one by one */
	  while (object_bond_ptr != NULL)
	    {
	      /* Get the current bond */
	      bond_ptr = object_bond_ptr->bond_ptr;

	      /* Get pointers to the two atoms at either end of the
		 bond */
	      atom1_ptr = bond_ptr->first_atom_ptr;
	      atom2_ptr = bond_ptr->second_atom_ptr;

	      /* Compute the PostScript coordinates of the two ends
		 of the bond */
	      ps_coord_ax = Plot_Centre_x
		+ scale_factor * (atom1_ptr->x - centre_x);
	      ps_coord_ay = Plot_Centre_y
		+ scale_factor * (atom1_ptr->y - centre_y);
	      ps_coord_bx = Plot_Centre_x
		+ scale_factor * (atom2_ptr->x - centre_x);
	      ps_coord_by = Plot_Centre_y
		+ scale_factor * (atom2_ptr->y - centre_y);

	      /* If double- and triple-bonds to be plotted, calculate
		 coordinate-shifts required */
	      if (Include->Double_Bonds == TRUE)
		double_bond_split(atom1_ptr->x,atom1_ptr->y,atom2_ptr->x,
				  atom2_ptr->y,line_width,scale_factor,
				  bond_ptr->bond_order,
				  &ps_dx,&ps_dy);

	      /* Plot the bond */

	      /* If bond of a single colour, then plot */
	      if ((object_type == LIGAND &&
		   Split_Colour_Ligand_Bonds == FALSE) ||
		  (object_type == HGROUP &&
		   Split_Colour_Nonligand_Bonds == FALSE))
		{
		  /* If this is a double- or triple-bond, plot a
		     double-line */
		  width_change = FALSE;
		  if (Include->Double_Bonds == TRUE &&
		      (bond_ptr->bond_order == 2 ||
		       bond_ptr->bond_order == 3))
		    {
		      /* If this is a ligand bond, then need to make thinner */
		      if (object_type == LIGAND)
			{
			  pslwid_(Size->Nonligand_Bonds);
			  width_change = TRUE;
			}

		      /* Draw the two lines representing the double bond */
		      psline_(ps_coord_ax - ps_dx,ps_coord_ay + ps_dy,
			      ps_coord_bx - ps_dx,ps_coord_by + ps_dy);
		      psline_(ps_coord_ax + ps_dx,ps_coord_ay - ps_dy,
			      ps_coord_bx + ps_dx,ps_coord_by - ps_dy);
		    }

		  /* For a single or triple-bond, print line between
		     atom centres */
		  if (Include->Double_Bonds == FALSE ||
		      (Include->Double_Bonds == TRUE &&
		       (bond_ptr->bond_order == 1 ||
			bond_ptr->bond_order == 3)))
		    psline_(ps_coord_ax,ps_coord_ay,ps_coord_bx,ps_coord_by);
	      
		  /* If ligand-line was made thin, revert to usual width */
		  if (width_change == TRUE)
		    pslwid_(Size->Ligand_Bonds);
		}

	      /* If split-colour bonds required, then process */
	      else
		{
		  /* Calculate mid-point coordinates */
		  ps_coord_midx = (ps_coord_ax + ps_coord_bx) / 2.0;
		  ps_coord_midy = (ps_coord_ay + ps_coord_by) / 2.0;

		  /* Get the colours of the two atoms at either end of
		     the bond */
		  get_atom_colour(atom1_ptr->atom_type,
				  atom1_ptr->residue_ptr->res_name,
				  atom1_colour);
		  get_atom_colour(atom2_ptr->atom_type,
				  atom2_ptr->residue_ptr->res_name,
				  atom2_colour);

		  /* Plot the two halves of the line */

		  /* If this is a double-bond, plot a double-line */
		  if (Include->Double_Bonds == TRUE &&
		      (bond_ptr->bond_order == 2 ||
		       bond_ptr->bond_order == 3))
		    {
		      /* If this is a ligand bond, then need to make thinner */
		      if (object_type == LIGAND)
			{
			  pslwid_(Size->Nonligand_Bonds);
			  width_change = TRUE;
			}

		      /* Draw the two lines representing the double bond */
		      pscolb_(atom1_colour);
		      psline_(ps_coord_ax - ps_dx,ps_coord_ay + ps_dy,
			      ps_coord_midx - ps_dx,ps_coord_midy + ps_dy);
		      psline_(ps_coord_ax + ps_dx,ps_coord_ay - ps_dy,
			      ps_coord_midx + ps_dx,ps_coord_midy - ps_dy);
		      pscolb_(atom2_colour);
		      psline_(ps_coord_midx - ps_dx,ps_coord_midy + ps_dy,
			      ps_coord_bx - ps_dx,ps_coord_by + ps_dy);
		      psline_(ps_coord_midx + ps_dx,ps_coord_midy - ps_dy,
			      ps_coord_bx + ps_dx,ps_coord_by - ps_dy);
		    }

		  /* For a single or triple-bond, print line between
		     atom centres */
		  if (Include->Double_Bonds == FALSE ||
		      (Include->Double_Bonds == TRUE &&
		       (bond_ptr->bond_order == 1 ||
			bond_ptr->bond_order == 3)))
		    {
		      pscolb_(atom1_colour);
		      psline_(ps_coord_ax,ps_coord_ay,ps_coord_midx,
			      ps_coord_midy);
		      pscolb_(atom2_colour);
		      psline_(ps_coord_midx,ps_coord_midy,ps_coord_bx,
			      ps_coord_by);
	      
		      /* If ligand-line was made thin, revert to usual width */
		      if (width_change == TRUE)
			pslwid_(Size->Ligand_Bonds);
		    }
		}

	      /* Get pointer to this object's next bond entry */
	      object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
	    }
	}

      /* Get pointer for next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

  /* Restore the graphics state */
    psrest_();
}
/***********************************************************************

print_hbond_lengths  -  Prints H-bond lengths on the H-bond lines

***********************************************************************/

void print_hbond_lengths(float bond_length,float x1,float y1,
			 float x2,float y2,float scale_factor,
			 float centre_x,float centre_y)
{
  char hbond_dist[5];
  float sizex, sizey;
  float ps_coord_x, ps_coord_y, xc, yc;

  /* Get the bond length */
  sprintf(hbond_dist,"%4.2f",bond_length);
  hbond_dist[4] = '\0';

  /* Get the PostScript position of the label */
  xc = (x1 + x2) / 2;
  yc = (y1 + y2) / 2;

  /* Save the coordinates of the label */
  if (n_labels < MAXLABELS)
    {
      label_coord[n_labels][0] = xc;
      label_coord[n_labels][1] = yc;
      n_labels++;
    }
      
  /* Save the coordinates of the quarter- and three-quarter
     points along the H-bond, too */
  if (n_labels < MAXLABELS)
    {
      label_coord[n_labels][0] = x1 + 0.25 * (x2 - x1);
      label_coord[n_labels][1] = y1 + 0.25 * (y2 - y1);
      n_labels++;
    }
  if (n_labels < MAXLABELS)
    {
      label_coord[n_labels][0] = x1 + 0.75 * (x2 - x1);
      label_coord[n_labels][1] = y1 + 0.75 * (y2 - y1);
      n_labels++;
    }

  /* Convert coordinates to PostScript coordinates */
  ps_coord_x = scale_factor * (xc - centre_x) + Plot_Centre_x;
  ps_coord_y = scale_factor * (yc - centre_y) + Plot_Centre_y;
      
  /* Calculate size of rectangle required to blank out a region
     over the line such that can write the H-bond length clearly */
  sizex = 4 * Text_Size_Val->Hbond_Lengths * CHAR_ASPECT / 2.0;
  sizey = 0.6 * Text_Size_Val->Hbond_Lengths;

  /* Blank out the rectangle */
  pssave_();
  pscolb_(Colour->Background);
  pssave_();
  pslwid_("Zero_linewidth");
  pscolb_(Colour->Background);
  psubox_(ps_coord_x - sizex,ps_coord_y + sizey,
	  ps_coord_x - sizex,ps_coord_y - sizey,
	  ps_coord_x + sizex,ps_coord_y - sizey,
	  ps_coord_x + sizex,ps_coord_y + sizey);

  /* Print the H-bond length */
  pscolb_(Colour->Hydrogen_Bonds);
  psctxt_(ps_coord_x,ps_coord_y,Text_Size->Hbond_Lengths,hbond_dist);

  /* Reset graphics state */
  psrest_();
}
/***********************************************************************

plot_hbond_lines  -  Print all the H-bond lines
  
***********************************************************************/

void plot_hbond_lines(float centre_x,float centre_y,float scale_factor)
{
  char atom_type[5];

  int bond_type, last_hbond;

  float x1, x2, y1, y2;
  float ps_coord_ax, ps_coord_ay, ps_coord_bx, ps_coord_by;
  
  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;

  /* Initialise print */
  atom_type[0] = '\0';
  pscomm_("Hydrogen bond lines");
  pssave_();
  pscolb_(Colour->Hydrogen_Bonds);
  pslwid_(Size->Hydrogen_Bonds);
  psdash_(3);
  last_hbond = TRUE;

  /* Initialise bond pointer to the first bond */
  printf("   Plotting H-bonds ...\n");
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds */
  while (bond_ptr != NULL)
    {
      /* Get this bond's type */
      bond_type = bond_ptr->bond_type;

      /* Process only if this is a hydrogen bond */
      if (bond_type == HBOND ||
	  (Include->Internal_Hbonds == TRUE && bond_type == HBOND) ||
	  (Include->Hydrophobic_Bonds == TRUE && bond_type == CONTACT))
	{
	  /* Get pointers to the two atoms at either end of the
	     bond */
	  atom1_ptr = bond_ptr->first_atom_ptr;
	  atom2_ptr = bond_ptr->second_atom_ptr;

	  /* Get the atom coordinates */
	  x1 = atom1_ptr->x;
	  y1 = atom1_ptr->y;
	  x2 = atom2_ptr->x;
	  y2 = atom2_ptr->y;

	  /* Compute the PostScript coordinates of the two ends
	     of the bond */
	  ps_coord_ax = Plot_Centre_x + scale_factor * (x1 - centre_x);
	  ps_coord_ay = Plot_Centre_y + scale_factor * (y1 - centre_y);
	  ps_coord_bx = Plot_Centre_x + scale_factor * (x2 - centre_x);
	  ps_coord_by = Plot_Centre_y + scale_factor * (y2 - centre_y);

	  /* Change colour, if necessary */
	  if (bond_type == HBOND && last_hbond == FALSE)
	    {
	      pscolb_(Colour->Hydrogen_Bonds);
	      last_hbond = TRUE;
	    }
	  else if (bond_type != HBOND && last_hbond == TRUE)
	    {
	      pscolb_(Colour->Hydrophobics);
	      last_hbond = FALSE;
	    }

	  /* Plot the H-bond */
	  psline_(ps_coord_ax,ps_coord_ay,ps_coord_bx,ps_coord_by);
	      
	  /* If Hbond-lengths are required, then print */
	  if (bond_type == HBOND && Include->Hbond_Lengths == TRUE)
	    print_hbond_lengths(bond_ptr->bond_length,x1,y1,x2,y2,
				scale_factor,centre_x,centre_y);
	}

      /* Get pointer for next object */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Restore the graphics state */
  psdash_(0);
  psrest_();
}
/***********************************************************************

plot_ext_bond_lines  -  Print any covalent bonds between objects
  
***********************************************************************/

void plot_ext_bond_lines(float centre_x,float centre_y,float scale_factor)
{
  char atom_type[5];

  int bond_type;

  float x1, x2, y1, y2;
  float ps_coord_ax, ps_coord_ay, ps_coord_bx, ps_coord_by;
  
  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;

  /* Initialise print */
  atom_type[0] = '\0';
  pscomm_("External bond lines");
  pssave_();
  pscolb_(Colour->External_Bonds);
  pslwid_(Size->External_Bonds);
/* v.3.2--> */
  if (Include->External_Bonds_Solid == FALSE)
/* <--v.3.2 */
    psdash_(3);

  /* Initialise bond pointer to the first bond */
  printf("   Plotting additional bonds ...\n");
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds */
  while (bond_ptr != NULL)
    {
      /* Get this bond's type */
      bond_type = bond_ptr->bond_type;

      /* Process only if this is a covalent, elastic bond */
      if (bond_type == COVALENT && bond_ptr->elastic == TRUE)
	{
	  /* Get pointers to the two atoms at either end of the
	     bond */
	  atom1_ptr = bond_ptr->first_atom_ptr;
	  atom2_ptr = bond_ptr->second_atom_ptr;

	  /* Get the atom coordinates */
	  x1 = atom1_ptr->x;
	  y1 = atom1_ptr->y;
	  x2 = atom2_ptr->x;
	  y2 = atom2_ptr->y;

	  /* Compute the PostScript coordinates of the two ends
	     of the bond */
	  ps_coord_ax = Plot_Centre_x + scale_factor * (x1 - centre_x);
	  ps_coord_ay = Plot_Centre_y + scale_factor * (y1 - centre_y);
	  ps_coord_bx = Plot_Centre_x + scale_factor * (x2 - centre_x);
	  ps_coord_by = Plot_Centre_y + scale_factor * (y2 - centre_y);

	  /* Plot the bond */
	  psline_(ps_coord_ax,ps_coord_ay,ps_coord_bx,ps_coord_by);
	}

      /* Get pointer for next object */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Restore the graphics state */
/* v.3.2--> */
  if (Include->External_Bonds_Solid == FALSE)
/* <--v.3.2 */
    psdash_(0);
  psrest_();
}
/***********************************************************************

get_segment  -  Determine which segment the given bond occupies

***********************************************************************/

void get_segment(float x1,float y1,float x2,float y2,int segment[8],
		 int first_time)
{
  int iseg;

  float angle, dx, dy;

  /* Initialise variables */
  if (first_time == TRUE)
    for (iseg = 0; iseg < 8; iseg++)
      segment[iseg] = FREE;

  /* Determine the direction of the current bond */
  dx =  x2 - x1;
  dy =  y2 - y1;
  if (fabs(dx) > 0.0001)
    angle = atan(fabs(dy) / fabs(dx));
  else
    angle = PI / 2.0;

  /* Convert angle into degrees */
  angle = angle * RADDEG;

  /* Adjust angle to get in right quadrant */
  if (dx < 0.0 && dy >= 0.0)
    angle = 180.0 - angle;
  else if (dx < 0.0 && dy < 0.0)
    angle = 180.0 + angle;
  else if (dx >= 0.0 && dy < 0.0)
    angle = 360.0 - angle;

  /* Determine which quadrant the bond lies in */
  if (angle > 10.0 && angle <= 80.0)
    segment[0] = NOT_FREE;
  else if (angle > 80.0 && angle <= 100.0)
    segment[1] = NOT_FREE;
  else if (angle > 100.0 && angle <= 170.0)
    segment[2] = NOT_FREE;
  else if (angle > 170.0 && angle <= 190.0)
    segment[3] = NOT_FREE;
  else if (angle > 190.0 && angle <= 260.0)
    segment[4] = NOT_FREE;
  else if (angle > 260.0 && angle <= 280.0)
    segment[5] = NOT_FREE;
  else if (angle > 280.0 && angle <= 350.0)
    segment[6] = NOT_FREE;
  else if ((angle > 350.0 && angle <= 360.0) ||
	   (angle >= 0.0 && angle <= 10.0))
    segment[7] = NOT_FREE;
}
/***********************************************************************

position_names  -  Check the segment-regions around the current atom to
                   see which does not contain a bond and hence is
                   a suitable site for placement of the atom name

***********************************************************************/

void position_names(struct coordinate *atom_ptr,float *x_pos,float *y_pos)
{
  int first_time;
  int segment[8];

  float x1, x2, y1, y2;

  struct bond *bond_ptr;
  struct coordinate *other_atom_ptr;

  /* Get the current atom's coordinates */
  x1 = atom_ptr->x;
  y1 = atom_ptr->y;
  first_time = TRUE;

  /* Get pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds, picking any involving the current atom */
  while (bond_ptr != NULL)
    {
      /* Process only if this is a covalent bond */
      if (bond_ptr->bond_type == COVALENT)
	{
	  /* Initialise pointer to other atom */
	  other_atom_ptr = NULL;

	  /* Test whether this bond links to the current atom */
	  if (bond_ptr->first_atom_ptr == atom_ptr)
	    other_atom_ptr = bond_ptr->second_atom_ptr;
	  if (bond_ptr->second_atom_ptr == atom_ptr)
	    other_atom_ptr = bond_ptr->first_atom_ptr;
	   
	  /* If bond contains the atom at either end, then process */
	  if (other_atom_ptr != NULL)
	    {
	      /* Get the other atom's coordinates */
	      x2 = other_atom_ptr->x;
	      y2 = other_atom_ptr->y;

	      /* Determine which segment the current bond lies in */
	      get_segment(x1,y1,x2,y2,segment,first_time);
	      first_time = FALSE;
	    }
	}

      /* Get pointer to next atom in linked-list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Check which segments contain one or more bonds and which are
     free for placement of the atom name */
  if (segment[0] == FREE)
    {
      *x_pos =  1.0;
      *y_pos =  1.0;
    }
  else if (segment[2] == FREE)
    {
      *x_pos = -1.0;
      *y_pos =  1.0;
    }
  else if (segment[4] == FREE)
    {
      *x_pos = -1.0;
      *y_pos = -1.0;
    }
  else if (segment[6] == FREE)
    {
      *x_pos =  1.0;
      *y_pos = -1.0;
    }
  else if (segment[1] == FREE)
    {
      *x_pos =  0.0;
      *y_pos =  1.0;
    }
  else if (segment[3] == FREE)
    {
      *x_pos = -1.0;
      *y_pos =  0.0;
    }
  else if (segment[5] == FREE)
    {
      *x_pos =  0.0;
      *y_pos = -1.0;
    }
  else if (segment[7] == FREE)
    {
      *x_pos =  1.0;
      *y_pos =  0.0;
    }
}
/***********************************************************************

print_current_atom  -  Print out the current atom as a sphere in the
                       PostScript file

***********************************************************************/

void print_current_atom(struct coordinate *atom_ptr,int in_ligand,
			float centre_x,float centre_y,float scale_factor)
{
  char atom_colour[COL_NAME_LEN + 1];

  int inligand, wanted;

  float ps_coord_ax, ps_coord_ay;

  /* Initialise flag */
  wanted = TRUE;
  inligand = atom_ptr->residue_ptr->inligand;

  /* Check whether this is a main-chain atom */
  if (atom_ptr->side_chain == FALSE)
    {
      /* If plotting the schematic plot, then don't want to plot
	 the mainchain atoms (unless carbonyl oxygen) */
      if (Include->Simple_Ligand_Residues == TRUE &&
	  inligand == TRUE &&
	  (strncmp(atom_ptr->atom_type," O  ",4) &&
	   strncmp(atom_ptr->atom_type," OXT",4)))
	wanted = FALSE;

      /* Don't plot atom if XSITE plot and occupancy is zero */
      if (xsite_file == TRUE && strncmp(xsite_probe,"  0.00",6) &&
	  !strncmp(atom_ptr->occupancy,"  0.00",6))
	wanted = FALSE;
    }

  /* Check whether atom previously set as non-plottable */
  if (atom_ptr->plot_atom == FALSE)
    wanted = FALSE;

  /* If atom is to be printed, proceed */
  if (wanted == TRUE)
    {
      /* Calculate PostScript coordinates */
      ps_coord_ax = scale_factor * (atom_ptr->x - centre_x)
	+ Plot_Centre_x;
      ps_coord_ay = scale_factor * (atom_ptr->y - centre_y)
	+ Plot_Centre_y;

      /* Get this atom's colour and set sphere colour accordingly */
      pscolb_(Colour->Atom_Edges);
      get_atom_colour(atom_ptr->atom_type,atom_ptr->residue_ptr->res_name,
		      atom_colour);
      psphcl_(atom_colour);

      /* If this is an XSITE plot and occupancy does not
	 correspond to the probe number read in at the start,
	 then plot reduced sphere size */
      if (xsite_file == TRUE)
	pspher_(ps_coord_ax,ps_coord_ay,Size->Nonligand_Atoms);
	    
      /* If water molecule, set its size */
      else if (!strncmp(atom_ptr->residue_ptr->res_name,"HOH",3))
	pspher_(ps_coord_ax,ps_coord_ay,Size->Waters);

      /* Otherwise, plot the standard sphere size according to 
	 whether this is a ligand or non-ligand atom */
      else
	{
	  if (in_ligand == TRUE)
	    pspher_(ps_coord_ax,ps_coord_ay,Size->Ligand_Atoms);
	  else
	    pspher_(ps_coord_ax,ps_coord_ay,Size->Nonligand_Atoms);
	}
    }
}
/***********************************************************************

get_atom_text_colour  -  Get the colour of the text for the name of the
                         current atom

***********************************************************************/

void get_atom_text_colour(int in_ligand,char res_name[4],
			  char text_colour[COL_NAME_LEN + 1])
{
  /* Determine the atom-type and its text colour */

  /* Water */
  if (!strncmp(res_name,"HOH",3))
    strcpy(text_colour,Text_Colour->Water_Names);

  /* Ligand atom */
  else if (in_ligand == TRUE)
      strcpy(text_colour,Text_Colour->Ligand_Atom_Names);

  /* Ligand atom */
  else if (in_ligand == FALSE)
      strcpy(text_colour,Text_Colour->Nonligand_Atom_Names);
}
/***********************************************************************

print_atom_name  -  Print the current atom's name

***********************************************************************/

void print_atom_name(struct coordinate *atom_ptr,int in_ligand,
		     float centre_x,float centre_y,float scale_factor)
{
  char text_size[COL_NAME_LEN + 1], text_colour[COL_NAME_LEN + 1];

  int wanted;

  float atom_radius, x_pos, y_pos, text_height, text_x, text_y;
  float ps_coord_x, ps_coord_y;

  /* Get the current atom's radius and text-size */
  atom_radius = atom_ptr->atom_size / scale_factor;
  if (in_ligand == TRUE)
    {
      strcpy(text_size,Text_Size->Ligand_Atom_Names);
      text_height = (Text_Size_Val->Ligand_Atom_Names / scale_factor) / 3.0;
    }
  else
    {
      strcpy(text_size,Text_Size->Nonligand_Atom_Names);
      text_height = (Text_Size_Val->Nonligand_Atom_Names / scale_factor) / 3.0;
    }

  /* Get the atom's radius */
  atom_radius = atom_ptr->atom_size;

  /* Set default position for atom label diagonally above and to the
     right of the atom */
  x_pos = y_pos = 1.0;

  /* Determine whether this atom is required for printing */
  wanted = TRUE;
  if (atom_ptr->print_name[0] == '\0')
    wanted = FALSE;

  /* If atom-name is required, then proceed to print it */
  if (wanted == TRUE)
    {
      /* Find the best position for the atom name */
      position_names(atom_ptr,&x_pos,&y_pos);
      text_x = (atom_radius + text_height) * x_pos;
      text_y = (atom_radius + text_height) * y_pos;

      /* Calculate the PostScript coordinates for name */
      ps_coord_x = scale_factor * (atom_ptr->x - centre_x + text_x)
	+ Plot_Centre_x;
      ps_coord_y = scale_factor * (atom_ptr->y - centre_y + text_y)
	+ Plot_Centre_y;

      /* Print the atom name */
      get_atom_text_colour(in_ligand,atom_ptr->residue_ptr->res_name,
			   text_colour);
      pscolb_(text_colour);
      psctxt_(ps_coord_x,ps_coord_y,text_size,atom_ptr->print_name);
    }
}
/***********************************************************************

plot_atoms  -  Plot all the atoms in the picture

***********************************************************************/

void plot_atoms(float centre_x,float centre_y,float scale_factor)
{
  int iatom, iobject, iresid, natoms, nbonds, nresid;
  int deleted, in_ligand;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr[2];
  struct object *object_ptr;
  struct object_bond *object_bond_ptr;
  struct residue *residue_ptr;
    
  /* Initialise */
  pscomm_("Atoms");
  pssave_();
  pslwid_("Default_linewidth");

  /* Initialise all the atom flags */
  initialise_atoms();

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* If object is not a hydrophobic group, and not a simplified
	 H-group, then plot all its atoms */
      if (object_ptr->object_type != HYDROPHOBIC &&
	  object_ptr->object_type != SIMPLE_HGROUP)
	{
	  /* Get pointer to first residue to determine whether this
	     object is a ligand or non-ligand one */
	  residue_ptr = object_ptr->first_residue_ptr;
	  in_ligand = residue_ptr->inligand;

	  /* If object is a water, then plot its single atom */
	  if (object_ptr->object_type == WATER)
	    {
	      /* Get its atom pointer */
	      atom_ptr[0] = residue_ptr->first_atom_ptr;

	      /* Plot the water atom */
/* v.4.0--> */
/*	      if (Include->Nonligand_Atoms == TRUE) */
	      if (Include->Water_Atoms == TRUE)
/* <--v.4.0 */
		print_current_atom(atom_ptr[0],FALSE,
				   centre_x,centre_y,scale_factor);
	    }

	  /* Otherwise, loop through the object's bonds and plot all
	     atoms on either end */
	  else
	    {
	      /* Initialise count of bonds encountered */
	      nbonds = 0;

	      /* Get pointer to the first of this object's bonds */
	      object_bond_ptr = object_ptr->first_object_bond_ptr;
      
	      /* Loop through all object's bonds */
	      while (object_bond_ptr != NULL)
		{
		  /* Get the bond pointer */
		  bond_ptr = object_bond_ptr->bond_ptr;

		  /* Get pointers to the bond's two atoms */
		  atom_ptr[0] = bond_ptr->first_atom_ptr;
		  atom_ptr[1] = bond_ptr->second_atom_ptr;

		  /* Check that neither atom has been deleted */
		  deleted = FALSE;
		  if (atom_ptr[0]->deleted == TRUE ||
		      atom_ptr[1]->deleted == TRUE)
		    deleted = TRUE;
		  else
		    nbonds++;

		  /* Loop to plot the bond's two atoms */
		  for (iatom = 0; iatom < 2 && deleted == FALSE; iatom++)
		    {
		      /* Check that this atom hasn't already been plotted */
		      if (atom_ptr[iatom]->checked == FALSE)
			{
			  /* Plot the current atom */
			  if ((object_ptr->object_type == LIGAND &&
			       Include->Ligand_Atoms == TRUE) ||
/* v.4.0--> */
			      (object_ptr->object_type == WATER &&
			       Include->Water_Atoms == TRUE) ||
/* <--v.4.0 */
/* v.4.0.1--> */
/*			      Include->Nonligand_Atoms == TRUE) */
			      (object_ptr->object_type != LIGAND &&
			       object_ptr->object_type != WATER &&
			       Include->Nonligand_Atoms == TRUE))
/* <--v.4.0.1 */
			    print_current_atom(atom_ptr[iatom],in_ligand,
					       centre_x,centre_y,
					       scale_factor);

			  /* If atom name is required, print it by the atom */
			  if (Include->Atom_Names == TRUE)
			    print_atom_name(atom_ptr[iatom],in_ligand,
					    centre_x,centre_y,scale_factor);

			  /* Update flag to indicate atom has already been
			     plotted */
			  atom_ptr[iatom]->checked = TRUE;
			}
		    }
      
		  /* Get pointer to this object's next bond entry */
		  object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
		}

	      /* If the object contains no bonds, loop through all its
		 atoms and print them one by one */
	      if (nbonds == 0)
		{
		  /* Get pointer to the first of this object's residues */
		  residue_ptr = object_ptr->first_residue_ptr;
		  nresid = object_ptr->nresidues;
		  iresid = 0;

		  /* Loop over all this object's residues */
		  while (iresid < nresid && residue_ptr != NULL)
		    {
		      /* Get pointer to the first of this residue's atoms */
		      atom_ptr[0] = residue_ptr->first_atom_ptr;
		      iatom = 0;
		      natoms = residue_ptr->natoms;

		      /* Loop over all the residue's atoms to plot */
		      while (iatom < natoms && atom_ptr != NULL)
			{
			  /* If atom not deleted, then plot */
			  if (atom_ptr[0]->deleted == FALSE)
			    {
			      /* Plot the current atom */
			      if ((object_ptr->object_type == LIGAND &&
				   Include->Ligand_Atoms == TRUE) ||
/* v.4.0--> */
				  (object_ptr->object_type == WATER &&
				   Include->Water_Atoms == TRUE) ||
/* <--v.4.0 */
				  Include->Nonligand_Atoms == TRUE)
				print_current_atom(atom_ptr[0],in_ligand,
						   centre_x,centre_y,
						   scale_factor);

			      /* If atom name is required, print it by
				 the atom */
			      if (Include->Atom_Names == TRUE)
				print_atom_name(atom_ptr[0],in_ligand,
						centre_x,centre_y,
						scale_factor);
			    }

			  /* Get pointer to the next atom */
			  atom_ptr[0] = atom_ptr[0]->next;
			  iatom++;
			}

		      /* Get pointer to the next residue in the list */
		      residue_ptr = residue_ptr->next_residue_ptr;
		      iresid++;
		    }
		}
	    }
	}

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

    /* Reset graphics state */
    psrest_();
}
/***********************************************************************

convert_residue_name  -  Convert residue name and number as it is to
                         appear on the plot

***********************************************************************/

void convert_residue_name(char res_name[4],char res_num[6],char chain,
			  char res_name_string[14],int *n_chars)
{
  int i, j, found, number;
  char converted_res_name[4];
  char upper[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  char lower[] = "abcdefghijklmnopqrstuvwxyz";
  char digit[] = "0123456789";

  /* Initialise variables */
  *n_chars = 0;

  /* Check whether any digit is a number */
  number = FALSE;
/* v.4.0--> */
  for (i = 0; i < 3; i++)
    {
/* <--v.4.0 */
      for (j = 0; j < 10 && number == FALSE; j++)
	{
/* v.4.0--> */
/*	  if (res_name[2] == digit[j]) */
	  if (res_name[i] == digit[j])
/* <--v.4.0 */
	    number = TRUE;
	}
/* v.4.0--> */
    }
/* <--v.4.0 */

  /* If residue name is a water, or last character is a number,
     or first character is a space, then keep in upper case */
/* v.4.0--> */
/*  if (!strncmp(res_name,"HOH",3) || !strncmp(res_name,"WAT",3) ||
      number == TRUE) */
  if (!strncmp(res_name,"HOH",3) || !strncmp(res_name,"WAT",3) ||
      res_name[0] == ' ' || number == TRUE)
/* <--v.4.0 */
    {
      strncpy(converted_res_name,res_name,3);
      converted_res_name[3] = '\0';
    }
  
  /* Otherwise, convert 2nd and 3rd letters of residue-name to lower-case */
  else
    {
      converted_res_name[0] = res_name[0];
      for (i = 1; i < 3; i++)
	{
	  converted_res_name[i] = res_name[i];
	  found = FALSE;
	  for (j = 0; j < 26 && found == FALSE; j++)
	    {
	      if (res_name[i] == upper[j])
		{
		  found = TRUE;
		  converted_res_name[i] = lower[j];
		}
	    }
	}
      converted_res_name[3] = '\0';
    }

  /* Transfer new residue name into string holding final residue
     details */
  strncpy(res_name_string,converted_res_name,3);
  res_name_string[3] = ' ';
  *n_chars = 4;

  /* Add the residue number */
  if (res_num[2] == ' ')
    {
      strncpy(res_name_string+4,res_num+3,2);
      *n_chars = 6;
    }
  else if (res_num[1] == ' ')
    {
      strncpy(res_name_string+4,res_num+2,3);
      *n_chars = 7;
    }
  else if (res_num[0] == ' ')
    {
      strncpy(res_name_string+4,res_num+1,4);
      *n_chars = 8;
    }
  else
    {
      strncpy(res_name_string+4,res_num,5);
      *n_chars = 9;
    }
  if (res_num[4] == ' ')
    (*n_chars)--;

  /* Add the chain-id */
  if (chain == ' ')
    res_name_string[*n_chars]  = '\0';
  else
    {
      res_name_string[(*n_chars)] = '(';
      res_name_string[(*n_chars) + 1] = chain;
      res_name_string[(*n_chars) + 2] = ')';
      res_name_string[(*n_chars) + 3]  = '\0';
      *n_chars = (*n_chars) + 3;
    }
}
/***********************************************************************

print_name  -  Print residue name

***********************************************************************/

void print_name(char res_name[4],char res_num[6],char chain,
		float x,float y,float centre_x,float centre_y,
		float text_size,char text_size_string[PS_STRING_LENGTH],
		char text_colour[PS_STRING_LENGTH],float scale_factor,
		int blank_out,int border,int shift_x,int shift_y,
		int stack)
{
  char res_name_string[14], rname[4];

  int n_chars;

  float ps_coord_x, ps_coord_y, ps_y, sizex, sizey;

  /* Calculate PostScript coordinates */
  ps_coord_x = scale_factor * (x - centre_x) + Plot_Centre_x;
  ps_coord_y = scale_factor * (y - centre_y) + Plot_Centre_y;

  /* Prepare string containing residue name, number and chain-id */
  convert_residue_name(res_name,res_num,chain,res_name_string,&n_chars);

  /* Get size of box just bounding residue-name label */
  sizex = 10 * text_size * CHAR_ASPECT / 2.0;
  sizey = 0.6 * text_size;

  /* If shift is required, then apply it */
  ps_coord_x = ps_coord_x - shift_x * sizex;
  ps_coord_y = ps_coord_y - shift_y * sizey;

  /* Set the text colour */
  pscolb_(text_colour);
	      
  /* Blank out region for residue name, if required */
  if (blank_out == TRUE)
    {
      pslwid_("Zero_linewidth");
      pshade_(Colour->Background);

      /* Plot a bounded or unbounded box, as required */
      if (border == TRUE)
	psbbox_(ps_coord_x - sizex,ps_coord_y + sizey,
		ps_coord_x - sizex,ps_coord_y - sizey,
		ps_coord_x + sizex,ps_coord_y - sizey,
		ps_coord_x + sizex,ps_coord_y + sizey);
      else
	psubox_(ps_coord_x - sizex,ps_coord_y + sizey,
		ps_coord_x - sizex,ps_coord_y - sizey,
		ps_coord_x + sizex,ps_coord_y - sizey,
		ps_coord_x + sizex,ps_coord_y + sizey);
    }

  /* Print the residue name, either with name stacked above number,
     or in one long string */
  if (stack == TRUE)
    {
      /* Extract just the residue name */
      strncpy(rname,res_name_string,3);
      rname[3] = '\0';

      /* Print the residue name above the number and chain-id */
      ps_y = ps_coord_y + text_size / 2.0;
      psctxt_(ps_coord_x,ps_y,text_size_string,rname);
      ps_y = ps_coord_y - text_size / 2.0;
      psctxt_(ps_coord_x,ps_y,text_size_string,res_name_string + 4);
    }
  else
    psctxt_(ps_coord_x,ps_coord_y,text_size_string,res_name_string);
}
/***********************************************************************

plot_hydrophobics  -  Plot all the symbols representing hydrophobic
                      interactions

***********************************************************************/

void plot_hydrophobics(float centre_x,float centre_y,float scale_factor)
{
/* v.4.0--> */
  char text_colour[PS_STRING_LENGTH];
/* <--v.4.0 */

  int iatom, iobject, jatom, spokes_only;
/* v.4.0--> */
  int interface;
/* <--v.4.0 */

  float radius, xc, x[2], yc, y[2];

  struct bond *bond_ptr;
  struct coordinate *atom_ptr[2];
  struct object *object_ptr;
  struct residue *residue_ptr;
    
  /* Initialise */
  pscomm_("Hydrophobic interactions");
  pssave_();
  pslwid_("Default_linewidth");
  printf("   Plotting hydrophobic contacts ...\n");

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* If object is a hydrophobic group, then print residue name */
      if (object_ptr->object_type == HYDROPHOBIC)
	{
	  /* Get the coordinates of the object */
	  xc = object_ptr->minx;
	  yc = object_ptr->miny;

	  /* Get pointer to residue so that can pick up residue
	     details */
	  residue_ptr = object_ptr->first_residue_ptr;

/* v.4.0--> */
	  /* Get the appropriate colour for the text of the residue name */
	  if (Interface_Plot == FALSE)
	    strcpy(text_colour,Text_Colour->Hydrophobic_Names);
	  else if (object_ptr->interface == 1)
	    strcpy(text_colour,Text_Colour->Ligand_Residue_Names);
	  else
	    strcpy(text_colour,Text_Colour->Nonligand_Residue_Names);
/* <--v.4.0 */

	  /* Plot the residue name at this point */
	  pscomm_("Hydrophobic contact");
	  print_name(residue_ptr->res_name,residue_ptr->res_num,
		     residue_ptr->chain,xc,yc,centre_x,centre_y,
		     Text_Size_Val->Hydrophobic_Names,
		     Text_Size->Hydrophobic_Names,
/* v.4.0--> */
/*		     Text_Colour->Hydrophobic_Names,scale_factor, */
		     text_colour,scale_factor,
/* <--v.4.0 */
		     TRUE,FALSE,0,0,FALSE);
	}

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

  /* Initialise pointer to the first stored bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds to pick out those representing
     hydrophobic contacts */
  while (bond_ptr != NULL)
    {
      /* If this represents a non-bonded hydrophobic contact, then process */
      if (bond_ptr->bond_type == CONTACT)
	{
	  /* Get the two atoms involved in the contact */
	  atom_ptr[0] = bond_ptr->first_atom_ptr;
	  atom_ptr[1] = bond_ptr->second_atom_ptr;

	  /* Get the atom coordinates */
	  x[0] = atom_ptr[0]->x;
	  y[0] = atom_ptr[0]->y;
	  x[1] = atom_ptr[1]->x;
	  y[1] = atom_ptr[1]->y;

	  /* Loop over the two atoms to plot their interactions */
	  for (iatom = 0; iatom < 2; iatom++)
	    {
	      /* Plot atom if not deleted */
	      if (atom_ptr[iatom]->deleted == FALSE)
		{
		  /* Get index of the other atom */
		  jatom = 1 - iatom;

		  /* Determine whether this atom belongs to a hydrophobic
		     group or to a ligand/H-bond group */
		  residue_ptr = atom_ptr[iatom]->residue_ptr;
		  object_ptr = residue_ptr->object_ptr;

/* v.4.0--> */
		  /* Initialise interface side */
		  interface = 0;
/* <--v.4.0 */

		  /* Determine radius of arc from which spokes will emanate */
		  spokes_only = TRUE;
		  if (object_ptr->object_type == HYDROPHOBIC)
		    {
		      radius = Size_Val->Hydrophobics;
		      spokes_only = FALSE;

/* v.4.0--> */
		      /* For interface plot, determine which surface this
			 hydrophobic group belongs to */
		      if (Interface_Plot == TRUE)
			interface = object_ptr->interface;
/* <--v.4.0 */
		    }
		  else
		    radius = scale_factor * atom_ptr[iatom]->atom_size;

/* Old ***	  else if (residue_ptr->residue_type == SIMPLE_LIGAND &&
			   !strncmp(atom_ptr[iatom]->atom_type," CA ",4))
		    radius = Size_Val->Simple_Residues;
		  else if (residue_ptr->inligand == TRUE)
		    radius = Size_Val->Ligand_Atoms;
		  else
		    radius = Size_Val->Nonligand_Atoms; */

		  /* Plot the spokes */
		  plot_contact_group(x[iatom],y[iatom],x[jatom],y[jatom],
				     centre_x,centre_y,scale_factor,
/* v.4.0--> */
/*				     radius,120.0,spokes_only); */
				     radius,120.0,spokes_only,interface);
/* <--v.4.0 */
		}
	    }
	}

      /* Go to the next bond in the linked list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Reset graphics state */
  psrest_();
}
/***********************************************************************

plot_simple_ligand_residues  -  Plot all the residues in the simplified
                                representation

***********************************************************************/

void plot_simple_ligand_residues(float centre_x,float centre_y,
				 float scale_factor)
{
  int iatom, iobject, iresid, natoms, nresid;

  float ps_coord_x, ps_coord_y, xc, yc;

  struct coordinate *atom_ptr, *ca_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;
  char radius[COL_NAME_LEN + 1];
    
  /* Initialise */
  pscomm_("Simplified ligand residues");
  pssave_();
  pslwid_("Default_linewidth");
  printf("   Plotting ligand residues ...\n");

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* If object is a ligand, then process its residues */
      if (object_ptr->object_type == LIGAND)
	{
	  /* Get pointer to the first of this object's residues */
	  residue_ptr = object_ptr->first_residue_ptr;
	  nresid = object_ptr->nresidues;
	  iresid = 0;

	  /* Loop over all this object's residues */
	  while (iresid < nresid && residue_ptr != NULL)
	    {
	      /* If this residue is to be plotted in the simplified
		 representation, then do so */
	      if (residue_ptr->residue_type == SIMPLE_LIGAND)
		{
		  /* Get pointer to the first of this residue's atoms */
		  atom_ptr = residue_ptr->first_atom_ptr;
		  iatom = 0;
		  natoms = residue_ptr->natoms;
		  ca_ptr = NULL;
		  
		  /* Loop over all the residue's atoms to locate
		     the CA */
		  while (iatom < natoms && atom_ptr != NULL &&
			 ca_ptr == NULL)
		    {
		      /* If this is the CA, then save its pointer */
		      if (!strncmp(atom_ptr->atom_type," CA ",4))
			ca_ptr = atom_ptr;

		      /* Get pointer to the next atom */
		      atom_ptr = atom_ptr->next;
		      iatom++;
		    }

		  /* If have found the CA atom, then plot the residue
		     name at that position */
		  if (ca_ptr != NULL)
		    {
		      /* Get the coordinates of the CA */
		      xc = ca_ptr->x;
		      yc = ca_ptr->y;
		      ps_coord_x = scale_factor * (xc - centre_x)
			+ Plot_Centre_x;
		      ps_coord_y = scale_factor * (yc - centre_y)
			+ Plot_Centre_y;

		      /* Plot the circle representing the residue */
		      strcpy(radius,Size->Simple_Residues);
		      psccol_(Colour->Background);
		      pscolb_(Colour->Simple_Residues);
		      pslwid_("Simple_width");
		      pscirc_(ps_coord_x,ps_coord_y,radius);

		      /* Plot the residue name at this point */
		      pscomm_("Ligand residue");
		      print_name(residue_ptr->res_name,residue_ptr->res_num,
				 residue_ptr->chain,xc,yc,centre_x,centre_y,
				 Text_Size_Val->Simple_Residue_Names,
				 Text_Size->Simple_Residue_Names,
				 Text_Colour->Ligand_Residue_Names,
				 scale_factor,FALSE,FALSE,0,0,TRUE);
		    }
		}

	      /* Get pointer to the next residue in the list */
	      residue_ptr = residue_ptr->next_residue_ptr;
	      iresid++;
	    }
	}

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

  /* Reset graphics state */
  psrest_();
}
/***********************************************************************

plot_simple_hgroups  -  Plot all the symbols representing the simplified
                       H-bonded residues

***********************************************************************/

void plot_simple_hgroups(float centre_x,float centre_y,float scale_factor)
{
  int iobject;
  int shift_x, shift_y;
  int first_time;
  int diff, iseg, jseg, max_score, segment[8], seg_score[8], sdist;

  float xc, yc;
  float x1, x2, y1, y2;

  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;

  static int score[] = { 10, 9, 7, 4, 0 };
    
  /* Initialise */
  pscomm_("Simplified H-groups");
  pssave_();
  pslwid_(Size->Nonligand_Bonds);
  printf("   Plotting simplified H-groups ...\n");

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* If object is a hydrophobic group, then plot it */
      if (object_ptr->object_type == SIMPLE_HGROUP)
	{
	  /* Get the coordinates of the object */
	  xc = object_ptr->minx;
	  yc = object_ptr->miny;

	  /* Get pointer to residue so that can pick up residue
	     details */
	  residue_ptr = object_ptr->first_residue_ptr;

	  /* Determine best position of residue label so that it
	     doesn't interfere with the bonds coming from it */
	  shift_x = 1;
	  shift_y = 0;
	  first_time = TRUE;

	  /* Get pointer to first bond */
	  bond_ptr = first_bond_ptr;

	  /* Loop through all the hydrogen bonds involving this H-group
	     and determine the direction in which each bond lies */
	  while (bond_ptr != NULL)
	    {
	      /* If this is a hydrogen bond, then get its two atoms */
	      if (bond_ptr->bond_type == HBOND)
		{
		  /* Get the bond's two atoms */
		  atom1_ptr = bond_ptr->first_atom_ptr;
		  atom2_ptr = bond_ptr->second_atom_ptr;

		  /* If either atom involves the current object, then
		     get direction of the bond */
		  if (atom1_ptr->residue_ptr->object_ptr == object_ptr ||
		      atom2_ptr->residue_ptr->object_ptr == object_ptr)
		    {
		      /* Get the two atoms' coordinates */
		      if (atom1_ptr->residue_ptr->object_ptr == object_ptr)
			{
			  x1 = atom1_ptr->x;
			  y1 = atom1_ptr->y;
			  x2 = atom2_ptr->x;
			  y2 = atom2_ptr->y;
			}
		      else
			{
			  x1 = atom2_ptr->x;
			  y1 = atom2_ptr->y;
			  x2 = atom1_ptr->x;
			  y2 = atom1_ptr->y;
			}

		      /* Determine which segment the current bond lies in */
		      get_segment(x1,y1,x2,y2,segment,first_time);
		      first_time = FALSE;
		    }
		}

	      /* Go to the next bond in the linked list */
	      bond_ptr = bond_ptr->next_bond_ptr;
	    }

	  /* Initialise the segment scores */
	  for (iseg = 0; iseg < 8; iseg++)
	    seg_score[iseg] = 0;

	  /* Loop through to see which segments are free and which are
	     taken, and accumulate the scores accordingly */
	  for (iseg = 0; iseg < 8; iseg++)
	    {
	      /* If this segment is taken, score the nearest segments
		 most highly */
	      if (segment[iseg] == NOT_FREE)
		{
		  /* Loop through all segment positions, scoring
		     each one according to its proximity to the
		     current segment */
		  for (jseg = 0; jseg < 8; jseg++)
		    {
		      /* Get the distance between the segments and
			 score accordingly */
		      diff = jseg - iseg;
		      if (diff > 4)
			sdist = 8 - diff;
		      else if (diff < -4)
			sdist = 8 + diff;
		      else if (diff < 0)
			sdist = - diff;
		      else
			sdist = diff;

		      /* Compute the score */
		      seg_score[jseg] = seg_score[jseg] + score[sdist];
		    }
		}
	    }

	  /* Determine which is the maximum-scoring position */
	  jseg = 0;
	  max_score = 0;
	  for (iseg = 0; iseg < 8; iseg++)
	    {
	      if (seg_score[iseg] > max_score)
		{
		  max_score = seg_score[iseg];
		  jseg = iseg;
		}
	    }

	  /* Convert the optimal segment position into x- and y-shifts */
	  if (jseg == 0)
	    {
	      shift_x = 1;
	      shift_y = 1;
	    }
	  else if (jseg == 1)
	    {
	      shift_x =  0;
	      shift_y =  1;
	    }
	  else if (jseg == 2)
	    {
	      shift_x = -1;
	      shift_y =  1;
	    }
	  else if (jseg == 3)
	    {
	      shift_x = -1;
	      shift_y =  0;
	    }
	  else if (jseg == 4)
	    {
	      shift_x = -1;
	      shift_y = -1;
	    }
	  else if (jseg == 5)
	    {
	      shift_x =  0;
	      shift_y = -1;
	    }
	  else if (jseg == 6)
	    {
	      shift_x =  1;
	      shift_y = -1;
	    }
	  else if (jseg == 7)
	    {
	      shift_x =  1;
	      shift_y =  0;
	    }
	  
	  /* Print the residue name, with a border around it */
	  pscomm_("H-group");
	  print_name(residue_ptr->res_name,residue_ptr->res_num,
		     residue_ptr->chain,xc,yc,centre_x,centre_y,
		     Text_Size_Val->Nonligand_Residue_Names,
		     Text_Size->Nonligand_Residue_Names,
		     Text_Colour->Nonligand_Residue_Names,scale_factor,
		     TRUE,TRUE,shift_x,shift_y,FALSE);
	}

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

  /* Reset graphics state */
  psrest_();
}
/***********************************************************************

extract_search_atoms  -  Extract all the atoms that are within range
                         of the given range-limits

***********************************************************************/

void extract_search_atoms(struct residue *current_residue_ptr,
			  struct coordinate **first_stack_ptr,
			  float search_min_x,float search_max_x,
			  float search_min_y,float search_max_y,
			  float scale_factor)
{
  int current, iatom, natoms, wanted;
/* v.4.0--> */
  int plotted_accessibility;
/* <--v.4.0 */

  float acc_size, atom_size, x, y;

  struct coordinate *atom_ptr, *last_atom_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  acc_size = 4.0 * Size_Val->Ligand_Atoms / scale_factor;

  /* Initialise all the atom flags */
  initialise_atoms();

  /* Initialise stack pointers */
  last_atom_ptr = NULL;
  *first_stack_ptr = NULL;

  /* Initialise pointer to the first stored residue */
  residue_ptr = first_residue_ptr;

  /* Loop through all residues */
  while (residue_ptr != NULL)
    {
      /* Determine whether this residue is in range of the search region */
      wanted = TRUE;
      if (residue_ptr->minx > search_max_x ||
	  residue_ptr->maxx < search_min_x ||
	  residue_ptr->miny > search_max_y ||
	  residue_ptr->maxy < search_min_y)
	wanted = FALSE;

      /* If atoms from this residue likely to be wanted, then loop through
	 them */
      if (wanted == TRUE)
	{
	  /* If this is the same residue as the current on, then set
	     flag */
	  current = FALSE;
	  if (residue_ptr == current_residue_ptr)
	    current = TRUE;

	  /* Get pointer to this residue's first atom */
	  atom_ptr = residue_ptr->first_atom_ptr;
	  iatom = 0;
	  natoms = residue_ptr->natoms;

	  /* Loop over all the residue's atoms, performing the
	     transformation */
	  while (iatom < natoms && atom_ptr != NULL)
	    {
	      /* Get this atom's coordinates and its size */
	      x = atom_ptr->x;
	      y = atom_ptr->y;
	      atom_size = atom_ptr->atom_size;

/* v.4.0--> */
	      /* Determine whether accessibility will have been plotted
		 for this atom */
	      plotted_accessibility = FALSE;
	      if (Include->Accessibilities == TRUE)
		{
		  if (residue_ptr->inligand == TRUE ||
		      Include->Ligand_Accessibilities_Only == FALSE)
		    plotted_accessibility = TRUE;
		}
/* <--v.4.0 */

	      /* If have plotted accessibility shading, need to increase
		 atom's size (if not already done) */
/* v.4.0--> */
/*	      if (Include->Accessibilities == TRUE) */
	      if (plotted_accessibility == TRUE)
/* <--v.4.0 */
		{
		  if (atom_size < (acc_size))
		    {
		      atom_size = acc_size;
		      atom_ptr->atom_size = atom_size;
		    }
		}

	      /* Determine whether this atom is within range */
	      wanted = TRUE;
	      if ((x - atom_size > search_max_x) ||
		  (x + atom_size < search_min_x) ||
		  (y - atom_size > search_max_y) ||
		  (y + atom_size < search_min_y))
		wanted = FALSE;

	      /* If atom is within range, then add it to the stack */
	      if (wanted == TRUE)
		{
		  /* If this is the first atom, get the first pointer
		     to point to the second atom */
		  if (*first_stack_ptr == NULL)
		    *first_stack_ptr = atom_ptr;

		  /* Otherwise, set the last atom's pointer to the
		     other atom */
		  else
		    last_atom_ptr->next_stack_ptr = atom_ptr;

		  /* Save pointer to the current atom */
		  last_atom_ptr = atom_ptr;

		  /* If this atom belongs to the current residue,
		     flag it as checked */
		  if (current == TRUE)
		    atom_ptr->checked = TRUE;
		}

	      /* Get pointer to the next atom */
	      atom_ptr = atom_ptr->next;
	      iatom++;
	    }
	}

      /* Get pointer to the next residue in the list */
      residue_ptr = residue_ptr->next_residue_ptr;
    }
}
/***********************************************************************

label_dist  -  Calculate (squared) distance of point from residue-label
               position

***********************************************************************/

float label_dist(float x, float y, float label_coordx, float label_coordy,
		 float zone_x, float zone_y)
{
  int i, j;
  float corner_x, corner_y, distance_sqrd, dist_min, dist_x, dist_y;

  /* Check if the atom is overlapping any part of the label area */
  if (x < label_coordx + zone_x && x > label_coordx - zone_x && 
      y < label_coordy + zone_y && y > label_coordy - zone_y)
      distance_sqrd = 0.0;

  /* If the atom's x-coordinate falls between the
     corner x-coordinates, then use the straight-line
     y-distance */
  else if (x < label_coordx + zone_x && x > label_coordx - zone_x)
    {
      distance_sqrd = fabs(y - (label_coordy + zone_y));
      dist_y = fabs(y - (label_coordy - zone_y));
      if (dist_y < distance_sqrd)
	distance_sqrd = dist_y;
      distance_sqrd = distance_sqrd * distance_sqrd;
    }
  
  /* If the atom's y-coordinate falls between the
     corner y-coordinates, then use the straight-line
     x-distance */
  else if (y < label_coordy + zone_y && y > label_coordy - zone_y)
    {
      distance_sqrd = fabs(x - (label_coordx + zone_x));
      dist_x = fabs(x - (label_coordx - zone_x));
      if (dist_x < distance_sqrd)
	distance_sqrd = dist_x;
      distance_sqrd = distance_sqrd * distance_sqrd;
    }
  
  /* Otherwise, calculate the distance of the atom
     from the nearest corner */
  else
    {
      dist_min = 100000.0;

      /* Loop across the four corners */
      for (i = 0; i < 2; i++)
	{
	  corner_x =  label_coordx + (i * 2 - 1) * zone_x;
	  for (j = 0; j < 2; j++)
	    {
	      corner_y =  label_coordy + (j * 2 - 1) * zone_y;
	      
	      /* Calculate the distance of the atom from
		 this corner */
	      dist_x = x - corner_x;
	      dist_y = y - corner_y;
	      distance_sqrd = dist_x * dist_x + dist_y * dist_y;
	      
	      /* Save this distance if it is the minimum
		 so far */
	      if (distance_sqrd < dist_min)
		dist_min = distance_sqrd;
	    }
	}
      distance_sqrd = dist_min;
    }
  
  /* Compute the effective energy from the closest approach */
  if (distance_sqrd < 0.00001)
    distance_sqrd = 0.00001;

  return(distance_sqrd);
}
/***********************************************************************

perform_search  -  Search for a suitable location for the placement of
                   the residue label, avoiding overlaps with atoms and
		   other labels as far as possible

***********************************************************************/

void perform_search(struct coordinate *first_stack_ptr,float grid_size,
		    float box_min_x,float box_max_x,
		    float box_min_y,float box_max_y,
		    float half_width,float half_height,
		    float *pos_x,float *pos_y)
{
  int ilabel;
  int npoint_x, npoint_y, point_x, point_y;
  int no_good;

  float x, y;
  float atom_size, length_x, length_y;
  float label_coordx, label_coordy;
  float energy, energy_diff, energy_store, label_energy;
  float lowest_energy, highest_energy, worst_energy, worst_other;
  float distance, distance_sqrd, closest_distance_sqrd;

  struct coordinate *atom_ptr, *closest_atom_ptr;

  /* Initialise variables */
  lowest_energy = 100000.0;
  highest_energy = 0.0;
  energy_store = 0.0;

  /* Work out the dimensions of the box in which the residue name can
     be placed */
  length_x = (box_max_x - box_min_x);
  length_y = (box_max_y - box_min_y);
  *pos_x = box_min_x;
  *pos_y = box_min_y;

  /* Work out number of grid-points required for trial placement of
     residue name */
  npoint_x = 2 + length_x / grid_size;
  npoint_y = 2 + length_y / grid_size;

  /* Put the residue name at each grid-point and calculate an "energy"
     for it, depending on how close it is to its residue and whether
     it clashes with any other atoms */

  /* Loop over the grid-points in the y-direction */
  for (point_y = 0; point_y < npoint_y; point_y++)
    {
      label_coordy = box_min_y + point_y * grid_size;

      /* Loop over the grid-points in the x-direction */
      for (point_x = 0; point_x < npoint_x; point_x++)
	{
	  label_coordx = box_min_x + point_x * grid_size; 
	  worst_energy = 0.0;
	  label_energy = 0.0;
	  worst_other = 0.0;
	  closest_atom_ptr = NULL;
	  closest_distance_sqrd = 1000000.0;
	  no_good = FALSE;

	  /* Get pointer to the first atom on the stack of stored atoms
	     within range */
	  atom_ptr = first_stack_ptr;

	  /* Loop through all the stored atoms in the vicinity of
	     the current residue */
	  while (atom_ptr != NULL && no_good == FALSE)
	    {
	      /* Get the coordinates and size of this atom */
	      x = atom_ptr->x;
	      y = atom_ptr->y;
	      atom_size = atom_ptr->atom_size;

	      /* Calculate the nearest distance of atom from
		 any part of the residue label */
	      distance_sqrd = label_dist(x,y,label_coordx,label_coordy,
					 half_width,half_height);

	      /* If atom overlaps the label, then discard this position */
	      if (distance_sqrd < (atom_size + 0.005) * (atom_size + 0.005))
		{
		  no_good = TRUE;
		  energy_store = 10000000.0;
		  distance_sqrd = 0.0;
		}

		/* Otherwise, compute the effective energy from the closest
		   approach */
		else
		  {
		    distance = sqrt((double)distance_sqrd) - atom_size;
		    distance_sqrd = distance * distance;
		    energy_store = (1.0 / distance_sqrd);
		  }

		/* If this is the closest atom so far, save it */
		if (distance_sqrd < closest_distance_sqrd)
		  {
		    closest_atom_ptr = atom_ptr;
		    closest_distance_sqrd = distance_sqrd;
		  }

		/* If the current atom belongs to the correct residue
		   save the energy only if it is the highest so far */
		if (atom_ptr->checked == TRUE)
		  {
		    /* Save the worst energy encountered so far */
		    energy = 10.0 * (IDEAL_ENERGY - energy_store);
		    if (energy < 0)
		      energy = 10000000.0;
		    worst_energy = worst_energy + energy;
		  }

		/* Otherwise, save the worst energy from close approaches
		   to other residues */
		else
		    worst_other = worst_other + energy_store;

	      /* Get pointer to the next atom in the stack */
	      atom_ptr = atom_ptr->next_stack_ptr;
	    }

	  /* Loop through all the labels written so far (ie H-bond
	     distances and residue labels done to date) */
	  for (ilabel = 0; ilabel < n_labels && no_good == FALSE; ilabel++)
	    {
	      /* Get the coordinates of this label */
	      x = label_coord[ilabel][0];
	      y = label_coord[ilabel][1];

	      /* Calculate the nearest distance of this point from
		 any part of the residue label */
	      distance = label_dist(x,y,label_coordx,label_coordy,
				    half_width,half_height);
	      if (distance < 0.0001)
		{
		  no_good = TRUE;
		  energy_store = 10000000.0;
		}

	      /* Compute the effective energy from the closest
		 aproach */
	      else
		energy_store = (1.0 / distance);

	      /* Save the worst label-energy encountered so far */
	      label_energy = label_energy + energy_store;
	    }

	  /* Check difference between the highest energy (ie coming from
	     the atom that's closest to the label) and the "ideal"
	     energy */
	  energy_diff = worst_energy + worst_other + 10.0 * label_energy;
	  if (energy_diff < 0.0)
	    energy_diff = 10000000.0;

	  /* If this is the very first loop, set current
	     energy as the lowest so far encountered */
	  if (point_x == 0 && point_y == 0)
	    lowest_energy = energy_diff + 1.0;

	  /* Now pick out the coordinates of the lowest energy */
	  if (energy_diff < lowest_energy)
	    {
	      *pos_x = label_coordx;
	      *pos_y = label_coordy;
	      lowest_energy = energy_diff;
	    }
	  else if (energy_diff < 1000000.0 &&
		   energy_diff > highest_energy)
	    highest_energy = energy_diff;
	}
    }

  /* Save the coords actually used */
  if (n_labels < MAXLABELS)
    {
      label_coord[n_labels][0] = *pos_x;
      label_coord[n_labels][1] = *pos_y;
      n_labels++;
    }
}
/***********************************************************************

resname_grid_search  -  Move the residue name so it doesn't overlap with
                        the atoms

***********************************************************************/

void resname_grid_search(struct residue *current_residue_ptr,
			 float text_size,int n_chars,
			 float *new_coordx,float *new_coordy,
			 float min_x,float min_y,
			 float max_x,float max_y,
			 float scale_factor)
{
  float grid_size;
  float box_min_x, box_max_x, box_min_y, box_max_y;
  float search_min_x, search_max_x, search_min_y, search_max_y;
  float x, y;
  float half_width, half_height, zone_x, zone_y;

  float upper_width;

  struct coordinate *first_stack_ptr, *last_atom_ptr;

  /* Initialise variables */
  grid_size = 0.3;
  first_stack_ptr = NULL;
  last_atom_ptr = NULL;

  /* Compute width of residue name label, given its size and number
     of characters */
  half_width = n_chars * text_size * CHAR_ASPECT / 2.0;
  half_height = text_size / 2.0;
  upper_width = 15 * text_size * CHAR_ASPECT / 2.0;

  /* Calculate the size of the exclusion zone for the residue label */
  zone_x = half_width + Max_object_size / scale_factor;
  zone_y = half_height + Max_object_size / scale_factor;

  /* Define the minimum and maximum coordinates of the box within which
     the residue name can be placed and of the search region for
     atoms */
  box_min_x = min_x - 2.0 * half_width;
  box_max_x = max_x + 2.0 * half_width;
  box_min_y = min_y - 4.0 * half_height;
  box_max_y = max_y + 4.0 * half_height;
  search_min_x = box_min_x - zone_x;
  search_max_x = box_max_x + zone_x;
  search_min_y = box_min_y - zone_y;
  search_max_y = box_max_y + zone_y;

  /* Initialise all the atom flags */
  initialise_atoms();

  /* Extract all the atoms that are within the search box */
  extract_search_atoms(current_residue_ptr,&first_stack_ptr,
		       search_min_x,search_max_x,search_min_y,
		       search_max_y,scale_factor);

  /* Perform the search for a suitable location for the residue name */
  perform_search(first_stack_ptr,grid_size,box_min_x,box_max_x,box_min_y,
		 box_max_y,half_width,half_height,&x,&y);

  /* Return the calculated coordinates for the residue label */
  *new_coordx = x;
  *new_coordy = y;
}
/***********************************************************************

print_residue_names  -  Print all the residue names for the ligand and
                        nonligand sidechains

***********************************************************************/

void print_residue_names(float centre_x,float centre_y,float scale_factor)
{
  char chain, res_name[4], res_num[6], res_name_string[14];
  char text_size[COL_NAME_LEN + 1];

  int iobject, iresid, n_chars, nresid;
  int object_type;

  float min_x, min_y, max_x, max_y, new_coordx, new_coordy;
  float ps_coord_x, ps_coord_y, text_size_val;
  
  struct object *object_ptr;
  struct residue *residue_ptr;

  /* Initialise print */
  pscomm_("Residue names");
  pssave_();

  /* Initialise pointer to the first stored object */
  printf("   Printing residue names ...\n");
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Get this object's type */
      object_type = object_ptr->object_type;

      /* Process only if not a hydrophobic group */
      if (object_type != HYDROPHOBIC && object_type != SIMPLE_HGROUP)
	{
	  /* Set default text colour and size according to the
	     object type */
	  if (object_type == LIGAND)
	    {
	      pscolb_(Text_Colour->Ligand_Residue_Names);
	      strcpy(text_size,Text_Size->Ligand_Residue_Names);
	      text_size_val = Text_Size_Val->Ligand_Residue_Names
		/ scale_factor;
	    }

	  /* Default text colour and size for waters */
	  else if (object_type == WATER)
	    {
/* v.4.0--> */
/*	      pscolb_(Text_Size->Water_Names); */
	      pscolb_(Text_Colour->Water_Names);
/* <--v.4.0 */
	      strcpy(text_size,Text_Size->Water_Names);
	      text_size_val = Text_Size_Val->Water_Names / scale_factor;
	    }

	  /* Default text colour and size for non-ligand bonds */
	  else
	    {
	      pscolb_(Text_Colour->Nonligand_Residue_Names);
	      strcpy(text_size,Text_Size->Nonligand_Residue_Names);
	      text_size_val = Text_Size_Val->Nonligand_Residue_Names
		/ scale_factor;
	    }

	  /* Get pointer to the first of this object's residues */
	  residue_ptr = object_ptr->first_residue_ptr;
	  nresid = object_ptr->nresidues;
	  iresid = 0;

	  /* Loop through all this object's residues, one by one */
	  while (iresid < nresid && residue_ptr != NULL)
	    {
	      /* Precess only if this is a standard residue (ie not a
		 schematic in the simplified ligand representation */
	      if (residue_ptr->residue_type == STANDARD)
		{
		  /* Get the residue name, number and chain ID */
		  strncpy(res_name,residue_ptr->res_name,3);
		  res_name[3] = '\0';
		  strncpy(res_num,residue_ptr->res_num,5);
		  res_num[5] = '\0';
		  chain = residue_ptr->chain;

		  /* Prepare string containing residue name, number
		     and chain-id */
		  convert_residue_name(res_name,res_num,chain,
				       res_name_string,&n_chars);

		  /* Get the minimum and maximum coordinates of the atoms
		     belonging to this residue */
		  min_x = residue_ptr->minx;
		  max_x = residue_ptr->maxx;
		  min_y = residue_ptr->miny;
		  max_y = residue_ptr->maxy;

		  /* Perform the grid-search for the optimal positioning
		     of the residue name */
		  resname_grid_search(residue_ptr,text_size_val,n_chars,
				      &new_coordx,&new_coordy,
				      min_x,min_y,max_x,max_y,
				      scale_factor);

		  /* Calculate PostScript coords for residue label and text
		     height */
		  ps_coord_x = scale_factor * (new_coordx - centre_x)
		    + Plot_Centre_x;
		  ps_coord_y = scale_factor * (new_coordy - centre_y)
		    + Plot_Centre_y;

		  /* Print out string containing residue name and number */
		  psctxt_(ps_coord_x,ps_coord_y,text_size,res_name_string);
		}

	      /* Get pointer to the next residue */
	      residue_ptr = residue_ptr->next_residue_ptr;
	      iresid++;
	    }
	}

      /* Get pointer for next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

  /* Reset graphics state */
  psrest_();
}
/***********************************************************************

plot_postscript_picture  -  Plot the PostScript picture

***********************************************************************/

void plot_postscript_picture(float centre_x, float centre_y,
			     float scale_factor)
{
    /* If accessibility shading required, do it */
    if (Include->Accessibilities == TRUE)
      shade_accessibilities(centre_x,centre_y,scale_factor);

    /* Draw all the H-bond lines */	    
    plot_hbond_lines(centre_x,centre_y,scale_factor);
    
    /* Draw any additional covalent bonds */
    if (Include->External_Bonds == TRUE)
      plot_ext_bond_lines(centre_x,centre_y,scale_factor);
    
    /* Draw all the covalent bond lines */
    plot_bond_lines(centre_x,centre_y,scale_factor);
    
    /* Plot the atoms */
    plot_atoms(centre_x,centre_y,scale_factor);

    /* Plot the hydrophobic symbols */
    plot_hydrophobics(centre_x,centre_y,scale_factor);

    /* If this is a schematic ligand plot, plot the individual components */
    if (Include->Simple_Ligand_Residues == TRUE)
      {
	/* Print the residue names for the ligand */
	plot_simple_ligand_residues(centre_x,centre_y,scale_factor);
      }

    /* If the non-ligand residues are schematic, plot them */
    if (Include->Simple_Nonligand_Residues == TRUE)
      {
	/* Print the residue names for the H-groups */
	plot_simple_hgroups(centre_x,centre_y,scale_factor);
      }

    /* Print the residue names */
    print_residue_names(centre_x,centre_y,scale_factor);
}




/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  int atom_counter, nhbonds, ncontacts, nlinks;
/* v.3.1--> */
/*  int lig_atom_start, lig_atom_end, lig_end, lig_start; */
/* <--v.3.1 */
  int ndeleted, nligands, nobjects;
  int iseed;
/* v.3.1--> */
  int nligand_objects;
/* <--v.3.1 */
/* v.3.2--> */
  int show_deletions;
/* <--v.3.2 */

  float scale_factor;
  float centre_x, centre_y, min_x, max_x, min_y, max_y, x;
/* v.3.2--> */
  float matrix[3][3];
/* <--v.3.2 */

  char bonds_name[FILENAME_LEN], nnb_name[FILENAME_LEN],
  hhb_name[FILENAME_LEN], pdb_name[FILENAME_LEN];
  char chain_identi, res_num1[6], res_num2[6];
  char special_res[MAXSPECIAL_RES][15];
/* v.4.0--> */
  char res_name1[4], res_name2[4];
/* <--v.4.0 */

  /* Initialise global variables */
  scale_factor = 1.0;
/* v.4.0--> */
  Interface_Plot = FALSE;
/* <--v.4.0 */
  Nwarnings = 0;
  Water_as_Ligand = FALSE;
/* v.3.1.2--> */
  Metal_as_Ligand = FALSE;
/* <--v.3.1.2 */
/* v.3.2--> */
  show_deletions = FALSE;
  Write_Res_File = FALSE;
/* <--v.3.2 */
  
  /* Initialise the plot parameters and global variables */
  initialise_parameters();

  /* Open the parameter file and read user-definable information in */
  read_in_parameters(special_res);

  /* Initialise random number generator */
  iseed = -1;
  x = get_random_number(&iseed,Random_Start);

  /* Determine the maximum object and text sizes in the plot */
  define_object_sizes(scale_factor,FALSE);
  define_text_sizes(scale_factor,FALSE);

  /* sort out the command line arguments*/
/* v.3.2--> */
/*  get_command_arguments(argv,argc - 1,res_num1,res_num2,&chain_identi,
			hhb_name,nnb_name,bonds_name,Print_Title); */
/* v.4.0--> */
/*  get_command_arguments(argv,argc - 1,pdb_name,res_num1,res_num2,
			&chain_identi,hhb_name,nnb_name,bonds_name,
			Print_Title); */
  get_command_arguments(argv,argc - 1,pdb_name,res_name1,res_num1,
			res_name2,res_num2,&chain_identi,hhb_name,
			nnb_name,bonds_name,Print_Title);
/* <--v.4.0 */
/* <--v.3.2 */

/* v.3.2--> */
/*  strcpy(pdb_name,argv[1]); */
/* <--v.3.2 */

  /* Locate the start- and end-atoms and residues of the ligand
     in the PDB file */
/* v.3.1--> */
/*  get_ligand_start_end(&lig_start,&lig_end,&lig_atom_start,
		       &lig_atom_end,pdb_name,res_num1,res_num2,
		       chain_identi); */
/* v.4.0--> */
  if (Interface_Plot == FALSE)
/* <--v.4.0 */
/* v.4.0--> */
/*    get_ligand_start_end(pdb_name,res_num1,res_num2,chain_identi); */
    get_ligand_start_end(pdb_name,res_name1,res_num1,res_name2,res_num2,
			 chain_identi);
/* <--v.4.0 */
/* <--v.3.1 */

  /* If printing molecule as it stand, read in all the bonds and bond-
     lengths from the ligplot.bonds file */
  if (Print_as_is == TRUE)
    {
      /* Read in all the atom coordinates from the PDB file */
      atom_counter = read_pdb_file(pdb_name);

      /* Read in the bonds from ligplot.bonds */
      get_bonds(bonds_name);

      /* Define the objects that will appear on the final LIGPLOT diagram */
/* v.3.1--> */
/*      define_objects(&nobjects); */
      define_objects(&nobjects,&nligand_objects);
/* <--v.3.1 */

      /* Update all the residue boundaries */
      update_all_boundaries();
    }

  /* Otherwise, define all the bonds using the CONEC records, distances,
     and the HBPLUS output */
  else
    {
      /* Use the CONEC records to find which protein residues the ligand
	 is covalently bonded to */
      nlinks = 0;
/* v.3.1--> */
/*      read_conec_records(lig_atom_start,lig_atom_end,pdb_name,&nlinks); */
      read_conec_records(pdb_name,&nlinks);
/* <--v.3.1 */

      /* Check the stored CONEC records by calculating distances between
	 the atoms involved */
      if (nlinks > 0)
	verify_conec_records(pdb_name);

      /* Read in the ligand's hydrogen bonds from the .hhb file */
      if (Include->Hbonds == TRUE)
	{
	  read_hbplus_file(hhb_name,HHB_FILE,&nhbonds);

	  /* Include special case residues into the H_bond list - ie 
	     residues that are connected to other residues connected to 
	     the ligand */
	  if (Include->Linked_Residues == TRUE)
	    special_resinc(hhb_name,HHB_FILE,special_res);
	}

      /* Read in the hydrophobic contacts from the .nnb file */
      if (Include->Hydrophobics == TRUE)
	read_hbplus_file(nnb_name,NNB_FILE,&ncontacts);

      /* Read in all the atom coordinates from the PDB file */
      atom_counter = read_pdb_file(pdb_name);
      printf("   Number of hydrogen bonds stored   = %7d\n",nhbonds);
      printf("   Number of contacts stored         = %7d\n",ncontacts);

      /* If there have been no atoms read in, then abort */
      if (atom_counter == 0)
	{
	  printf("*** No atoms stored. Program terminated\n");
	  return(-1);
	}

      /* Define the objects that will appear on the final LIGPLOT diagram */
/* v.3.1--> */
/*      define_objects(&nobjects); */
      define_objects(&nobjects,&nligand_objects);
/* <--v.3.1 */
  
      /* Calculate covalent connectivity */
      covalent_connectivity();

      /* Combine objects where necessary */
      combine_objects();

/* v.3.1--> */
      /* Get rid of any objects that claim to be ligands, but aren't
         (eg protein residues with the same residue numbers as the
         ligand residues) */
      if (nligand_objects > 1)
	delete_false_ligands();
/* <--v.3.1 */

      /* Pick up H-bonds from .hhb file HBPLUS output, and add H-bond
	 connectivity */
      add_hbplus_bonds();

      /* Open the output PDB file */
/* v.3.2--> */
/*      open_pdb_output(); */
/* v.4.0--> */
/*      open_pdb_output(res_num1,res_num2,chain_identi); */
      open_pdb_output(res_name1,res_num1,res_name2,res_num2,chain_identi);
/* <--v.4.0 */
/* <--v.3.2 */

      /* Open the output bonds files: ligplot.bonds, ligplot.hhb
         and ligplot.nnb */
      open_bonds_output();
    }

  /* For each bond, determine which other bonds sprout off its two
     ends */
  bond_linkage();

  /* Remove any atoms that are not connected to other atoms in the
     same residue */
/* v.4.0--> */
  if (Print_as_is == FALSE)
/* <--v.4.0 */
    check_connections();

/* v.3.1.1--> */
  /* Check whether any objects are cyclic */
  check_for_cyclics();
/* <--v.3.1.1 */

  /* Split objects if there are no covalent bonds spanning adjacent
     residues */
  split_objects(&nobjects,&nligands);

  /* Determine which bonds relate to each object */
/* v.3.2--> */
/*  assign_bonds_to_objects(&ndeleted); */
  assign_bonds_to_objects(&ndeleted,show_deletions);
/* <--v.3.2 */

/* v.3.1.1--> */
  /* If any cyclic peptides identified, need to delete any bond connections
     involving covalent bonds that have now been made elastic */
  delete_unwanted_bond_connections();
/* <--v.3.1.1 */

  /* Mark all atoms that are reachable from the ligand residue(s) via
     one or more bonds */
/* v.4.0--> */
  if (Interface_Plot == FALSE && Print_as_is == FALSE)
    {
/* <--v.4.0 */
      mark_reachable_atoms();

      /* Delete any residues that are unreachable */
/* v.3.2--> */
/*  delete_unreachable_residues(&ndeleted); */
      delete_unreachable_residues(&ndeleted,show_deletions);
/* <--v.3.2 */
/* v.4.0--> */
    }
/* <--v.4.0 */

  /* If any residues have been deleted, then show explanatory message */
/* v.3.2--> */
/*  if (ndeleted > 0) */
  if (show_deletions == TRUE && ndeleted > 0)
/* <--v.3.2 */
    show_deleted_residues_message(ndeleted);

  /* Determine each atom's size */
  get_atom_sizes();

  /* For each atom, create a linked list of the atoms covalently bonded
     to it */
  atom_linkage();

  /* Determine which bonds in the structure are "rotatable" */
  identify_rotatable_bonds();

  /* If structure not to be plotted just as it stands, then perform
     all the necessary processing */
  if (Print_as_is == FALSE)
    {
/* v.4.0--> */
      /* Write out any covalent bonds between the ligand and the protein
	 to the .res file */
      if (Write_Res_File == TRUE)
	write_lig_prot_bonds();
/* <--v.4.0 */

      /* Print out the ligand and H-bonded sidechains into ligplot.frm */
      write_frm_file(pdb_name);

      /* Mark for deletion unwanted atoms in hydrophobic groups */
      pare_down_hydrophobics();

      /* Delete any unwanted atoms */
      delete_atoms_and_bonds();

      /* Transform the original coordinates, with each object oriented
	 along its principal axes, for fitting of objects during flattening
	 to be as close to original conformation as possible, and for
	 laying the objects out in their approximate relative arrangements */
      calc_fit_coordinates(nligands);

      /* Flatten all the objects, one by one */
      flatten_all_objects(&iseed);

      /* Score all the objects and then sort them by their relative
	 importance */
      score_objects();

      /* Sort objects by their weight scores */
      sort_objects();

/* v.4.0--> */
      /* Calculate the mean coordinates of each residue, for use in the
	 energy minimization process */
      calc_residue_means();

      /* Check for an input .rcm file containing residue centres of mass
	 to which the current plot's residues are to be restrained */
      read_rcm_file(pdb_name);

      /* For interface plot, work out how to place the residues on either
	 side of an interface */
      if (Interface_Plot == TRUE)
	lay_interface_objects_out();
      else
/* <--v.4.0 */

	/* Otherwise, lay objects out on the page around the ligand */
	lay_objects_out();

      /* Minimise the energy of the whole plot to reduce the numbers
	 of clashes and overlaps to a minimum */
      energy_minimize(&iseed,nligands);

      /* Reorientate ligand so its principal axis lies up the y-axis*/
/* v.4.0--> */
      if (Interface_Plot == FALSE && Have_Anchors == FALSE)
/* <--v.4.0 */
	ligand_orientate();

/* v.3.2--> */
    }

  /* If a further rotation is required, apply it now */
  if (Picture->Rotation_Angle != 0.0)
    {
      /* Calculate the rotation matrix */
      calc_rotation_matrix(Picture->Rotation_Angle,matrix);
      
      /* Apply the rotation to the whole structure */
      rotate_whole_structure(matrix);

      /* Update all residue boundaries */
      update_all_boundaries();
    }

  if (Print_as_is == FALSE)
    {
/* <--v.3.2 */

      /* Print out the new atom coordinates as a PDB file */
      write_pdb_file();

      /* Write out the bond data to the ligplot.bonds, ligplot.hhb
	 and ligplot.nnb files */
      write_bonds();

/* v.4.0--> */
      /* Write the residue centres of mass to file ligplot.rcm */
      write_rcm_file();
/* <--v.4.0 */
    }

  /* Plot molecule as a PostScript file */
  
  /* Open the PostScript file */
  psopen_();

  /* If picture to be plotted in landscape mode, then rotate page */
  if (Picture->Portrait == FALSE)
    {
      adjust_for_landscape();
      psrot_(BBOXX2 + BBOXX1,0.0);
    }

  /* Print the title, if there is one */
  if (Print_Title[0] != '\0')
    print_title(Print_Title);

  /* Write out the key to symbols used, if required */
  if (Include->Key == TRUE)
    print_symbol_key();

  /* Calculate the minimum and maximum x-y extent of the atoms to
     be plotted */
  extract_minmax(&min_x,&max_x,&min_y,&max_y);

  /* Calculate the scaling factor to fit the coords on the PostScript page */
  scale(&scale_factor,min_x,max_x,min_y,max_y);

  /* Get the equivalent transformation of the centre-point of the
     picture to be plotted, to the centre-point of the PostScript page */
  centre_point(&centre_x,&centre_y,scale_factor,
	       min_x,max_x,min_y,max_y);

  /* Calculate scaled sizes of the different items in the plot */
  pscomm_(" ");
  pscomm_("MAIN PLOT");
  pscomm_("Scaled object and text sizes");
  define_object_sizes(scale_factor,TRUE);
  define_text_sizes(scale_factor,TRUE);

  /* Plot the PostScript picture */
  plot_postscript_picture(centre_x,centre_y,scale_factor);
    
  /* Close the PostScript file */
  psclos_();

  printf("\nProgram complete\n\n");

  if (Nwarnings == 1)
    printf("+++ There was 1 warning or advisory message\n\n");
  else if (Nwarnings > 1)
    printf("+++ There were %d warning or advisory messages\n\n",Nwarnings);

  return(0);
}
