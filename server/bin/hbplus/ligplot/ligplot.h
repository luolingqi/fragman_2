#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>

/* Flag-values */
#define FALSE           0
#define TRUE            1
#define FREE         TRUE
#define NOT_FREE    FALSE
#define OFF         FALSE
#define ON           TRUE
#define RESTORE     FALSE
#define SAVE         TRUE
#define PLUS            1
#define BOTH            0
#define MINUS          -1
#define ALL_ATOMS      -1

/* Object types */
#define LIGAND          1
#define HGROUP          2
#define HYDROPHOBIC     3
#define WATER           4
#define SIMPLE_HGROUP   5

/* Residue types */
#define STANDARD        1
#define SIMPLE_LIGAND   2

/* Bond types */
#define COVALENT        0
#define HBOND           1
#define CONTACT         2
#define INTERNAL        3
#define DELETED       -99

/* Bond sources */
#define CONECT          0
#define CALCULATED      1
#define HBPLUS          2

/* Hhb_info sources */
/* v.3.1--> */
/*#define EXTRA          -1 */
/* <--v.3.1 */
#define CONECT          0
#define HHB_FILE        1
#define NNB_FILE        2

/* Bond ends */
#define EITHER          0
#define FIRST           1
#define SECOND          2

/* Move types */
#define MOVE_OBJECT     0
#define SWIVEL_OBJECT   1
#define SWING_OBJECT    2
#define FLIP_OBJECT     3
#define FLIP_ABOUT_ATOM 4

/* v.3.2--> */
/* Non-bonded contact types */
#define HYDROPHOBIC_ONLY 0
#define HYDROPHOBIC_ANY  1
#define ANY_ANY          2
/* <--v.3.2 */

/* Loop maxima */
#define MAXCYCLE       50
#define START_SWIVELS  (MAXCYCLE / 3)

/* Array parameters */
#define COL_NAME_LEN   15
#define CONNECT        20
#define CONNECTIONS    30
#define FILENAME_LEN   80
#define HHB            24
#define LINELEN       100
/* v.4.0--> */
/* #define MAXATOMS     1200 */
#define MAXATOMS     5000
/* <--v.4.0 */
#define MAXBONDS      600
#define MAX_COLOURS    20
#define MAXHBONDS     800
#define MAXLIGRES    1000
#define MAXNATM      1000
#define MAXNRES       150
#define MAX_PARAMETERS  4
#define MAXRESIDUES 10000
#define MAXLABELS  (2 * MAXRESIDUES)
/* v.4.0--> */
/* #define MAXSIDESTACK (MAXATOMS * MAXBONDS) */
/* <--v.4.0 */
/* v.3.2--> */
/* v.4.0.2--> */
/* #define MAXTOKENS       8 */
#define MAXTOKENS       20
/* <--v.4.0.2 */
/* <--v.3.2 */
#define MAX_FLAT_LOOP  10
#define MAXELBONDS     50
#define MAXSPECIAL_RES 10
#define MAXSTACK      600
#define POINTER         3
#define TITLE_LEN      80

/* Minimization parameters */
/* v.4.0--> */
/* #define NTERMS          5 */
#define NTERMS          7
/* <--v.4.0 */
#define MOVE_TYPES      4
#define MIN_ANGLE    10.0
#define MIN_MOVE      1.5

/* Miscellaneous parameters */
/* v.4.0--> */
/* #define INT_MAX       32767 */
/* <--v.4.0 */
#define SWIVEL_ANGLE  10.0

/* Plot parameters */
#define BBOXX1         30.0
#define BBOXX2        550.0
#define BBOXY1         50.0
#define BBOXY2        780.0
#define BORDER_MARGIN   0.93
#define CHAR_ASPECT    0.484
#define SPOKE_EXTENT   1.4
#define SPREAD_FACTOR  1.0
#define TITLE_FRACTION 0.04
#define TITLE_SIZE     20.0

/* Parameter arrays */
/* v.3.2--> */
/* #define PICTURE_OPTIONS 2
   #define INCLUDE_PARAMETERS      19 */
#define PICTURE_OPTIONS          3
/* v.4.0--> */
/* #define INCLUDE_PARAMETERS      21 */
#define INCLUDE_PARAMETERS      23
/* <--v.4.0 */
/* <--v.3.2 */
#define OBJECT_COLOURS          18
#define OBJECT_SIZES             9
#define TEXT_COLOURS             9
#define TEXT_SIZES               8
/* v.4.0--> */
/* #define MINIMIZATION_PARAMETERS 11 */
#define MINIMIZATION_PARAMETERS 15
/* <--v.4.0 */

/* Parameters defining layout of Key to symbols */
#define KEY_SCALE      0.15

/* Hydrophobic symbol parameters */
#define HYDROPHOBIC_RADIUS  1.1
#define SPOKE_ANGLE  12.0

#define NON_FLAT      0.2
#define KEY_HEIGHT    0.14
#define MARGIN_ACC    0.04
#define MARGIN_NO_HPH 0.04
#define IDEAL_ENERGY 10.0
#define SWING_STEPS  45
#define MOVE_STEPS   30

#define ATOM_INTERACT_DIST      1.0
#define ATOM_DIST2      (ATOM_INTERACT_DIST * ATOM_INTERACT_DIST)
#define BOND_INTERACT_DIST      0.6
#define INTERACT_DIST2  (BOND_INTERACT_DIST * BOND_INTERACT_DIST)
#define INTERNAL_INTERACT_DIST  0.3
#define INTERNAL_DIST2  (INTERNAL_INTERACT_DIST * INTERNAL_INTERACT_DIST)
#define PS_STRING_LENGTH        132
#define MAXH         4
#define HSCALE_FACTOR 0.0

/* Angle constants */
#define  PI          3.141592654
#define  RADDEG      (180.0 / PI)
#define  PI_BY_2     (PI / 2.0)

/* File pointers */
FILE *fil_bonds;
FILE *fil_hbplus;
FILE *fil_pdb;
FILE *ligplot_bonds_out;
FILE *ligplot_frm;
FILE *ligplot_hhb_out;
FILE *ligplot_nnb_out;
FILE *ligplot_par;
FILE *ligplot_pdb;
/* v.3.2--> */
FILE *ligplot_res;
/* <--v.3.2 */
/* v.4.0--> */
FILE *ligplot_rcm;
/* <--v.4.0 */

FILE *ps_file;

/* Miscellaneous counters, flags and arrays */
int Nobjects;
int Water_as_Ligand;
/* v.3.1.2--> */
int Metal_as_Ligand;
/* <--v.3.1.2 */
/* v.3.2--> */
int Write_Res_File;
/* <--v.3.2 */
int xsite_file;
char xsite_probe[7];
/* v.4.0--> */
int Interface_Plot;
/* <--v.4.0 */

float Maximum_accessibility, label_coord[MAXLABELS][2];
int n_labels;

char root_name[TITLE_LEN];

char lig_res_name_store[MAXLIGRES][4];
char lig_res_num_store[MAXLIGRES][6];
char lig_chain_store[MAXLIGRES];
int ligand_residues;

/* Structure for information read in from HBPLUS files
   --------------------------------------------------- */
struct hhb_info
{
  int                source;
  int                atom_number1;
  char               atom_type1[5];
  char               res_name1[4];
  char               res_num1[6];
  char               chain1;
  float              x1, y1, z1;
  int                atom_number2;
  char               atom_type2[5];
  char               res_name2[4];
  char               res_num2[6];
  char               chain2;
  float              x2, y2, z2;
  float              bond_length;
  struct hhb_info    *next_hhb_info_ptr;
};
struct hhb_info *first_hhb_info_ptr, *last_hhb_info_ptr;

/* Structure for objects to be drawn on the final LIGPLOT picture
   -------------------------------------------------------------- */
struct object
{
  int                object_type;
  int                nresidues;
  int                nbonds;
  int                nrot_bonds;
  int                nmain_chain;
  int                n_ca;
  int                weight;
  float              minx, miny, maxx, maxy;
/* v.4.0--> */
  int                interface;
  int                have_anchors;
/* <--v.4.0 */
  float              max_atom_size;
  float              place_x, place_y;
  float              internal_energy;
  float              total_energy;
  struct residue     *first_residue_ptr;
  struct object_bond *first_object_bond_ptr;
  struct object      *next_object_ptr;
};
struct object *first_object_ptr;

/* Structure for residue data read in from PDB file
   ------------------------------------------------ */
struct residue
{
  char               res_name[4];
  char               res_num[6];
  char               chain;
  int                residue_type;
  int                natoms;
  int                inligand;
  float              minx, miny, maxx, maxy;
  float              max_atom_size;
/* v.3.2--> */
  int                deleted;
/* <--v.3.2 */
/* v.4.0--> */
  float              original_mean_x, original_mean_y, original_mean_z;
  float              flattened_mean_x, flattened_mean_y;
  int                have_anchor;
  float              anchor_pstn_x, anchor_pstn_y;
  int                nanchored_atoms;
/* <--v.4.0 */
  struct coordinate  *first_atom_ptr;
  struct object      *object_ptr;
  struct residue     *next_residue_ptr;
};
struct residue *first_residue_ptr;

/* Structure for atom data
   ----------------------- */
struct coordinate
{
  char               atom_number[6];
  char               atom_type[5];
  char               print_name[5];
  struct residue     *residue_ptr;
  float              original_x, original_y, original_z;
  float              x, y, z;
  float              fit_x, fit_y, fit_z;
  float              save_x, save_y;
  char               occupancy[7];
  char               bvalue[7];
  float              accessibility;
  float              atom_size;
  int                side_chain;
  int                checked;
  int                deleted;
  int                natom_links;
  int                plot_atom;
/* v.4.0--> */
  int                have_anchor;
  float              anchor_pstn_x, anchor_pstn_y;
/* <--v.4.0 */
  struct coordinate  *next_stack_ptr;
  struct atom_link   *first_atom_link_ptr;
  struct coordinate  *next;
};
struct coordinate *first_atom_ptr;

/* Structure containing info about each bond
   ---------------------------------------- */
struct bond
{
  struct coordinate  *first_atom_ptr;
  struct coordinate  *second_atom_ptr;
  int                bond_type;
  int                elastic;
  int                rotatable_bond;
  int                ring_bond;
  int                end_bond;
  int                bond_order;
  int                bond_source;
  float              bond_length;
  int                checked;
  int                flattened;
  int                nfirst_atom_links;
  int                nsecond_atom_links;
  int                nbonds_from_first_atom;
  int                nbonds_from_second_atom;
  struct bond_link   *first_bond_link_ptr;
  struct bond        *next_stack_ptr;
  struct bond        *next_link_ptr;
  struct bond        *next_bond_ptr;
};
struct bond *first_bond_ptr;


/* Structure giving atoms covalently bonded to current atom
   -------------------------------------------------------- */
struct atom_link
{
  struct coordinate  *atom_ptr;
  struct atom_link   *next_atom_link_ptr;
};

/* Structure giving bonds linked to either end of each bond
   -------------------------------------------------------- */
struct bond_link
{
  struct bond        *bond_ptr;
  int                bond_end;
  struct bond_link   *next_bond_link_ptr;
};

/* Structure giving bonds belonging to each object
   ----------------------------------------------- */
struct object_bond
{
  struct bond        *bond_ptr;
  struct object_bond *next_object_bond_ptr;
};




/* Plot parameters defined by user
   ------------------------------- */
struct Picture_Options
{
int In_Colour, Portrait;
/* v.3.2--> */
float Rotation_Angle;
/* <--v.3.2 */
};
struct Picture_Options *Picture;

struct Include_Parameters
{
int Hydrophobics, Waters, Mainchain_Atoms, Linked_Residues, Hbonds,
    Internal_Hbonds, External_Bonds, Hydrophobic_Bonds,
    Simple_Ligand_Residues, Simple_Nonligand_Residues, Accessibilities,
    Ligand_Atoms, Nonligand_Atoms, Double_Bonds, Key,
    Residue_Names, Atom_Names, Hbond_Lengths, Filename_for_Title;
/* v.3.2--> */
int External_Bonds_Solid, Contact_Type;
/* v.4.0--> */
int Ligand_Accessibilities_Only, Water_Atoms;
/* <--v.4.0 */
/* <--v.3.2 */
};
struct Include_Parameters *Include;

int Include_Hydrogens;

struct Sizes
{
char Ligand_Atoms[COL_NAME_LEN + 1],
     Nonligand_Atoms[COL_NAME_LEN + 1],
     Waters[COL_NAME_LEN + 1],
     Hydrophobics[COL_NAME_LEN + 1],
     Simple_Residues[COL_NAME_LEN + 1],
     Ligand_Bonds[COL_NAME_LEN + 1],
     Nonligand_Bonds[COL_NAME_LEN + 1],
     Hydrogen_Bonds[COL_NAME_LEN + 1],
     External_Bonds[COL_NAME_LEN + 1];
};
struct Sizes *Size;

struct Size_Vals
{
float Ligand_Atoms, Nonligand_Atoms, Waters, Hydrophobics, Simple_Residues,
      Ligand_Bonds, Nonligand_Bonds, Hydrogen_Bonds, External_Bonds;
};
struct Size_Vals *Size_Val;

float object_size[OBJECT_SIZES];

struct Text_Sizes
{
char Ligand_Residue_Names[COL_NAME_LEN + 1],
     Nonligand_Residue_Names[COL_NAME_LEN + 1],
     Water_Names[COL_NAME_LEN + 1],
     Hydrophobic_Names[COL_NAME_LEN + 1],
     Simple_Residue_Names[COL_NAME_LEN + 1],
     Ligand_Atom_Names[COL_NAME_LEN + 1],
     Nonligand_Atom_Names[COL_NAME_LEN + 1],
     Hbond_Lengths[COL_NAME_LEN + 1];
};
struct Text_Sizes *Text_Size;

struct Text_Size_Vals
{
float Ligand_Residue_Names, Nonligand_Residue_Names, Water_Names,
      Hydrophobic_Names, Simple_Residue_Names, Ligand_Atom_Names, 
      Nonligand_Atom_Names, Hbond_Lengths;
};
struct Text_Size_Vals *Text_Size_Val;

float text_size[TEXT_SIZES];

struct Colours
{
char Background[COL_NAME_LEN + 1], Ligand_Bonds[COL_NAME_LEN + 1],
    Nonligand_Bonds[COL_NAME_LEN + 1], Hydrogen_Bonds[COL_NAME_LEN + 1],
    External_Bonds[COL_NAME_LEN + 1], Hydrophobics[COL_NAME_LEN + 1],
    Accessibility_Min[COL_NAME_LEN + 1], Accessibility_Max[COL_NAME_LEN + 1],
    Nitrogen[COL_NAME_LEN + 1], Oxygen[COL_NAME_LEN + 1],
    Carbon[COL_NAME_LEN + 1], Sulphur[COL_NAME_LEN + 1],
    Water[COL_NAME_LEN + 1], Phosphorus[COL_NAME_LEN + 1],
    Iron[COL_NAME_LEN + 1], Other[COL_NAME_LEN + 1],
    Atom_Edges[COL_NAME_LEN + 1], Simple_Residues[COL_NAME_LEN + 1];
};
struct Colours *Colour;

char object_colour[OBJECT_COLOURS][COL_NAME_LEN + 1];

struct Text_Colours
{
char Title[COL_NAME_LEN + 1],
    Key_Text[COL_NAME_LEN + 1],
    Ligand_Residue_Names[COL_NAME_LEN + 1],
    Nonligand_Residue_Names[COL_NAME_LEN + 1],
    Water_Names[COL_NAME_LEN + 1],
    Hydrophobic_Names[COL_NAME_LEN + 1],
    Ligand_Atom_Names[COL_NAME_LEN + 1],
    Nonligand_Atom_Names[COL_NAME_LEN + 1],
    Hbond_Lengths[COL_NAME_LEN + 1];
};
struct Text_Colours *Text_Colour;

char text_colour[TEXT_COLOURS][COL_NAME_LEN + 1];

char Colour_Table_Name[MAX_COLOURS][COL_NAME_LEN + 1];
float Colour_Table[MAX_COLOURS][3];

char special_res[MAXSPECIAL_RES][15];

char Print_Title[TITLE_LEN + 1];

/* Minimization parameters */
int Max_Loops, Random_Start;
float Atom_Atom_Clash, Bond_Atom_Clash, Bond_Stretch, Energy_Drop_Maximum;
float HB_Weight, Internal_Energy_Weight, Max_HBdist, Nonbond_Weight;
float Overlap_Score;
/* v.4.0--> */
float Anchor_Weight, Boundary_Energy_Weight, Min_Boundary_Dist;
float Rel_Pstn_Energy_Weight;
/* <--v.4.0 */

/* Global Plot variables
   --------------------- */
int Nwarnings;
int Print_as_is, Split_Colour_Ligand_Bonds, Split_Colour_Nonligand_Bonds;
float Margin_y;
float Page_Max_x, Page_Min_x, Page_Max_y, Page_Min_y;
float Plot_Centre_x, Plot_Centre_y;
float Accessibility_Max[3], Accessibility_Min[3], Max_atom_radius,
      Max_object_size, Max_object_height, Max_object_width,
      Max_Textsize, Max_Textwidth;
/* v.4.0--> */
int Have_Anchors;
/* <--v.4.0 */

/*************************************************************************

Include variables for Andrew Martin's least-squares fitting routines,
matfit and qikfit

*************************************************************************/

#define SMALL  1.0e-20     /* Convergence cutoffs                       */
#define SMALSN 1.0e-10

typedef short  BOOL;
typedef double REAL;
typedef struct
{  REAL x, y, z;
}  VEC3F;

typedef VEC3F COOR;

#define ABS(x)   (((x)<0)   ? (-(x)) : (x))

/* Prototypes for functions defined in fit.c */
/* BOOL matfit(COOR *x1, COOR *x2, REAL rm[3][3], int n, REAL *wt1, 
            BOOL column);

static void qikfit(REAL umat[3][3], REAL rm[3][3], BOOL column);  */

