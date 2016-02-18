/**************************************************************************
 *
 * hbadd.h - Include file for program hbadd.c
 *
 *************************************************************************/

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>

/* Flag-values */
#define FALSE           0
#define TRUE            1

/* Array parameters */
#define FILENAME_LEN  130         /* Max filename-length */
#define LINE_LEN      130         /* Max line-length */
#define MAXCONNECT      9         /* Max no. of connections per atom */
#define MAXANGLES      10         /* Max no. of bond angles per atom */
/* v.3.2--> */
#define MAXSTACK     1000         /* Max no. of entries on the heap-stack */
/* <--v.3.2 */

/* Miscellaneous parameters */
#define COVAL_DIST    2.0
#define COVAL_DIST2 (COVAL_DIST * COVAL_DIST)
/* v.3.1--> */
#define HATOM_DIST    1.42
#define HATOM_DIST2 (HATOM_DIST * HATOM_DIST)
/* <--v.3.1 */
#define INT_MAX     32767

/* Angle constants */
#define  PI          3.141592654
#define  RADDEG      (180.0 / PI)
#define  PI_BY_2     (PI / 2.0)
#define  ANGLE_ERROR    4.0       /* Allowed bond-angle error */

/* Error flag */
int error_flag;

/* File pointers */
FILE *fil_out;
FILE *fil_dic;
FILE *fil_pdb;

/* Structure for residue data read in from PDB file
   ------------------------------------------------ */
struct residue
{
  char              res_name[4];
  char              res_num[6];
  char              chain;
  int               natoms;
/* v.3.2--> */
  int               non_hydrogen;
/* <--v.3.2 */
  struct atom       *first_atom_ptr;
  struct dictionary *dictionary_ptr;
  struct residue    *next_residue_ptr;
};

/* Structure for atom coordinates data read in from PDB file
   --------------------------------------------------------- */
struct atom
{
  int               atom_number;
  char              atom_name[5];
/* v.3.2--> */
  int               hydrogen;
/* <--v.3.2 */
  float             x, y, z;
  struct residue    *residue_ptr;
/* v.3.1 -> */
  struct dic_data   *dic_data_ptr;
/* <- v.3.1 */
/* v.3.2--> */
  int               checked;
  struct atom       *connect_atom_ptr[MAXCONNECT];
/* <--v.3.2 */
  struct atom       *next_atom_ptr;
};

/* Structure for atom connectivities, as read in from PDB file
   ----------------------------------------------------------- */
struct connect
{
  int               atom_number1;
  int               atom_number2;
  struct connect    *next_connect_ptr;
};

/* Structure for HET group entries from HET Group Dictionary
   --------------------------------------------------------- */
struct dictionary
{
  char              res_name[4];
  int               natoms;
  int               nmatch;
/* v.3.2--> */
  int               non_hydrogen;
  int               have_hydrogens;
  int               fake_entry;
/* <--v.3.2 */
  struct residue    *residue_ptr;
  struct dic_data   *first_dic_data_ptr;
  struct dictionary *next_dictionary_ptr;
};


/* Structure for atom data read in from HET Group Dictionary
   --------------------------------------------------------- */
struct dic_data
{
  char              atom_name[5];
/* v.3.2--> */
  char              output_atom_name[5];
  int               hydrogen;
/* <--v.3.2 */
  int               number;
  int               nbound;
  char              boundto[MAXCONNECT][5];   
  int               boundto_num[MAXCONNECT];
/* v.3.2--> */
  int               checked;
/* <--v.3.2 */
  struct atom       *atom_ptr;
  struct dic_data   *connect_atom_ptr[MAXCONNECT];
  struct dic_data   *next_dic_data_ptr;
};

/* v.3.2--> */
/* Structure for graph nodes
   ------------------------- */
struct node
{
  int               node_number;
  int               checked;
  int               in_best_clique;
  struct atom       *atom_ptr;
  struct dic_data   *dic_data_ptr;
  struct node       *next_node_ptr;
};

/* Structure for links between graph nodes
   --------------------------------------- */
struct node_link
{
  int               link_number;
  int               checked;
  struct node       *node1_ptr;
  struct node       *node2_ptr;
  struct node_link  *next_node_link_ptr;
};

/* <--v.3.2 */
