/**************************************************************************
 *
 * Dimer.h - Include file for program Dimer.c
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

/* Output file names */
#define OUTFILE_1 "dimplot.pdb"
#define OUTFILE_2 "dimplot.hhb"
#define OUTFILE_3 "dimplot.nnb"

/* Array parameters */
#define FILENAME_LEN  120         /* Max filename-length */
#define LINE_LEN      250         /* Max line-length */
#define NUMBER_LEN     20         /* Max length for number-string */
#define TOKEN_LEN      20         /* Max length for token-string */

/* HBPLUS file-types */
#define HBONDS          1
#define CONTACTS        2

/* Structure for atom data
   ----------------------- */
struct atom 
{
  char             label[7];
  int              atom_number;
  char             atom_name[5];
  float            x,y,z;
  char             occupancy[7];
  char             bvalue[7];
  struct residue   *residue_ptr;
  struct atom      *next_atom_ptr;                  
};
struct atom *first_atom_ptr, *last_atom_ptr;

/* Structure for residue data read in from PDB file
   ------------------------------------------------ */
struct residue
{
  char             res_name[4];
  char             res_num[6];
  char             chain;
  int              num_atoms;
  int              chain_domain;
  int              interacting;
  struct atom      *first_atom_ptr;
  struct residue   *next_residue_ptr;
};
struct residue *first_residue_ptr, *last_residue_ptr; 

/* Structure for information read in from HBPLUS files
   --------------------------------------------------- */
struct hhb_info
{
  int                type;
  char               atom_name1[5];
  char               res_name1[4];
  char               res_num1[6];
  char               chain1;
  int                chain_domain1;
  char               atom_name2[5];
  char               res_name2[4];
  char               res_num2[6];
  char               chain2;
  int                chain_domain2;
  float              bond_length;
  struct residue     *residue1_ptr;
  struct residue     *residue2_ptr;
  struct hhb_info    *next_hhb_info_ptr;
};
struct hhb_info *first_hhb_info_ptr, *last_hhb_info_ptr;

/* Structure for domain data read in from .dom file
   ------------------------------------------------ */
struct domain 
{
  int              domain_number;
  char             chain;
  char             start_residue_number[6];
  char             end_residue_number[6];
  struct domain    *next_domain_ptr;                  
};
struct domain *first_domain_ptr, *last_domain_ptr;


