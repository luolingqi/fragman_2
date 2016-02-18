/**************************************************************************
 *
 * hbadd.c - Program to prepare a set of hbplus.rc entries for a given PDB
 *           file. An hbplus.rc entry is made for every HET group in the
 *           PDB file that is correctly located in the HET group dictionary,
 *           het_dictionary.txt, obtained from the PDB (from the PDB
 *           documentation page at http://www.pdb.bnl.gov/doc_help.html).
 *
 *           The hbplus.rc entry defines the connectivities of the atoms
 *           in the HET group and each atom's numbers of available donors
 *           and acceptors. It is used by program HBPLUS for calculating
 *           H-bonds.
 *
 **************************************************************************
 *
 * Version:-      v.1.0  -  30 January 1997
 * Written by:-   Nicholas M Luscombe(1) and Roman A Laskowski(2)
 * 
 *            (1) Department of Biochemistry and Molecular Biology
 *                University College London
 *                Gower Street, London WC1E 6BT, UK.
 * 
 *            (2) Department of Crystallography
 *                Birkbeck College
 *                Malet Street,  London WC1E 7HX, UK.
 *
 * Requirements:- An input PDB file for which the hbplus.rc entries are
 *                required. The HET Group Dictionary, het_dictionary.txt,
 *                obtained from the PDB from:-
 *
 *                http://www.pdb.bnl.gov/doc_help.html
 *
 *
 *
 * ------------------------------------------------------------------------
 *
 * Program description
 * -------------------
 *
 * Program reads in the given PDB file and identifies all the HET groups
 * involved. Computes connectivities using CONECT records and distances
 * between atoms.
 *
 * Locates each HET group in the HET Group Dictionary. Connectivities are
 * recorded. Atoms are matched by name and then by connectivity. Where the
 * HET group in the PDB file successfully matches the dictionary definition
 * on both these criteria, the relevant bond angles are computed where
 * necessary.
 *
 * Rules for H-bond formation are applied to all O and N atoms according to 
 * the following rules (John Mitchell, personal communication):-
 * 
 * H-bond donors: 
 *    Any O, S or N potentially donates the number of Hs bound to it.
 *
 * H-bond acceptors:
 *    O/S: sp3 - can accept up to 2
 *         sp2 - can accept up to 2
 *
 *    N:   planar   - no acceptance
 *         sp3      - can accept 1
 *         sp2      - can accept 1
 *         aromatic - can accept (3 - no. of bonds)
 *         amide    - no acceptance
 *
 * ------------------------------------------------------------------------
 *
 * Original version was part of v.3.0 of LIGPLOT. Subsequent amendments
 * have been labelled by v.m.n--> and <--v.m.n where m.n is the
 * version number corresponding to the change
 *
 * v.3.1          Amendment to write out the hbplus.rc file for HET groups
 *                even if not in the Het Group Dictionary, using the
 *                connectivity information in the PDB file.
 *                                                        12 Apr 1997 (RAL)
 *                Bug-fix to write out 2-character residue-names in the
 *                format that HBPLUS requires.
 *                                                        15 Apr 1997 (RAL)
 * v.3.1.2        Amendment to stop reading PDB file ATOM data when have
 *                hit an ENDMDL record.
 *                Amendment for metals to fool HBPLUS in recognizing bonds
 *                to metal ions.
 *                Bug-fix in routine validate_pdb_hets from v.3.1.
 *                                                        23 Apr 1997 (RAL)
 * v.3.2          Addition of graph-matching routines to identify atom
 *                correspondences between the PDB het group and the entry
 *                in the Het Group Dictionary when the atom names don't
 *                match. Functions check_atoms_present and
 *                check_connectivities replaced by function graph_match and
 *                its subfunctions.
 *                Generation of fake entries in cases when no match with
 *                the Het Group Dictionary possible. Numbers of acceptors
 *                and donors set to maximum for given atom type and bond
 *                connectivity in these cases.
 *                                                         1 May 1997 (RAL)
 *                Amendement for atom names containing double-quotes (eg
 *                in 1gac). Double-quotes written out as @ so that HBPLUS
 *                doesn't fail on reading string in. Amendment made to
 *                HBPLUS to convert @'s back to double-quotes.
 *                                                        25 May 1997 (RAL)
 *                Amendement to open hbplus.rc file only the first time it
 *                is needed, so that an empty file is not created if no data
 *                needs to be written out (or this causes problems for
 *                HBPLUS)
 *                                                         8 Jun 1997 (RAL)
 * v.3.3          Amendment to include ATOM-record residues as well as 
 *                HETATM ones (specifically for dealing with non-standard
 *                amino acids and nucleic acids).
 *                                                         6 Oct 1997 (RAL)
 *
 * ------------------------------------------------------------------------
 *
 * Datafiles
 * ---------
 *
 * Actual name          File pointer    Description
 * -----------          ------------    -----------
 * <filename>.pdb       fil_pdb         Input PDB file holding coords of
 *                                      protein and for which hbplus.rc 
 *                                      entries are required.
 * het_dictionary.txt   fil_dic         The Het Group Dictionary.
 * hbplus.rc            fil_out         File containing the hbplus.rc
 *                                      entries.
 * 
 * ------------------------------------------------------------------------
 * 
 * Function calling-tree
 * ---------------------
 * 
 * main
 *    -> get_command_arguments
 *    -> read_pdbfile
 *          -> store_connections
 *    -> read_dictionary
 *    -> dic_connectivities
 *    -> validate_pdb_hets
 *          -> print_missing_hets
 *          -> generate_fake_entry
 *                -> generate_connectivities
 *                      -> check_stored_conecs
 *          -> graph_match
 *                -> calc_covalent_connectivities
 *                -> create_nodes
 *                -> calc_node_links
 *                -> clique_detection
 *                      -> initialise_atoms
 *                      -> initialise_dic_data
 *                      -> initialise_nodes
 *                      -> initialise_node_links
 *                      -> count_atom_matches
 *                      -> store_best_clique
 *                -> get_best_match
 *    -> generate_hbplusrc
 *          -> open_hbplusrc
 *          -> writename_fn
 *                -> convert_quotes
 *          -> writecon_fn
 *          -> writehbdon_fn
 *                -> check_for_metal
 *          -> writehbacc_fn
 *                -> check_atoms_present
 *
 *************************************************************************/


/* H E A D E R   F I L E */

#include "hbadd.h"
#include <math.h>


/***********************************************************************

get_command_arguments  -  Get the name of the input PDB file and Het Group
                          Dictionary from the command-line arguments

***********************************************************************/

void get_command_arguments(char *strptrs[],int num_args,
			   char pdb_name[FILENAME_LEN],
			   char dictionary_name[FILENAME_LEN])
{
  int namlen;

  /* If no command-line arguments entered, then report error */
  if (num_args < 2)
    {
      printf("\n");
      printf("*** ERROR. Program hbadd requires 2 arguments:\n");
      printf("\n");
      printf("        hbadd  pdb_filename  dictionary_filename\n");
      printf("\n");
      printf("    where  pdb_filename is the PDB file to be processed\n");
      printf("      and  dictionary_filename is the name of the Het");
      printf(" Group Dictionary\n");
      printf("\n");
      printf("*** Program terminated with error.\n");
      exit(1);
    }

  /* Get the PDB file name */

  /* Check that length of filename not too long */
  namlen = strlen(strptrs[1]);
  if (namlen > FILENAME_LEN - 1)
    {
      printf("\n");
      printf("*** ERROR. Length of PDB filename longer than %d characters\n",
	     FILENAME_LEN);
      printf("***        [%s]\n",strptrs[1]);
      printf("***        Program hbadd terminated with error.\n");
      printf("\n");
      exit(1);
    }

  /* Store the PDB filename */
  strcpy(pdb_name,strptrs[1]);

  /* Get the Het Group Dictionary file name */

  /* Check that length of filename not too long */
  namlen = strlen(strptrs[2]);
  if (namlen > FILENAME_LEN - 1)
    {
      printf("\n");
      printf("*** ERROR. Length of dictionary filename longer ");
      printf("than %d characters\n",FILENAME_LEN);
      printf("***        [%s]\n",strptrs[2]);
      printf("***        Program hbadd terminated with error.\n");
      printf("\n");
      exit(1);
    }

  /* Store the PDB filename */
  strcpy(dictionary_name,strptrs[2]);
}
/***********************************************************************

store_connections  -  Process the current CONECT record and store its
                      connections if they relate to any of the HET groups

***********************************************************************/

void store_connections(char line[LINE_LEN + 1],int first_atom_number, 
/* v.3.1.2--> */
		       int last_atom_number,
/* <--v.3.1.2 */
		       struct connect **first_connect_ptr,
		       struct connect **last_connect_ptr,int *nlinks)
{
  char atmnum[6];
  int done, inhetgroup;
  int iatom, icon, ipos, jatom, ncon;
  int inhetgroup_atom[8], store_atom[8];

  struct connect *connect_ptr;
  
  /* Initialise variables */
  inhetgroup = 0;
  ipos = 6;
  ncon = 0;
  done = FALSE;

  /* Loop through the character positions in the record to extract
     all the atom numbers in the CONECT record */
  while (done == FALSE && ncon < 5)
    {
      /* Extract this atom number and store */
      strncpy(atmnum,line+ipos,5);
      atmnum[5] = '\0';

      /* If have blank, then end of connections reached */
      if (!strncmp(atmnum,"     ",5))
	done = TRUE;
      
      /* Otherwise, store the number */
      else
	{
	  store_atom[ncon] = atoi(atmnum);
	  
	  /* Check whether the atom might belong to one of the HET groups */
/* v.3.1.2--> */
/*	  if (store_atom[ncon] >= first_atom_number) */
	  if (store_atom[ncon] >= first_atom_number &&
	      store_atom[ncon] <= last_atom_number)
/* <--v.3.1.2 */
	    {
	      inhetgroup++;
	      inhetgroup_atom[ncon] = TRUE;
	    }
	  else
	    {
	      inhetgroup_atom[ncon] = FALSE;
	    }

	  /* Increment variables */
	  ncon++;
	  ipos = ipos + 5;
	}
    }

  /* If have possible covalent bonds between HET group atoms, then store */
  if (inhetgroup > 1)
    {
      /* Consider each connection in turn */
      iatom = store_atom[0];
      for (icon = 1; icon < ncon; icon++)
	{
	  jatom = store_atom[icon];

	  /* If both atoms are in a HET groups then store this connection */
	  if ((inhetgroup_atom[0] == TRUE && inhetgroup_atom[icon] == TRUE))
	    {
	      /* Allocate memory for structure to hold current atom's
		 details */
	      connect_ptr
		= (struct connect*)malloc(sizeof(struct connect));
	      if (connect_ptr == NULL)
		{
		  printf("*** ERROR. Unable to allocate memory for struct");
		  printf(" connect\n");
		  printf("***        Program hbadd terminated with error.\n");
		  exit (1);
		}

	      /* If this is the first connection stored, then update
		 pointer to head of linked list */
	      if (*first_connect_ptr == NULL)
		*first_connect_ptr = connect_ptr;
	      
	      /* Otherwise, make previous connection point to the current
		 one */
	      else
		(*last_connect_ptr)->next_connect_ptr = connect_ptr;

	      /* Save the current record's pointer */
	      *last_connect_ptr = connect_ptr;

	      /* Store the atom numbers of the two atoms */
	      connect_ptr->atom_number1 = iatom;
	      connect_ptr->atom_number2 = jatom;
/* v.3.2--> */
	      connect_ptr->next_connect_ptr = NULL;
/* <--v.3.2 */

	      /* Increment number of connections stored */
	      (*nlinks)++;
	    }
	}
    }
}
/* v.3.3--> */
/***********************************************************************

check_residue_name  -  Check whether residue name is one of the standard
                       amino acids or a nucleic acid base

***********************************************************************/

int check_residue_name(char res_name[4])
{
  int iamino;
  static int namino = 0;
  int done, standard;

  static char *amino[] = {
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS",
    "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN",
    "ARG", "SER", "THR", "VAL", "TRP", "TYR",
    "  A", "  C", "  G", "  T", "  T", "XXX"
  };

  /* Initialise flag */
  standard = FALSE;

  /* If this is the first call to this routine, count the number of
     standard residue names */
  if (namino == 0)
    {
      /* Initialise flag */
      done = FALSE;

      /* Loop through the codes until reach the last one */
      while (done == FALSE)
	{
	  /* Check for the last code */
	  if (!strncmp(amino[namino],"XXX",3))
	    done = TRUE;

	  /* Otherwise, increment the count */
	  else
	    namino++;
	}
    }

  /* Check whether residue appears in the list */
  for (iamino = 0; iamino < namino && standard == FALSE; iamino++)
    if (!strncmp(res_name,amino[iamino],3))
      standard = TRUE;

  /* Return whether residue is one of the standard one */
  return(standard);
}
/* <--v.3.3 */
/***********************************************************************

read_pdbfile  -  Read through PDB file to get pick up all the HET groups

***********************************************************************/

void read_pdbfile(char pdb_name[FILENAME_LEN],
		  struct atom **first_atom_ptr,
		  struct residue **first_residue_ptr,
		  struct connect **first_connect_ptr,
/* v.3.1--> */
/*		  int *natoms,int *nhetgroups) */
		  int *natoms,int *nhetgroups,int *have_hydrogens)
/* <--v.3.1 */
{
  char line[LINE_LEN + 1];
  char chain, last_chain, last_resnum[6], resnum[6];
  char atmnum[6], resnam[4];
/* v.3.1.2--> */
  char last_resnam[4];
/* <--v.3.1.2 */

/* v.3.1.2--> */
/*  int first_atom_number, nconnect, nlinks; */
  int first_atom_number, last_atom_number, nconnect, nlinks;
/* <--v.3.1.2 */
  int atom_count, residue_count;
/* v.3.1--> */
  int nhydrogens;
/* <--v.3.1 */
/* v.3.1.2--> */
/* v.3.3--> */
/*  int keep_reading; */
  int keep_reading, standard, wanted;
/* <--v.3.3 */
/* <--v.3.1.2 */

  struct atom *atom_ptr, *last_atom_ptr;
  struct residue *residue_ptr, *last_residue_ptr;
  struct connect *last_connect_ptr;

/* v.3.1--> */
  /* Initialise variables */
  *have_hydrogens = FALSE;
  nhydrogens = 0;
/* <--v.3.1 */

  /* Open the pdb file */
  printf("  Opening PDB file [%s] ...\n",pdb_name);
  if ((fil_pdb = fopen(pdb_name,"r")) == NULL)
    {
      printf("\n*** Unable to open your PDB file [%s]\n",pdb_name);
      exit(1);
    }
    
  /* Prepare to start reading in the data from the PDB file */
  printf("\n");
  printf("    Searching for HET groups in PDB file ...\n");

  /* Initialise variables */
  atom_count = 0;
  residue_count = 0;
  *first_atom_ptr = NULL;
  *first_residue_ptr = NULL;
  *first_connect_ptr = NULL;
  first_atom_number = 0;
/* v.3.1.2--> */
  keep_reading = TRUE;
  last_atom_number = 999999;
/* <--v.3.1.2 */
  last_atom_ptr = NULL;
  last_residue_ptr = NULL;
  last_connect_ptr = NULL;
  *natoms = 0;
  *nhetgroups = 0;
  nconnect = 0;
  nlinks = 0;
/* v.3.1.2--> */
  last_resnam[0] = '\0';
/* <--v.3.1.2 */
  last_resnum[0] = '\0';
  last_chain = '\0';
/* v.3.3--> */
  wanted = FALSE;
/* <--v.3.3 */
  
  /* Search through the PDB file to find all the HET groups */
  while (fgets(line,LINE_LEN,fil_pdb) != NULL)
    {
      /* If this is an HETATM record, process it */
/* v.3.1.2--> */
/*      if (!strncmp(line,"HETATM",6)) */
/* v.3.3--> */
/*      if (!strncmp(line,"HETATM",6) && keep_reading == TRUE) */
      if ((!strncmp(line,"HETATM",6) || !strncmp(line,"ATOM  ",6)) &&
	  keep_reading == TRUE)
/* <--v.3.3 */
/* <--v.3.1.2 */
	{
	  printf("%s\n", line);
	  /* Get the residue name, number and chain-ID */
	  strncpy(resnam,line+17,3);
	  resnam[3] = '\0';
	  strncpy(resnum,line+22,5);
	  resnum[5] = '\0';
	  chain = line[21];
	  
/* v.3.3--> */
	  /* If this is a new residue type, check whether it is a standard
	     residue type */
	  if (strncmp(resnam,last_resnam,3))
	    {
	      /* Check whether this is a standard amino acid or nucleic acid
		 base */
	      standard = check_residue_name(resnam);

	      /* If standard residue, then don't want it */
	      if (standard == TRUE)
		wanted = FALSE;
	      else
		wanted = TRUE;
	    }
/* <--v.3.3 */

	  /* Check that this is not a water */
/*	  if (strncmp(resnam,"HOH",3) && strncmp(resnam,"WAT",3)) */
	  if (!strncmp(resnam,"HOH",3) || !strncmp(resnam,"WAT",3))
/* v.3.3--> */
	    wanted = FALSE;

	  if (wanted == TRUE)
/* <--v.3.3 */
	    {
	      /* Store the current atom's details */

	      /* Allocate memory for structure to hold current atom's
		 details */
	      atom_ptr
		= (struct atom*)malloc(sizeof(struct atom));
	      if (atom_ptr == NULL)
		{
		  printf("*** ERROR. Unable to allocate memory for struct");
		  printf(" atom\n");
		  printf("***        Program hbadd terminated with error.\n");
		  exit (1);
		}

	      /* If this is the first atom stored, then update
		 pointer to head of linked list */
	      if (*first_atom_ptr == NULL)
		*first_atom_ptr = atom_ptr;
	      
	      /* Otherwise, make previous residue point to the current
		 one */
	      else
		last_atom_ptr->next_atom_ptr = atom_ptr;

	      /* Save the current atom's pointer */
	      last_atom_ptr = atom_ptr;

	      /* Extract the relevant atom data from the record */
/* v.3.1.2--> */
/*	      sscanf(line+6,"%d",&atom_ptr->atom_number); */
	      strncpy(atmnum,line+6,5);
	      atmnum[5] = '\0';
	      atom_ptr->atom_number = atoi(atmnum);
/* <--v.3.1.2 */
	      strncpy(atom_ptr->atom_name,line+12,4);
	      atom_ptr->atom_name[4] = '\0';
	      sscanf(line+30," %f %f %f ",&atom_ptr->x,&atom_ptr->y,
		     &atom_ptr->z);

/* v.3.1--> */
	      /* If this is a hydrogen atom, increment count */
/* v.3.2--> */
/*	      if (atom_ptr->atom_name[1] == 'H')
		nhydrogens++; */
	      if (atom_ptr->atom_name[0] == 'H' ||
		  atom_ptr->atom_name[1] == 'H')
		{
		  nhydrogens++;
		  atom_ptr->hydrogen = TRUE;
		}
	      else
		atom_ptr->hydrogen = FALSE;
/* <--v.3.2 */
/* <--v.3.1 */

	      /* If this is the first stored atom, store its number */
	      if (first_atom_number == 0)
		first_atom_number = atom_ptr->atom_number;

	      /* Store number of last atom read in */
	      last_atom_number = atom_ptr->atom_number;
	      
	      /* Check whether this is a new residue */
/* v.3.1.2--> */
/*	      if (chain != last_chain || strncmp(resnum,last_resnum,5)) */
	      if (chain != last_chain || strncmp(resnum,last_resnum,5) ||
		  strncmp(resnam,last_resnam,3))
/* <--v.3.1.2 */
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
		      printf("***        Program hbadd terminated with");
		      printf(" error.\n");
		      exit (1);
		    }

		  /* Store the residue details */
		  strncpy(residue_ptr->res_name,resnam,3);
		  residue_ptr->res_name[3] = '\0';
		  strncpy(residue_ptr->res_num,resnum,5);
		  residue_ptr->res_num[5] = '\0';
		  residue_ptr->chain = chain;
		  residue_ptr->natoms = 0;
/* v.3.2--> */
		  residue_ptr->non_hydrogen = 0;
/* <--v.3.2 */
		  residue_ptr->first_atom_ptr = atom_ptr;
		  residue_ptr->dictionary_ptr = NULL; 
		  residue_ptr->next_residue_ptr = NULL;

		  /* If this is the first residue stored, then update
		     pointer to head of linked list */
		  if (*first_residue_ptr == NULL)
		    *first_residue_ptr = residue_ptr;

		  /* Otherwise, make previous residue point to the current
		     one */
		  else
		    last_residue_ptr->next_residue_ptr = residue_ptr;

		  /* Save the current residue's pointer */
		  last_residue_ptr = residue_ptr;

		  /* Save current residue */
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

	      /* Save pointer to residue in current atom's record */
	      atom_ptr->residue_ptr = residue_ptr;
/* v.3.1--> */
	      atom_ptr->dic_data_ptr = NULL;
/* <--v.3.1 */
	      atom_count++;

	      /* Increment count of atoms in this residue */
	      residue_ptr->natoms++;

/* v.3.2--> */
	      /* If this is not a hydrogen, then increment count of
		 non-hydrogen atoms */
	      if (atom_ptr->hydrogen == FALSE)
		residue_ptr->non_hydrogen++;
/* <--v.3.2 */
	    }
/* v.3.3--> */
	  /* Save current residue */
	  strncpy(last_resnam,resnam,3);
	  last_resnam[3] = '\0';
/* <--v.3.3 */
	}

/* v.3.1.2--> */
      /* If this is an ENDMDL record, then don't want to read in any more
         HET groups */
      else if (!strncmp(line,"ENDMDL",6))
	keep_reading = FALSE;
/* <--v.3.1.2 */

      /* If this is an CONECT record, process it */
      else if (!strncmp(line,"CONECT",6))
	{
	  /* If this is the first CONECT record encountered, then
	     show message */
	  if (nconnect == 0)
	    {
	      printf("\n");
	      printf("    Processing CONECT records ...\n");
	    }
	  nconnect++;

/* v.3.1.2--> */
	  /* Check whether first atom corresponds to any atom within
	     the range of atom numbers stored */
/* <--v.3.1.2 */

	  /* Process the current CONECT record */
/* v.3.1.2--> */
/*	  store_connections(line,first_atom_number,&(*first_connect_ptr),
			    &last_connect_ptr,&nlinks); */
	  store_connections(line,first_atom_number,last_atom_number,
			    &(*first_connect_ptr),&last_connect_ptr,&nlinks);
/* <--v.3.1.2 */
	}
    }

  /* Print counts of data read in */
  printf("\n");
  printf("      Number of HET residues read in        = %7d\n",residue_count);
  *nhetgroups = residue_count;
  printf("      Number of atom records stored         = %7d\n",atom_count);
  *natoms = atom_count;
  printf("      Number of CONECT details stored       = %7d\n",nlinks);
  printf("\n");

/* v.3.1--> */
  /* If have an adequate number of hydrogen atoms, then set flag
     accordingly */
  if (atom_count > 0)
    if ((float) nhydrogens / (float) atom_count > 0.2)
      *have_hydrogens = TRUE;
/* <--v.3.1 */

  /* Close the PDB file */
  fclose(fil_pdb);
}
/***********************************************************************

read_dictionary  -  Read through HET Group Dictionary to pick up all the
                    required HET group definitions

***********************************************************************/

void read_dictionary(char dictionary_name[],
		  struct residue *first_residue_ptr,
		  struct dictionary **first_dictionary_ptr,
/* v.3.1--> */
		  struct dictionary **last_dictionary_ptr,
/* <--v.3.1 */
		  struct dic_data **first_dic_data_ptr,
/* v.3.1--> */
		  struct dic_data **last_dic_data_ptr,
/* <--v.3.1 */
		  int *ndic_entries)
{
  int i, j, nmatch, wanted;
/* v.3.2--> */
  int hydrogen;
/* <--v.3.2 */

  char atom_name[5], res_name[4], line[LINE_LEN + 1], nbound_char[3];
/* v.4.0--> */
  char ctrl_M = 13;
/* <--v.4.0 */

/* v.3.1--> */
/*  struct dictionary *dictionary_ptr, *last_dictionary_ptr;
  struct dic_data *dic_data_ptr, *last_dic_data_ptr; */
  struct dictionary *dictionary_ptr;
  struct dic_data *dic_data_ptr;
/* <--v.3.1 */
  struct residue *residue_ptr;

  /* Initialise variables */
  dictionary_ptr = NULL;
  *first_dictionary_ptr = NULL;
/* v.3.1--> */
/*  last_dictionary_ptr = NULL; */
  *last_dictionary_ptr = NULL;
/* <--v.3.1 */
  dic_data_ptr = NULL;
  *first_dic_data_ptr = NULL;
/* v.3.1--> */
/*  last_dic_data_ptr = NULL; */
  *last_dic_data_ptr = NULL;
/* <--v.3.1 */
  *ndic_entries = 0;
  wanted = FALSE;
  
  /* Open the HET Group Dictionary */
  printf("  Opening HET Group Dictionary ...\n");
  if ((fil_dic = fopen(dictionary_name,"r")) == NULL)
    {
      printf("\n*** Unable to open HET Group Dictionary [%s]\n",
	     dictionary_name);
/* v.3.1--> */
/*      exit(1); */
      return;
/* <--v.3.1 */
    }
    
  /* Prepare to start reading in the data from the file */
  printf("\n");
  printf("    Scanning HET Group Dictionary ...\n");

  /* Read through the dictionary */
  while (fgets(line,LINE_LEN,fil_dic) != NULL)
    {
/* v.4.0--> */
      /* Strip out any ^M line-feed characters */
      for (i = 0; i < LINE_LEN - 1; i++)
	if (line[i] == ctrl_M)
	  {
	    line[i] = ' ';
	    line[i + 1] = '\0';
	  }
/* <--v.4.0 */

      /* If this is the start of a new residue, then check whether
	 it is one of the ones required */
      if (!strncmp(line,"RESIDUE",7))
	{
	  /* Extract residue name */ 
	  strncpy(res_name,line+10,3);
	  res_name[3] = '\0';

        /* Dave's Code */
          if( res_name[2] == ' ' && res_name[1] == ' ')
            {
                res_name[2] = res_name[0];
                res_name[1] = ' ';
                res_name[0] = ' ';
            }

          if( res_name[2] == ' ' )
            {
                res_name[2] = res_name[1];
                res_name[1] = res_name[0];
                res_name[0] = ' ';
            }







	  /* Initialise variables for search through stored residues */
	  residue_ptr = first_residue_ptr;
	  nmatch = 0;
	  wanted = FALSE;

	  /* Loop through all HET groups read in from the PDB file to
	     see if this is amongst them */
	  while (residue_ptr != NULL)
	    {
	      /* If the residue names match, then want to store this
		 dictionary entry */
	      if (!strncmp(res_name,residue_ptr->res_name,3))
		{
		  wanted = TRUE;

		  /* If this is the first match for this het group,
		     create an entry in the dictionary linked list */
		  if (nmatch == 0)
		    {
		      dictionary_ptr =
			(struct dictionary*)malloc(sizeof(struct dictionary));
		      if (dictionary_ptr == NULL)
			{
			  printf("*** ERROR. Unable to allocate memory for");
			  printf(" struct dictionary\n");
			  printf("***        Program hbadd terminated with");
			  printf(" error.\n");
			  exit (1);
			}

		      /* Store the details for this entry */
		      strncpy(dictionary_ptr->res_name,res_name,3);
		      dictionary_ptr->res_name[3] = '\0';
		      dictionary_ptr->natoms = 0;
		      dictionary_ptr->nmatch = 0;
/* v.3.2--> */
		      dictionary_ptr->non_hydrogen = 0;
		      dictionary_ptr->have_hydrogens = TRUE;
		      dictionary_ptr->fake_entry = FALSE;
/* <--v.3.2 */
		      dictionary_ptr->residue_ptr = NULL;
		      dictionary_ptr->first_dic_data_ptr = NULL;
		      dictionary_ptr->next_dictionary_ptr = NULL;

		      /* If this is the first entry stored, then update
			 pointer to head of linked list */
		      if (*first_dictionary_ptr == NULL)
			*first_dictionary_ptr = dictionary_ptr;
		      
		      /* Otherwise, make previous entry point to the current
			 one */
		      else
/* v.3.1--> */
/*			last_dictionary_ptr->next_dictionary_ptr
			  = dictionary_ptr; */
			(*last_dictionary_ptr)->next_dictionary_ptr
			  = dictionary_ptr;
/* <--v.3.1 */
		      
		      /* Save the current residue's pointer */
/* v.3.1--> */
/*		      last_dictionary_ptr = dictionary_ptr; */
		      *last_dictionary_ptr = dictionary_ptr;
/* <--v.3.1 */
		      (*ndic_entries)++;
		    }

		  /* Add link from residue to dictionary entry */
		  residue_ptr->dictionary_ptr = dictionary_ptr;

		  /* Increment count of matches */
		  nmatch++;
		}

	      /* Get pointer to the next residue in the list */
	      residue_ptr = residue_ptr->next_residue_ptr;
	    }

	  /* Show number of matches found */
	  if (wanted == TRUE)
	    {
	      printf("      Found match for: %s ",res_name);
	      if (nmatch > 1)
		printf("x %d\n",nmatch);
	      else
		printf("\n");
	    }
	}

      /* If this is a CONECT record for a wanted residue, then store
	 it */
      else if (!strncmp(line,"CONECT",6) && wanted == TRUE)
	{
	  /* Extract the atom name */ 
	  strncpy(atom_name,line+11,4);
	  atom_name[4] = '\0';

/* v.3.2--> */
	  /* Determine whether this is a hydrogen atom */
	  hydrogen = FALSE;
	  if (atom_name[0] == 'H' || atom_name[1] == 'H')
	    hydrogen = TRUE;
/* <--v.3.2 */

	  /* Create entry in the atom-data linked list */
	  dic_data_ptr
	    = (struct dic_data*)malloc(sizeof(struct dic_data));
	  if (dic_data_ptr == NULL)
	    {
	      printf("*** ERROR. Unable to allocate memory for");
	      printf(" struct dic_data\n");
	      printf("***        Program hbadd terminated with");
	      printf(" error.\n");
	      exit (1);
	    }

	  /* Remove double-quotes("), present in some atom names  */
	  for (i = 0; i < 4; i++)
	    if (atom_name[i] == '"' || atom_name[i] == '\0'
		|| atom_name[i] == '\n')
	      atom_name[i] = ' ';

	  /* Store the details for this entry */
	  strncpy(dic_data_ptr->atom_name,atom_name,4);
	  dic_data_ptr->atom_name[4] = '\0';
/* v.3.2--> */
	  dic_data_ptr->hydrogen = hydrogen;
	  strncpy(dic_data_ptr->output_atom_name,atom_name,4);
	  dic_data_ptr->output_atom_name[4] = '\0';
/* <--v.3.2 */
	  dic_data_ptr->number = dictionary_ptr->natoms;

	  /* Number of connections */ 
	  if ((int) strlen(line) >= 19)
	    {
	      strncpy(nbound_char, line+18, 2);
	      nbound_char[2] = '\0';
	      dic_data_ptr->nbound = atoi(nbound_char);
	    }
	  else
	    dic_data_ptr->nbound = 0;
	 
	  /* Connectivity information */
	  for (i = 0; i < dic_data_ptr->nbound; i++)
	    {
	      strncpy(dic_data_ptr->boundto[i], line+(20+(i*5)), 5);      

	      /* Remove carriage returns and other stray characters */
	      for (j = 0; j < 4; j++)
		if ((dic_data_ptr->boundto[i][j] == '\n') 
/* v.4.0--> */
		    || (dic_data_ptr->boundto[i][j] == ctrl_M)
/* <--v.4.0 */
		    || (dic_data_ptr->boundto[i][j] == '"')
		    || (dic_data_ptr->boundto[i][j] == '\0'))
		  dic_data_ptr->boundto[i][j] = ' ';
	      dic_data_ptr->boundto[i][4] = '\0';
	    }

	  /* Initialize other variables */
	  for (i = 0; i < MAXCONNECT; i++)
	    {
	      dic_data_ptr->connect_atom_ptr[i] = NULL;
	    }
	  
	  /* If this is the first entry stored, then update
	     pointer to head of linked list */
	  if (*first_dic_data_ptr == NULL)
	    *first_dic_data_ptr = dic_data_ptr;
	  
	  /* Otherwise, make previous entry point to the current
	     one */
	  else
/* v.3.1--> */
/*	    last_dic_data_ptr->next_dic_data_ptr = dic_data_ptr; */
	    (*last_dic_data_ptr)->next_dic_data_ptr = dic_data_ptr;
/* <--v.3.1 */
	  
	  /* Save the current atom's pointer */
/* v.3.1--> */
/*	  last_dic_data_ptr = dic_data_ptr; */
	  *last_dic_data_ptr = dic_data_ptr;
/* <--v.3.1 */

	  /* If this is the first CONECT record for this entry, then
	     store pointer to the first atom */
	  if (dictionary_ptr->natoms == 0)
	    dictionary_ptr->first_dic_data_ptr = dic_data_ptr;

	  /* Increment atom-count for this dictionary entry */
	  dictionary_ptr->natoms++;
/* v.3.2--> */
	  if (hydrogen == FALSE)
	    dictionary_ptr->non_hydrogen++;
/* <--v.3.2 */
	}
    }

  /* Print counts of data read in */
  printf("\n");
  printf("      Number of dictionary entries read in  = %7d\n",*ndic_entries);
  printf("\n");

  /* Close the PDB file */
  fclose(fil_dic);
}
/***********************************************************************

dic_connectivities  -  Scan through the connectivity info for each stored
                       CONEC record from the HET Groups Dictionary and
                       match up the other atoms that each atom is bonded to

***********************************************************************/

void dic_connectivities(struct dictionary *first_dictionary_ptr)
{
  int found, iatom, iconnect, jatom, natoms;

  struct dictionary *dictionary_ptr;
  struct dic_data *dic_data_ptr, *other_dic_data_ptr;

  /* Initialise HET group pointer to start of linked list */
  dictionary_ptr = first_dictionary_ptr;

  /* Loop through all the HET group residues read in from the dictionary */
  while (dictionary_ptr != NULL)
    {
      /* Get the number of atoms making up this HET group and the pointer
	 to the first one */
      dic_data_ptr = dictionary_ptr->first_dic_data_ptr;
      natoms = dictionary_ptr->natoms;

      /* Loop through each of the atoms in turn and process */
      for (iatom = 0; iatom < natoms && dic_data_ptr != NULL; iatom++)
	{
	  /* Loop through the list of atoms bonded to this one */
	  for (iconnect = 0; iconnect < dic_data_ptr->nbound; iconnect++)
	    {
	      /* Initialise search list to first of the atoms in this
		 HET group */
	      other_dic_data_ptr = dictionary_ptr->first_dic_data_ptr;
	      found = FALSE;

	      /* Loop through all the atoms in the HET group until find
		 the atom we're after */
	      for (jatom = 0; jatom < natoms && found == FALSE; jatom++)
		{
		  /* Check whether this is the atom we want */
		  if (!strncmp(dic_data_ptr->boundto[iconnect],
			       other_dic_data_ptr->atom_name, 4))
		    {
		      /* Have a match, so store the pointer matched atom */
		      dic_data_ptr->connect_atom_ptr[iconnect]
			= other_dic_data_ptr;
		      dic_data_ptr->boundto_num[iconnect]
			= other_dic_data_ptr->number;
		      found = TRUE;
		    }
		  /* Get pointer to the next atom in the list */
		  other_dic_data_ptr = other_dic_data_ptr->next_dic_data_ptr;
		}

	      /* If the atom not found, then show warning */
	      if (found == FALSE)
		{
		  printf("*** Warning. Error in HET Group Dictionary");
		  printf(" for residue-type [%s]\n",
			 dictionary_ptr->res_name);
		  printf("***          Connect atom %d [%s] undefined",
			 iconnect + 1,dic_data_ptr->boundto[iconnect]);
		  printf(" for atom type [%s]\n",dic_data_ptr->atom_name);
		}
	    }

	  /* Get pointer to the next HET group atom in the list */
	  dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
	}

      /* Get pointer to the next HET group entry in the list */
      dictionary_ptr = dictionary_ptr->next_dictionary_ptr;
    }
}
/***********************************************************************

print_missing_hets  -  Print this HET group as missing from the HET
                       Group Dictionary, unless a warning has already
		       been printed for this residue type

***********************************************************************/

void print_missing_hets(struct residue *residue_ptr,
			struct residue *first_residue_ptr)
{
  int print_warning;

  struct residue *residue_comp_ptr;

  /* Initialise search of residues already done to check that
     we haven't already shown a warning message for this type */
  residue_comp_ptr = first_residue_ptr;
  print_warning = TRUE;
	  
  /* Loop through all the HET group residues already processed */
  while (residue_comp_ptr != NULL && print_warning == TRUE)
    {
      /* If this is the current residue, then end search */
      if (residue_comp_ptr == residue_ptr)
	residue_comp_ptr = NULL;

      /* Otherwise, check whether the residue name is the same */
      else
	{
	  /* If residue name the same, then don't need to print
	     the warning again */
	  if (!strncmp(residue_comp_ptr->res_name,
		       residue_ptr->res_name,3))
	    print_warning = FALSE;

	  /* Get pointer to the next residue in the list */
	  residue_comp_ptr = residue_comp_ptr->next_residue_ptr;
	}
    }
  
  /* Print the warning message, if required */
  if (print_warning == TRUE)
    {
      printf("*** Warning: Residue-type [%s] not found in HET Group",
	     residue_ptr->res_name);
      printf(" Dictionary\n");
    }
}
/***********************************************************************

check_stored_conecs  -  Loop through all the stored connections from the
                        CONECT records in the PDB file and see if there
			is one for the two atoms in question

***********************************************************************/

int check_stored_conecs(struct connect *first_connect_ptr,
			int atom_number1,int atom_number2)
{
  int connected;

  struct connect *connect_ptr;

  /* Initialise variables */
  connected = FALSE;

  /* Initialise connect-pointer to start of linked list */
  connect_ptr = first_connect_ptr;

  /* Loop through all the HET group connects */
  while (connect_ptr != NULL && connected == FALSE)
    {
      /* Check whether this bond corresponds to the two atoms in question */
      if ((atom_number1 == connect_ptr->atom_number1 &&
	   atom_number2 == connect_ptr->atom_number2) ||
	  (atom_number1 == connect_ptr->atom_number2 &&
	   atom_number2 == connect_ptr->atom_number1))
	connected = TRUE;

      /* Get pointer to the next connect in the list */
      connect_ptr = connect_ptr->next_connect_ptr;
    }

  /* Return the answer */
  return(connected);
}
/* v.3.1--> */
/***********************************************************************

generate_connectivities  -  Generate all the atom-atom connectivities
                            for the fake dictionary entry using the
                            atom-atom distances and CONECT records in
                            the PDB 

***********************************************************************/

void generate_connectivities(struct residue *residue_ptr,
                             struct dictionary *dictionary_ptr,
                             struct connect *first_connect_ptr)
{
  int iatom, iconnect, jatom, natoms;
  int connected, hydrogen1, hydrogen2;

  float connect_dist2, distance2, x1, x2, y1, y2, z1, z2;

  struct atom *atom_ptr, *other_atom_ptr;
  struct dic_data *dic_data_ptr, *other_dic_data_ptr;

  /* Initialise pointer to first atom in the dictionary list of
     atoms making up this HET group */
  dic_data_ptr = dictionary_ptr->first_dic_data_ptr;
  natoms = dictionary_ptr->natoms;

  /* Loop through each of the atoms in turn and check each of its bonds */
  for (iatom = 0; iatom < natoms && dic_data_ptr != NULL; iatom++)
    {
      /* Get the pointer to the atom coordinates */
      atom_ptr = dic_data_ptr->atom_ptr;
      if (atom_ptr != NULL)
        {
	  /* Get the atom's coordinates */
	  x1 = atom_ptr->x;
	  y1 = atom_ptr->y;
	  z1 = atom_ptr->z;

	  /* Check if this atom is a hydrogen */
	  hydrogen1 = FALSE;
	  if (atom_ptr->atom_name[0] == 'H' ||
	      atom_ptr->atom_name[1] == 'H')
	    hydrogen1 = TRUE;

          /* Loop through all this residue's atoms to find all atoms
             bonded to this one */

          /* Get pointer to first of residue's atoms */
          other_atom_ptr = residue_ptr->first_atom_ptr;

          /* Loop through the residue's atoms */
          for (jatom = 0; jatom < natoms && other_atom_ptr != NULL;
               jatom++)
            {
	      /* Process if not the same atom */
	      if (other_atom_ptr != atom_ptr)
		{
		  /* Check whether the two atoms are covalently
		     bonded */

		  /* Get the second atom's coordinates */
		  x2 = other_atom_ptr->x;
		  y2 = other_atom_ptr->y;
		  z2 = other_atom_ptr->z;
		  
		  /* Check if this atom is a hydrogen */
		  hydrogen2 = FALSE;
		  if (other_atom_ptr->atom_name[0] == 'H' ||
		      other_atom_ptr->atom_name[1] == 'H')
		    hydrogen2 = TRUE;

		  /* Calculate the distance between them */
		  distance2 = (x1 - x2) * (x1 - x2)
		    + (y1 - y2) * (y1 - y2)
		    + (z1 - z2) * (z1 - z2);

		  /* Get the other atom's dictionary entry */
		  other_dic_data_ptr = other_atom_ptr->dic_data_ptr;

		  /* If either atom is a hydrogen, then reduce
		     connect distance */
		  if (hydrogen1 == TRUE || hydrogen2 == TRUE)
		    connect_dist2 = HATOM_DIST2;
		  else
		    connect_dist2 = COVAL_DIST2;

		  /* If bond distance is not within allowed tolerances,
		     then check the connectivities as read in from
		     the CONECT records in the PDB file */
		  connected = FALSE;
		  if (distance2 > connect_dist2)
		    {
		      /* Check for a CONECT record between these two
			 atoms */
		      connected
			= check_stored_conecs(first_connect_ptr,
					      atom_ptr->atom_number,
					      other_atom_ptr->atom_number);
		    }
		  else
		    connected = TRUE;

		  /* If atoms are connected, then add this bond to the
		     dictionary entry */
		  if (connected == TRUE)
		    {
		      /* Get next available slot for connectivity data */
		      iconnect = dic_data_ptr->nbound;

		      /* Check that maximum number of possible connections
			 not exceeded */
		      if (iconnect > MAXCONNECT - 1)
			{
			  printf("*** Warning. Maximum number of bonds");
			  printf(" exceeded for atom [%s %s %s %c] %d",
				 atom_ptr->atom_name,
				 atom_ptr->residue_ptr->res_name,
				 atom_ptr->residue_ptr->res_num,
				 atom_ptr->residue_ptr->chain,
				 iconnect + 1);
			}

		      /* Otherwise, store this bond */
		      else
			{
			  /* Store the details of the other atom */
			  dic_data_ptr->connect_atom_ptr[iconnect]
			    = other_dic_data_ptr;
			  dic_data_ptr->boundto_num[iconnect]
			    = other_atom_ptr->atom_number;
			  strncpy(dic_data_ptr->boundto[iconnect],
				  other_atom_ptr->atom_name,4);
			  dic_data_ptr->boundto[iconnect][4] = '\0';
		      
			  /* Increment count of bonded connections */
			  dic_data_ptr->nbound++;
			}
		    }
		}

              /* Get pointer to next atom */
              other_atom_ptr = other_atom_ptr->next_atom_ptr;
            }
	}
      
      /* Get pointer to the next HET group atom in the list */
      dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
    }
}
/***********************************************************************

generate_fake_entry  -  Check that the atoms in the HET group from the
                        PDB file match those in the dictionary definition

***********************************************************************/

void generate_fake_entry(struct residue *residue_ptr,
                         struct residue *first_residue_ptr,
                         struct dictionary **first_dictionary_ptr,
                         struct dictionary **last_dictionary_ptr,
                         int *ndic_entries,
/* v.3.2--> */
/*                         struct dic_data **first_dic_data_ptr, */
/* <--v.3.2 */
                         struct dic_data **last_dic_data_ptr,
                         struct connect *first_connect_ptr)
{
  int iatom, iconnect, natoms;
  int already_got_it;

  struct atom *atom_ptr;
/* v.3.2--> */
/*  struct dictionary *dictionary_ptr; */
  struct dictionary *dictionary_ptr, *save_dictionary_ptr;
/* <--v.3.2 */
  struct dic_data *dic_data_ptr;
  struct residue *check_residue_ptr;

  /* Show message */
  printf("Generating fake HBPLUS entry from PDB file for residue type [%s]\n",
	 residue_ptr->res_name);

  /* Check that don't already have a dictionary entry for this residue
     type */
  already_got_it = FALSE;
  check_residue_ptr = first_residue_ptr;
/* v.3.2--> */
  save_dictionary_ptr = NULL;
/* <--v.3.2 */

  /* Loop through all the residues */
/* v.3.2--> */
/*  while (check_residue_ptr != NULL) */
  while (check_residue_ptr != NULL && already_got_it == FALSE)
/* <--v.3.2 */
    {
      /* If residue-type matches and have a dictionary pointer, then
         must already have an entry in the dictionary */
      if (!strncmp(check_residue_ptr->res_name,residue_ptr->res_name,3) &&
	  check_residue_ptr->dictionary_ptr != NULL)
/* v.3.2--> */
	{
	  /* Store the pointer */
	  save_dictionary_ptr = check_residue_ptr->dictionary_ptr;
/* <--v.3.2 */
	  already_got_it = TRUE;
/* v.3.2--> */
	}
/* <--v.3.2 */

      /* Get pointer to the next residue */
      check_residue_ptr = check_residue_ptr->next_residue_ptr;
    }

  /* If already have this residue type, then don't need to create another
     entry */
  if (already_got_it == TRUE)
/* v.3.2--> */
/*    return; */
    {
      /* Create the appropriate link */
      residue_ptr->dictionary_ptr = save_dictionary_ptr;
      return;
    }
/* <--v.3.2 */

  /* Create a new (fake) dictionary entry corresponding to the unmatched
     residue */
  dictionary_ptr
    =(struct dictionary*)malloc(sizeof(struct dictionary));
  if (dictionary_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct dictionary\n");
      printf("***        Program hbadd terminated with");
      printf(" error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  strncpy(dictionary_ptr->res_name,residue_ptr->res_name,3);
  dictionary_ptr->res_name[3] = '\0';
  dictionary_ptr->natoms = residue_ptr->natoms;
  dictionary_ptr->nmatch = 0;
/* v.3.2--> */
  dictionary_ptr->non_hydrogen = residue_ptr->non_hydrogen;
/* <--v.3.2 */
  dictionary_ptr->have_hydrogens = FALSE;
  dictionary_ptr->fake_entry = TRUE;
  dictionary_ptr->residue_ptr = residue_ptr;
  dictionary_ptr->first_dic_data_ptr = NULL;
  dictionary_ptr->next_dictionary_ptr = NULL;

  /* If this is the first entry stored, then update pointer to head 
     of linked list */
  if (*first_dictionary_ptr == NULL)
    *first_dictionary_ptr = dictionary_ptr;
		      
  /* Otherwise, make previous entry point to the current one */
  else
    (*last_dictionary_ptr)->next_dictionary_ptr = dictionary_ptr;
		      
  /* Save the current residue's pointer */
  *last_dictionary_ptr = dictionary_ptr;
  (*ndic_entries)++;

  /* Add link from residue to dictionary entry */
  residue_ptr->dictionary_ptr = dictionary_ptr;

  /* Initialise search to first PDB atom in the current residue */
  atom_ptr = residue_ptr->first_atom_ptr;
  natoms = residue_ptr->natoms;

  /* Loop through all the atoms, creating a dictionary entry for 
     each one */
  for (iatom = 0; iatom < natoms && atom_ptr != NULL; iatom++)
    {
      /* Create entry in the atom-data linked list for this atom */
      dic_data_ptr
	= (struct dic_data*)malloc(sizeof(struct dic_data));
      if (dic_data_ptr == NULL)
	{
	  printf("*** ERROR. Unable to allocate memory for");
	  printf(" struct dic_data\n");
	  printf("***        Program hbadd terminated with");
	  printf(" error.\n");
	  exit (1);
	}

      /* Store the details for this entry */
      strncpy(dic_data_ptr->atom_name,atom_ptr->atom_name,4);
      dic_data_ptr->atom_name[4] = '\0';
/* v.3.2--> */
      strncpy(dic_data_ptr->output_atom_name,atom_ptr->atom_name,4);
      dic_data_ptr->output_atom_name[4] = '\0';
/* <--v.3.2 */
      dic_data_ptr->number = atom_ptr->atom_number;
      dic_data_ptr->atom_ptr = atom_ptr;
      dic_data_ptr->next_dic_data_ptr = NULL;   

/* v.3.2--> */
      /* If atom is a hydrogen, note that have hydrogens in the residue */
      if (atom_ptr->atom_name[0] == 'H' || atom_ptr->atom_name[1] == 'H')
	dictionary_ptr->have_hydrogens = FALSE;
/* <--v.3.2 */

      /* Store pointer the other way */
      atom_ptr->dic_data_ptr = dic_data_ptr;

      /* Initialise the connectivity information */
      for (iconnect = 0; iconnect < MAXCONNECT; iconnect++)
	{
	  dic_data_ptr->boundto[iconnect][0] = '\0';
	  dic_data_ptr->boundto_num[iconnect] = 0;
	  dic_data_ptr->connect_atom_ptr[iconnect] = NULL;
	}

      /* Initialise number of connections */ 
      dic_data_ptr->nbound = 0;

      /* If this is the first atom, then store pointer to it
	 in the dictionary entry */
      if (dictionary_ptr->first_dic_data_ptr == NULL)
	dictionary_ptr->first_dic_data_ptr = dic_data_ptr;

      /* Otherwise, make previous entry point to the current
	 one */
      else
	(*last_dic_data_ptr)->next_dic_data_ptr = dic_data_ptr;

      /* Save pointer to the current entry */
      *last_dic_data_ptr = dic_data_ptr;

      /* Get pointer to the next coordinate record in the list */
      atom_ptr = atom_ptr->next_atom_ptr;
    }

  /* Calculate all the connectivities using distances and CONECT data */
  generate_connectivities(residue_ptr,dictionary_ptr,first_connect_ptr);
}
/* <--v.3.1 */
/* v.3.2--> */
/***********************************************************************

calc_covalent_connectivities  -  Calculate the covalent connectivities
                                 of all the atoms in the given residue

***********************************************************************/

void calc_covalent_connectivities(struct residue *residue_ptr,
				  struct connect *first_connect_ptr)
{
  int connected, iatom, iconnect, jatom, natoms;

  float distance2, x1, x2, y1, y2, z1, z2;

  struct atom *atom_ptr, *other_atom_ptr;

  /* Initialise pointer to first atom in the residue */
  atom_ptr = residue_ptr->first_atom_ptr;
  natoms = residue_ptr->natoms;

  /* Loop through each of the atoms in turn and check each of its bonds */
  for (iatom = 0; iatom < natoms && atom_ptr != NULL; iatom++)
    {
      /* Get this atom's coordinates */
      x1 = atom_ptr->x;
      y1 = atom_ptr->y;
      z1 = atom_ptr->z;

      /* Initialise all its connectivity pointers */
      for (iconnect = 0; iconnect < MAXCONNECT; iconnect++)
	atom_ptr->connect_atom_ptr[iconnect] = NULL;

      /* Re-initialise connections counter */
      iconnect = 0;

      /* Loop over all the other atoms */
      other_atom_ptr = residue_ptr->first_atom_ptr;

      /* Loop through each of the atoms in turn and check each of
	 its bonds */
      for (jatom = 0; jatom < natoms && other_atom_ptr != NULL; jatom++)
	{
	  /* Proceed if not the same as the first atom */
	  if (other_atom_ptr != atom_ptr)
	    {
	      /* Get this atom's coordinates */
	      x2 = other_atom_ptr->x;
	      y2 = other_atom_ptr->y;
	      z2 = other_atom_ptr->z;

	      /* Check whether the two atoms are covalently bonded in
		 the structure */

	      /* Calculate the distance between them */
	      distance2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)
		+ (z1 - z2) * (z1 - z2);

	      /* If bond distance is not within allowed tolerances, then
		 check the connectivities as read in from the CONECT
		 records in the PDB file */
	      connected = FALSE;
	      if (distance2 > COVAL_DIST2)
		{
		  /* Check for a CONECT record between these two
		     atoms */
		  connected
		    = check_stored_conecs(first_connect_ptr,
					  atom_ptr->atom_number,
					  other_atom_ptr->atom_number);
		}
	      else
		connected = TRUE;

	      /* If atoms are bonded, add to the first atom's list of
		 connections */
	      if (connected == TRUE)
		{
		  if (iconnect > MAXCONNECT - 1)
		    {
		      printf("*** Warning. Maximum number of bonds");
		      printf(" exceeded for atom [%s %s %s %c] %d\n",
			     atom_ptr->atom_name,
			     atom_ptr->residue_ptr->res_name,
			     atom_ptr->residue_ptr->res_num,
			     atom_ptr->residue_ptr->chain,
			     iconnect + 1);
		    }

		  /* Otherwise, store this bond */
		  else
		    {
		      /* Store pointer to bonded atom */
		      atom_ptr->connect_atom_ptr[iconnect] = other_atom_ptr;
		      iconnect++;
		    }
		}
	    }

	  /* Get pointer to the next atom in the list */
	  other_atom_ptr = other_atom_ptr->next_atom_ptr;
	}
      
      /* Get pointer to the next atom in the list */
      atom_ptr = atom_ptr->next_atom_ptr;
    }
}
/***********************************************************************

create_nodes  -  Create nodes of the graph, each corresponding to a
                 pairing between similar atoms of the residue and from
		 the Het Group Dictionary entry

***********************************************************************/

void create_nodes(struct residue *residue_ptr,
		  struct node **first_node_ptr,int *nnodes)
{
  char atom_name[5], dic_atom_name[5];

  int iatom, jatom, ndic_atoms, natoms;
  int hydrogen;

  struct atom *atom_ptr;
  struct dictionary *dictionary_ptr;
  struct dic_data *dic_data_ptr;
  struct node *last_node_ptr, *node_ptr;

  /* Initialise variables */
  last_node_ptr = NULL;
  *first_node_ptr = NULL;
  *nnodes = 0;

  /* Get dictionary entry corresponding to this residue */
  dictionary_ptr = residue_ptr->dictionary_ptr;

  /* Initialise pointer to first atom in the residue */
  atom_ptr = residue_ptr->first_atom_ptr;
  natoms = residue_ptr->natoms;

  /* Loop through each of the atoms in turn */
  for (iatom = 0; iatom < natoms && atom_ptr != NULL; iatom++)
    {
      /* Get the atom name */
      strncpy(atom_name,atom_ptr->atom_name,4);
      atom_name[4] = '\0';

      /* Check if this atom is a hydrogen */
      hydrogen = FALSE;
      if (atom_name[0] == 'H' || atom_name[1] == 'H')
	hydrogen = TRUE;

      /* Skip atom if it is a hydrogen */
      if (hydrogen == FALSE)
	{
	  /* Get the first of the atoms corresponding to the dictionary
	     entry */
	  dic_data_ptr = dictionary_ptr->first_dic_data_ptr;
	  ndic_atoms = dictionary_ptr->natoms;

	  /* Loop through each of the atoms in the dictionary entry */
	  for (jatom = 0; jatom < ndic_atoms && dic_data_ptr != NULL;
	       jatom++)
	    {
	      /* Get this atom's name */
	      strncpy(dic_atom_name,dic_data_ptr->atom_name,4);
	      dic_atom_name[4] = '\0';

	      /* If the two atoms are of the same type, create a
		 graph-node to correspond to their potentially
		 being equivalent */
	      if (atom_name[1] == dic_atom_name[1])
		{
		  /* Create a node equivalencing these two atoms */
		  node_ptr
		    = (struct node*)malloc(sizeof(struct node));
		  if (node_ptr == NULL)
		    {
		      printf("*** ERROR. Unable to allocate memory");
		      printf(" for struct node\n");
		      printf("***        Program hbadd terminated with");
		      printf(" error.\n");
		      exit (1);
		    }

		  /* If this is the first node, then update pointer
		     to head of linked list */
		  if (*first_node_ptr == NULL)
		    *first_node_ptr = node_ptr;

		  /* Otherwise, make previous node point to the current
		     one */
		  else
		    last_node_ptr->next_node_ptr = node_ptr;

		  /* Save the current record's pointer */
		  last_node_ptr = node_ptr;

		  /* Store the two potentially equivalent atoms */
		  node_ptr->atom_ptr = atom_ptr;
		  node_ptr->dic_data_ptr = dic_data_ptr;

		  /* Fill in the remaining data items */
		  node_ptr->node_number = *nnodes;
		  node_ptr->checked = FALSE;
		  node_ptr->in_best_clique = FALSE;
		  node_ptr->next_node_ptr = NULL;

		  /* Increment count of graph nodes created */
		  (*nnodes)++;
		}

	      /* Get pointer to the next HET group atom in the list */
	      dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
	    }
	}
      
      /* Get pointer to the next atom in the list */
      atom_ptr = atom_ptr->next_atom_ptr;
    }
}
/***********************************************************************

calc_node_links  -  Calculate all the links between nodes - each link,
                    between node (Si,Pj) and node (Sk,Pl), represents
                    possibly-equivalent atoms that also share possibly
		    equivalent bonds (ie bonds Si-Sk and Pj-Pl exist)

***********************************************************************/

void calc_node_links(struct node *first_node_ptr,
		     struct node_link **first_node_link_ptr,
		     int *nnode_links)
{
  int done, have_bond;
  int iconnect;

  struct atom *atom_ptr, *other_atom_ptr;
  struct dic_data *dic_data_ptr, *other_dic_data_ptr;
  struct node *node_ptr, *other_node_ptr;
  struct node_link *last_node_link_ptr, *node_link_ptr;


  /* Initialise variables */
  last_node_link_ptr = NULL;
  *first_node_link_ptr = NULL;
  *nnode_links = 0;

  /* Get pointer to first node */
  node_ptr = first_node_ptr;

  /* Loop through all the nodes */
  while (node_ptr != NULL)
    {
      /* Get the atom pointers to the corresponding atoms */
      atom_ptr = node_ptr->atom_ptr;
      dic_data_ptr = node_ptr->dic_data_ptr;

      /* Get pointer to first node */
      other_node_ptr = node_ptr->next_node_ptr;

      /* Loop through all other nodes */
      while (other_node_ptr != NULL)
	{
	  /* Get the atom pointers to the corresponding atoms */
	  other_atom_ptr = other_node_ptr->atom_ptr;
	  other_dic_data_ptr = other_node_ptr->dic_data_ptr;

	  /* Skip if any of the atoms in the two nodes are the
	     same */
	  if (other_atom_ptr != atom_ptr &&
	      other_dic_data_ptr != dic_data_ptr)
	    {
	      /* Check whether have a bond between the two atoms in
		 the residue */
	      have_bond = FALSE;

	      /* Loop through first atom's connections to see if
		 other atom is amongst them */
	      done = FALSE;
	      for (iconnect = 0; iconnect < MAXCONNECT && done == FALSE;
		   iconnect++)
		{
		  /* Check for end of the connections */
		  if (atom_ptr->connect_atom_ptr[iconnect] == NULL)
		    done = TRUE;

		  /* Otherwise, check for the second atom */
		  else if (atom_ptr->connect_atom_ptr[iconnect]
			   == other_atom_ptr)
		    {
		      have_bond = TRUE;
		      done = TRUE;
		    }
		}

	      /* If have a bond, then repeat for the two atoms from
		 the dictionary */
	      if (have_bond == TRUE)
		{
		  /* Initialise flag */
		  have_bond = FALSE;

		  /* Loop through first dictionary atom's connections
		     to see if other ictionary atom is amongst them */
		  done = FALSE;
		  for (iconnect = 0; iconnect < MAXCONNECT && done == FALSE;
		       iconnect++)
		    {
		      /* Check for end of the connections */
		      if (dic_data_ptr->connect_atom_ptr[iconnect] == NULL)
			done = TRUE;

		      /* Otherwise, check for the second atom */
		      else if (dic_data_ptr->connect_atom_ptr[iconnect]
			       == other_dic_data_ptr)
			{
			  have_bond = TRUE;
			  done = TRUE;
			}
		    }

		  /* If have a bond here, too, then can connect the
		     two nodes */
		  if (have_bond == TRUE)
		    {
		      /* Create a connection between the two nodes */
		      node_link_ptr
			= (struct node_link*)malloc(sizeof(struct node_link));
		      if (node_link_ptr == NULL)
			{
			  printf("*** ERROR. Unable to allocate memory");
			  printf(" for struct node_link\n");
			  printf("***        Program hbadd terminated with");
			  printf(" error.\n");
			  exit (1);
			}

		      /* If this is the first link, then update pointer
			 to head of linked list */
		      if (*first_node_link_ptr == NULL)
			*first_node_link_ptr = node_link_ptr;

		      /* Otherwise, make previous node_link point to
			 the current one */
		      else
			last_node_link_ptr->next_node_link_ptr
			  = node_link_ptr;

		      /* Save the current record's pointer */
		      last_node_link_ptr = node_link_ptr;

		      /* Store the link details */
		      node_link_ptr->node1_ptr = node_ptr;
		      node_link_ptr->node2_ptr = other_node_ptr;
		      node_link_ptr->link_number = *nnode_links;
		      node_link_ptr->next_node_link_ptr = NULL;
		      (*nnode_links)++;
		    }
		}
	    }

	  /* Get pointer to the next node */
	  other_node_ptr = other_node_ptr->next_node_ptr;
	}

      /* Get pointer to the next node */
      node_ptr = node_ptr->next_node_ptr;
    }
}
/***********************************************************************

initialise_atoms  -  Initialise all the atom flags for the given residue

***********************************************************************/

void initialise_atoms(struct residue *residue_ptr)
{
  int iatom, natoms;

  struct atom *atom_ptr;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = residue_ptr->first_atom_ptr;
  natoms = residue_ptr->natoms;

  /* Loop through all the atoms to initialise the checked flag */
  for (iatom = 0; iatom < natoms && atom_ptr != NULL; iatom++)
    {
      /* Initialise the flag */
      atom_ptr->checked = FALSE;

      /* Blank out the pointer to the dictionary entry */
      atom_ptr->dic_data_ptr = NULL;
      
      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next_atom_ptr;
    }	       
}
/***********************************************************************

initialise_dic_data  -  Initialise all the dic_data flags for the given
                        dictionary entry

***********************************************************************/

void initialise_dic_data(struct dictionary *dictionary_ptr)
{
  int iatom, natoms;

  struct dic_data *dic_data_ptr;

  /* Get pointer to the start of linked list of dic_data coords */
  dic_data_ptr = dictionary_ptr->first_dic_data_ptr;
  natoms = dictionary_ptr->natoms;

  /* Loop through all the dic_datas to initialise the checked flag */
  for (iatom = 0; iatom < natoms && dic_data_ptr != NULL; iatom++)
    {
      /* Initialise the flag */
      dic_data_ptr->checked = FALSE;
      
      /* Blank out the pointer to the dictionary entry */
      dic_data_ptr->atom_ptr = NULL;

      /* Get pointer to next dic_data in linked-list */
      dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
    }	       
}
/***********************************************************************

initialise_nodes  -  Initialise all the node flags

***********************************************************************/

void initialise_nodes(struct node *first_node_ptr)
{
  struct node *node_ptr;

  /* Get pointer to the start of linked list of node coords */
  node_ptr = first_node_ptr;

  /* Loop through all the nodes to initialise the checked flag */
  while (node_ptr != NULL)
    {
      /* Initialise the flag */
      node_ptr->checked = FALSE;
      
      /* Get pointer to next node in linked-list */
      node_ptr = node_ptr->next_node_ptr;
    }	       
}
/***********************************************************************

initialise_node_links  -  Initialise all the node_link flags

***********************************************************************/

void initialise_node_links(struct node_link *first_node_link_ptr)
{
  struct node_link *node_link_ptr;

  /* Get pointer to the start of linked list of node_link coords */
  node_link_ptr = first_node_link_ptr;

  /* Loop through all the node_links to initialise the checked flag */
  while (node_link_ptr != NULL)
    {
      /* Initialise the flag */
      node_link_ptr->checked = FALSE;
      
      /* Get pointer to next node_link in linked-list */
      node_link_ptr = node_link_ptr->next_node_link_ptr;
    }	       
}
/***********************************************************************

count_atom_matches  -  Loop through all the nodes in the current clique
                       and count the number of times that an atom-atom
                       mapping, represented by a node, gives two atoms
		       with identical names

***********************************************************************/

void count_atom_matches(struct node *first_node_ptr,int *nmatch_atoms)
{
  struct node *node_ptr;

  /* Get pointer to the start of linked list of node coords */
  node_ptr = first_node_ptr;

  /* Loop through all the nodes to check the mappings */
  while (node_ptr != NULL)
    {
      /* Check only if node belongs to the current clique */
      if (node_ptr->checked == TRUE)
	{
	  /* If the atom-names are identical, increase the score */
	  if (!strncmp(node_ptr->atom_ptr->atom_name,
		       node_ptr->dic_data_ptr->atom_name,4))
	    (*nmatch_atoms)++;
	}

      /* Get pointer to next node in linked-list */
      node_ptr = node_ptr->next_node_ptr;
    }	       
}
/***********************************************************************

store_best_clique  -  Mark all the nodes according to whether they belong
	              to the current clique or not

***********************************************************************/

void store_best_clique(struct node *first_node_ptr)
{
  struct node *node_ptr;

  /* Get pointer to the start of linked list of nodes */
  node_ptr = first_node_ptr;

  /* Loop through all the nodes to initialise the checked flag */
  while (node_ptr != NULL)
    {
      /* Set the clique flag according to whether node belongs to the
         current clique or not */
      node_ptr->in_best_clique = node_ptr->checked;

      /* Get pointer to next node in linked-list */
      node_ptr = node_ptr->next_node_ptr;
    }	       
}
/***********************************************************************

clique_detection  -  Determine which mapping of residue atoms onto
                     dictionary atoms gives the best match

***********************************************************************/

void clique_detection(struct residue *residue_ptr,
		      struct node *first_node_ptr,
		      struct node_link *first_node_link_ptr)
{
  int add_to_stack, new_node;
  int nmatch_atoms, nmatch_links;
  int nstack;
  int best_score, link_score, match_score, score;

  struct atom *atom_ptr, *atom1_ptr, *atom2_ptr;
  struct dic_data *dic_data_ptr, *dic_data1_ptr, *dic_data2_ptr;
  struct dictionary *dictionary_ptr;
  struct node *node1_ptr, *node2_ptr;
  struct node *other_node1_ptr, *other_node2_ptr;
  struct node_link *trial_node_link_ptr, *node_link_ptr;
  struct node_link *node_link_stack_ptr[MAXSTACK], *other_node_link_ptr;

/* debug
  int inode;

  inode = 0;
debug */

  /* Initialise variables */
  best_score = 0;
  link_score = 2;
  match_score = 1;

  /* Get the corresponding dictionary entry */
  dictionary_ptr = residue_ptr->dictionary_ptr;

  /* Get pointer to first node link */
  trial_node_link_ptr = first_node_link_ptr;

  /* Loop through all the node links, using each one as the start-point */
  while (trial_node_link_ptr != NULL)
    {
      /* Initialise score for this matching */
      score = 0;
      nmatch_atoms = 0;
      nmatch_links = 0;
      nstack = 0;

      /* Initialise all nodes, links and atoms */
      initialise_nodes(first_node_ptr);
      initialise_node_links(first_node_link_ptr);
      initialise_atoms(residue_ptr);
      initialise_dic_data(dictionary_ptr);

      /* Get the two nodes joining this node, and their corresponding
         atoms */
      node1_ptr = trial_node_link_ptr->node1_ptr;
      node2_ptr = trial_node_link_ptr->node2_ptr;
      atom1_ptr = node1_ptr->atom_ptr;
      atom2_ptr = node2_ptr->atom_ptr;
      dic_data1_ptr = node1_ptr->dic_data_ptr;
      dic_data2_ptr = node2_ptr->dic_data_ptr;
/* debug
      inode++;
      printf("%3d. Node:  [%s] -> [%s] to [%s] -> [%s]\n",inode,
	     atom1_ptr->atom_name,dic_data1_ptr->atom_name,
	     atom2_ptr->atom_name,dic_data2_ptr->atom_name);
debug */

      /* Set current node-link and its constituent nodes and atoms as
	 checked */
      trial_node_link_ptr->checked = TRUE;
      node1_ptr->checked = TRUE;
      node2_ptr->checked = TRUE;
      atom1_ptr->checked = TRUE;
      atom2_ptr->checked = TRUE;
      dic_data1_ptr->checked = TRUE;
      dic_data2_ptr->checked = TRUE;

      /* Add node-link to the stack */
      node_link_stack_ptr[nstack] = trial_node_link_ptr;
      nstack++;

      /* Increment count of matched links */
      nmatch_links++;

      /* Loop until all stack entries have been exhausted */
      while (nstack > 0)
	{
	  /* Get next entry from the top of the stack */
	  nstack--;
	  node_link_ptr = node_link_stack_ptr[nstack];

	  /* Get the two nodes at either end of this node link */
	  node1_ptr = node_link_ptr->node1_ptr;
	  node2_ptr = node_link_ptr->node2_ptr;

	  /* Get first node-link again to loop through all links */
	  other_node_link_ptr = first_node_link_ptr;

	  /* Loop through all the node links adding any attached to
	     the current one onto the stack */
	  while (other_node_link_ptr != NULL)
	    {
	      /* Proceed only if not already marked as checked */
	      if (other_node_link_ptr->checked == FALSE)
		{
		  /* Get the two nodes joined by this node-link */
		  other_node1_ptr = other_node_link_ptr->node1_ptr;
		  other_node2_ptr = other_node_link_ptr->node2_ptr;

		  /* Initialise flag */
		  add_to_stack = FALSE;

		  /* Check whether either node corresponds to one of
		     the nodes of the link pulled off the stack */
		  if (other_node1_ptr == node1_ptr ||
		      other_node1_ptr == node2_ptr)
		    {
		      /* Node links connected via first node */
		      add_to_stack = TRUE;

		      /* Get the two atoms equivalenced by the other node */
		      atom_ptr = other_node2_ptr->atom_ptr;
		      dic_data_ptr = other_node2_ptr->dic_data_ptr;

		      /* Determine whether the other node is a new node */
		      new_node = !(other_node2_ptr->checked);
		    }

		  else if (other_node2_ptr == node1_ptr ||
			   other_node2_ptr == node2_ptr)
		    {
		      /* Node links connected via second node */
		      add_to_stack = TRUE;

		      /* Get the two atoms equivalenced by the other node */
		      atom_ptr = other_node1_ptr->atom_ptr;
		      dic_data_ptr = other_node1_ptr->dic_data_ptr;

		      /* Determine whether the other node is a new node */
		      new_node = !(other_node1_ptr->checked);
		    }

		  /* If both nodes have been checked already, then
		     don't need to add to stack, but do need to increment
		     score */
		  if (add_to_stack == TRUE && new_node == FALSE)
		    {
		      /* Increment the score */
		      nmatch_links++;
		      add_to_stack = FALSE;

		      /* Mark the link as processed */
		      other_node_link_ptr->checked = TRUE;
		    }

		  /* If this node-link a possibility for adding to
		     stack, perform final checks */
		  if (add_to_stack == TRUE)
		    {
		      /* Check that the atoms in the other node
			 have not already been included */
		      if (atom_ptr->checked == FALSE &&
			  dic_data_ptr->checked == FALSE)
			{
			  /* Increment the score */
			  nmatch_links++;

			  /* Mark the link as processed */
			  other_node_link_ptr->checked = TRUE;

			  /* Add this node-link to the stack */
			  node_link_stack_ptr[nstack] = other_node_link_ptr;
			  nstack++;

			  /* Mark its two nodes and its atoms as checked */
			  other_node1_ptr->checked = TRUE;
			  other_node2_ptr->checked = TRUE;
			  other_node1_ptr->atom_ptr->checked = TRUE;
			  other_node2_ptr->atom_ptr->checked = TRUE;
			  other_node1_ptr->dic_data_ptr->checked = TRUE;
			  other_node2_ptr->dic_data_ptr->checked = TRUE;
			}
		    }
		}
	      
	      /* Get pointer to the next node link */
	      other_node_link_ptr = other_node_link_ptr->next_node_link_ptr;
	    }
	}
      
      /* Count number of identical atom names in the current mapping */
      count_atom_matches(first_node_ptr,&nmatch_atoms);
      
      /* Calculate a score for the current clique from the numbers of
	 matched links and matched atom names */
      score = nmatch_links * link_score + nmatch_atoms * match_score;

      /* If this is the best score so far, store it as the best-matching
	 clique */
      if (score > best_score)
	{
	  /* Mark all the nodes according to whether they belong to the
	     clique or not */
	  store_best_clique(first_node_ptr);
	  best_score = score;
	}

      /* Get pointer to the next node link */
      trial_node_link_ptr = trial_node_link_ptr->next_node_link_ptr;
    }
}
/***********************************************************************

get_best_match  -  Use the best atom-atom match between the residue and
                   the corresponding Het Group Dictionary entry, linking
		   the dictionary atoms to the corresponding residue atoms

***********************************************************************/

void get_best_match(struct residue *residue_ptr,
		    struct node *first_node_ptr)
{
  int iatom, natoms;
  int nmapped, nmatch, nmissing;

  struct atom *atom_ptr;
  struct dic_data *dic_data_ptr;
  struct dictionary *dictionary_ptr;
  struct node *node_ptr;

  /* Get the dictionary pointer corresponding to this residue */
  dictionary_ptr = residue_ptr->dictionary_ptr;

  /* Print match details */
  printf("Residue %s %s %c  matched with Het Group Dictionary entry %s\n",
	 residue_ptr->res_name,
	 residue_ptr->res_num,
	 residue_ptr->chain,
	 dictionary_ptr->res_name);

  /* If numbers of atoms don't agree, then show the numbers */
  if (residue_ptr->non_hydrogen != dictionary_ptr->non_hydrogen)
    {
      printf("   Number of non-hydrogen atoms in residue = %d.",
	     residue_ptr->non_hydrogen);
      printf("  In Het Dic entry = %d\n",dictionary_ptr->non_hydrogen);
      printf("\n");
    }

  /* Initialise count of matching atoms */
  nmatch = 0;
  nmapped = 0;

  /* Initialise indicators of which atoms have been processed */
  initialise_atoms(residue_ptr);
  initialise_dic_data(dictionary_ptr);

  /* Get pointer to the start of linked list of nodes */
  node_ptr = first_node_ptr;

  /* Loop through all the nodes to pick out all the activated ones */
  while (node_ptr != NULL)
    {
      /* Check if node is an activated one */
      if (node_ptr->in_best_clique == TRUE)
	{
	  /* Increment count of matches */
	  nmatch++;

	  /* Get the mapping represented by this node */
	  atom_ptr = node_ptr->atom_ptr;
	  dic_data_ptr = node_ptr->dic_data_ptr;

	  /* Mark these as accounted for */
	  atom_ptr->checked = TRUE;
	  dic_data_ptr->checked = TRUE;

	  /* Establish pointers between the dictionary atom and the
	     residue atom */
	  dic_data_ptr->atom_ptr = atom_ptr;
	  atom_ptr->dic_data_ptr = dic_data_ptr;

	  /* If atom names not identical, then show mapping applied */
	  if (strncmp(atom_ptr->atom_name,dic_data_ptr->atom_name,4))
	    {
	      /* If this is the first mis-match, show the residue name */
	      if (nmapped == 0)
		printf("*** Residue:-  %s %s %c\n",
		       residue_ptr->res_name,
		       residue_ptr->res_num,
		       residue_ptr->chain);

	      /* Increment count of mismatched mappings */
	      nmapped++;
	      printf("*** Warning: Name mismatch.");
	      printf(" Residue atom [%s] maps to dictionary atom [%s]\n",
		     node_ptr->atom_ptr->atom_name,
		     node_ptr->dic_data_ptr->atom_name);
	    }

	  /* Copy atom-name from residue into the output name to be
	     written to the hbplus.rc file */
	  strncpy(dic_data_ptr->output_atom_name,atom_ptr->atom_name,4);
	  dic_data_ptr->output_atom_name[4] = '\0';
	}

      /* Get pointer to next node in linked-list */
      node_ptr = node_ptr->next_node_ptr;
    }   

  /* Show number of atoms mapped onto atoms with a different name */
  if (nmapped > 0)
    {
      printf("\n");
      printf("   Number of atoms mapped to atoms with a different ");
      printf("name:  %5d\n",nmapped);
      printf("   Number of identical matches:                     ");
      printf("       %5d\n",nmatch - nmapped);
      printf("\n");
    }

  /* If too few atoms matched, then take this to be an unsuccessful
     match, and abort */
  if ((residue_ptr->non_hydrogen < dictionary_ptr->non_hydrogen &&
       nmatch < residue_ptr->non_hydrogen - 1) ||
      (residue_ptr->non_hydrogen >= dictionary_ptr->non_hydrogen &&
       nmatch < dictionary_ptr->non_hydrogen - 1))
    {
      printf("*** Warning. Can't map residue to corresponding ");
      printf("dictionary entry\n");
      printf("***          Only mapped %d of %d atoms\n",
	     nmatch,residue_ptr->non_hydrogen);
      printf("\n");
      return;
    }

  /* If have atoms present in the dictionary but not in the residue,
     list them */
  if (residue_ptr->non_hydrogen < dictionary_ptr->non_hydrogen)
    {
      /* Initialise count of missing atoms */
      nmissing = 0;

      /* Get pointer to dictionary's first atom entry */
      dic_data_ptr = dictionary_ptr->first_dic_data_ptr;
      natoms = dictionary_ptr->natoms;

      /* Loop through all the dic_data atoms to pick up the missed
	 ones */
      for (iatom = 0; iatom < natoms && dic_data_ptr != NULL; iatom++)
	{
	  /* If entry has no pointer to a residue atom, then is one
	     that is missing from the residue */
	  if (dic_data_ptr->atom_ptr == NULL &&
	      dic_data_ptr->hydrogen == FALSE)
	    {
	      printf("*** Warning. Atom [%s] missing from residue %s\n",
		     dic_data_ptr->atom_name,dictionary_ptr->res_name);
	      nmissing++;
	    }

	  /* Get pointer to next dic_data in linked-list */
	  dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
	}

      /* Print count of missing atoms */
      printf("\n");
      printf("Number of missing atoms:  %3d\n",nmissing);
      printf("\n");
    }

  /* Otherwise, is have atoms present in the residue but not in the
     dictionary, list these */
  else if (residue_ptr->non_hydrogen > dictionary_ptr->non_hydrogen)
    {
      /* Initialise count of missing atoms */
      nmissing = 0;

      /* Get pointer to residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
      natoms = residue_ptr->natoms;

      /* Loop through all the atoms to pick up the missed ones */
      for (iatom = 0; iatom < natoms && atom_ptr != NULL; iatom++)
	{
	  /* If atom has no pointer to a dictionary entry, then is one
	     that is extra to the residue */
	  if (atom_ptr->dic_data_ptr == NULL &&
	      atom_ptr->hydrogen == FALSE)
	    {
	      printf("*** Warning. Extra atom, not present in");
	      printf(" dictionary [%s]\n",atom_ptr->atom_name);
	      nmissing++;
	    }
      
	  /* Get pointer to next atom in linked-list */
	  atom_ptr = atom_ptr->next_atom_ptr;
	}	       

      /* Print count of missing atoms */
      printf("\n");
      printf("Number of atoms in residue but not in dictionary:  %3d\n",
	     nmissing);
      printf("\n");
    }

  /* If satisfied that the residue matches the dictionary entry
     sufficiently well, then link the dictionary entry to this residue */
  dictionary_ptr->residue_ptr = residue_ptr;
}
/***********************************************************************

graph_match  -  Perform a simple graph-match between the atoms of the
                residue and those of the entry from the Het Group
		Dictionary to get the best atom-atom correspondences

***********************************************************************/

void graph_match(struct residue *residue_ptr,
		 struct connect *first_connect_ptr)
{
  int nnodes, nnode_links;

  struct node *first_node_ptr;
  struct node_link *first_node_link_ptr;

  /* Calculate all the covalent connectivities for the atoms of the
     residue */
  calc_covalent_connectivities(residue_ptr,first_connect_ptr);

  /* Create the nodes of the graph by pairing of like-minded atoms
     between the residue and the entry from the Het Group Dictionary */
  create_nodes(residue_ptr,&first_node_ptr,&nnodes);

  /* Calculate links between nodes */
  calc_node_links(first_node_ptr,&first_node_link_ptr,&nnode_links);

  /* Perform simple clique detection to find the best mapping between
     the atoms of the residue and the atom names in the entry from
     the Het Group Dictionary */
  clique_detection(residue_ptr,first_node_ptr,first_node_link_ptr);

  /* Create the links between the dictionary atoms and the residue
     atoms corresponding to the best match found in the clique
     detection routine */
  get_best_match(residue_ptr,first_node_ptr);
}
/***********************************************************************

validate_pdb_hets  -  Validate atom names and connectivities of the HET
                      groups in the PDB file according to the data in the
                      HET Group Dictionary, listing any groups not found
		      in the dictionary

***********************************************************************/

void validate_pdb_hets(struct residue *first_residue_ptr,
                       struct dictionary **first_dictionary_ptr,
		       struct dictionary **last_dictionary_ptr,
		       int *ndic_entries,
		       struct dic_data **last_dic_data_ptr,
		       struct connect *first_connect_ptr)
{
  int nmatched_bonds, printed_warning;

  struct residue *residue_ptr;
  struct dictionary *dictionary_ptr;

  /* Initialise variables */
  printed_warning = FALSE;
  nmatched_bonds = 0;

  /* Initialise residue-pointer to start of linked list */
  residue_ptr = first_residue_ptr;

  /* Loop through all the HET group residues */
  while (residue_ptr != NULL)
    {
      /* If have a matching entry in the PDB, then validate all the atom
	 names and connectivities */
      if (residue_ptr->dictionary_ptr != NULL)
	{
	  /* Get the dictionary entry corresponding to this residue */
	  dictionary_ptr = residue_ptr->dictionary_ptr;

	  /* If this dictionary entry has not yet been matched against
	     a residue, then perform the matching */
	  if (dictionary_ptr->residue_ptr == NULL)
	    {
	      /* Apply a simple graph-matching technique to match
		 atoms between the residue and the entry in the Het
		 Group Dictionary */
	      graph_match(residue_ptr,first_connect_ptr);
	    }

	}

      /* Get pointer to the next residue in the list */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Re-initialise residue-pointer to start of linked list */
  residue_ptr = first_residue_ptr;

  /* Loop through all the HET group residues to pick up any for which
     no dictionary entry has been written out */
  while (residue_ptr != NULL)
    {
      /* If this het group has no pointer to a dictionary entry, then
	 must be missing from the dictionary */
      if (residue_ptr->dictionary_ptr == NULL)
	{
	  /* Print warning message about this HET group */
	  print_missing_hets(residue_ptr,first_residue_ptr);
	  printed_warning = TRUE;

          /* Do the best that can be done with the information
	     available - ie create a dictionary entry from the
             atom-names and their connectivities, as in the PDB file */
	  generate_fake_entry(residue_ptr,first_residue_ptr,
			      &(*first_dictionary_ptr),
			      &(*last_dictionary_ptr),
			      &(*ndic_entries),
			      &(*last_dic_data_ptr),first_connect_ptr);
	}

      /* Get pointer to the next residue in the list */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* If any warning messages printed, then throw a blank line */
  if (printed_warning == TRUE)
    printf("\n");
}
/***********************************************************************

convert_quotes - Convert any double-quotes in the atom name to @'s

***********************************************************************/

void convert_quotes(char atom_name[5])
{
  int ipos;

  /* Loop through all the characters, making the conversion if necessary */
  for (ipos = 0; ipos < 4; ipos++)
    if (atom_name[ipos] == '"')
      atom_name[ipos] = '@';
}
/* <--v.3.2 */
/***********************************************************************

writename_fn - Function to write atom and residue names to hbplus.rc file

***********************************************************************/

void writename_fn(char res_name[4],int natoms,
		  struct dic_data *first_dic_data_ptr)
{
   int i, iatom;
   struct dic_data *dic_data_ptr;

   /* Initialize variables */
   i = iatom = 0;

   /* Write out the residue name */
/* v.3.1--> */
/*   fprintf(fil_out, "-u %s\n",res_name); */
   fprintf(fil_out, "-u \'%s\'\n",res_name);
/* <--v.3.1 */

   /* Initialise pointer to the first of the atom names */
   dic_data_ptr = first_dic_data_ptr;

   /* Loop through all the atoms to write out to file */
   while (iatom < natoms && dic_data_ptr != NULL)
     {
/* v.3.2--> */
       /* Convert any double-quotes in the atom name to @'s */
       convert_quotes(dic_data_ptr->output_atom_name);
/* <--v.3.2 */

       /* Write out only if a non-hydrogen atom */
       if (dic_data_ptr->atom_name[1] != 'H')
	 {
	   /* Start new line every 10 atoms */
	   if (i >= 10)
	     {
	       fprintf(fil_out,"\"\n");
	       i = 0;
	     }
	   if (i == 0)
/* v.3.1--> */
/*	     fprintf(fil_out,"-M %s \"",res_name); */
	     fprintf(fil_out,"-M \'%s\' \"",res_name);
/* <--v.3.1 */
	   
	   /* Write out the current atom name */
/* v.3.2--> */
/*	   fprintf(fil_out, "%s",dic_data_ptr->atom_name); */
	   fprintf(fil_out, "%s",dic_data_ptr->output_atom_name);
/* <--v.3.2 */
	   i++;
	 }

       /* Set pointer to the next atom and increment counter */
       dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
       iatom++;
     }
   /* Write out the end part of the current line */
   fprintf(fil_out,"\"\n");
}
/***********************************************************************

writecon_fn - Function to write the connectivities to the hbplus.rc file

***********************************************************************/

void writecon_fn(char res_name[4],int natoms,
		 struct dic_data *first_dic_data_ptr)
{
  int i, iatom, jatom;
  int flag;
  int flagb;
/* v.3.2--> */
/*  struct dic_data *dic_data_ptr; */
  struct dic_data *dic_data_ptr, *other_dic_data_ptr;
/* <--v.3.2 */

  /* Initialize variables */
  i = iatom = 0;
  flag  = 0;
  flagb = 0;

  /* Initialise pointer to the first of the atom names */
  dic_data_ptr = first_dic_data_ptr;
  
  /* Loop through all the atoms to write out to file */
  while (iatom < natoms && dic_data_ptr != NULL)
    {
      /* Loop through the connected atoms */
      for (jatom = 0; jatom < dic_data_ptr->nbound; jatom++)
	{
	  /* Write out only if both atoms are non-hydrogens */
	  if (dic_data_ptr->atom_name[1] != 'H' &&
	      dic_data_ptr->boundto[jatom][1] != 'H')
	    {
	      flagb = 1;
	      /* Start new line every 5 pairs */
	  
	      /* Flag is to prevent writing "-T" twice on same line */
	      if (i >= 5)
		{
		  fprintf(fil_out,"\"\n");
		  i    = 0;
		  flag = 0;
		}
	      if (i == 0 && flag == 0)
		{
/* v.3.1--> */
/*		  fprintf(fil_out,"-T %s \"",res_name); */
		  fprintf(fil_out,"-T \'%s\' \"",res_name);
/* <--v.3.1 */
		  flag = 1;
		}
	      
	      /* Write connectivity information */
	      if (dic_data_ptr->boundto_num[jatom] >= dic_data_ptr->number)
		{
/* v.3.2--> */
/*		  fprintf(fil_out,"%s%s:",dic_data_ptr->atom_name,
			  dic_data_ptr->boundto[jatom]); */
		  other_dic_data_ptr = dic_data_ptr->connect_atom_ptr[jatom];
		  fprintf(fil_out,"%s%s:",dic_data_ptr->output_atom_name,
			  other_dic_data_ptr->output_atom_name);
/* <--v.3.2 */
		  i++;
		}
	    }
	}

      /* Set pointer to next atom */
      dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
      iatom++;
    }
  /* Close line */
  if (flagb == 1)
    fprintf(fil_out,"\"\n");
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

writehbdon_fn - Function to write Hbond donor information to hbplus.rc 

***********************************************************************/

void writehbdon_fn(char res_name[4],int natoms,
/* v.3.2--> */
/*		   struct dic_data *first_dic_data_ptr) */
		   struct dic_data *first_dic_data_ptr,int fake_entry,
		   int have_hydrogens)
/* <--v.3.2 */
{
  int i, iatom, jatom;
  int nhyd;                              /* Count no. of H atoms bonded */   
/* v.3.1.2--> */
  int ncoval;
  int metal;
/* <--v.3.1.2 */

  struct dic_data *dic_data_ptr;

  /* Initialize variables */
  i = iatom = 0;
   
  /* Initialise pointer to the first of the atom names */
  dic_data_ptr = first_dic_data_ptr;

  /* Loop through all the atoms to write out to file */
  while (iatom < natoms && dic_data_ptr != NULL)
    {
/* v.3.1.2--> */
      /* If residue consists of a single atom, check whether it is a metal */
      metal = check_for_metal(dic_data_ptr->atom_name);
/* <--v.3.1.2 */

      /* Check if this atom is O, S or N - these potentially donate
         the same number of hydrogens as H atoms bonded to them */
      if (dic_data_ptr->atom_name[1] == 'O' ||
	  dic_data_ptr->atom_name[1] == 'S' ||
/* v.3.1.2--> */
/*	  dic_data_ptr->atom_name[1] == 'N') */
	  dic_data_ptr->atom_name[1] == 'N' || metal == TRUE)
/* <--v.3.1.2 */
	{
	  /* Initialise count of hydrogens attached to current atom */
	  nhyd = 0;
/* v.3.1.2--> */
	  ncoval = 0;
/* <--v.3.1.2 */

	  /* Count number of H atoms bonded to */
	  for (jatom = 0; jatom < dic_data_ptr->nbound; jatom++)
	    {
/* v.3.2--> */
/*	      if (dic_data_ptr->boundto[jatom][1] == 'H') */
	      if (dic_data_ptr->boundto[jatom][0] == 'H' ||
		  dic_data_ptr->boundto[jatom][1] == 'H')
/* <--v.3.2 */
		nhyd++;
/* v.3.1.2--> */
	      else
		ncoval++;
/* <--v.3.1.2 */
	    }

/* v.3.1.2--> */
/*	  if (nhyd > 0) */

	  /* If this is a metal, then write out up to 6 potential donors
	     so that HBPLUS picks up polar atoms coordinating this metal */
	  if (metal == TRUE)
	    {
	      nhyd = 6 - ncoval;
/* v.3.2--> */
/*	      fprintf(fil_out,"-E \'%s\' \"%s\" %d\n",
		      res_name,dic_data_ptr->atom_name,nhyd); */
	      fprintf(fil_out,"-E \'%s\' \"%s\" %d\n",
		      res_name,dic_data_ptr->output_atom_name,nhyd);
/* <--v.3.2 */
	    }

	  /* Write out number of Hbond donors = number of hydrogens found */
	  else if (nhyd > 0)
/* <--v.3.1.2 */
/* v.3.1--> */
/*	    fprintf(fil_out,"-E %s \"%s\" %d\n", */
	    fprintf(fil_out,"-E \'%s\' \"%s\" %d\n",
/* <--v.3.1 */
/* v.3.2--> */
/*		    res_name,dic_data_ptr->atom_name,nhyd); */
		    res_name,dic_data_ptr->output_atom_name,nhyd);

	  /* Otherwise, if no hydrogens found, and this is a fake entry,
	     then write out the maximum number of donors for this type
	     of atom */
	  else if (fake_entry == TRUE && have_hydrogens == FALSE)
	    {
	      /* Assign maximum number of hydrogens */
	      if (dic_data_ptr->atom_name[1] == 'O')
		nhyd = 1;
	      else if (dic_data_ptr->atom_name[1] == 'S')
		nhyd = 1;
	      else if (dic_data_ptr->atom_name[1] == 'N')
		nhyd = 2;

	      /* Write out the maximum number of hydrogen donors */
	      if (nhyd - ncoval > 0)
		fprintf(fil_out,"-E \'%s\' \"%s\" %d\n",
			res_name,dic_data_ptr->output_atom_name,
			nhyd - ncoval);
/* <--v.3.1.2 */
	    }
/* <--v.3.2 */
	}
      /* Set pointer to next atom */
      dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
      iatom++;
    }
}
/***********************************************************************

calc_angle  -  Calculate angles between the 3 sets of coordinates
               supplied

***********************************************************************/

float calc_angle(float coords1[3], float coords2[3], float coords3[3])
{
  int icoord;
  float angle, cosang, diff1, diff2, diff3, dist1, dist2, dist3;

  /* Initialise variables */
  angle = 0.0;
  dist1 = 0;
  dist2 = 0;
  dist3 = 0;

  /* Calculate distances between points 1->2, 2->3, 3->1 */
  for (icoord = 0; icoord < 3; icoord++)
    {
      diff1 = coords2[icoord] - coords1[icoord];
      diff2 = coords3[icoord] - coords2[icoord];
      diff3 = coords1[icoord] - coords3[icoord];
      dist1 = dist1 + diff1 * diff1;
      dist2 = dist2 + diff2 * diff2;
      dist3 = dist3 + diff3 * diff3;
    }

  /* Calculate angle */
  if (dist1 > 0.0 && dist2 > 0.0 && dist3 > 0.0)
    {
      dist1 = sqrt((double)dist1);
      dist2 = sqrt((double)dist2);
      dist3 = sqrt((double)dist3);
      cosang = (dist3 * dist3 - dist1 * dist1 - dist2 * dist2)
	/ (-2.0 * dist1 * dist2);
      if (abs(cosang) <= 1.0)
	angle = RADDEG * acos(cosang);
      else
	angle = 0.0;
    }

  /* Return the angle */
  return(angle);
}
/***********************************************************************

writehbacc_fn - Function to write Hbond acceptor information to hbplus.rc 

***********************************************************************/

void writehbacc_fn(char res_name[4],int natoms,
		   struct dictionary *dictionary_ptr,
		   struct dic_data *first_dic_data_ptr)
{
  int i, iatom, istore, jatom, jstore, naccept, nangles;
  int nplanar, nstore, nsp3, printed_warning;

  float angle, coords_N[3], coords_other[MAXCONNECT][3];

  struct atom *atom_ptr;
  struct dic_data *dic_data_ptr, *other_dic_data_ptr;
  struct residue *residue_ptr;

  /* Initialize variables */
  i = iatom = 0;
  printed_warning = FALSE;

  /* Get the HET group from the PDB file that matched all the
     atom names and connectivities */
  residue_ptr = dictionary_ptr->residue_ptr;

  /* If a valid pointer, then map all the atoms in the HET group to the
     corresponding coordinates in the PDB file */
/* v.3.2--> */
/*  if (residue_ptr != NULL)
    check_atoms_present(residue_ptr,&printed_warning,&nmatch); */
/* <--v.3.2 */
  
  /* Re-initialise pointer to the first of the atom names */
  dic_data_ptr = first_dic_data_ptr;
  iatom = 0;  

  /* Loop through all the atoms to write out to file */
  while (iatom < natoms && dic_data_ptr != NULL)
    {
      /* Check if this atom is O or S - these potentially accept
         2 hydrogens each */
      if (dic_data_ptr->atom_name[1] == 'O' ||
	  dic_data_ptr->atom_name[1] == 'S')
	{
	  /* Write out 2 Hbond acceptors */
/* v.3.1--> */
/*	  fprintf(fil_out,"-e %s \"%s\" 2\n", */
	  fprintf(fil_out,"-e \'%s\' \"%s\" 2\n",
/* <--v.3.1 */
/* v.3.2--> */
/*		  res_name,dic_data_ptr->atom_name); */
		  res_name,dic_data_ptr->output_atom_name);
/* <--v.3.2 */
	}

      /* If this atom is N, then more complicated rules apply according
         to connectivities */
      else if (dic_data_ptr->atom_name[1] == 'N')
	{
	  /* Case 1: sp2 or aromatic - always accepts 1 */
	  if (dic_data_ptr->nbound == 2)
/* v.3.1--> */
/*	    fprintf(fil_out, "-e %s \"%s\" 1\n", */
	    fprintf(fil_out, "-e \'%s\' \"%s\" 1\n",
/* <--v.3.1 */
/* v.3.2--> */
/*		    res_name,dic_data_ptr->atom_name); */
		    res_name,dic_data_ptr->output_atom_name);
/* <--v.3.2 */

	  /* Case 2: sp3 or amide - need to check bond angles */
	  else if (dic_data_ptr->nbound == 3)
	    {
	      /* Calculate all the bond-angles at this nitrogen position,
		 using the corresponding coordinates from the atoms in the
		 PDB file */

	      /* Get the coordinates of the N atom */
	      atom_ptr = dic_data_ptr->atom_ptr;

	      /* If have the N coordinates, then retrieve them */
	      if (atom_ptr != NULL)
		{
		  coords_N[0] = atom_ptr->x;
		  coords_N[1] = atom_ptr->y;
		  coords_N[2] = atom_ptr->z;

		  /* Initialise count of stored coordinates of atoms
		     bonded to this one */
		  nstore = 0;

		  /* Loop through the connected atoms and store all
		     their coordinates */
		  for (jatom = 0; jatom < dic_data_ptr->nbound; jatom++)
		    {
		      /* Get the pointer to the bonded atom */
		      other_dic_data_ptr
			= dic_data_ptr->connect_atom_ptr[jatom];
		      if (other_dic_data_ptr != NULL)
			atom_ptr = other_dic_data_ptr->atom_ptr;

		      /* If have a valid coordinate record, then store
			 the coordinates */
		      if (atom_ptr != NULL)
			{
			  coords_other[nstore][0] = atom_ptr->x;
			  coords_other[nstore][1] = atom_ptr->y;
			  coords_other[nstore][2] = atom_ptr->z;
			  nstore++;
			}
		    }

		  /* Calculate all possible angles and store */
		  nangles = 0;
		  nplanar = 0;
		  nsp3 = 0;

		  /* Loop over all pairs of stored coordinates */
		  for (istore = 0; istore < nstore - 1; istore++)
		    {
		      for (jstore = istore + 1; jstore < nstore; jstore++)
			{
			  /* Calculate this angle */
			  angle = calc_angle(coords_other[istore],coords_N,
					     coords_other[jstore]);
			  nangles++;

			  /* Angle close to 120 => planar, amide bond */
			  if (angle > (120 - ANGLE_ERROR))
			    nplanar++;
			  
			  /* Angle close to 110 => sp3 */
			  else if (angle < (110 + ANGLE_ERROR))
			    nsp3++;
			}
		    }

		  /* Now use the computed angles to determine the
		     hybridisation state of the nitrogen */

		  /* Determine number of H-bond acceptors according to
		     the angles just calculated */
		  naccept = 0;

		  /* No bond angles within limits */
		  if (nangles == 0)
		    {
		      printf("*** Warning. Cannot determine hybridization");
		      printf(" of %s atom in residue-type %s\n",
/* v.3.2--> */
/*			     dic_data_ptr->atom_name, */
			     dic_data_ptr->output_atom_name,
/* <--v.3.2 */
			     dictionary_ptr->res_name);
		      printf("***          H-bond acceptance set to 1\n");
		      naccept = 1;
		    }

		  /* More angles at 110 than at 120, so assume sp3 */
		  else if (nsp3 > nplanar)
		    naccept = 1;

		  /* Hybridization state ambiguous */
		  else if (nsp3 == nplanar)
		    {
		      printf("*** Warning. Hybridization ambiguous");
		      printf(" for %s atom in residue-type %s\n",
/* v.3.2--> */
/*			     dic_data_ptr->atom_name, */
			     dic_data_ptr->output_atom_name,
/* <--v.3.2 */
			     dictionary_ptr->res_name);
		      printf("***          H-bond acceptance set to 1\n");
		      naccept = 1;
		    }

		  /* Write out the number of H-bond acceptors */
		  if (naccept > 0)
/* v.3.1--> */
/*		    fprintf(fil_out,"-e %s \"%s\" %d\n", */
		    fprintf(fil_out,"-e \'%s\' \"%s\" %d\n",
/* <--v.3.1 */
/* v.3.2--> */
/*			    res_name,dic_data_ptr->atom_name,naccept); */
			    res_name,dic_data_ptr->output_atom_name,
			    naccept);
/* <--v.3.2 */

		}
	      /* Otherwise, print an error message */
	      else
		printf("*** Warning. Coords for N missing\n");
	    }
	}

      /* Set pointer to next atom */
      dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
      iatom++;
    }
}
/* v.3.2--> */
/***********************************************************************

open_hbplusrc  -  Open the hbplus.rc file

***********************************************************************/

void open_hbplusrc(char out_name[])
{
  /* Open the output hbplus.rc file */
  printf("  Opening HBPLUS input file, hbplus.rc ...\n");
  if ((fil_out = fopen(out_name,"w")) == NULL)
    {
      printf("\n*** Unable to open HBPLUS file [%s] for output\n",out_name);
      exit(1);
    }
}    
/* <--v.3.2 */
/***********************************************************************

generate_hbplusrc  -  Process each correctly-matched dictionary HET
                      group entry and write out a set of hbplus.rc
		      records for it

***********************************************************************/

void generate_hbplusrc(char out_name[],
		       struct dictionary *first_dictionary_ptr)
{
/* v.3.2--> */
  static int file_open = FALSE;
/* <--v.3.2 */
  struct dictionary *dictionary_ptr;

  /* Open the output hbplus.rc file */
/* v.3.2-->
  printf("  Opening HBPLUS input file, hbplus.rc ...\n");
  if ((fil_out = fopen(out_name,"w")) == NULL)
    {
      printf("\n*** Unable to open HBPLUS file [%s] for output\n",out_name);
      exit(1);
    }
<--v.3.2 */
    
  /* Initialise residue-pointer to start of linked list */
  dictionary_ptr = first_dictionary_ptr;

  /* Loop through all the HET group residues */
  while (dictionary_ptr != NULL)
    {
      /* If this het group has no pointer to a residue from the PDB, then
	 PDB entry must in some way differ (either in terms of mismatched
	 atom names or connectivities) */
      if (dictionary_ptr->residue_ptr == NULL)
	{
	  /* Print warning message about this HET group */
	  printf("\n");
	  printf("*** Warning. Could not match residue type [%s] ",
		 dictionary_ptr->res_name);
	  printf("to entry in dictionary.\n");
	  printf("***          Will not be written to hbplus.rc file\n");
	}

      /* If have a matching entry in the PDB, then process */
      else
	{
/* v.3.2--> */
	  /* If this is the first entry to be written out, open the
	     hbplus.rc file */
	  if (file_open == FALSE)
	    {
	      open_hbplusrc(out_name);
	      file_open = TRUE;
	    }
/* <--v.3.2 */

	  /* Write out the residue and atom names to the hbplus.rc file */
	  printf("    Writing H-bonding definitions for residue-type");
	  printf(" [%s] ...\n",dictionary_ptr->res_name);
	  writename_fn(dictionary_ptr->res_name,dictionary_ptr->natoms,
		       dictionary_ptr->first_dic_data_ptr);

	  /* Write out all the connectivities to the hbplus.rc file */
	  writecon_fn(dictionary_ptr->res_name,dictionary_ptr->natoms,
		      dictionary_ptr->first_dic_data_ptr);

	  /* Calculate numbers of hydrogen donors on each atom */
	  writehbdon_fn(dictionary_ptr->res_name,dictionary_ptr->natoms,
/* v.3.2--> */
/*			dictionary_ptr->first_dic_data_ptr); */
			dictionary_ptr->first_dic_data_ptr,
			dictionary_ptr->fake_entry,
			dictionary_ptr->have_hydrogens);
/* <--v.3.2 */

	  /* Calculate numbers of hydrogen acceptors on each atom */
	  writehbacc_fn(dictionary_ptr->res_name,dictionary_ptr->natoms,
 			dictionary_ptr,dictionary_ptr->first_dic_data_ptr);
	}

      /* Get pointer to the next residue in the list */
      dictionary_ptr = dictionary_ptr->next_dictionary_ptr;
    }

/* v.3.2--> */
  /* Close the hbplus.rc file */
  if (file_open == TRUE)
    fclose(fil_out);
/* <--v.3.2 */
}
/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  char dictionary_name[FILENAME_LEN];
  char out_name[] = "hbplus.rc";
  char pdb_name[FILENAME_LEN];

  int natoms, nhetgroups, ndic_entries;
/* v.3.1--> */
  int have_hydrogens;
/* <--v.3.1 */

  struct atom *first_atom_ptr;
  struct residue *first_residue_ptr;
  struct connect *first_connect_ptr;
  struct dictionary *first_dictionary_ptr;
  struct dic_data *first_dic_data_ptr;
/* v.3.1--> */
  struct dictionary *last_dictionary_ptr;
  struct dic_data *last_dic_data_ptr;
/* <--v.3.1 */

  /* Get the PDB file name from the command-line arguments*/
  get_command_arguments(argv,argc - 1,pdb_name,dictionary_name);
  
  /* Read in all the HET groups in the PDB file */
  read_pdbfile(pdb_name,&first_atom_ptr,&first_residue_ptr,
/* v.3.1--> */
/*	       &first_connect_ptr,&natoms,&nhetgroups); */
	       &first_connect_ptr,&natoms,&nhetgroups,&have_hydrogens);
/* <--v.3.1 */

  /* Read in any matching HET group definitions from the HET Group
     Dictionary */
/* v.3.1--> */
/*  read_dictionary(dictionary_name,first_residue_ptr,&first_dictionary_ptr,
                  &first_dic_data_ptr,&ndic_entries); */
  read_dictionary(dictionary_name,first_residue_ptr,&first_dictionary_ptr,
		  &last_dictionary_ptr,&first_dic_data_ptr,
		  &last_dic_data_ptr,&ndic_entries);
/* <--v.3.1 */

  /* Match up all the connected atoms from the data read in from the 
     HET group definitions */
  dic_connectivities(first_dictionary_ptr);

  /* Print any HET groups in the PDB file that weren't found in the 
     HET Group Dictionary, or don't quite match */
/* v.3.1--> */
/*  validate_pdb_hets(first_residue_ptr,first_connect_ptr);  */
  validate_pdb_hets(first_residue_ptr,&first_dictionary_ptr,
/* v.3.2--> */
/*		    &last_dictionary_ptr,&ndic_entries,&first_dic_data_ptr,
		    &last_dic_data_ptr,first_connect_ptr,have_hydrogens); */
		    &last_dictionary_ptr,&ndic_entries,
		    &last_dic_data_ptr,first_connect_ptr);
/* <--v.3.2 */
/* <--v.3.1 */

  /* Process each correctly-matched dictionary HET group entry and write
     out a set of hbplus.rc records for it */
  generate_hbplusrc(out_name,first_dictionary_ptr);

  /* Finish */
  return(0);
}
