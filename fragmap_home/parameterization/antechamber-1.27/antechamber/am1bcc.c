/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    bcc                                                        *
*  Version: version 1.0                                                *
*  Author:  Junmei Wang                                                *
*                                                                      *
*  Department of Pharmaceutical Chemistry                              *
*  School of Pharmacy                                                  *
*  University of California                                            *
*  San Francisco   CA 94143                                            *
*  Octomber, 2001                                                      *
************************************************************************
*/

# include "common.h"
# include "define.h"
# include "atom.h"
# include "utility.c"
# include "common.c"
# include "rotate.c"
# include "ac.c"
# include "pdb.c"
# define MAX_BCCATOMTYPENUM 110
# define MAX_BCCBONDTYPENUM 15
# define debug 0

ATOM *atom;
BOND *bond;
int atomnum = 0;
int bondnum = 0;
int *bcctype;
int overflow_flag = 0;
MOLINFO minfo;
CONTROLINFO cinfo;

FILE *fpparm;
FILE *fpin;
FILE *fpout;

char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char pfilename[MAXCHAR] = "BCCPARM.DAT";
char *system_env;
char line[MAXCHAR];
int i, j, k, l;
int intstatus =1;
double bccparm[MAX_BCCATOMTYPENUM][MAX_BCCATOMTYPENUM][MAX_BCCBONDTYPENUM];
double totalcharge = 0.0;


void rparm(char *filename)
{
	int i, j, k;
	int tmpint1, tmpint2, tmpint3, tmpint4;
	double tmpfloat;
	char line[LINELEN_MAX];
	FILE *fpin;

	if ((fpin = fopen(filename, "r")) == NULL) {
		printf("\n Cannot open %s , exit", filename);
		exit(0);
	}

	for (i = 0; i < MAX_BCCATOMTYPENUM; i++)
		for (j = 0; j < MAX_BCCATOMTYPENUM; j++)
			for (k = 0; k < MAX_BCCBONDTYPENUM; k++)
				bccparm[i][j][k] = 0.0;

	for (;;) {
		if (fgets(line, LINELEN_MAX, fpin) == NULL)
			break;
		sscanf(line, "%d%d%d%d%lf", &tmpint1, &tmpint2, &tmpint3, &tmpint4,
			   &tmpfloat);
		bccparm[tmpint2][tmpint3][tmpint4] = tmpfloat;
	}
	if (debug == 1) {
		for (i = 0; i < MAX_BCCATOMTYPENUM; i++)
			for (j = i; j < MAX_BCCATOMTYPENUM; j++)
				for (k = 0; k < MAX_BCCBONDTYPENUM; k++)
					if (fabs(bccparm[i][j][k]) < 0.00000000001)
						continue;
					else
						printf("\n%5d%5d%5d%10.4lf", i, j, k,
							   bccparm[i][j][k]);
	}
}

void charge(void)
{
	int i, j, k, l, m, n, p, q, r, s;
	for (l = 0; l < bondnum; l++) {
		i = bond[l].bondi;
		j = bond[l].bondj;
		k = bond[l].type;
		m = bcctype[i];
		n = bcctype[j];
		if (m > n) {
			p = j;
			q = i;
			r = n;
			s = m;
		} else {
			p = i;
			q = j;
			r = m;
			s = n;
		}

		if (debug == 1) {
			printf("\n%5d %5d %5d %5d %5s %5s", l, i, j, k, atom[i].name,
				   atom[j].name);
			printf("\n%5d %5d %5d %9.4lf %9.4lf %9.4lf", l, p, q,
				   atom[p].charge, atom[q].charge, bccparm[r][s][k]);
		}
		atom[p].charge += bccparm[r][s][k];
		atom[q].charge -= bccparm[r][s][k];
	}

	if (debug == 1)
		for (i = 0; i < atomnum; i++)
			printf("\n%5d%9.4lf", i + 1, atom[i].charge);
}


int main(int argc, char *argv[])
{
	int i;
	int index;
	int format = 0;
	int judge_flag = 0;
	int status = 0;
	char command_at[2 * MAXCHAR];
	char command_bt[2 * MAXCHAR];

	default_minfo(&minfo);
	default_cinfo(&cinfo);

	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf
				("[31mUsage: am1bcc -i[0m input file name in ac format \n"
				 "[31m              -o[0m output file name \n"
				 "[31m              -f[0m output file format(pdb or ac, optional, default is ac)\n"
				 "[31m              -p[0m bcc parm file name (optional))\n"
				 "[31m              -s[0m status information, can be 0 (brief), 1 (the default) and 2 (verbose)\n"
				 "[31m              -j[0m atom and bond type judge option, default is 0)\n"
				 "[32m                 0[0m: No judgement\n"
				 "[32m                 1[0m: Atom type\n"
				 "[32m                 2[0m: Full bond type\n"
				 "[32m                 3[0m: Partial bond type\n"
				 "[32m                 4[0m: Atom and full bond type\n"
				 "[32m                 5[0m: Atom and partial bond type\n");
			exit(0);
		}
		if (argc != 7 && argc != 9 && argc != 11 &argc != 13) {
			printf
				("[31mUsage: am1bcc -i[0m input file name in ac format \n"
				 "[31m              -o[0m output file name \n"
				 "[31m              -f[0m output file format(pdb or ac, optional, default is ac)\n"
				 "[31m              -p[0m bcc parm file name (optional))\n"
				 "[31m              -s[0m status information, can be 0 (brief), 1 (the default) and 2 (verbose)\n"
				 "[31m              -j[0m atom and bond type judge option, default is 0)\n"
				 "[32m                 0[0m: No judgement\n"
				 "[32m                 1[0m: Atom type\n"
				 "[32m                 2[0m: Full bond type\n"
				 "[32m                 3[0m: Partial bond type\n"
				 "[32m                 4[0m: Atom and full bond type\n"
				 "[32m                 5[0m: Atom and partial bond type\n");
			exit(0);
		}
	} else {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("Usage: am1bcc -i input file name in ac format \n"
				   "              -o output file name \n"
				   "              -f output file format(pdb or ac, optional, default is ac)\n"
				   "              -p bcc parm file name (optional))\n"
				   "              -s status information, can be 0 (brief), 1 (the default) and 2 (verbose)\n"
				   "              -j atom and bond type judge option, default is 0)\n"
				   "                 0: No judgement\n"
				   "                 1: Atom type\n"
				   "                 2: Full bond type\n"
				   "                 3: Partial bond type\n"
				   "                 4: Atom and full bond type\n"
				   "                 5: Atom and partial bond type\n");
			exit(0);
		}
		if (argc != 7 && argc != 9 && argc != 11 && argc != 13) {
			printf("Usage: am1bcc -i input file name in ac format \n"
				   "              -o output file name \n"
				   "              -f output file format(pdb or ac, optional, default is ac)\n"
				   "              -p bcc parm file name (optional))\n"
				   "              -s status information, can be 0 (brief), 1 (the default) and 2 (verbose)\n"
				   "              -j atom and bond type judge option, default is 0)\n"
				   "                 0: No judgement\n"
				   "                 1: Atom type\n"
				   "                 2: Full bond type\n"
				   "                 3: Partial bond type\n"
				   "                 4: Atom and full bond type\n"
				   "                 5: Atom and partial bond type\n");
			exit(0);
		}

	}

	index = 0;
	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-i") == 0)
			strcpy(ifilename, argv[i + 1]);
		if (strcmp(argv[i], "-o") == 0)
			strcpy(ofilename, argv[i + 1]);
		if (strcmp(argv[i], "-j") == 0)
			judge_flag = atoi(argv[i + 1]);
		if (strcmp(argv[i], "-s") == 0)
			intstatus = atoi(argv[i + 1]);
		if (strcmp(argv[i], "-f") == 0)
			if (strcmp(argv[i + 1], "PDB") == 0
				|| strcmp(argv[i + 1], "pdb") == 0)
				format = 1;
			else
				format = 0;
		if (strcmp(argv[i], "-p") == 0) {
			strcpy(pfilename, argv[i + 1]);
			index = 1;
		}
	}
	system_env = (char *) getenv("ACHOME");
	if (system_env == NULL) 
		system_env = (char *) getenv("AMBERHOME");
	if (system_env != NULL) {
		if (index == 0) {
			pfilename[0] = '\0';
			strcpy(pfilename, system_env);
			strcat(pfilename, "/dat/antechamber/BCCPARM.DAT");
		}
		command_at[0] = '\0';
		command_bt[0] = '\0';
		strcpy(command_at, system_env);
		strcpy(command_bt, system_env);
		strcat(command_at,
			   "/exe/atomtype -f ac -p bcc -o ANTECHAMBER_AM1BCC.AC -i ");
		strcat(command_bt,
			   "/exe/bondtype -f ac -o ANTECHAMBER_AM1BCC.AC -i ");
	} else {
		strcpy(pfilename, "BCCPARM.DAT");
		strcpy(command_at,
			   "atomtype -f ac -p bcc -o ANTECHAMBER_AM1BCC.AC -i ");
		strcpy(command_bt, "bondtype -f ac -o ANTECHAMBER_AM1BCC.AC -i ");
	}

	if (judge_flag == 1) {
		strcat(command_at, ifilename);
		if (intstatus == 2)
			fprintf(stderr, "\nRunning: %s\n", command_at);
		status = system(command_at);
		if(status != 0) {
                	fprintf(stderr, "Error: cannot run \"%s\" in main() of am1bcc.c properly, exit\n", command_at);
                	exit(0);
		}
	}
	if (judge_flag == 2) {
		strcat(command_bt, ifilename);
		strcat(command_bt, " -j full");
		if (intstatus == 2)
			fprintf(stderr, "\nRunning: %s\n", command_bt);
		status = system(command_bt);
	        if(status != 0) {
                        fprintf(stderr, "Error: cannot run \"%s\" in main() of am1bcc.c properly, exit\n", command_bt);
                        exit(0);
                }
	}
	if (judge_flag == 3) {
		strcat(command_bt, ifilename);
		if (intstatus == 2)
			fprintf(stderr, "\nRunning: %s\n", command_bt);
		status = system(command_bt);
                if(status != 0) {
                        fprintf(stderr, "Error: cannot run \"%s\" in main() of am1bcc.c properly, exit\n", command_bt);
                        exit(0);
                }
	}
	if (judge_flag == 4) {
		strcat(command_bt, ifilename);
		strcat(command_bt, " -j full");
		if (intstatus == 2)
			fprintf(stderr, "\nRunning: %s\n", command_bt);
		status = system(command_bt);
                if(status != 0) {
                        fprintf(stderr, "Error: cannot run \"%s\" in main() of am1bcc.c properly, exit\n", command_bt);
                        exit(0);
                }
		strcat(command_at, "ANTECHAMBER_AM1BCC.AC");
		if (intstatus == 2)
			fprintf(stderr, "\nRunning: %s\n", command_at);
		status = system(command_at);
                if(status != 0) {
                        fprintf(stderr, "Error: cannot run \"%s\" in main() of am1bcc.c properly, exit\n", command_at);
                        exit(0);
                }
	}
	if (judge_flag == 5) {
		strcat(command_bt, ifilename);
		if (intstatus == 2)
			fprintf(stderr, "\nRunning: %s\n", command_bt);
		status = system(command_bt);
                if(status != 0) {
                        fprintf(stderr, "Error: cannot run \"%s\" in main() of am1bcc.c properly, exit\n", command_bt);
                        exit(0);
                }
		strcat(command_at, "ANTECHAMBER_AM1BCC.AC");
		if (intstatus == 2)
			fprintf(stderr, "\nRunning: %s\n", command_at);
		status = system(command_at);
                if(status != 0) {
                        fprintf(stderr, "Error: cannot run \"%s\" in main() of am1bcc.c properly, exit\n", command_at);
                        exit(0);
                }
	}

	if (judge_flag != 0)
		strcpy(ifilename, "ANTECHAMBER_AM1BCC.AC");
/*allocate memory*/
	atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
	if (atom == NULL) {
		fprintf(stderr, "memory allocation error for *atom\n");
		exit(0);
	}
	bond = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
	if (bond == NULL) {
		fprintf(stderr, "memory allocation error for *bond\n");
		exit(0);
	}
	bcctype = (int *) malloc(sizeof(int) * cinfo.maxatom);
	if (bcctype == NULL) {
		fprintf(stderr, "memory allocation error for *bcctype\n");
		exit(0);
	}
/*read ac file*/
	overflow_flag =
		rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
	if (overflow_flag) {
		cinfo.maxatom = atomnum + 10;
		cinfo.maxbond = bondnum + 10;
		free(atom);
		free(bond);
		free(bcctype);
		atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
		if (atom == NULL) {
			fprintf(stderr, "memory allocation error for *atom\n");
			exit(0);
		}
		bond = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
		if (bond == NULL) {
			fprintf(stderr, "memory allocation error for *bond\n");
			exit(0);
		}
		bcctype = (int *) malloc(sizeof(int) * cinfo.maxatom);
		if (bcctype == NULL) {
			fprintf(stderr, "memory allocation error for *bcctype\n");
			exit(0);
		}
		overflow_flag =
			rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
	}
	atomicnum(atomnum, atom);
	adjustatomname(atomnum, atom, 1);

	for (i = 0; i < atomnum; i++)
		bcctype[i] = atoi(atom[i].ambername);
	rparm(pfilename);
	charge();
	if (format == 0)
		wac(ofilename, atomnum, atom, bondnum, bond, cinfo, minfo);
	if (format == 1)
		wpdb(ofilename, atomnum, atom);
/*
	 free(atom);
	 free(bond);
	 free(bcctype);
*/
	return (0);
}
