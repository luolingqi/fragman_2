# include <stdio.h>
# include <math.h>
# include <ctype.h>
# include <stdlib.h>
# include <string.h>
# define COLORTEXT "YES"
# define MAXCHAR 256
# define MAXPARM 500
# define MAXATOM 1024
# define MAXSET 10
# define debug 1

typedef struct {
	char name[20];
	double radius;
	double coefficient;
	double value;
	int corrid; 
} PARM;

typedef struct {
	char name[20];
	char type[20];
	double radius;
	double coefficient[MAXSET];
	double area;
	double x;
	double y;
	double z;
	int id;
} ATOM;

PARM parm[MAXPARM];
ATOM atom[MAXATOM];
int atomnum;
int parmnum;
int corrnum;
int setnum=1; 
double totalarea = 0;
double value = 0.0;

char line[MAXCHAR];
char ifilename[MAXCHAR];
char afilename[MAXCHAR];
char ofilename[MAXCHAR];
char pfilename[MAXCHAR];

char compname[MAXCHAR];

FILE *fpin;
FILE *fpa;
FILE *fpout;
FILE *fpparm;

int input_flag ;
int output_flag;
int parm_flag;
int iformat;
int command =1;
int colid = 1;
int colid_flag = 0;
char probe_radius[MAXCHAR] = "1.4";

void read_ac(char* filename) {
        int i;
        int num;
        int tmpint;
        double tmpfloat1, tmpfloat2;
        FILE *fpin;
        if ((fpin = fopen(filename, "r")) == NULL) {
                printf("\n Cannot open file %s, exit", filename);
                exit(0);
        }
        num = 0;
        for (;;) {
                if (fgets(line, MAXCHAR, fpin) == NULL) break;
                if (strncmp(line, "ATOM", 4) == 0) {
                        sscanf(&line[13], "%s",  atom[num].name);
                        sscanf(&line[26], "%lf%lf%lf",  &atom[num].x, &atom[num].y, &atom[num].z);
                        sscanf(&line[64], "%s",  atom[num].type); 
                        num++;
                        if(num > MAXATOM) {
                                printf("\nERROR: the atom number (%d) exceeds MAXATOM (%d), increase MAXATOM and recompile the program", num, MAXATOM);
                                exit(0);
                        }
                }
        }
        atomnum = num;
	if(debug == 1) 
		for(i=0;i<atomnum;i++)
			printf("\nATOM  %5d %5s %5s %9.3lf %9.3lf %9.3lf", i+1, atom[i].name, atom[i].type, atom[i].x, atom[i].y, atom[i].z);			
}
void read_mol2(char* filename) {
        int i;
        int num;
        int tmpint;
	int flag;
        double tmpfloat1, tmpfloat2;
	char tmpchar[MAXCHAR];
        FILE *fpin;
        if ((fpin = fopen(filename, "r")) == NULL) {
                printf("\n Cannot open file %s, exit", filename);
                exit(0);
        }
        num = 0;
	flag = 0;
        for (;;) {
                if (fgets(line, MAXCHAR, fpin) == NULL) break;
		sscanf(line, "%s", tmpchar);
		if(flag == 0 && strcmp(tmpchar, "@<TRIPOS>ATOM") == 0) {
			flag = 1;
			continue;
		}
		if(flag == 1 && strcmp(tmpchar, "@<TRIPOS>BOND") == 0) {
			flag = -1;
			break;
		}
		if(flag == 1) {
                        sscanf(line, "%d%s%lf%lf%lf%s", &tmpint, atom[num].name, &atom[num].x, &atom[num].y, &atom[num].z, atom[num].type);
                        num++;
                        if(num > MAXATOM) {
                                printf("\nERROR: the atom number (%d) exceeds MAXATOM (%d), increase MAXATOM and recompile the program", num, MAXATOM);
                                exit(0);
                        }
                }
        }
        atomnum = num;
	if(debug == 1) 
		for(i=0;i<atomnum;i++)
			printf("\nATOM  %5d %5s %5s %9.3lf %9.3lf %9.3lf", i+1, atom[i].name, atom[i].type, atom[i].x, atom[i].y, atom[i].z);			
}
void read_atc(char* filename) {
	int i,j;
	int num;
	int count;
	int colnum;
	int tmpint;
	double tmpfloat1, tmpfloat2;
	FILE *fpin;
	if ((fpin = fopen(filename, "r")) == NULL) {
		printf("\n Cannot open file %s, exit", filename);
		exit(0);
	}
	num = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) break;
		if (strncmp(line, "TOTAL_PARAMETER_SETS", 20) == 0) {
			sscanf(&line[20], "%d", &setnum);
			if(setnum >= MAXSET) {
				printf("\nThe number of coefficient parameter sets (%d) exceeds MAXSET (%d) , increase MAXSET and recompile", setnum, MAXSET); 
				exit(0);
			}
		}
		if (strncmp(line, "ATOM", 4) == 0) {
			sscanf(&line[4], "%d%s%s", &tmpint, atom[num].name, atom[num].type);
			count = 0;
			for(i = 0; i< strlen(line) -1; i++) 
				if(line[i]==' ' && line[i+1] != ' ') {
					count ++;
					if(count < 4) continue;
					if(count > 4 + setnum) break;
					sscanf(&line[i+1], "%lf", &atom[num].coefficient[count-4]);
				}
			num++;
			if(num > MAXATOM) {
				printf("\nERROR: the atom number (%d) exceeds MAXATOM (%d), increase MAXATOM and recompile the program", num, MAXATOM); 
				exit(0);
			}
		}
	}
	atomnum = num;
	if(debug == 1) {
		printf("\nThere are %d set(s) of coefficient parameters", setnum); 
		for(i=0;i<atomnum;i++) {
			printf("\nATOM  %5d %5s %5s", i+1, atom[i].name, atom[i].type);			
			for(j=0;j<setnum;j++)
				printf("%9.4lf", atom[i].coefficient[j]);			
		}
	}
}
void read_sac(char* filename){
	int i,j;
	int count;
	int colnum;
	int num;
	int tmpint;
	double tmpfloat1, tmpfloat2;
	FILE *fpin;
	if ((fpin = fopen(filename, "r")) == NULL) {
		printf("\n Cannot open file %s, exit", filename);
		exit(0);
	}
	num = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) break;
		if (strncmp(line, "TOTAL_PARAMETER_SETS", 20) == 0) {
			sscanf(&line[20], "%d", &setnum);
			if(setnum >= MAXSET) {
				printf("\nThe number of coefficient parameter sets (%d) exceeds MAXSET (%d) , increase MAXSET and recompile", setnum, MAXSET); 
				exit(0);
			}
		}
		if (strncmp(line, "ATOM", 4) == 0) {
			sscanf(&line[4], "%d%s%s%lf", &tmpint, atom[num].name, atom[num].type, &atom[num].area);
			count = 0;
			for(i = 0; i< strlen(line) -1; i++) 
				if(line[i]==' ' && line[i+1] != ' ') {
					count ++;
					if(count < 4) continue;
					if(count > 4 + setnum) break;
					sscanf(&line[i+1], "%lf", &atom[num].coefficient[count-4]);
				}
			num++;
			if(num > MAXATOM) {
				printf("\nERROR: the atom number (%d) exceeds MAXATOM (%d), increase MAXATOM and recompile the program", num, MAXATOM); 
				exit(0);
			}
		}
	}
	atomnum = num;
	if(debug == 1) {
		printf("\nThere are %d set(s) of coefficient parameters", setnum); 
		for(i=0;i<atomnum;i++) {
			printf("\nATOM  %5d %5s %5s %9.4lf", i+1, atom[i].name, atom[i].type, atom[i].area);			
			for(j=0;j<setnum;j++)
				printf("%9.2lf", atom[i].coefficient[j]);			
		}
	}
}
void read_parm(char* filename) {
	int i, j;
	int num1 = 0;
	int num2 = 0;
	int error_flag;
	char tmpchar1[10];
	char tmpchar2[10];
	double tmpfloat1, tmpfloat2;
	FILE *fpin;
	if ((fpin = fopen(filename, "r")) == NULL) {
		printf("\n Cannot open file %s, exit", filename);
		exit(0);
	}
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) break;
		if (strncmp(line, "RADIUS", 6) == 0) {
			parm[num1].coefficient = 0;
			sscanf(&line[7], "%s%lf%lf", tmpchar1, &tmpfloat1, &tmpfloat2);
			strcpy(parm[num1].name, tmpchar1);
			parm[num1].radius = tmpfloat1;
			parm[num1].coefficient = tmpfloat2;
			parm[num1].corrid = -1; 
			num1++;
			if(num1 > MAXPARM) {
				printf("\nERROR: the parameter number (%d) exceeds MAXPARM (%d), increase MAXPARM and recompile the program", num1, MAXPARM); 
				exit(0);
			}
		}
	}
	rewind(fpin);
	num2 = num1;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) break;
		if (strncmp(line, "CORR", 4) == 0) {
			parm[num2].coefficient = 0;
			sscanf(&line[4], "%s%s", tmpchar1, tmpchar2) ;
			strcpy(parm[num2].name, tmpchar1);
			error_flag = 1;		
			for(i = 0; i<num1; i++)  
				if(strcmp(parm[i].name, tmpchar2) == 0) {
					parm[num2].radius = parm[i].radius;
					parm[num2].coefficient = parm[i].coefficient;
					parm[num2].corrid = i;
					error_flag = 0;	
					num2++;
					if(num2 > MAXPARM) {
						printf("\nERROR: the parameter number (%d) exceeds MAXPARM (%d), increase MAXPARM and recompile the program", num2, MAXPARM); 
						exit(0);
					}
					break;
				}
			if(error_flag == 0) {
				fprintf(stderr, "\nERROR: no corresponding atom type parameter (%s) is avaliable for %s", tmpchar2, tmpchar1);
				exit(0);
			}	
		}
	}
	fclose(fpin);
	parmnum = num2;
	corrnum = num2 - num1;

	if(debug == 1) {
		printf("\nThe following atom type have their own parameters\n");
		for(i=0;i<num1;i++)
			printf("\nRADIUS\t%5s\t%9.4lf\t%9.4lf", parm[i].name, parm[i].radius, parm[i].coefficient);
		printf("\nThe following atom types borrow parameters from other atom types\n");
		for(i=num1;i<num2;i++)
			printf("\nRADIUS\t%5s\t%9.4lf\t%9.4lf\t%d(%s)", parm[i].name, parm[i].radius, 
			parm[i].coefficient, parm[i].corrid + 1, parm[parm[i].corrid].name);
	}
	for (i = 0; i < atomnum; i++) {
		error_flag = 1;
		for (j = 0; j < parmnum; j++)
			if (strcmp(parm[j].name, atom[i].type) == 0) {
				atom[i].radius = parm[j].radius;
				atom[i].coefficient[colid-1] = parm[j].coefficient;
				if(parm[j].corrid == -1)
					atom[i].id = j; 
				else
					atom[i].id = parm[j].corrid; 
				error_flag = 0;
				break;
			}
		if(error_flag == 1) {
			printf("\nERROR: no parameters for %s (%s)", atom[i].name, atom[i].type);
			exit(0);
		}
	} 
}
void write_atc(char* filename){
	FILE *fpout;
	int i,j;
	if ((fpout = fopen(ofilename, "w")) == NULL) {
        	printf("\n Cannot open file %s, exit", filename);
        	exit(0);
	}
	fprintf(fpout, "FILE TYPE: ATC (Atom Type Coefficient)\n");
	fprintf(fpout, "FORMAT: ATOM ID NAME TYPE COEFFICIENT)\n");
	fprintf(fpout, "PARAMETER: %s\n", pfilename);
	fprintf(fpout, "TOTAL_PARAMETER_SETS  %d\n", setnum);
	for(i=0;i<atomnum; i++) {
		fprintf(fpout, "ATOM %5d %8s %8s \n", i+1, atom[i].name, atom[i].type);
		for(j=0;j<setnum-1; j++) 
			fprintf(fpout, "%9.4lf ", atom[i].coefficient[j]);
		fprintf(fpout, "%9.4lf\n", atom[i].coefficient[setnum-1]);
	}
	fclose(fpout);
}
void write_sac(char* filename){
	FILE *fpout;
	int i,j;
	if ((fpout = fopen(ofilename, "w")) == NULL) {
        	printf("\n Cannot open file %s, exit", filename);
        	exit(0);
	}
	fprintf(fpout, "FILE TYPE: SAC (Surface Area Coefficient)\n");
	fprintf(fpout, "FORMAT: ATOM ID NAME TYPE AREA COEFFICIENT)\n");
	fprintf(fpout, "PARAMETER: %s\n", pfilename);
	fprintf(fpout, "PROBE RADIUS: %10s\n", probe_radius);
	fprintf(fpout, "TOTAL_PARAMETER_SETS  %d\n", setnum);
	for(i=0;i<atomnum; i++) {
		fprintf(fpout, "ATOM %5d %8s %8s %9.4lf\n", i+1, atom[i].name, atom[i].type, atom[i].area);
		for(j=0;j<setnum-1; j++) 
			fprintf(fpout, "%9.4lf ", atom[i].coefficient[j]);
		fprintf(fpout, "%9.4lf\n", atom[i].coefficient[setnum-1]);
	}
	fclose(fpout);
}
void write_des(char* filename){
	FILE *fpout;
	int i,j;
	if ((fpout = fopen(ofilename, "w")) == NULL) {
        	printf("\n Cannot open file %s, exit", filename);
        	exit(0);
	}
	for(i=0;i<parmnum-corrnum-1;i++) 
		parm[i].value = 0.0;
	for(j=0;j<atomnum;j++) {
		if(iformat == 2) parm[atom[j].id].value ++ ;
		if(iformat == 3) parm[atom[j].id].value += atom[j].area;
	}
	fprintf(fpout, "%s, ", compname);
	for(i=0;i<parmnum-corrnum-1;i++) {
		if(iformat == 2) fprintf(fpout, "%9.0lf, ", parm[i].value);
		if(iformat == 3) fprintf(fpout, "%9.4lf, ", parm[i].value);
	}
	if(iformat == 2) fprintf(fpout, "%9.0lf, ", parm[parmnum - corrnum -1].value);
	if(iformat == 3) fprintf(fpout, "%9.4lf, ", parm[parmnum - corrnum -1].value);
	fclose(fpout);
}
void write_sas(char *filename) {
	int i;
	int index;
	int id; 
	double radius = 0.0, coefficient = 0.0;
	char line[MAXCHAR];
	FILE *fpout;

	if ((fpout = fopen(filename, "w")) == NULL) {
		printf("\n Cannot open file %s, exit", filename);
		return;
	}
	for (i = 0; i < atomnum; i++) 
        	fprintf(fpout, "%8.3lf%8.3lf%8.3lf%8.4lf%8.4lf\n", 
		atom[i].x,atom[i].y,atom[i].z, atom[i].radius, atom[i].coefficient);
	fclose(fpout);
}
void read_area(char* filename) {
	FILE *fpsas; 	
	int num;
	int tmpint;
	double tmpf1, tmpf2;
	if ((fpsas = fopen(filename, "r")) == NULL) {
        	printf("\n Cannot open expt file to read%s, exit", filename);
        	exit(0);
	}
	if (fgets(line, MAXCHAR, fpsas) == NULL) return;
	num = 0;
	for(;;) {
        	if (fgets(line, MAXCHAR, fpsas) == NULL) break;
        	sscanf(line, "%ld%lf%lf", &tmpint, &tmpf1, &tmpf2);
		atom[num].area = tmpf2;
		num++;
	}	
	fclose(fpsas);
}

int main(int argc, char *argv[]) {
	int i;
	char tmpchar[MAXCHAR];
	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-H") == 0)) {
			printf("[31mUsage: molprop -i[0m  input file name (ac or mol2, sac, atc)\n"  
				   "[31m               -f[0m  input file type (ac or mol2, atc, sac)\n"  
				   "[31m               -o[0m  output file name (atc, sac, des)\n"  
				   "[31m               -p[0m  parameter file name (par)\n" 
				   "[31m               -r[0m  probe radius for SAS calculation for \"-c 3\"\n"  
				   "[31m               -ci[0m column id for command 1 or 5\n"  
				   "[31m               -co[0m command type\n"  
				   "                   1: calculate value(s) for one (specified by -ci) or all sets (the default) of coefficients (need a sac or atc file, and/or -ci):\n"  
				   "                   2: generate atc file (need ac or mol2, par) \n"  
				   "                   3: generate sac file (need ac or mol2, par, do SAS calculation)\n"  
				   "                   4: append coefficients to an atc or a sac file (need atc/sac and par)\n"  
				   "                   5: replace coefficients in an atc or a sac file (need atc/sac, par, and -ci)\n"  
				   "                   6: generate des file (need sac or atc and par\n") ;
			exit(0);
		}
		if (argc != 7 && argc != 9 && argc != 11 && argc != 13) {
			printf("[31mUsage: molprop -i[0m  input file name (ac or mol2, sac, atc)\n"  
				   "[31m               -f[0m  input file type (ac or mol2, atc, sac)\n"  
				   "[31m               -o[0m  output file name (atc, sac, des)\n"  
				   "[31m               -p[0m  parameter file name (par)\n" 
				   "[31m               -r[0m  probe radius for SAS calculation for \"-c 3\"\n"  
				   "[31m               -ci[0m column id for command 1 or 5\n"  
				   "[31m               -co[0m command type\n"  
				   "                   1: calculate value(s) for one (specified by -ci) or all sets (the default) of coefficients (need a sac or atc file, and/or -ci):\n"  
				   "                   2: generate atc file (need ac or mol2, par) \n"  
				   "                   3: generate sac file (need ac or mol2, par, do SAS calculation)\n"  
				   "                   4: append coefficients to an atc or a sac file (need atc/sac and par)\n"  
				   "                   5: replace coefficients in an atc or a sac file (need atc/sac, par, and -ci)\n"  
				   "                   6: generate des file (need sac or atc and par\n") ;
			exit(0);
		}
	} 
	else {
		if (argc == 2)
			if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-H") == 0) {
				printf("Usage: molprop -i  input file name (ac or mol2, sac, atc)\n"
				   "               -f  input file type (ac or mol2, atc, sac)\n"  
				   "               -o  output file name (atc, sac, des)\n"
				   "               -p  parameter file name (par)\n"
				   "               -r  probe radius for SAS calculation for \"-c 3\"\n"  
				   "               -ci column id for command 1 or 5\n"  
				   "               -co command type\n"  
				   "                   1: calculate value(s) for one (specified by -ci) or all sets (the default) of coefficients (need a sac or atc file, and/or -ci):\n"  
				   "                   2: generate atc file (need ac or mol2, par) \n"  
				   "                   3: generate sac file (need ac or mol2, par, do SAS calculation)\n"  
				   "                   4: append coefficients to an atc or a sac file (need atc/sac and par)\n"  
				   "                   5: replace coefficients in an atc or a sac file (need atc/sac, par, and -ci)\n"  
				   "                   6: generate des file (need sac or atc and par\n") ;
				exit(0);
			}
		if (argc != 7 && argc != 9 && argc != 11 && argc != 13) {
				printf("Usage: molprop -i  input file name (ac or mol2, sac, atc)\n"
				   "               -f  input file type (ac or mol2, atc, sac)\n"  
				   "               -o  output file name (atc, sac, des)\n"
				   "               -p  parameter file name (par)\n"
				   "               -r  probe radius for SAS calculation for \"-c 3\"\n"  
				   "               -ci column id for command 1 or 5\n"  
				   "               -co command type\n"  
				   "                   1: calculate value(s) for one (specified by -ci) or all sets (the default) of coefficients (need a sac or atc file, and/or -ci):\n"  
				   "                   2: generate atc file (need ac or mol2, par) \n"  
				   "                   3: generate sac file (need ac or mol2, par, do SAS calculation)\n"  
				   "                   4: append coefficients to an atc or a sac file (need atc/sac and par)\n"  
				   "                   5: replace coefficients in an atc or a sac file (need atc/sac, par, and -ci)\n"  
				   "                   6: generate des file (need sac or atc and par\n") ;
				exit(0);
			}
	}
	input_flag = 0;
	output_flag = 0;
	parm_flag = 0;
	strcpy(probe_radius, "1.4");
	command = 1;
	iformat = 0;
	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-i") == 0) {
			strcpy(ifilename, argv[i + 1]);
			input_flag = 1;
		}
		if (strcmp(argv[i], "-o") == 0) {
			strcpy(ofilename, argv[i + 1]);
			output_flag = 1;
		}
		if (strcmp(argv[i], "-p") == 0) {
			strcpy(pfilename, argv[i + 1]);
			parm_flag = 1;
		}
		if (strcmp(argv[i], "-r") == 0)
			strcpy(probe_radius, argv[i+1]);
		if (strcmp(argv[i], "-f") == 0)  
			if(strcmp(argv[i+1], "ac") == 0)
				iformat = 0;
			else if(strcmp(argv[i+1], "mol2") == 0)
				iformat = 1;
			else if(strcmp(argv[i+1], "atc") == 0)
				iformat = 2;
			else if(strcmp(argv[i+1], "sac") == 0)
				iformat = 3;
		if (strcmp(argv[i], "-ci") == 0)  {
			colid = atoi(argv[i+1]);
			colid_flag = 1;
			if(colid < 1) {
				printf("\nParameter Set ID sepcified by \"-ci\" cannot smaller than 1, exit");
				exit(0);
			}
		}
		if (strcmp(argv[i], "-co") == 0) 
			command = atoi(argv[i+1]);
	}
	if(command == 1) { 
		if(input_flag == 0 || (iformat != 2 && iformat != 3)) {
			printf("\nCommand = 1 requires a sac or atc file as the input");	
			exit(0);
		}
		parm_flag = 0;
		output_flag = 0;
	}
	if(command == 2)  {
		if(input_flag == 0 || (iformat != 0 && iformat != 1)) {
			printf("\nCommand = 2 requires an ac or mol2 file as the input");	
			exit(0);
		}
		if(parm_flag == 0) {
			printf("\nCommand = 2 requires to read in a parameter file");	
			exit(0);
		}
		if(output_flag == 0) {
			printf("\nCommand = 2 requires to specify an output file");	
			exit(0);
		}
		colid = 1;
	}
	if(command == 3)  {
		if(input_flag == 0 || (iformat != 0 && iformat != 1)) {
			printf("\nCommand = 3 requires an ac or mol2 file as the input");	
			exit(0);
		}
		if(parm_flag == 0) {
			printf("\nCommand = 3 requires to read in a parameter file");	
			exit(0);
		}
		if(output_flag == 0) {
			printf("\nCommand = 3 requires to specify an output file");	
			exit(0);
		}
		colid = 1;
	}
	if(command == 4)  {
		if(input_flag == 0 || (iformat != 2 && iformat != 3)) {
			printf("\nCommand = 4 requires an atc/sac file as the input");	
			exit(0);
		}
		if(parm_flag == 0) {
			printf("\nCommand = 4 requires to read in a parameter file");	
			exit(0);
		}
		if(output_flag == 0) {
			printf("\nCommand = 4 requires to specify an output file");	
			exit(0);
		}
	}
	if(command == 5)  {
		if(input_flag == 0 || (iformat != 2 && iformat != 3)) {
			printf("\nCommand = 5 requires an atc/sac file as the input");	
			exit(0);
		}
		if(colid_flag == 0) {
			printf("\nCommand = 5 requires to specify the # set of coefficients to work on with -ci");	
			exit(0);
		}
		if(parm_flag == 0) {
			printf("\nCommand = 5 requires to read in a parameter file");	
			exit(0);
		}
		if(output_flag == 0) {
			printf("\nCommand = 5 requires to specify an output file");	
			exit(0);
		}
	}
	if(command == 6)  {
		if(input_flag == 0 || (iformat != 2 && iformat != 3)) {
			printf("\nCommand = 5 requires a sac or atc file as the input");	
			exit(0);
		}
		if(parm_flag == 0) {
			printf("\nCommand = 5 requires to read in a parameter file");	
			exit(0);
		}
		if(output_flag == 0) {
			printf("\nCommand = 5 requires to specify an output file");	
			exit(0);
		}
	}

	if(input_flag == 1) {
		strcpy(compname, "");
		compname[0]='\0';
		for(i=0;i<strlen(ifilename);i++) {	
			if(ifilename[i]=='.') {
				compname[i]=='\0';
				break;
			}
			compname[i] = ifilename[i];
		}
	}
/* read in files */
	if(input_flag == 1 && iformat == 0) read_ac(ifilename);
	if(input_flag == 1 && iformat == 1) read_mol2(ifilename);
	if(input_flag == 1 && iformat == 2) read_atc(ifilename);
	if(input_flag == 1 && iformat == 3) read_sac(ifilename);
	if(parm_flag == 1) read_parm(pfilename);
	
/* run command*/
	if(command == 1) {
		value = 0;
		if(iformat == 2) read_atc(ifilename);
		if(iformat == 3) read_sac(ifilename);
		if(colid_flag == 1) {
			for(i=0;i<atomnum;i++) {
				if(iformat == 2) value+= atom[i].coefficient[colid-1]; 
				if(iformat == 3) value+= atom[i].coefficient[colid-1] * atom[i].area; 
			}
			printf("\nCalculated property is %9.4lf", value);
		}
		if(colid_flag == 0) {
			printf("\nCalculated property are");
			for(i=0;i<setnum;i++) {
				value = 0.0;
				for(i=0;i<atomnum;i++) {
					if(iformat == 2) value+= atom[i].coefficient[i]; 
					if(iformat == 3) value+= atom[i].coefficient[i] * atom[i].area; 
				}
			 	printf("%9.4lf", value);
			}
		}
		printf("\n");
		exit(0);
	}

	if(command == 2) {
		if(iformat == 0) read_ac(ifilename);
		if(iformat == 1) read_mol2(ifilename);
		read_parm(pfilename);
		write_atc(ofilename);
		exit(0);
	}
	if(command == 3) {
		if(iformat == 0) read_ac(ifilename);
		if(iformat == 1) read_mol2(ifilename);
		read_parm(pfilename);
		write_sas("ms.crd");
		strcpy(tmpchar, "ms -i ms.crd -o ms.out -r ");
		strcat(tmpchar, probe_radius);
		system(tmpchar);
		read_area("ms.out");
		write_sac(ofilename);
		if(debug == 1) {
			totalarea = 0.0;
			for(i=0;i<atomnum;i++) {
				printf("\nATOM  %5d %5s %5s %9.3lf %9.3lf %9.3lf %9.4lf %9.4lf", i+1, 
				atom[i].name, atom[i].type, atom[i].x, atom[i].y, atom[i].z, atom[i].area, atom[i].coefficient);			
				totalarea += atom[i].area ;
			}
			printf("\nTotal SAS is %9.4lf", totalarea);
		}
		exit(0);
	}
	if(command == 4) {
		if(iformat == 2) read_atc(ifilename);
		if(iformat == 3) read_sac(ifilename);
		setnum++;
                if(setnum >= MAXSET) {
                	printf("\nThe number of coefficient parameter sets including the appended one (%d) exceeds MAXSET (%d) , increase MAXSET and recompile", setnum, MAXSET);
                        exit(0);
                }
		colid = setnum ;
		read_parm(pfilename);
		if(iformat == 2) write_atc(ofilename);
		if(iformat == 3) write_sac(ofilename);
		exit(0);
	}
	if(command == 5) {
		if(iformat == 2) read_atc(ifilename);
		if(iformat == 3) read_sac(ifilename);
		if(colid < setnum) {
			printf("\nParameter Set ID sepcified by \"-ci\" (%d) is larger than the total set (%d) in the input file" , colid, setnum);
			exit(0);
		}
		read_parm(pfilename);
		if(iformat == 2) write_atc(ofilename);
		if(iformat == 3) write_sac(ofilename);
		exit(0);
	}
	if(command == 6) {
		if(iformat == 2) read_atc(ifilename);
		if(iformat == 3) read_sac(ifilename);
		read_parm(pfilename);
		write_des(ofilename);
		if(debug == 2) {
			totalarea = 0.0;
			for(i=0;i<parmnum - corrnum;i++) {
				printf("\nPARM %5d %5s  %9.4lf %9.4lf %9.4lf", i+1, parm[i].name, parm[i].radius,
				parm[i].coefficient, parm[i].value);
				if(iformat == 3) totalarea += parm[i].value ;
			}
			if(iformat == 3) printf("\nTotal SAS is %9.4lf", totalarea);
		}
		exit(0);
	}
	return (0);
}
