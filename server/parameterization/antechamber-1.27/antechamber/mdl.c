/* MDL */
int rmdl(char *filename, int *atomnum, ATOM * atom, int *bondnum,
		 BOND * bond, CONTROLINFO cinfo, MOLINFO minfo)
{
	int i, j;
	int index = 0;
	int overflow_flag = 0;
	int tmpint1;
	int tmpint2;
	int tmpint3;
	int tmpint4;
	int itype;
	char type[20];
	char line[MAXCHAR];
	FILE *fpin;

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open file %s, exit\n", filename);
		exit(0);
	}
	initial(cinfo.maxatom, atom, minfo.resname);
	i = 0;
	j = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) {
/*       printf("\nFinished reading %s file.", cinfo.ifilename); */
			break;
		}
		index++;
		if (index == 4) {
			sscanf(line, "%d%d", &tmpint1, &tmpint2);
			continue;
		}
		if (index > 4 && index < tmpint1 + 5) {
			if (overflow_flag == 0)
				sscanf(line, "%lf%lf%lf%s", &atom[i].x, &atom[i].y,
					   &atom[i].z, atom[i].name);
			i++;
			if (i >= cinfo.maxatom && overflow_flag == 0) {
				printf
					("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
				overflow_flag = 1;
			}
		}
		if (index >= tmpint1 + 5 && index < tmpint1 + tmpint2 + 5) {
			sscanf(line, "%d%d%s", &tmpint3, &tmpint4, type);
			if (overflow_flag == 0) {
				atom[tmpint3 - 1].con[atom[tmpint3 - 1].connum++] =
					tmpint4 - 1;
				atom[tmpint4 - 1].con[atom[tmpint4 - 1].connum++] =
					tmpint3 - 1;
				itype = 0;
				if (strcmp(type, "1") == 0)
					itype = 1;
				if (strcmp(type, "2") == 0)
					itype = 2;
				if (strcmp(type, "3") == 0)
					itype = 3;
				if (strcmp(type, "am") == 0)
					itype = 1;
				if (strcmp(type, "ar") == 0)
					itype = 10;
				if (strcmp(type, "SINGLE") == 0)
					itype = 1;
				if (strcmp(type, "DOUBLE") == 0)
					itype = 2;
				if (strcmp(type, "TRIPLE") == 0)
					itype = 3;
				bond[j].bondi = tmpint3 - 1;
				bond[j].bondj = tmpint4 - 1;
				bond[j].type = itype;
			}
			j++;
			if (j >= cinfo.maxbond && overflow_flag == 0) {
				printf
					("\nInfo: the bond number exceeds the MAXBOND, reallocate memory automatically");
				overflow_flag = 1;
			}
		}
	}
	*atomnum = i;
	*bondnum = j;
/*  printf("\n Info: the atomic number is %5d%5d%5d", *atomnum,tmpint1, tmpint2); */
	fclose(fpin);
/*   for(i=0;i<*atomnum;i++)
   printf("\n%5d%5s%9.3lf%9.3lf%9.3lf%5d%5d%5d%5d%5d%5d", i+1, atom[i].name,
   atom[i].x, atom[i].y, atom[i].z,atom[i].con[0], atom[i].con[1], atom[i].con[2], atom[i].con[3], atom[i].con[4],atom[i].con[5]);
*/
	return overflow_flag;
}

void wmdl(char *filename, int atomnum, ATOM * atom, int bondnum,
		  BOND * bond, CONTROLINFO cinfo)
{
        typedef struct {
                char name[6];
        } ATOMNAME;
        ATOMNAME *name;
	
	int i, j, k;
	/* int breakindex; */
	char tmpchar[MAXCHAR];
	char type[5];
	FILE *fpout;
	if ((fpout = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Cannot open file %s, exit\n", filename);
		return;
	}
        name = (ATOMNAME *) malloc(sizeof(ATOMNAME) * (atomnum +10));
        if (name == NULL) {
                fprintf(stderr, "memory allocation error for *name in wmdl()\n");
                exit(0);
        }

	for (i = 0; i < strlen(filename); i++) {
		if (filename[i] == '.')
			break;
		tmpchar[i] = filename[i];
	}
	tmpchar[i] = '\0';
	element(atomnum, atom);
	for (i = 0; i < atomnum; i++)
		strcpy(name[i].name, atom[i].element);

	fprintf(fpout, "%s\n", tmpchar);
	fprintf(fpout, "  -ISIS-            3D\n\n");
	fprintf(fpout, "%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d\n", atomnum, bondnum,
			0, 0, 0, 0, 0, 0, 0, 0, 0);
	for (i = 0; i < atomnum; i++)
		fprintf(fpout, "%10.4lf%10.4lf%10.4lf %-2s%3d%3d%3d%3d%3d\n",
				atom[i].x, atom[i].y, atom[i].z, name[i].name, 0, 0, 0, 0, 0);

	for (i = 0; i < atomnum; i++)
		for (j = i + 1; j < atomnum; j++)
			for (k = 0; k < bondnum; k++)
				if ((bond[k].bondi == i && bond[k].bondj == j) ||
					(bond[k].bondi == j && bond[k].bondj == i)) {
					fprintf(fpout, "%3d%3d", bond[k].bondi + 1,
							bond[k].bondj + 1);
					strcpy(type, "0");
					if (bond[k].type == 1)
						strcpy(type, "1");
					if (bond[k].type == 2)
						strcpy(type, "2");
					if (bond[k].type == 3)
						strcpy(type, "3");
					if (bond[k].type == 7)
						strcpy(type, "1");
					if (bond[k].type == 8)
						strcpy(type, "2");
					if (bond[k].type == 9)
						strcpy(type, "2");
					if (bond[k].type == 10)
						strcpy(type, "1");
					fprintf(fpout, "%3s%3d%3d%3d\n", type, 0, 0, 0);
				}
	fprintf(fpout, "M  END\n");
	free(name);
	fclose(fpout);


}
