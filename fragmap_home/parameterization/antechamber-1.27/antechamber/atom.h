# define LINELEN_MAX 160
typedef struct {
	char name[10];
	char aa[20];
	char ambername[10];
	char bcctype[10];
	char element[10];
	char resid[10];
        char chain[5];
	double x;
	double y;
	double z;
	double radius;
	double charge;
	double eps;
	double lja;
	double ljb;
	double bond;
	double angle;
	double twist;
        int ter;
	int resno;
	int hdonor;
	int haccept;
	int atomtype;
	int con[6];
	int bondatom;
	int angleatom;
	int twistatom;
	int connum;
	int mainatom;
	int atomicnum;
	int select;
	int ewd;		/* electrwithdraw or not: 1 ewd; -1 edonor; 0 others */
	int arom;		/* 5: 5 mem.ring aromatic atom; 6: 6 mem.ring aromatic atom;
       
								   -1: aliphatic; 0: others */
	int aliph;
	int saturate;				/* 1: saturated; -1 unsaturated; */
	int improper;
} ATOM;
typedef struct {
	char name[10];
	double x;
	double y;
	double z;
} ATOMSML;
typedef struct {
	char name[10];
	int atomtype;
	int con[6];
	double x;
	double y;
	double z;
	double charge;
} ATOMDBS;
typedef struct {
	char name[10];
	double x;
	double y;
	double z;
	double radius;
	double charge;
	double eps;
	double lja;
	double ljb;
	int hdonor;
	int haccept;
	int atomtype;
	int con[6];
} LIGATOM;
typedef struct {
	int atm1;
	int atm2;
	int atm3;
	int atm4;
	char type;
	double k;
	double kexi;
	int multi;
} TWIST;
typedef struct {
	double x;
	double y;
	double z;
} POINT;
typedef struct {
	int atomno;
	double x;
	double y;
	double z;
	double esp;
	double dist;
} ESP;
typedef struct {
	double x;
	double y;
	double z;
	double total;
} DM;
typedef struct {
	double xx;
	double yy;
	double zz;
	double xy;
	double xz;
	double yz;
} QM;
typedef struct {
	int atconnum;
        int apindex;
        char ap[512];
	char atname[10];
	char cesname[10];
} CHEMENV;

typedef struct {
        int bondi;
        int bondj;
        int type;
	double bcc;
} BOND;

typedef struct {
	char name1[10];
	char name2[10];
	double length;
	double force;
} BOND_FF;

typedef struct {
        int num;
        int atomno[12];
} RING;

typedef struct {
        int nr;
        int ar1;
        int ar2;
        int ar3;
        int ar4;
        int ar5;
        int rg[12];
} AROM;

typedef struct {
        char gkeyword[MAXCHAR] ;
        char mkeyword[MAXCHAR];
        char resname[20];
	char chkfile[20];
        char atom_type_def[10];
        char resfilename[MAXCHAR];
        char gfilename[MAXCHAR];
        char connect_file[MAXCHAR];
        char radius_file[MAXCHAR];
        char longresname[MAXCHAR];
        char inf_filename[MAXCHAR];
        int multiplicity; /*molecular multiplicity*/
        int icharge ; /*integer charge */
        int usercharge; /*user read in charge*/
	int envtype; /*environment type: 1-$ACHOME; 2-$AMBERHOME; 3-both set; 0-not set*/
        double dcharge; /*float charge*/
} MOLINFO;

typedef struct {
        char intype[10];       
        char outtype[10];       
        char atype[10];
        char chargetype[10];
        int rnindex; /*read in residue index ?*/ 
        int intstatus; /*information statius*/  
        int pfindex; /*purify intermediate files*/
        int prediction_index; /*the index of performing atomtype */
        int bpindex; /*bond type prediciton index*/
        int  maxatom; 
        int  maxbond;
        int  maxring;
} CONTROLINFO;

typedef struct {
        char name[10];
        double a;
        double b;
        double c;
	double d;
	double charge;
} GASTEIGER;

typedef struct {
        char name[10];
} NAME;

typedef struct {
        int aps[10];
} AV;
