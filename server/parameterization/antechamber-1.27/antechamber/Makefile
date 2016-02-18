# makefile for antechamber

include ../config.h

PROGS = \
	am1bcc$(SFX) \
	antechamber$(SFX) \
	atomtype$(SFX) \
	bondtype$(SFX) \
	charmmgen$(SFX) \
	database$(SFX) \
	espgen$(SFX) \
	parmcal$(SFX) \
	parmchk$(SFX) \
	prepgen$(SFX) \
	respgen$(SFX) \
	translate$(SFX) \
	top2mol2$(SFX)

#############################################################################

all: $(PROGS)

install: all
	cp -f mopac.sh ../exe
	cp -f dacdif ../exe
	mv -f $(PROGS) ../exe
        
clean:
	-rm -f *.o $(PROGS)

#############################################################################

#_____________________________________________________________________________
# Dependencies, for ./makedepend:

am1bcc.o: common.h
am1bcc.o: define.h
am1bcc.o: atom.h
am1bcc.0: utility.c
am1bcc.o: common.c
am1bcc.o: rotate.c
am1bcc.o: ac.c
am1bcc.o: pdb.c
antechamber.o: common.h
antechamber.o: define.h
antechamber.o: atom.h
antechamber.o: utility.c
antechamber.o: common.c
antechamber.o: ring.c
antechamber.o: rotate.c
antechamber.o: ac.c
antechamber.o: mol2.c
antechamber.o: mopcrt.c
antechamber.o: mopint.c
antechamber.o: mopout.c
antechamber.o: gcrt.c
antechamber.o: gzmat.c
antechamber.o: gout.c
antechamber.o: pdb.c
antechamber.o: csd.c
antechamber.o: mdl.c
antechamber.o: alc.c
antechamber.o: hin.c
antechamber.o: prep.c
antechamber.o: rst.c
antechamber.o: charge.c
antechamber.o: jcrt.c
antechamber.o: jzmat.c
antechamber.o: jout.c
atomtype.o: common.h
atomtype.o: define.h
atomtype.o: atom.h
atomtype.o: utility.c
atomtype.o: common.c
atomtype.o: ring.c
atomtype.o: rotate.c
atomtype.o: ac.c
atomtype.o: mol2.c
charmmgen.o: common.h
charmmgen.o: define.h
charmmgen.o: atom.h
charmmgen.o: utility.c
charmmgen.o: common.c
charmmgen.o: rotate.c
charmmgen.o: ac.c
charmmgen.o: mol2.c
charmmgen.o: prep.c
crdgrow.o: common.h
crdgrow.o: define.h
crdgrow.o: atom.h
crdgrow.o: common.c
crdgrow.o: rotate.c
crdgrow.o: pdb.c
crdgrow.o: ac.c
crdgrow.o: prep.c
database.o: common.h
delphigen.o: common.h
delphigen.o: define.h
delphigen.o: atom.h
delphigen.o: common.c
delphigen.o: ac.c
delphigen.o: pdb.c
espgen.o: common.h
espgen.o: define.h
espgen.o: atom.h
parmcal.o: common.h
parmcal.o: define.h
parmcal.o: utility.c
parmchk.o: common.h
parmchk.o: define.h
parmchk.o: atom.h
parmchk.o: utility.c
parmchk.o: common.c
parmchk.o: rotate.c
parmchk.o: ac.c
parmchk.o: mol2.c
parmchk.o: prep.c
parmjoin.o: common.h
parmjoin.o: define.h
parmjoin.o: utility.c
prepgen.o: common.h
prepgen.o: define.h
prepgen.o: atom.h
prepgen.o: utility.c
prepgen.o: common.c
prepgen.o: rotate.c
prepgen.o: ring.c
prepgen.o: ac.c
respgen.o: common.h
respgen.o: define.h
respgen.o: atom.h
respgen.o: utility.c
respgen.o: common.c
respgen.o: ac.c
translate.o: common.h 
translate.o: define.h 
translate.o: atom.h 
translate.o: utility.c 
translate.o: common.c 
translate.o: rotate.c 
translate.o: ac.c 
translate.o: mol2.c
translate.o: prep.c
translate.o: lsfit.c
rotate.o: common.h
smiles.o: common.h
smiles.o: define.h
smiles.o: atom.h
smiles.o: globals.h

#_____________________________________________________________________________
# Programs

am1bcc$(SFX): am1bcc.o
	$(CC) $(CFLAGS) am1bcc.o $(LM) -o am1bcc$(SFX)  

antechamber$(SFX): antechamber.o
	$(CC) $(CFLAGS) antechamber.o $(LM) -o antechamber$(SFX)  

atomtype$(SFX): atomtype.o
	$(CC) $(CFLAGS) atomtype.o $(LM) -o atomtype$(SFX)  

bondtype.o: bondtype.C
	$(CXX) -c $(CPPFLAGS) -o bondtype.o bondtype.C

bondtype$(SFX): bondtype.o
	$(CXX) bondtype.o $(LM) -o bondtype$(SFX)  

charmmgen$(SFX): charmmgen.o
	$(CC) $(CFLAGS) charmmgen.o $(LM) -o charmmgen$(SFX)

crdgrow$(SFX): crdgrow.o
	$(CC) $(CFLAGS) crdgrow.o $(LM) -o crdgrow$(SFX)

database$(SFX): database.o
	$(CC) $(CFLAGS) database.o $(LM) -o database$(SFX)

delphigen$(SFX): delphigen.o
	$(CC) $(CFLAGS) delphigen.o -o delphigen$(SFX)

espgen$(SFX): espgen.o
	$(CC) $(CFLAGS) espgen.o -o espgen$(SFX)

parmcal$(SFX): parmcal.o 
	$(CC) $(CFLAGS) parmcal.o $(LM) -o parmcal$(SFX)

parmchk$(SFX): parmchk.o 
	$(CC) $(CFLAGS) parmchk.o $(LM) -o parmchk$(SFX)

parmjoin$(SFX): parmjoin.o 
	$(CC) $(CFLAGS) parmjoin.o -o parmjoin$(SFX)

prepgen$(SFX): prepgen.o 
	$(CC) $(CFLAGS) prepgen.o $(LM) -o prepgen$(SFX)

respgen$(SFX): respgen.o 
	$(CC) $(CFLAGS) respgen.o $(LM) -o respgen$(SFX)

translate$(SFX): translate.o 
	$(CC) $(CFLAGS) translate.o $(LM) -o translate$(SFX)

top2mol2$(SFX): top2mol2.o
	$(CC) $(CFLAGS) top2mol2.o $(LM) -o top2mol2$(SFX)

