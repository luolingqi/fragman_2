install:
	cd tleap; make install
	cd antechamber; make install
	cd mopac6/src; make install
	cd resp; make install

clean:
	cd tleap; make clean
	cd antechamber; make clean
	cd mopac6/src; make clean
	cd resp; make clean

test::
	cd test/tp; Run.tp
	cd test/sustiva; Run.sustiva
	cd test/nad; Run.nad_bcc

RELEASE=1.25
TARFILE=antechamber/antechamber-$(RELEASE).tar
release:
#   make a tar file for "release"
#   Be sure to do "make clean" before running this
#
	( cd ..; \
	ln -s antechamber antechamber-$(RELEASE); \
	echo making $(TARFILE); \
	tar cvfX $(TARFILE) antechamber/excludes ./antechamber-$(RELEASE)/*; \
	gzip -v $(TARFILE); \
	rm antechamber-$(RELEASE) )
