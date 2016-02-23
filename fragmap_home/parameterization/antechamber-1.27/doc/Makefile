#
#		make doc		! make an ascii version
#		make x.doc		! single section .doc file
#		make ps			! all documentation in postscript form., with 
#                         table of contents and index
#       make x.ps		! single section .ps file
#

SHELL = /bin/sh
.SUFFIXES: .me .doc .ps .html
#
#  set site-specific ROFF stuff
#
NROFF  = troff -Tascii
TROFF  = troff -U
REFER  = refer
PIC    = pic
TBL    = tbl
EQN    = eqn
NEQN   = eqn -Tascii
GROTTY = grotty -u -b -f -o
GROPS  = grops

SRC = 	start.me antechamber.me tleap.me

ps:		AC.ps 
		/bin/rm -f AC.me

AC.me:	$(SRC)
		cat $(SRC) > AC.me
		./mk_index AC.me
		cat AC.me.index >> AC.me
		cat toc.me >> AC.me
		/bin/rm AC.me.index

AC.html:    $(SRC)
		cat $(SRC) > AC.me
		troff2html -split 2 -toc bottom -eqn ps -tbl ps AC.me
		perl -pi -e "s/AC Users' Manual/\<H1\>AC Users' Manual\<\/H1\>/;" AC.html
#		insertsearch
		perl -pi -e 's/^Table of Contents/AC Documentation/;s?^</ul>?</ul>\n<p>?' $@
		-mkdir www
		/bin/mv AC*html www
		/bin/mv AC*gif www

.me.html:    $<
		troff2html -split 2 -toc bottom -eqn ps -tbl nroff $<
		perl -pi -e 's/<body>/<body bgcolor="#A0F0F0" link="#EE0000" alink="#0000CC" vlink="#0000CC">/;s/_\[/</g;s/_\]/>/g;'  *.html
		perl -pi -e 's/^Table of Contents/AC Documentation/;s?^</ul>?</ul>\n<p>?' $@

.me.doc:	$<
		cat nroff.dot ac.mac $< | $(REFER) | $(PIC) | $(TBL) | $(NEQN) | \
			$(NROFF) -me 2> $@.errs | $(GROTTY) > $@ 
		/bin/rm $@.errs

.me.ps:		$<
		cat ac.mac $< | $(PIC) | $(REFER) | $(TBL) | $(EQN) -Tps | \
			$(TROFF) -me -Tps 2> $@.errs | $(GROPS) | fix_bookmarks > $@ 
		/bin/rm $@.errs

clean::
	-/bin/rm www/AC-* www/AC.html
