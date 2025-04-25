xtincludedir = $(includedir)/xt

nodist_include_HEADERS = yaxt.h

BUILT_SOURCES = yaxt.h

SUBDIRS = f90

../src/xtinclude_Makefile.inc:
	(cd ../src ; if test -r Makefile ; then \
	  $(MAKE) xtinclude_Makefile.inc ; fi)

@AMDEP_TRUE@@am__include@ @am__quote@../src/xtinclude_Makefile.inc@am__quote@

yaxt.h: Makefile $(xtinclude_HEADERS)
	( echo '/* license */' ; echo '#ifndef YAXT_H' ; echo '#define YAXT_H' ; \
	for header in $(xtinclude_HEADERS) ; \
	  do echo "#include <$$header>" ; done ; \
	echo '#endif' \
	) >$@.$$$$.tmp ; \
	if test ! -f $@ || ! cmp $@ $@.$$$$.tmp >/dev/null ; then \
	mv $@.$$$$.tmp $@ ; else rm $@.$$$$.tmp ; fi

$(xtinclude_HEADERS): ../src/$@ ../include/config.h \
	$(top_srcdir)/scripts/header2installheader.pl
	$(MKDIR_P) xt
	headerbasedir=`test -f '../src/$@' || echo '$(srcdir)/'` ; \
	  $(PERL) -I$(top_srcdir)/scripts \
	  $(top_srcdir)/scripts/header2installheader.pl \
	  --config-header=../include/config.h \
	  --srcdir=$$headerbasedir../src \
	  --`test $(C_LONG_IS_LARGER_THAN_INT) = 1 || echo \
	  no`c-long-is-larger-than-int \
	  --`test $(C_LONG_LONG_IS_LARGER_THAN_LONG) = 1 || echo \
	  no`c-long-long-is-larger-than-long \
	  --c-char-is-`test $(C_CHAR_IS_UNSIGNED) = 0 || echo un`signed \
	  --dstdir=. \
	  $$headerbasedir../src/$@

clean-local:
	rm -rf xt yaxt.h

install-data-local: yaxt.h $(xtinclude_HEADERS)
	@echo " $(MKDIR_P) '$(DESTDIR)$(xtincludedir)'"; \
          $(MKDIR_P) "$(DESTDIR)$(xtincludedir)" || exit 1; \
	for p in $(xtinclude_HEADERS); do \
          echo " $(INSTALL_HEADER) $$p '$(DESTDIR)$(xtincludedir)'"; \
          $(INSTALL_HEADER) $$p "$(DESTDIR)$(xtincludedir)" || exit $$?; \
        done

uninstall-local:
	rm -rf "$(DESTDIR)$(xtincludedir)"

dist-hook:
	-rm $(distdir)/xt/*.h
	-rmdir $(distdir)/xt
