SUBDIRS = . src inst_headers tests examples perf
DIST_SUBDIRS = src inst_headers tests examples perf .

ACLOCAL_AMFLAGS=-I m4

EXTRA_DIST = util/mpi_launch_poe util/sunf95preproc-wrapper \
	util/sxpreproc-wrapper util/xlfpreproc-wrapper \
	util/crayftnpreproc-wrapper \
	util/pbs_aprun \
	config/checksrc/openmpi_struct_dt.c \
	config/checkdoc/openmpi_struct_dt.txt \
	config/checksrc/openmpi_dup.c \
	config/checkdoc/openmpi_dup.txt \
	config/checksrc/pe5228_iv32584.c \
	config/checkdoc/pe5228_iv32584.txt \
	contrib/00nagfor-libtool-patch/README.txt \
	contrib/00nagfor-libtool-patch/nagfor-libtool.patch \
	contrib/00nagfor-libtool-patch/nagfor-libtool-2.4.2.patch \
	contrib/00nagfor-libtool-patch/nagfor-libtool-2.4.6.patch \
	contrib/01aix-deplib-rpath-patch/README \
	contrib/01aix-deplib-rpath-patch/aix-deplib-libtool.patch \
	contrib/02nagfor53-shared-patch/README.txt \
	contrib/02nagfor53-shared-patch/nagfor53-shared.patch \
	doc/Makefile \
	doc/Doxyfile.in doc/README.build doc/cfortran.doc \
	doc/cfortran.html doc/cfortran.txt doc/styleguide.txt \
	doc/images/decomp_idxlist_collection.png \
	doc/images/decomp_idxlist_collection.svg \
	doc/images/decomp_idxmod.png \
	doc/images/decomp_idxsection.png \
	doc/images/decomp_idxsection.svg \
	doc/images/decomp_idxstripes.png \
	doc/images/decomp_idxstripes.svg \
	doc/images/decomp_idxvec.png \
	doc/images/decomp_idxvec.svg \
	doc/images/neg_size_idxsection_nnnn.png \
	doc/images/neg_size_idxsection_nnnp.png \
	doc/images/neg_size_idxsection_nnpn.png \
	doc/images/neg_size_idxsection_nnpp.png \
	doc/images/neg_size_idxsection_npnn.png \
	doc/images/neg_size_idxsection_npnp.png \
	doc/images/neg_size_idxsection_nppn.png \
	doc/images/neg_size_idxsection_nppp.png \
	doc/images/neg_size_idxsection_pnnn.png \
	doc/images/neg_size_idxsection_pnnp.png \
	doc/images/neg_size_idxsection_pnpn.png \
	doc/images/neg_size_idxsection_pnpp.png \
	doc/images/neg_size_idxsection_ppnn.png \
	doc/images/neg_size_idxsection_ppnp.png \
	doc/images/neg_size_idxsection_pppn.png \
	doc/images/neg_size_idxsection_pppp.png \
	doc/images/stride_idxsection_x.png \
	doc/images/stride_idxsection_xy.png \
	doc/src/idxlists.dox \
	doc/src/idxsection.dox \
	doc/src/mainpage.dox \
	doc/src/redist_colls.dox \
	doc/src/rrobin.dox \
	doc/src/xmaps.dox \
	doc/src/yaxt_concept.dox \
	include/cfortran.h \
	m4/tj_find_type.m4 \
	scripts/cpp.pm \
	scripts/header2installheader.pl \
	scripts/generate_doc.pl \
	scripts/makef90depends


all-local:
	@for d in $(SUBDIRS) ; do \
          test $$d = . -o $(srcdir)/$$d/Makefile.in -ot $$d/Makefile \
	    || ./config.status $$d/Makefile ; done

clean-local:
	@for d in $(SUBDIRS) ; do \
          test $$d = . -o $(srcdir)/$$d/Makefile.in -ot $$d/Makefile \
	    || ./config.status $$d/Makefile ; done

distclean-local:
	for d in $(SUBDIRS) ; do \
		test $$d = . || rm -f $$d/Makefile ; done

install-exec-local:
	@for d in $(SUBDIRS); do \
	  test $$d = . -o $(srcdir)/$$d/Makefile.in -ot $$d/Makefile || \
	  ./config.status $$d/Makefile ; done

install-data-local:
	@for d in $(SUBDIRS); do \
	  test $$d = . -o $(srcdir)/$$d/Makefile.in -ot $$d/Makefile || \
	  ./config.status $$d/Makefile ; done

check-local:
	@for d in $(SUBDIRS); do \
	  test $$d = . -o $(srcdir)/$$d/Makefile.in -ot $$d/Makefile || \
	  ./config.status $$d/Makefile ; done
