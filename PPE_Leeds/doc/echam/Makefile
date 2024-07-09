SHELL = /bin/ksh

BIBLIOGRAPHY = echam6.bib echam.bst babelbst.tex englishbst.tex
BIBLIOGRAPHY_GUIDE =  guide/echam6_userguide_lit.bib \
	              guide/cr/cr2014_07_15_rjs/cr2014_07_15_rjs.bib \
	              guide/cr/cr2014_05_09_rjs/cr2014_05_09_rjs.bib

SCIENCE_SOURCES = echam6_scidoc.tex atmosphere/cloud.tex		\
	atmosphere/gwdrag.tex atmosphere/sso.tex			\
	atmosphere/convection.tex atmosphere/methox.tex			\
	atmosphere/surface.tex atmosphere/dynamic.tex			\
	atmosphere/radiation.tex atmosphere/transport.tex		\
	atmosphere/hdiff.tex                                            \
	atmosphere/ECHAM6_resolutions.tex atmosphere/slab_ocean.tex	\
	atmosphere/vdiff.tex atmosphere/seaice.tex                      \
	externaldata/aerosols.tex                                       \
	externaldata/ozone.tex externaldata/introduction.tex		\
	externaldata/solar.tex externaldata/orbits.tex			\
	externaldata/orbtab.tex atmosphere/baseqn.tex                   \
	atmosphere/erratum.tex

SCIENCE_FIGURES = externaldata/tsi_1850-2100.pdf		\
	science/dec_1870-2150.pdf science/product.pdf		\
	science/obliquity.pdf science/ra_1870-2150.pdf		\
	science/dec_1978-1996.pdf science/eccentricity.pdf	\
	science/ga_1870-2150.pdf science/pentagonal.pdf		\
	science/ra_1978-1996.pdf science/dis_1870-2150.pdf	\
	science/ga_1978-1996.pdf science/precession.pdf		\
	science/dis_1978-1996.pdf science/triangular.pdf        \
	science/jgrd14967-fig-0001.png science/jgrd14967-fig-0002.png \
        science/jgrd14967-fig-0003.png

GUIDE_SOURCES = echam6_userguide.tex guide/ECHAM.guide.compiling.tex	\
	guide/ECHAM.guide.running.tex                                   \
	guide/ECHAM.guide.namelistother.tex				\
	guide/ECHAM.guide.inputfiles.tex				\
	guide/ECHAM.guide.outputfiles.tex				\
	guide/ECHAM.guide.introduction.tex				\
	guide/ECHAM.guide.parallel.tex guide/ECHAM.guide.memory.tex	\
	guide/ECHAM.guide.postprocessing.tex				\
	guide/ECHAM.guide.namelistecham.tex				\
	guide/ECHAM.guide.runscripts.tex				\
	guide/ECHAM.guide.submodel.tex					\
	guide/ECHAM.guide.timemanager.tex                               \
	guide/ECHAM.guide.appendix.tex                                  \
	guide/cr/cr2009_03_31_rjs/cr2009_03_31_rjs.tex                  \
	guide/cr/cr2014_05_09_rjs/cr2014_05_09_rjs.tex                  \
	guide/cr/cr2014_07_15_rjs/cr2014_07_15_rjs.tex                  \
        ../jsbach/namelists/namelists_jsbach.tex

GUIDE_FIGURES = guide/transpose.pdf guide/cr/*/*.pdf

all: scidoc guide

scidoc: echam6_scidoc.pdf
guide: echam6_userguide.pdf

echam6_scidoc.pdf: $(SCIENCE_SOURCES) $(SCIENCE_FIGURES) $(BIBLIOGRAPHY)
	pdflatex echam6_scidoc
	bibtex echam6_scidoc
	pdflatex echam6_scidoc
	pdflatex echam6_scidoc
	pdflatex echam6_scidoc

echam6_userguide.pdf:  $(GUIDE_SOURCES) $(GUIDE_FIGURES) $(BIBLIOGRAPHY_GUIDE)
	./update_version.sh
	pdflatex echam6_userguide
	makeindex echam6_userguide
	bibtex echam6_userguide
	pdflatex echam6_userguide
	pdflatex echam6_userguide

docclean:
	rm -f echam6_scidoc.aux echam6_scidoc.bbl echam6_scidoc.blg echam6_scidoc.idx \
	      echam6_scidoc.log echam6_scidoc.toc atmosphere/*aux externaldata/*aux 

guideclean:	
	rm -f echam6_userguide.aux echam6_userguide.bbl echam6_userguide.blg echam6_userguide.idx \
	      echam6_userguide.log echam6_userguide.lot echam6_userguide.out echam6_userguide.toc echam6_userguide.pdf

clean:	
	rm -f *.aux *.bbl *.blg *.brf *.idx *.log *.lot *.out *.pdf *.toc \
	      atmosphere/*.aux externaldata/*.aux \
