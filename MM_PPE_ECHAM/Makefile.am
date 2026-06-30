SUBDIRS =
DIST_SUBDIRS =

if BUILD_INTERNAL_YAXT
SUBDIRS += yaxt
DIST_SUBDIRS += yaxt
endif

if BUILD_INTERNAL_CDI
SUBDIRS += cdi
DIST_SUBDIRS += cdi
endif

if BUILD_INTERNAL_LAPACK
SUBDIRS += blas lapack
DIST_SUBDIRS += blas lapack
endif

SUBDIRS += support src util
DIST_SUBDIRS += support src util

ACLOCAL_AMFLAGS = -I m4

program_transform_name = @program_transform_name@

index:
	-rm -f html/[a-z]*
	util/f2html.pl -f util/fgenrc -d html src include
