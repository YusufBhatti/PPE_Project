#
# Keywords:
# Maintainer: Jörg Behrens <behrens@dkrz.de>
#             Moritz Hanke <hanke@dkrz.de>
#             Thomas Jahns <jahns@dkrz.de>
# URL: https://doc.redmine.dkrz.de/yaxt/html/
#
# Redistribution and use in source and binary forms, with or without
# modification, are  permitted provided that the following conditions are
# met:
#
# Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# Neither the name of the DKRZ GmbH nor the names of its contributors
# may be used to endorse or promote products derived from this software
# without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

AM_CPPFLAGS = -I../include -I$(srcdir)/../include $(MPI_C_INCLUDE)
AM_FCFLAGS  = $(MPI_FC_INCLUDE) $(MPI_FC_MOD)
lib_LTLIBRARIES = libyaxt.la

# headers from subdir xt will be installed (after substituting
# definitions from config.h). These must be listed in a sequence in
# which they can be combined in yaxt.h.
libyaxt_la_SOURCES =						\
	xt/xt_core.h						\
	xt/xt_handles.h						\
	xt/sort_common.h					\
	xt/mergesort.h						\
	xt/quicksort.h						\
	xt/xt_sort.h						\
	xt_sort_internal.h					\
	xt/xt_idxlist.h						\
	xt/xt_idxempty.h					\
	xt/xt_idxlist_collection.h				\
	xt/xt_idxsection.h					\
	xt/xt_idxstripes.h					\
	xt/xt_idxvec.h 						\
	xt/xt_mpi.h						\
	xt/xt_redist.h						\
	xt/xt_redist_collection.h				\
	xt/xt_redist_collection_static.h			\
	xt/xt_redist_p2p.h					\
	xt/xt_redist_repeat.h					\
	xt/xt_stripe.h						\
	xt/xt_xmap.h						\
	xt/xt_xmap_all2all.h					\
	xt/xt_xmap_dist_dir.h					\
	xt/xt_xmap_intersection.h				\
	xt/xt_idxmod.h						\
	xt_arithmetic_util.h					\
	xt_cover.h						\
	xt_cover.c						\
	xt_handles.c						\
	xt_idxempty_internal.h					\
	xt_idxempty.c						\
	xt_idxlist_internal.h					\
	xt_idxlist.c						\
	xt_idxlist_intersection.c				\
	xt_idxlist_collection.c					\
	xt_idxlist_collection_f.f90				\
	xt_idxlist_collection_internal.h			\
	xt_idxlist_unpack.h					\
	xt_idxlist_unpack.c 					\
	xt_idxmod.c						\
	xt_idxsection_internal.h				\
	xt_idxsection.c						\
	xt_idxstripes_internal.h				\
	xt_idxstripes.c						\
	xt_idxvec.c 						\
	xt_idxvec_internal.h					\
	xt_mpi_internal.h					\
	xt_mpi.c						\
	xt_redist_internal.h					\
	xt_redist.c						\
	xt_redist_collection.c					\
	xt_redist_collection_static.c				\
	xt_redist_p2p.c						\
	xt_redist_repeat.c					\
	xt_redist_single_array_base.c				\
	xt_redist_single_array_base.h				\
	xt_stripe.c						\
	xt_stripe_util.h					\
	xt_xmap_internal.h					\
	xt_xmap.c						\
	xt_xmap_all2all.c					\
	xt_xmap_dist_dir.c					\
	xt_xmap_intersection.c					\
	xt_xmap_intersection_ext.c				\
	xt_xmap_intersection_f.f90				\
	xt_init.c						\
	xt_ut.f90						\
	xt_ut_c.c						\
	xt_ut_c.h						\
	core/core.c						\
	core/core.h						\
	core/ppm_xfuncs.h					\
	core/symprefix.h					\
	core/xmalloc.c						\
	core/xstdio.c						\
	ensure_array_size.c					\
	ensure_array_size.h					\
	quicksort.c						\
	mergesort.c						\
	instr.h							\
	xt_sort.c						\
	yaxt_f2c.c						\
	xt_mpi_f.f90						\
	xt_core_f.f90						\
	xt_idxlist_f.f90					\
	xt_idxvec_f.f90						\
	xt_idxstripes_f.f90					\
	xt_idxsection_f.f90					\
	xt_xmap_f.f90						\
	xt_redist_f.f90						\
	yaxt.f90						\
	xt_exchanger.c						\
	xt_exchanger.h						\
	xt_exchanger_simple_base.c				\
	xt_exchanger_simple_base.h				\
	xt_exchanger_irecv_send.c				\
	xt_exchanger_irecv_send.h				\
	xt_exchanger_irecv_isend.c				\
	xt_exchanger_irecv_isend.h				\
	xt_exchanger_irecv_isend_packed.c			\
	xt_exchanger_irecv_isend_packed.h			\
	xt_exchanger_mix_isend_irecv.c				\
	xt_exchanger_mix_isend_irecv.h

libyaxt_la_LIBADD = $(MPI_FC_LIB)

xtinclude_Makefile.inc: Makefile
	headers=`echo $(libyaxt_la_SOURCES) | sed 's/[ 	][ 	]*/ /g' | tr ' ' '\n' | grep 'xt/.*\.h$$'` ; \
	  echo "xtinclude_HEADERS = \\" >$@ ; \
	  for header in $$headers ; do echo "	$$header \\" ; done \
	    | sed '$$s/\\$$//' >>$@

all-local: xtinclude_Makefile.inc

clean-local:
	-rm -f *.$(FCMODEXT) i.*.L

distclean-local:
	-rm -f xtinclude_Makefile.inc

install-exec-local:
	$(mkinstalldirs) "$(DESTDIR)$(libdir)/pkgconfig"
	$(install_sh_DATA) pkgconfig/yaxt.pc \
		"$(DESTDIR)$(libdir)/pkgconfig/yaxt.pc"

uninstall-local:
	rm -f "$(DESTDIR)$(libdir)/pkgconfig/yaxt.pc"
	-rmdir "$(DESTDIR)$(libdir)/pkgconfig"

./$(DEPDIR)/FC.deps: $(SOURCES) Makefile
	@unique=`for i in $(SOURCES); do \
	  if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | grep '\.f90$$' ` ; \
	  FC="$(FC)" FPP="$(FPP)" \
	  $(PERL) -I$(top_srcdir)/scripts \
	  $(top_srcdir)/scripts/makef90depends --ignore-module iso_c_binding \
	  --fc-mod-flag=$(FC_MOD_FLAG) --fc-def-opt=$(FPP_DEFOPT) \
	  --fc-inc-flag=$(FPP_INCOPT) --fc-mod-suf='.$(FCMODEXT)' \
	  --fc-mod-case=$(FCMODCASE) --prog-fpp='$(FPP)' \
	  --`test $(C_LONG_IS_LARGER_THAN_INT) = 1 || echo \
	  no`fpp-long-is-larger-than-int \
	  --`test $(C_LONG_LONG_IS_LARGER_THAN_LONG) = 1 || echo \
	  no`fpp-long-long-is-larger-than-long \
	  --fpp-char-is-`test $(C_CHAR_IS_UNSIGNED) = 0 || echo un`signed \
	  --strip-obj-dirname --parse-cache=./$(DEPDIR)/FC.depdb \
	  --obj-suf='.lo' -- \
	  $(AM_FCFLAGS) $(FCFLAGS) $(FPPFLAGS) -- $$unique \
	  >./$(DEPDIR)/FC.deps.temp
	cmp ./$(DEPDIR)/FC.deps.temp ./$(DEPDIR)/FC.deps >/dev/null 2>&1 \
	  || ln -f ./$(DEPDIR)/FC.deps.temp ./$(DEPDIR)/FC.deps

@AMDEP_TRUE@@am__include@ @am__quote@./$(DEPDIR)/FC.deps@am__quote@
