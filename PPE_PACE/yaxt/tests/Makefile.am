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

noinst_PROGRAMS = 						\
	test_mpi_generate_datatype				\
	test_mpi_smartdedup					\
	test_exchanger_parallel					\
	test_handles						\
	test_idxempty						\
	test_idxempty_f						\
	test_idxvec						\
	test_idxvec_f						\
	test_idxlist_collection					\
	test_idxlist_collection_f				\
	test_idxmod						\
	test_idxmod_f						\
	test_idxsection						\
	test_idxsection_f					\
	test_idxstripes						\
	test_idxstripes_f					\
	test_redist_collection					\
	test_redist_collection_f				\
	test_redist_collection_parallel				\
	test_redist_collection_parallel_f			\
	test_redist_collection_static				\
	test_redist_collection_static_f				\
	test_redist_collection_static_parallel			\
	test_redist_collection_static_parallel_f		\
	test_redist_p2p						\
	test_redist_p2p_f					\
	test_redist_p2p_parallel				\
	test_redist_p2p_parallel_f				\
	test_redist_repeat					\
	test_redist_repeat_f					\
	test_redist_repeat_parallel				\
	test_redist_repeat_parallel_f				\
	test_xmap_all2all					\
	test_xmap_all2all_f					\
	test_xmap_all2all_fail					\
	test_xmap_all2all_fail_f				\
	test_xmap_all2all_parallel				\
	test_xmap_all2all_parallel_f				\
	test_xmap_dist_dir					\
	test_xmap_dist_dir_f					\
	test_xmap_dist_dir_parallel				\
	test_xmap_dist_dir_parallel_f				\
	test_xmap_intersection_parallel				\
	test_xmap_intersection_parallel_f			\
	test_initialized_finalized				\
	test_initialized_finalized_f				\
	test_ut							\
	test_perf						\
	test_perf_stripes					\
	test_sort						\
	test_uid						\
	test_yaxt


noinst_LIBRARIES = libtestutil.a

libtestutil_a_SOURCES = tests.c ftest_common.f90 \
	test_idxlist_utils.c test_idxlist_utils_f.f90 \
	test_xmap_common.c test_xmap_common.h test_xmap_common_f.f90 \
	test_xmap_common_parallel.c test_xmap_common_parallel_f.f90 \
	test_redist_common_f.f90

test_exchanger_parallel_SOURCES = test_exchanger_parallel.c tests.h
test_idxempty_SOURCES = test_idxempty.c tests.h test_idxlist_utils.h
test_idxempty_f_SOURCES = test_idxempty_f.f90
test_idxempty_f_LDADD = $(XT_FC_LDADD)
test_idxvec_SOURCES = test_idxvec.c tests.h test_idxlist_utils.h
test_idxvec_f_SOURCES = test_idxvec_f.f90
test_idxvec_f_LDADD = $(XT_FC_LDADD)
test_idxlist_collection_SOURCES = test_idxlist_collection.c tests.h \
	test_idxlist_utils.h
test_idxlist_collection_f_SOURCES = test_idxlist_collection_f.f90
test_idxlist_collection_f_LDADD = $(XT_FC_LDADD)
test_idxmod_SOURCES = test_idxmod.c tests.h
test_idxmod_f_SOURCES = test_idxmod_f.f90
test_idxmod_f_LDADD = $(XT_FC_LDADD)
test_idxsection_SOURCES = test_idxsection.c tests.h test_idxlist_utils.h
test_idxsection_f_SOURCES = test_idxsection_f.f90
test_idxsection_f_LDADD = $(XT_FC_LDADD)
test_idxstripes_SOURCES = test_idxstripes.c tests.h test_idxlist_utils.h
test_idxstripes_f_SOURCES = test_idxstripes_f.f90
test_idxstripes_f_LDADD = $(XT_FC_LDADD)
test_redist_collection_SOURCES = test_redist_collection.c tests.h
test_redist_collection_f_SOURCES = test_redist_collection_f.f90
test_redist_collection_f_LDADD = $(XT_FC_LDADD)
test_redist_collection_parallel_SOURCES = test_redist_collection_parallel.c \
	tests.h
test_redist_collection_parallel_f_SOURCES = \
	test_redist_collection_parallel_f.f90
test_redist_collection_parallel_f_LDADD = $(XT_FC_LDADD)
test_redist_collection_static_SOURCES = test_redist_collection_static.c tests.h
test_redist_collection_static_f_SOURCES = test_redist_collection_static_f.f90
test_redist_collection_static_f_LDADD = $(XT_FC_LDADD)
test_redist_collection_static_parallel_SOURCES = \
	test_redist_collection_static_parallel.c tests.h
test_redist_collection_static_parallel_f_SOURCES = \
	test_redist_collection_static_parallel_f.f90
test_redist_collection_static_parallel_f_LDADD = $(XT_FC_LDADD)
test_redist_p2p_SOURCES = test_redist_p2p.c tests.h
test_redist_p2p_f_SOURCES = test_redist_p2p_f.f90
test_redist_p2p_f_LDADD = $(XT_FC_LDADD)
test_redist_p2p_parallel_SOURCES = test_redist_p2p_parallel.c tests.h
test_redist_p2p_parallel_f_SOURCES = test_redist_p2p_parallel_f.f90
test_redist_p2p_parallel_f_LDADD = $(XT_FC_LDADD)
test_redist_repeat_SOURCES = test_redist_repeat.c tests.h
test_redist_repeat_f_SOURCES = test_redist_repeat_f.f90
test_redist_repeat_f_LDADD = $(XT_FC_LDADD)
test_redist_repeat_parallel_SOURCES = test_redist_repeat_parallel.c tests.h
test_redist_repeat_parallel_f_SOURCES = test_redist_repeat_parallel_f.f90
test_redist_repeat_parallel_f_LDADD = $(XT_FC_LDADD)
test_xmap_all2all_SOURCES = test_xmap_all2all.c test_xmap_common.h
test_xmap_all2all_f_SOURCES = test_xmap_all2all_f.f90
test_xmap_all2all_f_LDADD = $(XT_FC_LDADD)
test_xmap_all2all_fail_SOURCES = test_xmap_all2all_fail.c tests.h
test_xmap_all2all_fail_f_SOURCES = test_xmap_all2all_fail_f.f90
test_xmap_all2all_fail_f_LDADD = $(XT_FC_LDADD)
test_xmap_all2all_parallel_SOURCES = test_xmap_all2all_parallel.c tests.h
test_xmap_all2all_parallel_f_SOURCES = test_xmap_all2all_parallel_f.f90
test_xmap_all2all_parallel_f_LDADD = $(XT_FC_LDADD)
test_xmap_dist_dir_SOURCES = test_xmap_dist_dir.c test_xmap_common.h
test_xmap_dist_dir_f_SOURCES = test_xmap_dist_dir_f.f90
test_xmap_dist_dir_f_LDADD = $(XT_FC_LDADD)
test_xmap_dist_dir_parallel_SOURCES = test_xmap_dist_dir_parallel.c tests.h
test_xmap_dist_dir_parallel_f_SOURCES = test_xmap_dist_dir_parallel_f.f90
test_xmap_dist_dir_parallel_f_LDADD = $(XT_FC_LDADD)
test_xmap_intersection_parallel_SOURCES = test_xmap_intersection_parallel.c tests.h
test_xmap_intersection_parallel_f_SOURCES = test_xmap_intersection_parallel_f.f90
test_xmap_intersection_parallel_f_LDADD = $(XT_FC_LDADD)
test_initialized_finalized_SOURCES = test_initialized_finalized.c tests.h
test_initialized_finalized_f_SOURCES = test_initialized_finalized_f.f90
test_initialized_finalized_f_LDADD = $(XT_FC_LDADD)
test_handles_SOURCES = test_handles.c tests.h
test_ut_SOURCES = test_ut.f90
test_ut_LDADD = $(XT_FC_LDADD)
test_perf_SOURCES = test_perf.f90
test_perf_LDADD = $(XT_FC_LDADD)
test_perf_stripes_SOURCES = test_perf_stripes.f90 ftest_common.f90
test_perf_stripes_LDADD = $(XT_FC_LDADD)
test_sort_SOURCES = test_sort.c tests.h
test_uid_SOURCES = test_uid.c
test_yaxt_SOURCES = test_yaxt.f90 ftest_common.f90
test_yaxt_LDADD = $(XT_FC_LDADD)
test_mpi_generate_datatype_SOURCES = test_mpi_generate_datatype.c tests.h
test_mpi_smartdedup_SOURCES = test_mpi_smartdedup.c tests.h

AM_CFLAGS = -I../inst_headers $(MPI_C_INCLUDE)
AM_FCFLAGS = $(FC_MOD_FLAG)../inst_headers/f90 $(MPI_FC_INCLUDE) $(MPI_FC_MOD)
LDADD = libtestutil.a ../src/libyaxt.la $(MPI_C_LIB)
XT_FC_LDADD = libtestutil.a ../src/libyaxt.la $(MPI_FC_LIB)

TESTS = 							\
	test_mpi_generate_datatype_run				\
	test_mpi_smartdedup_run					\
	test_exchanger_parallel_run				\
	test_handles_run					\
	test_idxempty_run					\
	test_idxvec_run						\
	test_idxlist_collection_run				\
	test_idxmod_run						\
	test_idxsection_run					\
	test_idxstripes_run					\
	test_redist_p2p_run					\
	test_redist_collection_run				\
	test_redist_collection_parallel_run			\
	test_redist_collection_static_run			\
	test_redist_collection_static_parallel_run		\
	test_redist_p2p_parallel_run				\
	test_redist_repeat_run					\
	test_redist_repeat_parallel_run				\
	test_xmap_all2all_run					\
	test_xmap_all2all_fail_run				\
	test_xmap_all2all_parallel_run				\
	test_xmap_dist_dir_run					\
	test_xmap_dist_dir_parallel_run				\
	test_xmap_intersection_parallel_run			\
	test_initialized_finalized_run				\
	test_exported_symbols					\
	test_ut_run						\
	test_perf_run						\
	test_perf_stripes_run					\
	test_sort_run						\
	test_uid_run						\
	test_yaxt_run

XFAIL_TESTS =                                                   \
	test_xmap_all2all_fail_run


AUTOMAKE_OPTIONS = color-tests

EXTRA_DIST = core

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
	  --obj-suf='.$(OBJEXT)' -- \
	  $(AM_FCFLAGS) $(FCFLAGS) $(FPPFLAGS) -- $$unique \
	  >./$(DEPDIR)/FC.deps.temp
	cmp ./$(DEPDIR)/FC.deps.temp ./$(DEPDIR)/FC.deps >/dev/null 2>&1 \
	  || ln -f ./$(DEPDIR)/FC.deps.temp ./$(DEPDIR)/FC.deps

@AMDEP_TRUE@@am__include@ @am__quote@./$(DEPDIR)/FC.deps@am__quote@

clean-local:
	-rm -f *.$(FCMODEXT)
