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
	rrobin							\
	idxlist_examples					\
	idxsection_examples					\
	row2col							\
	row2col_f						\
	tr_gp2fs						\
	unstruct_halo_exchange_f

rrobin_SOURCES = rrobin.c

idxlist_examples_SOURCES = idxlist_examples.c print_index_list.h

idxsection_examples_SOURCES = idxsection_examples.c print_index_list.h

row2col_SOURCES = row2col.c

row2col_f_SOURCES = row2col_f.f90
row2col_f_LDADD = ../src/libyaxt.la $(MPI_FC_LIB)

tr_gp2fs_SOURCES = tr_gp2fs.c
tr_gp2fs_LDADD = $(LDADD) -lm

unstruct_halo_exchange_f_SOURCES = unstruct_halo_exchange_f.f90

AM_CFLAGS = -I../inst_headers $(MPI_C_INCLUDE)
AM_FCFLAGS = $(FC_MOD_FLAG)../inst_headers/f90 $(MPI_FC_INCLUDE) $(MPI_FC_MOD)
LDADD = ../src/libyaxt.la $(MPI_C_LIB)
