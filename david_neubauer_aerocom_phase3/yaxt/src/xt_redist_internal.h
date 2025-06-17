/**
 * @file xt_redist_internal.h
 * @brief redistribution of data, non-public declarations
 *
 * contains declaration the redistribution data structure, which
 * is derived from one or more xt_xmaps
 *
 * @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://doc.redmine.dkrz.de/yaxt/html/
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef XT_REDIST_INTERNAL_H
#define XT_REDIST_INTERNAL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>

#include "xt/xt_redist.h"

struct Xt_redist_msg {

  int rank;
  MPI_Datatype datatype;
};

struct xt_redist_vtable {

  void (*delete)(Xt_redist);
  void (*s_exchange)(Xt_redist, int, const void **, void **);
  void (*s_exchange1)(Xt_redist, const void *, void *);
  MPI_Datatype (*get_send_MPI_Datatype)(Xt_redist, int);
  MPI_Datatype (*get_recv_MPI_Datatype)(Xt_redist, int);
  MPI_Comm (*get_MPI_Comm)(Xt_redist);
};

struct Xt_redist_ {
  struct xt_redist_vtable *vtable;
};

extern void
(*xt_redist_s_exchange_internal)(const void *src_data, void *dst_data,
                                 int nsend, int nrecv,
                                 struct Xt_redist_msg * send_msgs,
                                 struct Xt_redist_msg * recv_msgs,
                                 MPI_Comm comm);

void
xt_redist_msgs_strided_copy(size_t n,
                            struct Xt_redist_msg *restrict src,
                            size_t src_stride,
                            struct Xt_redist_msg *restrict dst,
                            size_t dst_stride,
                            MPI_Comm comm);

void xt_redist_msgs_strided_destruct(size_t n, struct Xt_redist_msg *msgs,
                                     MPI_Comm comm, size_t ofs_stride);

static inline void
xt_redist_msgs_free(size_t n, struct Xt_redist_msg *msgs, MPI_Comm comm)
{
  xt_redist_msgs_strided_destruct(n, msgs, comm, sizeof (*msgs));
  free(msgs);
}


#endif
