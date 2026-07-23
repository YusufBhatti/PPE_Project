/**
 * @file xt_exchanger_irecv_send.c
 *
 * @copyright Copyright  (C)  2014 Jörg Behrens <behrens@dkrz.de>
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt_redist_internal.h"
#include "xt_exchanger_irecv_send.h"
#include "xt_exchanger_simple_base.h"

static void
xt_exchanger_irecv_send_s_exchange(const void *src_data, void *dst_data,
                                   int nsend, int nrecv,
                                   struct Xt_redist_msg * send_msgs,
                                   struct Xt_redist_msg * recv_msgs,
                                   int tag_offset, MPI_Comm comm) {

  MPI_Request *recv_request
    = xmalloc((size_t)nrecv * sizeof (*recv_request));

  for (int i = 0; i < nrecv; ++i)
    xt_mpi_call(MPI_Irecv(dst_data, 1, recv_msgs[i].datatype,
                          recv_msgs[i].rank,
                          tag_offset + xt_mpi_tag_exchange_msg, comm,
                          recv_request+i), comm);

  for (int i = 0; i < nsend; ++i)
    xt_mpi_call(MPI_Send((void *)src_data, 1, send_msgs[i].datatype,
                         send_msgs[i].rank,
                         tag_offset + xt_mpi_tag_exchange_msg, comm),
                comm);

  xt_mpi_call(MPI_Waitall(nrecv, recv_request, MPI_STATUSES_IGNORE), comm);

  free(recv_request);
}

Xt_exchanger
xt_exchanger_irecv_send_new(int nsend, int nrecv,
                            struct Xt_redist_msg * send_msgs,
                            struct Xt_redist_msg * recv_msgs,
                            MPI_Comm comm, int tag_offset) {

  return xt_exchanger_simple_base_new(nsend, nrecv, send_msgs, recv_msgs,
                                      comm, tag_offset,
                                      xt_exchanger_irecv_send_s_exchange);
}
