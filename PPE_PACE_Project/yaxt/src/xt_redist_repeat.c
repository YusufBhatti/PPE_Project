/**
 * @file xt_redist_repeat.c
 *
 * @copyright Copyright  (C)  2012 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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

#include <assert.h>
#include <limits.h>
#include <stdlib.h>

#include <mpi.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt/xt_redist_repeat.h"
#include "xt_redist_single_array_base.h"
#include "ensure_array_size.h"
#include "xt/xt_redist.h"
#include "xt_redist_internal.h"

static void
generate_msg_infos(struct Xt_redist_msg ** msgs, int * nmsgs,
                   MPI_Aint extent, int * displacements, Xt_redist redist,
                   int num_repetitions, MPI_Comm comm,
                   MPI_Datatype (*get_MPI_datatype)(Xt_redist,int)) {

  size_t msgs_array_size = 0;

  int comm_size;
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  MPI_Datatype datatype;

  assert(*nmsgs >= 0);
  size_t num_messages = (size_t)*nmsgs;
  struct Xt_redist_msg *p = *msgs;
  for (int i = 0; i < comm_size; ++i) {

    datatype = get_MPI_datatype(redist, i);

    if (datatype != MPI_DATATYPE_NULL) {

      MPI_Aint curr_lb, curr_extent, curr_true_lb, curr_true_extent;
      MPI_Datatype datatype_with_extent;

      // adjust extent of datatype to match the displacements
      xt_mpi_call(MPI_Type_get_extent(datatype, &curr_lb, &curr_extent), comm);
      xt_mpi_call(MPI_Type_get_true_extent(datatype, &curr_true_lb,
                                           &curr_true_extent), comm);
      if (curr_true_extent > extent)
        Xt_abort(comm, "ERROR: new datatype extent is too small (Xt_redist_repeat)\n",
                 __FILE__, __LINE__);
      xt_mpi_call(MPI_Type_create_resized(datatype, curr_lb, extent,
                                          &datatype_with_extent), comm);

      ENSURE_ARRAY_SIZE(p, msgs_array_size, num_messages+1);

      p[num_messages].rank = i;
      p[num_messages].datatype =
        xt_mpi_generate_datatype(displacements, num_repetitions,
                                 datatype_with_extent, comm);
      ++num_messages;

      MPI_Type_free(&datatype_with_extent);
      MPI_Type_free(&datatype);
    }
  }

  if (num_messages > 0)
    p = xrealloc(p, num_messages * sizeof(*p));
  *msgs = p;
  *nmsgs = (int)num_messages;
}

Xt_redist xt_redist_repeat_new(Xt_redist redist, MPI_Aint src_extent,
                               MPI_Aint dst_extent, int num_repetitions,
                               int displacements[num_repetitions]) {

  int nsend, nrecv;
  struct Xt_redist_msg * send_msgs, * recv_msgs;

  nsend = 0;
  nrecv = 0;
  send_msgs = NULL;
  recv_msgs = NULL;

  MPI_Comm comm;
  int tag_offset;
  comm = xt_mpi_comm_smart_dup(xt_redist_get_MPI_Comm(redist), &tag_offset);

  if (num_repetitions < 1)
    Xt_abort(comm, "ERROR: invalid number of repetitions (Xt_redist_repeat)\n",
             __FILE__, __LINE__);


  generate_msg_infos(&send_msgs, &nsend, src_extent,
                     displacements, redist, num_repetitions, comm,
                     xt_redist_get_send_MPI_Datatype);

  generate_msg_infos(&recv_msgs, &nrecv, dst_extent, 
                     displacements, redist, num_repetitions, comm,
                     xt_redist_get_recv_MPI_Datatype);

  Xt_redist result
    = xt_redist_single_array_base_new(nsend, nrecv, send_msgs, recv_msgs, comm);
  xt_mpi_comm_smart_dedup(&comm, tag_offset);
  return result;
}
