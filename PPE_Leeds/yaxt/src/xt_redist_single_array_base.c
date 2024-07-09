/**
 * @file xt_redist_single_array_base.c
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <mpi.h>

#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt_redist_single_array_base.h"
#include "xt_redist_internal.h"
#include "xt/xt_xmap.h"
#include "xt/xt_idxlist.h"
#include "core/ppm_xfuncs.h"
#include "core/core.h"
#include "xt_exchanger.h"

static void
redist_sab_delete(Xt_redist redist);

static void
redist_sab_s_exchange(Xt_redist redist, int num_arrays,
                      const void **src_data, void **dst_data);

static void
redist_sab_s_exchange1(Xt_redist redist, const void *src_data, void *dst_data);

static MPI_Datatype
redist_sab_get_send_MPI_Datatype(Xt_redist redist, int rank);

static MPI_Datatype
redist_sab_get_recv_MPI_Datatype(Xt_redist redist, int rank);

static MPI_Comm
redist_sab_get_MPI_Comm(Xt_redist redist);

static const struct xt_redist_vtable redist_sab_vtable = {
  .delete                = redist_sab_delete,
  .s_exchange            = redist_sab_s_exchange,
  .s_exchange1           = redist_sab_s_exchange1,
  .get_send_MPI_Datatype = redist_sab_get_send_MPI_Datatype,
  .get_recv_MPI_Datatype = redist_sab_get_recv_MPI_Datatype,
  .get_MPI_Comm          = redist_sab_get_MPI_Comm
};

typedef struct Xt_redist_sab_ *Xt_redist_sab;

struct Xt_redist_sab_ {

  const struct xt_redist_vtable *vtable;

  Xt_exchanger exchanger;

  int nsend, nrecv;

  struct Xt_redist_msg * send_msgs;
  struct Xt_redist_msg * recv_msgs;

  MPI_Comm comm;
  int tag_offset;
};

Xt_redist xt_redist_single_array_base_new(int nsend, int nrecv,
                                          struct Xt_redist_msg * send_msgs,
                                          struct Xt_redist_msg * recv_msgs,
                                          MPI_Comm comm) {

  Xt_redist_sab redist = xmalloc(sizeof (*redist));

  redist->comm = xt_mpi_comm_smart_dup(comm, &redist->tag_offset);
  redist->exchanger
    = xt_exchanger_default_constructor(nsend, nrecv, send_msgs,
                                       recv_msgs, redist->comm,
                                       redist->tag_offset);

  redist->vtable = &redist_sab_vtable;
  redist->nsend = nsend;
  redist->nrecv = nrecv;
  redist->send_msgs = send_msgs;
  redist->recv_msgs = recv_msgs;

  return (Xt_redist)redist;
}

static inline Xt_redist_sab
xrsab(void *redist)
{
  return (Xt_redist_sab)redist;
}

static void
redist_sab_delete(Xt_redist redist) {

  Xt_redist_sab redist_sab = xrsab(redist);

  xt_exchanger_delete(redist_sab->exchanger);

  xt_redist_msgs_free((size_t)redist_sab->nsend, redist_sab->send_msgs,
                      redist_sab->comm);
  xt_redist_msgs_free((size_t)redist_sab->nrecv, redist_sab->recv_msgs,
                      redist_sab->comm);

  xt_mpi_comm_smart_dedup(&redist_sab->comm, redist_sab->tag_offset);

  free(redist_sab);
}

static void
redist_sab_s_exchange(Xt_redist redist, int num_arrays,
                      const void **src_data, void **dst_data)
{
  Xt_redist_sab redist_rep = xrsab(redist);
  if (num_arrays == 1)
    redist_sab_s_exchange1(redist, src_data[0], dst_data[0]);
  else
    Xt_abort(redist_rep->comm, "ERROR: multi-array s_exchange is not"
             " implemented for this xt_redist type "
             "(Xt_redist_single_array_base)", __FILE__, __LINE__);
}

static void
redist_sab_s_exchange1(Xt_redist redist, const void *src_data, void *dst_data) {

  Xt_redist_sab redist_sab = xrsab(redist);

  xt_exchanger_s_exchange(redist_sab->exchanger, src_data, dst_data);
}

static MPI_Datatype
copy_msg_dt(int nmsg, const struct Xt_redist_msg msg[nmsg],
            MPI_Comm comm, int rank)
{
  MPI_Datatype datatype_copy = MPI_DATATYPE_NULL;

  for (int i = 0; i < nmsg; ++i)
    if (msg[i].rank == rank) {
      xt_mpi_call(MPI_Type_dup(msg[i].datatype, &datatype_copy), comm);
      break;
    }

  return datatype_copy;
}

static MPI_Datatype
redist_sab_get_send_MPI_Datatype(Xt_redist redist, int rank) {

  Xt_redist_sab redist_sab = xrsab(redist);
  return copy_msg_dt(redist_sab->nsend, redist_sab->send_msgs,
                     redist_sab->comm, rank);
}

static MPI_Datatype
redist_sab_get_recv_MPI_Datatype(Xt_redist redist, int rank) {

  Xt_redist_sab redist_sab = xrsab(redist);
  return copy_msg_dt(redist_sab->nrecv, redist_sab->recv_msgs,
                     redist_sab->comm, rank);
}

static MPI_Comm
redist_sab_get_MPI_Comm(Xt_redist redist) {

  Xt_redist_sab redist_sab = xrsab(redist);

  return redist_sab->comm;
}
