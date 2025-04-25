/**
 * @file xt_redist_collection.c
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
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt/xt_redist_collection.h"
#include "ensure_array_size.h"
#include "xt/xt_redist.h"
#include "xt_redist_internal.h"
#include "xt_exchanger.h"

#define DEFFAULT_DATATYPE_CACHE_SIZE (16)

static void
redist_collection_delete(Xt_redist redist);

static void
redist_collection_s_exchange(Xt_redist redist, int num_src_arrays,
                             const void **src_data, void **dst_data);

static void
redist_collection_s_exchange1(Xt_redist redist,
                              const void *src_data, void *dst_data);

static MPI_Datatype
redist_collection_get_send_MPI_Datatype(Xt_redist redist, int rank);

static MPI_Datatype
redist_collection_get_recv_MPI_Datatype(Xt_redist redist, int rank);

static MPI_Comm
redist_collection_get_MPI_Comm(Xt_redist redist);

static const struct xt_redist_vtable redist_collection_vtable = {
  .delete                = redist_collection_delete,
  .s_exchange            = redist_collection_s_exchange,
  .s_exchange1           = redist_collection_s_exchange1,
  .get_send_MPI_Datatype = redist_collection_get_send_MPI_Datatype,
  .get_recv_MPI_Datatype = redist_collection_get_recv_MPI_Datatype,
  .get_MPI_Comm          = redist_collection_get_MPI_Comm
};

struct redist_collection_msg {

  int rank;
  MPI_Datatype *component_dt; // datatypes of the redists (size == num_redists)
};

struct exchanger_cache
{
  size_t token;
  MPI_Aint *src_displacements, *dst_displacements;
  Xt_exchanger * exchangers;
  struct Xt_redist_msg * msgs;
};

typedef struct Xt_redist_collection_ *Xt_redist_collection;

struct Xt_redist_collection_ {

  const struct xt_redist_vtable *vtable;

  int num_redists;

  struct exchanger_cache cache;

  int ndst, nsrc;
  struct redist_collection_msg * send_msgs;
  struct redist_collection_msg * recv_msgs;

  size_t cache_size;

  MPI_Comm comm;
  int tag_offset;
};

static void copy_component_dt(struct redist_collection_msg ** msgs, int * nmsgs,
                              Xt_redist * redists, int num_redists,
                              MPI_Comm comm,
                              MPI_Datatype (*get_MPI_datatype)(Xt_redist,int))
{
  size_t msgs_array_size = 0;

  int comm_size;
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  MPI_Datatype datatypes[num_redists];
  assert(*nmsgs >= 0);
  size_t num_messages = (size_t)*nmsgs;
  struct redist_collection_msg *p = *msgs;

  for (int i = 0; i < comm_size; ++i) {

    int flag = 0;

    for (int j = 0; j < num_redists; ++j)
      flag |= ((datatypes[j] = get_MPI_datatype(redists[j], i))
               != MPI_DATATYPE_NULL);

    if (flag) {

        ENSURE_ARRAY_SIZE(p, msgs_array_size, num_messages+1);

        p[num_messages].rank = i;
        p[num_messages].component_dt = xmalloc((size_t)num_redists *
          sizeof(*(p[num_messages].component_dt)));
        memcpy(p[num_messages].component_dt, datatypes,
               (size_t)num_redists * sizeof(*datatypes));

        ++num_messages;
    }
  }

  if (num_messages > 0)
    p = xrealloc(p, num_messages * sizeof(**msgs));
  *msgs = p;
  *nmsgs = (int)num_messages;
}

/* not yet used cache entries are marked with -1 as first displacement,
 * which becomes 0 later on through use */
static inline void
init_cache(struct exchanger_cache *cache, size_t cache_size, size_t ntx,
           int num_redists)
{
  cache->exchangers = xcalloc(cache_size, sizeof(*(cache->exchangers)));
  size_t num_displ = cache_size * (size_t)num_redists;
  struct Xt_redist_msg *msgs = cache->msgs = xmalloc(ntx * sizeof (*msgs));
  for (size_t i = 0; i < ntx; ++i) msgs[i].datatype = MPI_DATATYPE_NULL;
  MPI_Aint
    *q = cache->src_displacements = xmalloc(num_displ * sizeof (*q)),
    *p = cache->dst_displacements = xmalloc(num_displ * sizeof (*p));
  for (size_t i = 0; i < num_displ; i += (size_t)num_redists)
    q[i] = p[i] = (MPI_Aint)-1;
  cache->token = 0;
}

static inline void
destruct_cache(struct exchanger_cache *cache,
               size_t cache_size, size_t ntx, MPI_Comm comm)
{
  for (size_t i = 0; i < cache_size; ++i)
    if (cache->exchangers[i] != NULL)
      xt_exchanger_delete(cache->exchangers[i]);
  free(cache->exchangers);

  xt_redist_msgs_free(ntx, cache->msgs, comm);
  free(cache->src_displacements);
  free(cache->dst_displacements);
}

static void check_comms(Xt_redist * redists, int num_redists,
                        MPI_Comm comm) {

  int result;

  for (int i = 0; i < num_redists; ++i) {

    xt_mpi_call(MPI_Comm_compare(xt_redist_get_MPI_Comm(redists[i]),
                                 comm, &result), comm);

    if ((result != MPI_IDENT) && (result != MPI_CONGRUENT))
      Xt_abort(comm, "ERROR: MPI communicators do not match; cannot build "
               "redist_collection_static\n", __FILE__, __LINE__);
  }
}

Xt_redist xt_redist_collection_new(Xt_redist * redists, int num_redists,
                                   int cache_size, MPI_Comm comm) {

  Xt_redist_collection redist_coll = xmalloc(sizeof (*redist_coll));

  redist_coll->vtable = &redist_collection_vtable;
  redist_coll->num_redists = num_redists;
  redist_coll->ndst = 0;
  redist_coll->nsrc = 0;
  redist_coll->send_msgs = NULL;
  redist_coll->recv_msgs = NULL;
  if (cache_size < -1)
    Xt_abort(comm, "ERROR: invalid cache size in xt_redist_collection_new",
             __FILE__, __LINE__);
  redist_coll->cache_size
    = (cache_size == -1)?(DEFFAULT_DATATYPE_CACHE_SIZE):(size_t)cache_size;

  redist_coll->comm = xt_mpi_comm_smart_dup(comm, &redist_coll->tag_offset);

  check_comms(redists, num_redists, comm);

  copy_component_dt(&(redist_coll->send_msgs), &(redist_coll->nsrc), redists,
                    num_redists, redist_coll->comm,
                    xt_redist_get_send_MPI_Datatype);
  copy_component_dt(&(redist_coll->recv_msgs), &(redist_coll->ndst), redists,
                    num_redists, redist_coll->comm,
                    xt_redist_get_recv_MPI_Datatype);

  init_cache(&redist_coll->cache, (size_t)redist_coll->cache_size,
             (size_t)redist_coll->nsrc + (size_t)redist_coll->ndst,
             num_redists);

  return (Xt_redist)redist_coll;
}

static MPI_Datatype
create_compound_dt(MPI_Aint *displacements, int *block_lengths,
                   struct redist_collection_msg *msg,
                   int num_redists, MPI_Comm comm)
{
  MPI_Datatype datatype;

  int num_datatypes = 0;

  for (int i = 0; i < num_redists; ++i)
    num_datatypes += (msg->component_dt[i] != MPI_DATATYPE_NULL);

  MPI_Datatype * datatypes;
  MPI_Aint * displacements_;

  if (num_datatypes != num_redists) {

    datatypes = xmalloc((size_t)num_datatypes * sizeof(*datatypes));
    displacements_ = xmalloc((size_t)num_datatypes * sizeof(*displacements));

    num_datatypes = 0;

    for (int i = 0; i < num_redists; ++i) {
      if (msg->component_dt[i] != MPI_DATATYPE_NULL) {
        datatypes[num_datatypes] = msg->component_dt[i];
        displacements_[num_datatypes] = displacements[i];
        ++num_datatypes;
      }
    }
  } else {
    datatypes = msg->component_dt;
    displacements_ = displacements;
  }

  assert(num_datatypes <= INT_MAX);
  xt_mpi_call(MPI_Type_create_struct((int)num_datatypes, block_lengths,
                                     displacements_, datatypes, &datatype),
              comm);

  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  if (num_datatypes != num_redists) {
    free(datatypes);
    free(displacements_);
  }
  return datatype;
}

static void
create_all_dt_for_dir(struct redist_collection_msg *msgs,
                      int num_messages, int num_redists,
                      MPI_Aint displacements[num_redists],
                      struct Xt_redist_msg redist_msgs[num_messages],
                      MPI_Comm comm)
{
  int block_lengths[num_redists];

  for (int i = 0; i < num_redists; ++i)
    block_lengths[i] = 1;
  for (int i = 0; i < num_messages; ++i)
  {
    if (redist_msgs[i].datatype != MPI_DATATYPE_NULL)
      xt_mpi_call(MPI_Type_free(&(redist_msgs[i].datatype)), comm);
    redist_msgs[i].datatype = create_compound_dt(displacements, block_lengths,
                                                 msgs+i, num_redists, comm);
    redist_msgs[i].rank = msgs[i].rank;
  }
}

static void
compute_displ(const void **data, int num_redists,
              MPI_Aint displacements[num_redists],
              MPI_Comm comm)
{
  if (num_redists)
  {
    MPI_Aint base_addr, offset;
    xt_mpi_call(MPI_Get_address((void *)data[0], &base_addr), comm);
    displacements[0] = 0;
    for (int i = 1; i < num_redists; ++i)
    {
      xt_mpi_call(MPI_Get_address((void *)data[i], &offset), comm);
      displacements[i] = offset - base_addr;
    }
  }
}

static size_t
lookup_cache_index(int num_redists,
                   MPI_Aint src_displacements[num_redists],
                   MPI_Aint dst_displacements[num_redists],
                   MPI_Aint (*cached_src_displacements)[num_redists],
                   MPI_Aint (*cached_dst_displacements)[num_redists],
                   size_t cache_size)
{
  for (size_t i = 0; i < cache_size &&
       cached_src_displacements[i][0] == (MPI_Aint)0 &&
       cached_dst_displacements[i][0] == (MPI_Aint)0; ++i) {
    int mismatch = 0;
    for (int j = 0; j < num_redists; ++j)
      mismatch |= (src_displacements[j] != cached_src_displacements[i][j]) ||
                  (dst_displacements[j] != cached_dst_displacements[i][j]);
    if (!mismatch) return i;
  }
  return cache_size;
}

static Xt_exchanger
get_exchanger(const void ** src_data, void ** dst_data,
              struct redist_collection_msg * send_msgs, int num_send_messages,
              struct redist_collection_msg * recv_msgs, int num_recv_messages,
              int num_redists, struct exchanger_cache *cache, size_t cache_size,
              MPI_Comm comm, int tag_offset)
{
  MPI_Aint displacements[2][num_redists];
  compute_displ(src_data, num_redists, displacements[0], comm);
  compute_displ((const void **)dst_data, num_redists, displacements[1], comm);

  Xt_exchanger exchanger;

  if (cache_size > 0)
  {
    size_t cache_index
      = lookup_cache_index(num_redists, displacements[0], displacements[1],
                           (MPI_Aint (*)[num_redists])cache->src_displacements,
                           (MPI_Aint (*)[num_redists])cache->dst_displacements,
                           cache_size);

    if (cache_index == cache_size)
    {
      cache_index = cache->token;
      create_all_dt_for_dir(send_msgs, num_send_messages, num_redists,
                            displacements[0], cache->msgs, comm);
      create_all_dt_for_dir(recv_msgs, num_recv_messages, num_redists,
                            displacements[1], cache->msgs +
                            (size_t)num_send_messages, comm);
      memcpy(cache->src_displacements + cache_index * (size_t)num_redists,
             displacements[0], sizeof (displacements[0]));
      memcpy(cache->dst_displacements + cache_index * (size_t)num_redists,
             displacements[1], sizeof (displacements[1]));

      if (cache->exchangers[cache_index] != NULL)
        xt_exchanger_delete(cache->exchangers[cache_index]);

      exchanger = cache->exchangers[cache_index] =
        xt_exchanger_default_constructor(num_send_messages, num_recv_messages,
                                         cache->msgs, cache->msgs +
                                         (size_t)num_send_messages,
                                         comm, tag_offset);

      cache->token = (cache->token + 1) % cache_size;
    }
    else
      exchanger = cache->exchangers[cache_index];
  }
  else
  {
    size_t nmsg = (size_t)num_send_messages + (size_t)num_recv_messages;
    struct Xt_redist_msg *p = xmalloc(nmsg * sizeof (*p));
    for (size_t i = 0; i < nmsg; ++i)
      p[i].datatype = MPI_DATATYPE_NULL;

    create_all_dt_for_dir(send_msgs, num_send_messages, num_redists,
                          displacements[0], p, comm);
    create_all_dt_for_dir(recv_msgs, num_recv_messages, num_redists,
                          displacements[1], p + num_send_messages, comm);

    exchanger =
      xt_exchanger_default_constructor(num_send_messages, num_recv_messages,
                                       p, p + (size_t)num_send_messages,
                                       comm, tag_offset);

    xt_redist_msgs_free(nmsg, p, comm);
  }

  return exchanger;
}

static inline Xt_redist_collection
xrc(void *redist)
{
  return (Xt_redist_collection)redist;
}

static void
redist_collection_s_exchange(Xt_redist redist, int num_arrays,
                             const void **src_data, void **dst_data) {

  Xt_redist_collection redist_coll = xrc(redist);

  if (num_arrays != redist_coll->num_redists)
    Xt_abort(redist_coll->comm, "ERROR: wrong number of arrays in "
             "redist_collection_s_exchange", __FILE__, __LINE__);


  Xt_exchanger exchanger = get_exchanger(src_data, dst_data,
                                         redist_coll->send_msgs,
                                         redist_coll->nsrc,
                                         redist_coll->recv_msgs,
                                         redist_coll->ndst,
                                         redist_coll->num_redists,
                                         &(redist_coll->cache),
                                         redist_coll->cache_size,
                                         redist_coll->comm,
                                         redist_coll->tag_offset);

  xt_exchanger_s_exchange(exchanger, src_data[0], dst_data[0]);

  if (redist_coll->cache_size == 0)
    xt_exchanger_delete(exchanger);
}

static void
free_redist_collection_msgs(struct redist_collection_msg * msgs,
                            int nmsgs, int num_redists,
                            MPI_Comm comm) {

  for (int i = 0; i < nmsgs; ++i) {

    for (int j = 0; j < num_redists; ++j)
      if (msgs[i].component_dt[j] != MPI_DATATYPE_NULL)
        xt_mpi_call(MPI_Type_free(msgs[i].component_dt+j), comm);
    free(msgs[i].component_dt);
  }
  free(msgs);
}

static void
redist_collection_delete(Xt_redist redist) {

  Xt_redist_collection redist_coll = xrc(redist);

  free_redist_collection_msgs(redist_coll->send_msgs, redist_coll->nsrc,
                              redist_coll->num_redists,
                              redist_coll->comm);

  free_redist_collection_msgs(redist_coll->recv_msgs, redist_coll->ndst,
                              redist_coll->num_redists,
                              redist_coll->comm);

  destruct_cache(&redist_coll->cache, redist_coll->cache_size,
                 (size_t)redist_coll->nsrc + (size_t)redist_coll->ndst,
                 redist_coll->comm);

  xt_mpi_comm_smart_dedup(&(redist_coll->comm), redist_coll->tag_offset);

  free(redist_coll);
}

static MPI_Datatype
redist_collection_get_send_MPI_Datatype(Xt_redist redist, int XT_UNUSED(rank))
{
  Xt_redist_collection redist_coll = xrc(redist);

  Xt_abort(redist_coll->comm, "ERROR: get_send_MPI_Datatype is not"
           " supported for this xt_redist type (Xt_redist_collection)",
           __FILE__, __LINE__);

  return MPI_DATATYPE_NULL;
}

static MPI_Datatype
redist_collection_get_recv_MPI_Datatype(Xt_redist redist, int XT_UNUSED(rank)) {

  Xt_redist_collection redist_coll = xrc(redist);

  Xt_abort(redist_coll->comm, "ERROR: get_recv_MPI_Datatype is not"
           " supported for this xt_redist type (Xt_redist_collection)",
           __FILE__, __LINE__);

  return MPI_DATATYPE_NULL;
}

static void
redist_collection_s_exchange1(Xt_redist redist,
                              const void *src_data, void *dst_data)
{

  Xt_redist_collection redist_coll = xrc(redist);
  if (redist_coll->num_redists == 1)
    redist_collection_s_exchange(redist, 1, &src_data, &dst_data);
  else
    Xt_abort(redist_coll->comm, "ERROR: s_exchange1 is not implemented for"
             " this xt_redist type (Xt_redist_collection)", __FILE__, __LINE__);
}

static MPI_Comm
redist_collection_get_MPI_Comm(Xt_redist redist) {

  Xt_redist_collection redist_coll = xrc(redist);

  return redist_coll->comm;
}
