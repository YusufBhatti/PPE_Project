/**
 * @file xt_exchanger.c
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

#include <string.h>
#include <stdlib.h>
#if ! defined TLS && defined HAVE_PTHREAD
#include <stdio.h>
#include <pthread.h>
#include "core/core.h"
#endif

#include "xt_exchanger.h"
#include "xt_exchanger_irecv_isend.h"
#include "xt_exchanger_mix_isend_irecv.h"
#include "xt_exchanger_irecv_isend_packed.h"

Xt_exchanger (*xt_exchanger_default_constructor)
  (int nsend, int nrecv, struct Xt_redist_msg * send_msgs,
   struct Xt_redist_msg * recv_msgs, MPI_Comm comm, int tag_offset) =
    xt_exchanger_mix_isend_irecv_new;
    // xt_exchanger_irecv_isend_new;

void xt_exchanger_delete(Xt_exchanger exchanger) {

  exchanger->vtable->delete(exchanger);
}

void xt_exchanger_s_exchange(Xt_exchanger exchanger, const void *src_data, void *dst_data) {

  exchanger->vtable->s_exchange(exchanger, src_data, dst_data);
}

#if ! defined TLS && defined HAVE_PTHREAD

static pthread_key_t tls_qsort_data_key;

void xt_exchanger_init(void) {
  int ierror = pthread_key_create(&tls_qsort_data_key, NULL);
  if (ierror) {
    fprintf(stderr, "%s: error creating pthread key: %s\n",
            __func__, strerror(ierror));
    die("Failed pthread_key_create in initialization!");
  }
}

void xt_exchanger_finalize(void) {
  int ierror = pthread_key_delete(tls_qsort_data_key);
  if (ierror) {
    fprintf(stderr, "%s: error deleting pthread key: %s\n",
            __func__, strerror(ierror));
    die("Failed pthread_key_delete in finalization!");
  }
}

#else

static TLS int comm_rank, comm_size;

#endif

static int compare_exchange_messages (const void * msg_a, const void * msg_b)
{
  int rank_a = *(const int *)msg_a;
  int rank_b = *(const int *)msg_b;
#if ! defined TLS && defined HAVE_PTHREAD
  int comm_rank, comm_size;
  {
    int *qsort_data = pthread_getspecific(tls_qsort_data_key);
    comm_size = qsort_data[0];
    comm_rank = qsort_data[1];
  }
#endif
  if (rank_a <= comm_rank) rank_a += comm_size;
  if (rank_b <= comm_rank) rank_b += comm_size;

  return ( rank_a - rank_b );
}

void xt_exchanger_internal_optimize(size_t n, void * msgs, size_t msg_type_size,
                                    MPI_Comm comm) {

#if ! defined TLS && defined HAVE_PTHREAD
  int qsort_data[2];
  MPI_Comm_rank(comm, qsort_data + 1);
  MPI_Comm_size(comm, qsort_data + 0);
  pthread_setspecific(tls_qsort_data_key, qsort_data);
#else
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);
#endif
  /* In order to avoid congestion of messages, the order of send and receive
   * messages is changed. This is done by sorting the messages according to the
   * rank of the respective message partner. Before the sorting to ranks that
   * are smaller or equal to the local rank the size of the communicator is
   * added.
   * example: process 5 is supposed to communicate with processes: 9, 5, 2, 6, 1
   * 1. add comm_size(10): 9, 15, 12, 6, 11
   * 2. sort: 6, 9, 11, 12, 15 -> final order: 6, 9, 1, 2, 5
   */
  qsort(msgs, n, msg_type_size, compare_exchange_messages);
}
