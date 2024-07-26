/**
 * @file test_xmap_common.c
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

#include <mpi.h>

#include <yaxt.h>

#define VERBOSE
#include "tests.h"
#include "test_xmap_common.h"

static void
test_xmap1(xmap_constructor new_xmap, MPI_Comm comm);
static void
test_xmap2(xmap_constructor new_xmap, MPI_Comm comm);

static int my_rank;

int
xt_xmap_self_test_main(xmap_constructor new_xmap)
{
  // init mpi
  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  MPI_Comm comms[2];
  comms[0] = MPI_COMM_WORLD;
  MPI_Comm_dup(MPI_COMM_WORLD, &comms[1]);
  xt_mpi_comm_mark_exclusive(comms[1]);

  for (size_t i = 0; i < 2; ++i) {
    test_xmap1(new_xmap, comms[i]);
    test_xmap2(new_xmap, comms[i]);
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void
test_xmap(const Xt_int *src_index_list, int num_src_indices,
          const Xt_int *dst_index_list, int num_dst_indices,
          xmap_constructor new_xmap, MPI_Comm comm)
{
  Xt_idxlist src_idxlist = xt_idxvec_new(src_index_list, num_src_indices);
  Xt_idxlist dst_idxlist = xt_idxvec_new(dst_index_list, num_dst_indices);

  // test of exchange map
  Xt_xmap xmap = new_xmap(src_idxlist, dst_idxlist, comm);
  xt_idxlist_delete(src_idxlist);
  xt_idxlist_delete(dst_idxlist);

  // test results
  if (xt_xmap_get_num_destinations(xmap) != 1)
    PUT_ERR("error in xt_xmap_get_num_destinations\n");

  if (xt_xmap_get_num_sources(xmap) != 1)
    PUT_ERR("error in xt_xmap_get_num_sources\n");

  int rank;

  xt_xmap_get_destination_ranks(xmap, &rank);
  if (rank != my_rank)
    PUT_ERR("error in xt_xmap_get_destination_ranks\n");

  xt_xmap_get_source_ranks(xmap, &rank);
  if (rank != my_rank)
    PUT_ERR("error in xt_xmap_get_source_ranks\n");
  // clean up
  xt_xmap_delete(xmap);
}

static inline void shift_idx(Xt_int idx[], int num, int offset)  {
  for (int i=0; i<num; i++) {
    idx[i] = (Xt_int)(idx[i] + my_rank * offset);
  }
}

static void
test_xmap1(xmap_constructor new_xmap, MPI_Comm comm)
{
  // source index list
  Xt_int src_index_list[] = {1,2,3,4,5,6,7};
  int num_src_indices
    = sizeof(src_index_list) / sizeof(src_index_list[0]);
  shift_idx(src_index_list, num_src_indices, 7);

  // destination index list
  Xt_int dst_index_list[] = {7,6,5,4,3,2,1};
  int num_dst_indices
    = sizeof(dst_index_list) / sizeof(dst_index_list[0]);
  shift_idx(dst_index_list, num_dst_indices, 7);

  test_xmap(src_index_list, num_src_indices, dst_index_list,
            num_dst_indices, new_xmap, comm);
}

static void
test_xmap2(xmap_constructor new_xmap, MPI_Comm comm)
{
  // source index list
  Xt_int src_index_list[] = {5,67,4,5,13,9,2,1,0,96,13,12,1,3};
  int num_src_indices
    = sizeof(src_index_list) / sizeof(src_index_list[0]);
  shift_idx(src_index_list, num_src_indices, 100);

  // destination index list
  Xt_int dst_index_list[] = {5,4,3,96,1,5,4,5,4,3,13,2,1};
  int num_dst_indices
    = sizeof(dst_index_list) / sizeof(dst_index_list[0]);
  shift_idx(dst_index_list, num_dst_indices, 100);

  test_xmap(src_index_list, num_src_indices,
            dst_index_list, num_dst_indices,
            new_xmap, comm);
}
