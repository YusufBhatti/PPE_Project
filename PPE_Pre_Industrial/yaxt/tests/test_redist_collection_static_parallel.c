/**
 * @file test_redist_collection_static_parallel.c
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
#include <stdlib.h>

#include <mpi.h>

#include <yaxt.h>

#include "tests.h"

static int test_communicator(MPI_Comm comm1, MPI_Comm comm2);

int main(void) {

  // init mpi

  int rank, size;

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  xt_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &rank), MPI_COMM_WORLD);
  xt_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &size), MPI_COMM_WORLD);

  if (size > 1) {

    { // redist test with four different redists

      Xt_idxlist indices_a, indices_b, indices_all;

      {
        Xt_idxlist indices_a_[2];

        Xt_int start = 0;
        assert(size <= XT_INT_MAX / size);
        Xt_int global_size[2] = {(Xt_int)(2*size), (Xt_int)(size*size)};
        int local_size[2] = {size,size};
        Xt_int local_start[2][2]
          = {{0, (Xt_int)(rank*size)},
             {(Xt_int)size, (Xt_int)(size*size-(rank+1)*size)}};

        indices_a_[0] = xt_idxsection_new(start, 2, global_size, local_size,
                                          local_start[0]);
        indices_a_[1] = xt_idxsection_new(start, 2, global_size, local_size,
                                          local_start[1]);

        indices_a = xt_idxlist_collection_new(indices_a_, 2);

        xt_idxlist_delete(indices_a_[0]);
        xt_idxlist_delete(indices_a_[1]);
      }

      {
        assert(size - 1 <= INT_MAX / 2 / size / size);
        struct Xt_stripe stripe = {.start = (Xt_int)(rank*2*size*size),
                                   .stride = 1, .nstrides = 2*size*size};

        indices_b = xt_idxstripes_new(&stripe, 1);
      }

      {
        assert(size <= INT_MAX / 2 / size / size);
        struct Xt_stripe stripe = {.start = 0, .stride = 1,
                                   .nstrides = 2*size*size*size};

        indices_all = xt_idxstripes_new(&stripe, 1);
      }

      Xt_int index_vector_a[2*size*size], index_vector_b[2*size*size],
             index_vector_all[2*size*size*size];

      xt_idxlist_get_indices(indices_a, index_vector_a);
      xt_idxlist_get_indices(indices_b, index_vector_b);
      xt_idxlist_get_indices(indices_all, index_vector_all);

      Xt_xmap xmaps[4] = {xt_xmap_all2all_new(indices_a, indices_b,
                                              MPI_COMM_WORLD),
                          xt_xmap_all2all_new(indices_b, indices_a,
                                              MPI_COMM_WORLD),
                          xt_xmap_all2all_new(indices_a, indices_all,
                                              MPI_COMM_WORLD),
                          xt_xmap_all2all_new(indices_b, indices_all,
                                              MPI_COMM_WORLD)};

      xt_idxlist_delete(indices_a);
      xt_idxlist_delete(indices_b);
      xt_idxlist_delete(indices_all);

      Xt_redist redists[4] = {xt_redist_p2p_new(xmaps[0], Xt_int_dt),
                              xt_redist_p2p_new(xmaps[1], Xt_int_dt),
                              xt_redist_p2p_new(xmaps[2], Xt_int_dt),
                              xt_redist_p2p_new(xmaps[3], Xt_int_dt)};

      xt_xmap_delete(xmaps[0]), xt_xmap_delete(xmaps[1]);
      xt_xmap_delete(xmaps[2]), xt_xmap_delete(xmaps[3]);

      Xt_int results_1[2*size*size], results_2[2*size*size];
      Xt_int results_3[2*size*size*size], results_4[2*size*size*size];

      MPI_Aint src_displacements[4]
        = {0, (MPI_Aint)((size_t)(index_vector_b-index_vector_a)
                         * sizeof(Xt_int)),
           0, (MPI_Aint)((size_t)(index_vector_b-index_vector_a)
                         * sizeof(Xt_int))};
      MPI_Aint dst_displacements[4]
        = {0, (MPI_Aint)((size_t)(results_2-results_1)*sizeof(Xt_int)),
           (MPI_Aint)((size_t)(results_3-results_1)*sizeof(Xt_int)),
           (MPI_Aint)((size_t)(results_4-results_1)*sizeof(Xt_int))};

      Xt_redist redist
        = xt_redist_collection_static_new(redists, 4, src_displacements,
                                          dst_displacements, MPI_COMM_WORLD);

      // test communicator of redist

      if (!test_communicator(xt_redist_get_MPI_Comm(redist), MPI_COMM_WORLD))
        PUT_ERR("error in xt_redist_get_MPI_Comm\n");

      xt_redist_delete(redists[0]), xt_redist_delete(redists[1]);
      xt_redist_delete(redists[2]), xt_redist_delete(redists[3]);

      xt_redist_s_exchange1(redist, (void*)index_vector_a, (void*)results_1);

      // check results

      for (int i = 0; i < 2*size*size; ++i)
        if (results_1[i] != index_vector_b[i])
          PUT_ERR("error on xt_redist_s_exchange\n");

      for (int i = 0; i < 2*size*size; ++i)
        if (results_2[i] != index_vector_a[i])
          PUT_ERR("error on xt_redist_s_exchange\n");

      for (int i = 0; i < 2*size*size*size; ++i)
        if (results_3[i] != i)
          PUT_ERR("error on xt_redist_s_exchange\n");

      for (int i = 0; i < 2*size*size*size; ++i)
        if (results_4[i] != i)
          PUT_ERR("error on xt_redist_s_exchange\n");

      // clean up

      xt_redist_delete(redist);
    }

    { // redist test with two redists that do a round robin exchange in
      // different directions

      Xt_idxlist src_indices, dst_indices[2];

      Xt_int src_indices_[5], dst_indices_[2][5];

      for (Xt_int i = 0; i < 5; ++i) {
        src_indices_[i] = (Xt_int)(rank * 5 + i);
        dst_indices_[0][i] = (Xt_int)((src_indices_[i] + 1) % (size * 5));
        Xt_int temp = (Xt_int)(src_indices_[i] - 1);
        dst_indices_[1][i] = (Xt_int)((temp < 0)?(size * 5 - 1):temp);
      }

      src_indices = xt_idxvec_new(src_indices_, 5);
      dst_indices[0] = xt_idxvec_new(dst_indices_[0], 5);
      dst_indices[1] = xt_idxvec_new(dst_indices_[1], 5);

      Xt_xmap xmaps[2] = {xt_xmap_all2all_new(src_indices, dst_indices[0],
                                              MPI_COMM_WORLD),
                          xt_xmap_all2all_new(src_indices, dst_indices[1],
                                              MPI_COMM_WORLD)};

      xt_idxlist_delete(src_indices);
      xt_idxlist_delete(dst_indices[0]);
      xt_idxlist_delete(dst_indices[1]);

      Xt_redist redists[2] = {xt_redist_p2p_new(xmaps[0], Xt_int_dt),
                              xt_redist_p2p_new(xmaps[1], Xt_int_dt)};

      xt_xmap_delete(xmaps[0]), xt_xmap_delete(xmaps[1]);

      Xt_int results_1[5] = {-1,-1,-1,-1,-1}, results_2[5] = {-1,-1,-1,-1,-1};

      MPI_Aint src_displacements[2] = {0, 0};
      MPI_Aint dst_displacements[2]
        = {0, (MPI_Aint)((size_t)(results_2-results_1)*sizeof(Xt_int))};

      Xt_redist redist
        = xt_redist_collection_static_new(redists, 2, src_displacements,
                                          dst_displacements, MPI_COMM_WORLD);

      // test communicator of redist

      if (!test_communicator(xt_redist_get_MPI_Comm(redist), MPI_COMM_WORLD))
        PUT_ERR("error in xt_redist_get_MPI_Comm\n");

      xt_redist_delete(redists[0]), xt_redist_delete(redists[1]);

      xt_redist_s_exchange1(redist, (void*)src_indices_, (void*)results_1);

      // check results

      for (int i = 0; i < 5; ++i) {
        if (results_1[i] != dst_indices_[0][i])
          PUT_ERR("error on xt_redist_s_exchange\n");
        if (results_2[i] != dst_indices_[1][i])
          PUT_ERR("error on xt_redist_s_exchange\n");
      }

      // clean up

      xt_redist_delete(redist);
    }
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static int test_communicator(MPI_Comm comm1, MPI_Comm comm2) {

  int result;

  xt_mpi_call(MPI_Comm_compare(comm1, comm2, &result), MPI_COMM_WORLD);

  return ((result == MPI_IDENT) || (result == MPI_CONGRUENT));
}
