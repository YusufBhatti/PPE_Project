/**
 * @file test_redist_collection_parallel.c
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
#include <yaxt.h>
#include "core/ppm_xfuncs.h"

#include "tests.h"
#include "test_idxlist_utils.h"

static int test_communicator(MPI_Comm comm1, MPI_Comm comm2);

static void
check_4redist_result(int size, void *results[4],
                     const Xt_int *index_vector_a,
                     const Xt_int *index_vector_b);

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
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 5)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wtype-limits"
#endif
        assert((long long)size * size < XT_INT_MAX
               && (long long)size * size < INT_MAX);
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 5)
#pragma GCC diagnostic pop
#endif
        Xt_int global_size[2]
          = { (Xt_int)((Xt_int)2 * size), (Xt_int)((Xt_int)size * size) };
        int local_size[2] = { size, size };
        Xt_int local_start[2][2]
          = {{0, (Xt_int)((Xt_int)rank * size)},
             { (Xt_int)size,
               (Xt_int)((Xt_int)size * size - (Xt_int)(rank+1) * size) } };

        indices_a_[0] = xt_idxsection_new(start, 2, global_size, local_size,
                                            local_start[0]);
        indices_a_[1] = xt_idxsection_new(start, 2, global_size, local_size,
                                            local_start[1]);

        indices_a = xt_idxlist_collection_new(indices_a_, 2);

        xt_idxlist_delete(indices_a_[0]);
        xt_idxlist_delete(indices_a_[1]);
      }

      {
        struct Xt_stripe stripe = {.start = (Xt_int)((Xt_int)rank*2*size*size),
                                   .stride = 1, .nstrides = 2*size*size};

        indices_b = xt_idxstripes_new(&stripe, 1);
      }

      {
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 5)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wtype-limits"
#endif
        assert((long long)2 * size * size * size <= XT_INT_MAX
               && (long long)2 * size * size * size <= INT_MAX);
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 5)
#pragma GCC diagnostic pop
#endif
        struct Xt_stripe stripe = { .start = 0, .stride = 1,
                                    .nstrides = 2*size*size*size };

        indices_all = xt_idxstripes_new(&stripe, 1);
      }

      enum
      {
        list_a = 0,
        list_b = 1,
        list_all = 2,
      };
      const int list_sizes[3]
        = { 2*size*size, 2*size*size, 2 * size * size * size };

      Xt_int *index_vector[2];
      for (size_t i = 0; i < 2; ++i)
        index_vector[i] = xmalloc((size_t)list_sizes[i] * sizeof (Xt_int));

      xt_idxlist_get_indices(indices_a, index_vector[list_a]);
      xt_idxlist_get_indices(indices_b, index_vector[list_b]);

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

      for (size_t i = 0; i < 4; ++i)
        xt_xmap_delete(xmaps[i]);

      Xt_redist redist = xt_redist_collection_new(redists, 4, -1,
                                                  MPI_COMM_WORLD);

      // test communicator of redist

      if (!test_communicator(xt_redist_get_MPI_Comm(redist), MPI_COMM_WORLD))
        PUT_ERR("error in xt_redist_get_MPI_Comm\n");

      for (size_t i = 0; i < 4; ++i)
        xt_redist_delete(redists[i]);

      enum { num_results = 4 };

      const size_t result_sizes[num_results] =
        { 2 * (size_t)size * (size_t)size,
          2 * (size_t)size * (size_t)size,
          2 * (size_t)size * (size_t)size * (size_t)size,
          2 * (size_t)size * (size_t)size * (size_t)size };

      void *results[num_results];
      for (size_t i = 0; i < num_results; ++i)
        results[i] = xmalloc(result_sizes[i] * sizeof (Xt_int));

      const void *input[num_results]
        = { index_vector[list_a], index_vector[list_b],
            index_vector[list_a], index_vector[list_b] };

      xt_redist_s_exchange(redist, num_results, input, results);

      check_4redist_result(size, results, index_vector[list_a],
                           index_vector[list_b]);

      /*
       * create another first buffer, to test successful
       * adaptation to different addresses...
       */
      if (rank == 0)
      {
        void *temp = results[0];
        results[0] = xmalloc(result_sizes[0] * sizeof(Xt_int));
        free(temp);
      }
      /* ...and repeat exchange */
      xt_redist_s_exchange(redist, num_results, input, results);

      check_4redist_result(size, results, index_vector[list_a],
                           index_vector[list_b]);

      // clean up
      for (size_t i = 0; i < num_results; ++i)
        free(results[i]);
      for (size_t i = 0; i < 2; ++i)
        free(index_vector[i]);
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

      Xt_redist redist = xt_redist_collection_new(redists, 2, -1,
                                                  MPI_COMM_WORLD);

      // test communicator of redist

      if (!test_communicator(xt_redist_get_MPI_Comm(redist), MPI_COMM_WORLD))
        PUT_ERR("error in xt_redist_get_MPI_Comm\n");

      xt_redist_delete(redists[0]), xt_redist_delete(redists[1]);

      Xt_int results_1[5] = {-1,-1,-1,-1,-1}, results_2[5] = {-1,-1,-1,-1,-1};

      void *results[2] = {results_1, results_2};

      const void *input[2] = {src_indices_, src_indices_};

      xt_redist_s_exchange(redist, 2, input, results);

      // check results

      for (int i = 0; i < 5; ++i) {
        if (((Xt_int *)results[0])[i] != dst_indices_[0][i])
          PUT_ERR("error on xt_redist_s_exchange\n");
        if (((Xt_int *)results[1])[i] != dst_indices_[1][i])
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

static void
check_4redist_result(int size, void *results[4],
                     const Xt_int *index_vector_a,
                     const Xt_int *index_vector_b)
{
  if (cmp_idx_arrays(2 * (size_t)size * (size_t)size,
                     (Xt_int *)results[0], index_vector_b))
    PUT_ERR("error on xt_redist_s_exchange\n");

  if (cmp_idx_arrays(2 * (size_t)size * (size_t)size,
                     (Xt_int *)results[1], index_vector_a))
    PUT_ERR("error on xt_redist_s_exchange\n");

  for (int i = 0; i < 2*size*size*size; ++i)
    if (((Xt_int *)results[2])[i] != i)
      PUT_ERR("error on xt_redist_s_exchange\n");

  for (int i = 0; i < 2*size*size*size; ++i)
    if (((Xt_int *)results[3])[i] != i)
      PUT_ERR("error on xt_redist_s_exchange\n");
}
