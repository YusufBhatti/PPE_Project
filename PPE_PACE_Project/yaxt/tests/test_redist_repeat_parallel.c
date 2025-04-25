/**
 * @file test_redist_repeat_parallel.c
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

    { // redist test with one redist and 4 levels

      Xt_idxlist indices_a, indices_b;

      {
        Xt_idxlist indices_a_[2];

        Xt_int start = 0;
        assert(size <= XT_INT_MAX / size);
        Xt_int global_size[2] = {(Xt_int)(2*size), (Xt_int)(size*size)};
        int local_size[2] = {size,size};
        Xt_int local_start[2][2]
          = {{ 0, (Xt_int)(rank*size)},
             { (Xt_int)size, (Xt_int)(size*size-(rank+1)*size) }};

        indices_a_[0] = xt_idxsection_new(start, 2, global_size, local_size,
                                          local_start[0]);
        indices_a_[1] = xt_idxsection_new(start, 2, global_size, local_size,
                                          local_start[1]);

        indices_a = xt_idxlist_collection_new(indices_a_, 2);

        xt_idxlist_delete(indices_a_[0]);
        xt_idxlist_delete(indices_a_[1]);
      }

      {
        assert(rank <= XT_INT_MAX / 2 / size / size);
        struct Xt_stripe stripe = {.start = (Xt_int)(rank*2*size*size),
                                   .stride = 1, .nstrides = 2*size*size};

        indices_b = xt_idxstripes_new(&stripe, 1);
      }

      Xt_int index_vector_a[2*size*size], index_vector_b[2*size*size];

      xt_idxlist_get_indices(indices_a, index_vector_a);
      xt_idxlist_get_indices(indices_b, index_vector_b);

      Xt_xmap xmap = xt_xmap_all2all_new(indices_a, indices_b, MPI_COMM_WORLD);

      xt_idxlist_delete(indices_a);
      xt_idxlist_delete(indices_b);

      Xt_redist redist_p2p = xt_redist_p2p_new(xmap, Xt_int_dt);

      xt_xmap_delete(xmap);

      Xt_int input[9][2*size*size];
      Xt_int result[4][2*size*size];
      Xt_int result_2[9][2*size*size];

      for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 2*size*size; ++j)
          input[i][j] = (Xt_int)(index_vector_a[j] + i * 4*size*size);

      int displacements[4] = {0, 1, 2, 3};
      int displacements_2[4] = {1, 2, 4, 8};
      MPI_Aint extent
        = (MPI_Aint)(2*(size_t)size*(size_t)size * sizeof(Xt_int));

      Xt_redist redist_repeat
        = xt_redist_repeat_new(redist_p2p, extent, extent, 4, displacements);
      Xt_redist redist_repeat_2
        = xt_redist_repeat_new(redist_p2p, extent, extent, 4, displacements_2);

      // test communicator of redist_repeat

      if (!test_communicator(xt_redist_get_MPI_Comm(redist_repeat), MPI_COMM_WORLD))
        PUT_ERR("error in xt_redist_get_MPI_Comm\n");
      if (!test_communicator(xt_redist_get_MPI_Comm(redist_repeat_2), MPI_COMM_WORLD))
        PUT_ERR("error in xt_redist_get_MPI_Comm\n");

      xt_redist_delete(redist_p2p);

      xt_redist_s_exchange1(redist_repeat, (void*)input, (void*)result);
      xt_redist_s_exchange1(redist_repeat_2, (void*)input, (void*)result_2);

      // check results

      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 2*size*size; ++j)
          if (result[i][j] != index_vector_b[j] + i * 4*size*size)
            PUT_ERR("ERROR: in first xt_redist_s_exchange1\n");

      for (int i = 1; i <= 8; i<<=1)
        for (int j = 0; j < 2*size*size; ++j)
          if (result_2[i][j] != index_vector_b[j] + i * 4*size*size)
            PUT_ERR("ERROR: in second xt_redist_s_exchange1\n");

      // clean up

      xt_redist_delete(redist_repeat_2);
      xt_redist_delete(redist_repeat);
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
