/**
 * @file test_redist_collection_static.c
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

#include <string.h>

#include <mpi.h>

#include <yaxt.h>

#include "tests.h"

static int test_communicator(MPI_Comm comm1, MPI_Comm comm2);

int main(void) {

  // init mpi

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  { // general test with one redist
    // set up data
    Xt_int src_index_list[] = {1,2,3,4,5};
    int src_num_indices
      = sizeof(src_index_list) / sizeof(src_index_list[0]);
    Xt_idxlist src_idxlist = xt_idxvec_new(src_index_list, src_num_indices);

    Xt_int dst_index_list[] = {1,3,5};
    int dst_num_indices
      = sizeof(dst_index_list) / sizeof(dst_index_list[0]);
    Xt_idxlist dst_idxlist = xt_idxvec_new(dst_index_list, dst_num_indices);

    Xt_xmap xmap;
    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_idxlist);

    Xt_redist redist;
    redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

    xt_xmap_delete(xmap);

    // generate redist_collection

    Xt_redist redist_coll;
    MPI_Aint src_displacements[1] = {0};
    MPI_Aint dst_displacements[1] = {0};

    redist_coll = xt_redist_collection_static_new(&redist, 1, src_displacements,
                                                  dst_displacements,
                                                  MPI_COMM_WORLD);

    // test communicator of redist

    if (!test_communicator(xt_redist_get_MPI_Comm(redist_coll), MPI_COMM_WORLD))
      PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    xt_redist_delete(redist);

    // test exchange

    double src_data[] = {1,2,3,4,5};
    double dst_data[] = {-1,-1,-1};

    xt_redist_s_exchange1(redist_coll, (void*)src_data, (void*)dst_data);

    double ref_dst_data[] = {1,3,5};

    for (size_t i = 0; i < sizeof(dst_data) / sizeof(dst_data[0]); ++i)
      if (ref_dst_data[i] != dst_data[i])
        PUT_ERR("error in xt_redist_s_exchange\n");

    // clean up

    xt_redist_delete(redist_coll);
  }

  { // test with one redist used three times (two exchanges)
    // set up data
    Xt_int src_index_list[] = {1,2,3,4,5};
    int src_num_indices = sizeof(src_index_list) / sizeof(src_index_list[0]);
    Xt_idxlist src_idxlist;
    src_idxlist = xt_idxvec_new(src_index_list, src_num_indices);

    Xt_int dst_index_list[] = {1,3,5};
    int dst_num_indices = sizeof(dst_index_list) / sizeof(dst_index_list[0]);
    Xt_idxlist dst_idxlist = xt_idxvec_new(dst_index_list, dst_num_indices);

    Xt_xmap xmap;
    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_idxlist);

    Xt_redist redist;
    redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

    xt_xmap_delete(xmap);

    // generate redist_collection

    Xt_redist redist_coll;
    Xt_redist redists[3] = {redist, redist, redist};
    static const double src_data[3][5]
      = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
    double dst_data[3][3];
    MPI_Aint src_displacements[3]
      = {0, (MPI_Aint)((size_t)(src_data[1] - src_data[0]) * sizeof (double)),
         (MPI_Aint)((size_t)(src_data[2] - src_data[0]) * sizeof(double)) };
    MPI_Aint dst_displacements[3]
      = {0, (MPI_Aint)((size_t)(dst_data[1] - dst_data[0]) * sizeof(double)),
         (MPI_Aint)((size_t)(dst_data[2] - dst_data[0]) * sizeof(double)) };

    redist_coll = xt_redist_collection_static_new(redists, 3, src_displacements,
                                                  dst_displacements,
                                                  MPI_COMM_WORLD);

    // test communicator of redist

    if (!test_communicator(xt_redist_get_MPI_Comm(redist_coll), MPI_COMM_WORLD))
      PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    xt_redist_delete(redist);

    // test exchange
    {
      for (size_t j = 0; j < 3; ++j)
        for (size_t i = 0; i < 3; ++i)
          dst_data[j][i] = -1.0;

      xt_redist_s_exchange1(redist_coll, (const void*)src_data, (void*)dst_data);

      static const double ref_dst_data[3][3] = {{1,3,5},{6,8,10},{11,13,15}};

      int i, j;

      for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
          if (ref_dst_data[i][j] != dst_data[i][j])
            PUT_ERR("error in xt_redist_s_exchange\n");
    }

    {
      for (size_t j = 0; j < 3; ++j)
        for (size_t i = 0; i < 3; ++i)
          dst_data[j][i] = -1.0;

      xt_redist_s_exchange1(redist_coll, (const void*)src_data, (void*)dst_data);

      static const double ref_dst_data[3][3] = {{1,3,5},{6,8,10},{11,13,15}};

      for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
          if (ref_dst_data[i][j] != dst_data[i][j])
            PUT_ERR("error in xt_redist_s_exchange\n");
    }

    // clean up

    xt_redist_delete(redist_coll);
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
