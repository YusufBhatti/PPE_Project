/**
 * @file test_redist_collection.c
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
#include "core/ppm_xfuncs.h"

static void
test_repeated_redist(int cache_size);

static int
test_communicator(MPI_Comm comm1, MPI_Comm comm2);

static Xt_xmap
build_odd_selection_xmap(int src_num_indices);


int main(void) {

  // init mpi

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  { // general test with one redist
    // set up data
    Xt_xmap xmap = build_odd_selection_xmap(5);

    Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

    xt_xmap_delete(xmap);

    // generate redist_collection

    Xt_redist redist_coll
      = xt_redist_collection_new(&redist, 1, -1, MPI_COMM_WORLD);

    // test communicator of redist

    if (!test_communicator(xt_redist_get_MPI_Comm(redist_coll), MPI_COMM_WORLD))
      PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    xt_redist_delete(redist);

    // test exchange

    double src_data[] = {1,2,3,4,5};
    double dst_data[] = {-1,-1,-1};

    double * src_data_p = &src_data[0];
    double * dst_data_p = &dst_data[0];

    xt_redist_s_exchange1(redist_coll, src_data_p, dst_data_p);

    double ref_dst_data[] = {1,3,5};

    for (size_t i = 0; i < sizeof(dst_data) / sizeof(dst_data[0]); ++i)
      if (ref_dst_data[i] != dst_data[i])
        PUT_ERR("error in xt_redist_s_exchange\n");

    // clean up

    xt_redist_delete(redist_coll);
  }

  // test with one redist used three times (with two different input data
  // displacements -> test of cache) (with default cache size)
  // set up data
  test_repeated_redist(-1);

  // test with one redist used three times (with two different input data
  // displacements -> test of cache) (with cache size == 0)
  // set up data
  test_repeated_redist(0);

  { // test with one redist used three times (with different input
    // data displacements until the cache is full)
    // set up data
    Xt_xmap xmap = build_odd_selection_xmap(5);

    Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

    xt_xmap_delete(xmap);

    // generate redist_collection

    Xt_redist redist_coll;
    Xt_redist redists[3] = {redist, redist, redist};

    redist_coll = xt_redist_collection_new(redists, 3, -1, MPI_COMM_WORLD);

    // test communicator of redist

    if (!test_communicator(xt_redist_get_MPI_Comm(redist_coll), MPI_COMM_WORLD))
      PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    xt_redist_delete(redist);

    double src_data[3][5] = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
    double dst_data[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

    const size_t cache_size = 16;

    double * src_data_, * dst_data_;

    src_data_ = malloc((5 + cache_size + 2) * sizeof(*src_data_));
    dst_data_ = malloc((3 + cache_size + 2) * sizeof(*dst_data_));

    const void *src_data_p[3] = {src_data[0],src_data[1],NULL};
    void *dst_data_p[3] = {dst_data[0],dst_data[1],NULL};

    // test exchange
    for (size_t k = 0; k < cache_size + 2; ++k) {

      memcpy(src_data_+k, src_data[2], 5 * sizeof(*src_data_));
      memcpy(dst_data_+k, dst_data[2], 3 * sizeof(*dst_data_));

      src_data_p[2] = src_data_+k;
      dst_data_p[2] = dst_data_+k;

      xt_redist_s_exchange(redist_coll, 3, src_data_p, dst_data_p);

      double ref_dst_data[3][3] = {{1,3,5},{6,8,10},{11,13,15}};

      for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
          if (ref_dst_data[i][j] != ((double *)dst_data_p[i])[j])
            PUT_ERR("error in xt_redist_s_exchange\n");
    }

    // clean up

    free(src_data_);
    free(dst_data_);
    xt_redist_delete(redist_coll);
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

enum {
  selection_stride = 2,
};

/*
 * build xmap for destination list containing all odd elements of
 * source list dimensioned 1 to src_slice_len
 */
static Xt_xmap
build_odd_selection_xmap(int src_num_indices)
{
  if (src_num_indices < 0)
    PUT_ERR("error: src_num_indices < 0");
  Xt_int *index_list = xmalloc((size_t)src_num_indices
                               * sizeof (index_list[0]));
  for (Xt_int i = 0; i < src_num_indices; ++i)
    index_list[i] = (Xt_int)(i + 1);
  Xt_idxlist src_idxlist = xt_idxvec_new(index_list, src_num_indices);
  int dst_num_indices
    = (int)((src_num_indices + selection_stride - 1) / selection_stride);
  for (int i = 0; i < dst_num_indices; ++i)
    index_list[i] = (Xt_int)(i * selection_stride + 1);
  Xt_idxlist dst_idxlist = xt_idxvec_new(index_list, dst_num_indices);
  free(index_list);
  Xt_xmap xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);
  xt_idxlist_delete(src_idxlist);
  xt_idxlist_delete(dst_idxlist);
  return xmap;
}


static void
test_repeated_redist(int cache_size)
{
  Xt_xmap xmap = build_odd_selection_xmap(5);

  Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

  xt_xmap_delete(xmap);

  // generate redist_collection

  Xt_redist redist_coll;
  Xt_redist redists[3] = {redist, redist, redist};

  redist_coll = xt_redist_collection_new(redists, 3, cache_size, MPI_COMM_WORLD);

  // test communicator of redist

  if (!test_communicator(xt_redist_get_MPI_Comm(redist_coll), MPI_COMM_WORLD))
    PUT_ERR("error in xt_redist_get_MPI_Comm\n");

  xt_redist_delete(redist);

  // test exchange
  {
    double src_data[3][5] = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
    double dst_data[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

    const void *src_data_p[3] = {src_data[0],src_data[1],src_data[2]};
    void *dst_data_p[3] = {dst_data[0],dst_data[1],dst_data[2]};

    xt_redist_s_exchange(redist_coll, 3, src_data_p, dst_data_p);

    double ref_dst_data[3][3] = {{1,3,5},{6,8,10},{11,13,15}};

    int i, j;

    for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j)
        if (ref_dst_data[i][j] != dst_data[i][j])
          PUT_ERR("error in xt_redist_s_exchange\n");
  }

  // test exchange with changed displacements
  {
    static const double src_data[3][5]
      = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
    double dst_data[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

    const void *src_data_p[3] = {src_data[1],src_data[0],src_data[2]};
    void *dst_data_p[3] = {dst_data[1],dst_data[0],dst_data[2]};

    xt_redist_s_exchange(redist_coll, 3, src_data_p, dst_data_p);

    double ref_dst_data[3][3] = {{1,3,5},{6,8,10},{11,13,15}};

    int i, j;

    for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j)
        if (ref_dst_data[i][j] != dst_data[i][j])
          PUT_ERR("error in xt_redist_s_exchange\n");
  }

  // test exchange with original displacements
  {
    static const double src_data[3][5]
      = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
    double dst_data[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

    const void *src_data_p[3] = {src_data[0],src_data[1],src_data[2]};
    void *dst_data_p[3] = {dst_data[0],dst_data[1],dst_data[2]};

    xt_redist_s_exchange(redist_coll, 3, src_data_p, dst_data_p);

    double ref_dst_data[3][3] = {{1,3,5},{6,8,10},{11,13,15}};

    int i, j;

    for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j)
        if (ref_dst_data[i][j] != dst_data[i][j])
          PUT_ERR("error in xt_redist_s_exchange\n");
  }

  // clean up

  xt_redist_delete(redist_coll);
}

static int test_communicator(MPI_Comm comm1, MPI_Comm comm2) {

  int result;

  xt_mpi_call(MPI_Comm_compare(comm1, comm2, &result), MPI_COMM_WORLD);

  return ((result == MPI_IDENT) || (result == MPI_CONGRUENT));
}

