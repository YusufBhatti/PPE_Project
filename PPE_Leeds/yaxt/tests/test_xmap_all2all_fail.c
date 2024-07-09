/**
 * @file test_xmap_all2all_fail.c
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

#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>

#include <yaxt.h>

#define VERBOSE
#include "tests.h"

static enum {
  SMALL,
  BIG,
} index_list_size = SMALL;

static void
parse_options(int *argc, char ***argv);

int main(int argc, char **argv) {

   // init mpi
   xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

   xt_initialize(MPI_COMM_WORLD);

   int my_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

   parse_options(&argc, &argv);

   {
      // source index list
      struct Xt_stripe src_stripe;
      src_stripe.nstrides = (index_list_size == SMALL)?7:1023;
      src_stripe.start = (Xt_int)(1 + (Xt_int)my_rank * src_stripe.nstrides);
      src_stripe.stride = 1;
      Xt_idxlist src_idxlist = xt_idxstripes_new(&src_stripe, 1);

      // destination index list
      struct Xt_stripe dst_stripe;
      dst_stripe.nstrides = src_stripe.nstrides;
      dst_stripe.start = (Xt_int)(src_stripe.start + src_stripe.nstrides);
      dst_stripe.stride = -1;
      Xt_idxlist dst_idxlist = xt_idxstripes_new(&dst_stripe, 1);

      // test of exchange map
      // NOTE: this should fail
      Xt_xmap xmap
        = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

      // test results

      if (xt_xmap_get_num_destinations(xmap) != 1)
        PUT_ERR("error in xt_xmap_get_num_destinations\n");

      if (xt_xmap_get_num_sources(xmap) != 1)
        PUT_ERR("error in xt_xmap_get_num_destinations\n");

      int rank;

      xt_xmap_get_destination_ranks(xmap, &rank);
      if (rank != my_rank)
        PUT_ERR("error in xt_xmap_get_destination_ranks\n");

      xt_xmap_get_source_ranks(xmap, &rank);
      if (rank != my_rank)
        PUT_ERR("error in xt_xmap_get_source_ranks\n");

      // clean up
      xt_xmap_delete(xmap);
      xt_idxlist_delete(src_idxlist);
      xt_idxlist_delete(dst_idxlist);
   }

   xt_finalize();
   xt_mpi_call(MPI_Finalize(), MPI_COMM_WORLD);

   return TEST_EXIT_CODE;
}

static void
parse_options(int *argc, char ***argv)
{
  int opt;
  while ((opt = getopt(*argc, *argv, "s:")) != -1) {
    switch (opt) {
    case 's':
      if (!strcmp(optarg, "small"))
        index_list_size = SMALL;
      else if (!strcmp(optarg, "big"))
        index_list_size = BIG;
      else
      {
        fprintf(stderr, "Unknown data size \"%s\"\n", optarg);
        exit(EXIT_FAILURE);
      }
    }
  }
}
