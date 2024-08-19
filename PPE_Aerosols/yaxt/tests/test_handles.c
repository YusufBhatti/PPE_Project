/**
 * @file test_handles.c
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

#include "tests.h"

int main(void) {
  enum {
    hsize = 32,
  };
  int data[hsize];
  int handle[hsize];
  int state[hsize];

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);
  xt_initialize(MPI_COMM_WORLD);

  for (int i = 0; i<hsize; i++) {
    state[i] = 0;
  }

  Xt_handle_set_type hset = xt_handle_set_new(4);

  // new handles
  for (int i = 0; i<hsize; i++) {
    if (i % 2 == 0) {
      data[i]=i;
      handle[i] = xt_handle_new(hset, &data[i]);
      if (!xt_handle_is_valid(hset,handle[i]))
        PUT_ERR("unexpected invalid handle\n");
      state[i] = 1;
    }
  }

  // delete some handles
  for (int i = 1; i<hsize; i++) {
    if (state[i]) {
      if (i % 3 == 0) {
        if (!xt_handle_is_valid(hset,handle[i]))
          PUT_ERR("unexpected invalid handle\n");
        xt_handle_delete(hset,handle[i]);
        if (xt_handle_is_valid(hset,handle[i]))
          PUT_ERR("unexpected valid handle\n");
        state[i] = 0;
      }
    }
  }

  // more new handles
  for (int i = 0; i<hsize; i++) {
    if (!state[i]) {
      data[i]=i;
      handle[i] = xt_handle_new(hset, &data[i]);
      if (!xt_handle_is_valid(hset,handle[i]))
        PUT_ERR("unexpected invalid handle\n");
      state[i] = 1;
    }
  }

  // delete some handles
  for (int i = 0; i<hsize; i++) {
    if (state[i]) {
      if (i % 5 == 0) {
        int h = handle[i];
        if (!xt_handle_is_valid(hset,h)) PUT_ERR("unexpected invalid handle\n");
        xt_handle_delete(hset,h);
        if (xt_handle_is_valid(hset,h)) PUT_ERR("unexpected valid handle\n");
        state[i] = 0;
      }
    }
  }

  // check content
  for (int i = 0; i<hsize; i++) {
    if (state[i]) {
      int h = handle[i];
      if (!xt_handle_is_valid(hset,h)) PUT_ERR("unexpected invalid handle\n");
      int *p = xt_handle2pointer(hset, h);
      if (*p != data[i]) PUT_ERR("data lookup failed\n");
    }
  }

  xt_handle_set_delete(hset);

  xt_finalize();
  MPI_Finalize();
  return TEST_EXIT_CODE;
}
