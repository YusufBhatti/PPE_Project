/**
 * @file test_idxmod.c
 *
 * @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
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

#include <stdlib.h>
#include <assert.h>

#include <mpi.h>

#include <yaxt.h>

#include "tests.h"
#include "test_idxlist_utils.h"
#include "core/ppm_xfuncs.h"

int main(void) {

  // init:
  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);
  xt_initialize(MPI_COMM_WORLD);

  { // idxvec modifier

    Xt_int g_src_idx[] = {1,2,3,4,5,6,7,8,9};
    int g_src_num = sizeof(g_src_idx) / sizeof(g_src_idx[0]);
    Xt_idxlist g_src_idxlist = xt_idxvec_new(g_src_idx, g_src_num);

    Xt_int g_dst_idx[] = {9,8,7,6,5,4,3,2,1};
    int g_dst_num = sizeof(g_dst_idx) / sizeof(g_dst_idx[0]);
    Xt_idxlist g_dst_idxlist = xt_idxvec_new(g_dst_idx, g_dst_num);

    struct Xt_modifier mod
      = {.extract=g_src_idxlist, .subst=g_dst_idxlist, .mask=0};

    Xt_int patch_idx[] = {3,4,4,4,7,7,8};
    int patch_num = sizeof(patch_idx) / sizeof(patch_idx[0]);
    Xt_idxlist patch_idxlist = xt_idxvec_new(patch_idx, patch_num);

    // idx:{3,4,4,4,7,7,8} -> pos:{2,3,3,3,6,6,7} => idx:{7,6,6,6,3,3,2}
    Xt_int ref_mpatch_idx[] = {7,6,6,6,3,3,2};

    Xt_idxlist mpatch_idxlist = xt_idxmod_new(patch_idxlist, &mod, 1, NULL);

    check_idxlist(mpatch_idxlist, ref_mpatch_idx, patch_num);

    xt_idxlist_delete(mpatch_idxlist);
    xt_idxlist_delete(patch_idxlist);
    xt_idxlist_delete(g_dst_idxlist);
    xt_idxlist_delete(g_src_idxlist);
  }

  { // idxstripes modifier

    struct Xt_stripe g_src_stripe[] = { {.start=1, .nstrides=20, .stride=1} };
    int g_src_num = sizeof(g_src_stripe) / sizeof(g_src_stripe[0]);
    Xt_idxlist g_src_idxlist = xt_idxstripes_new(g_src_stripe, g_src_num);

    struct Xt_stripe g_dst_stripe[]
      = { {.start=100, .nstrides=20, .stride=-1} };
    int g_dst_num = sizeof(g_dst_stripe) / sizeof(g_dst_stripe[0]);
    Xt_idxlist g_dst_idxlist = xt_idxstripes_new(g_dst_stripe, g_dst_num);

    struct Xt_modifier mod
      = {.extract=g_src_idxlist, .subst=g_dst_idxlist, .mask=32};

    Xt_int patch_idx[] = {0,1,3,3,5,50,100,150};
    int patch_num = sizeof(patch_idx) / sizeof(patch_idx[0]);
    Xt_idxlist patch_idxlist = xt_idxvec_new(patch_idx, patch_num);

    // inter:{1,3,3,5} => extract_pos:{0,2,2,4} => subst_idx:{100,98,98,96},
    // patch_pos:{1,2,3,4} = > mpatch:{0,100,98,98,96,50,100,150}
    Xt_int ref_mpatch_idx[] = {0,100,98,98,96,50,100,150};

    int mstate[]   = {1, 2, 3, 4, 5, 6, 7, 8};
    int ref_mstate[] = {1, 2|32, 3|32, 4|32, 5|32, 6, 7, 8};

    Xt_idxlist mpatch_idxlist = xt_idxmod_new(patch_idxlist, &mod, 1, mstate);

    check_idxlist(mpatch_idxlist, ref_mpatch_idx, patch_num);

    // check mstate:
    for(int i=0; i< patch_num; i++) {
      if (mstate[i] != ref_mstate[i])
        PUT_ERR("(mstate[i] != ref_mstate[i])\n");
    }

    xt_idxlist_delete(mpatch_idxlist);
    xt_idxlist_delete(patch_idxlist);
    xt_idxlist_delete(g_dst_idxlist);
    xt_idxlist_delete(g_src_idxlist);

  }

  { // multiple modifier, track modifier usage

    Xt_int g1_src_idx[] = {1,2,3,4,5,6,7,8,9};
    int g1_src_num = sizeof(g1_src_idx) / sizeof(g1_src_idx[0]);
    Xt_idxlist g1_src_idxlist = xt_idxvec_new(g1_src_idx, g1_src_num);

    Xt_int g1_dst_idx[] = {9,8,7,6,5,4,3,2,1};
    int g1_dst_num = sizeof(g1_dst_idx) / sizeof(g1_dst_idx[0]);
    Xt_idxlist g1_dst_idxlist = xt_idxvec_new(g1_dst_idx, g1_dst_num);

    Xt_int g2_src_idx[] = {1,2,8,9,10};
    int g2_src_num = sizeof(g2_src_idx) / sizeof(g2_src_idx[0]);
    Xt_idxlist g2_src_idxlist = xt_idxvec_new(g2_src_idx, g2_src_num);

    Xt_int g2_dst_idx[] = {8,2,8,2,5};
    int g2_dst_num = sizeof(g2_dst_idx) / sizeof(g2_dst_idx[0]);
    Xt_idxlist g2_dst_idxlist = xt_idxvec_new(g2_dst_idx, g2_dst_num);

    struct Xt_modifier mod[2]
      = { {.extract=g1_src_idxlist, .subst=g1_dst_idxlist, .mask=1},
          {.extract=g2_src_idxlist, .subst=g2_dst_idxlist, .mask=2} };

    Xt_int patch_idx[] = {6,7,25,8,9,10};
    int patch_num = sizeof(patch_idx) / sizeof(patch_idx[0]);
    Xt_idxlist patch_idxlist = xt_idxvec_new(patch_idx, patch_num);
    int ref_mstate[] = {1|0, 1|0, 0|0, 1|2, 1|2, 0|2};
    int mstate[patch_num];
    // reset mstate:
    for (int i=0; i<patch_num; i++) {
      mstate[i] = 0;
    }

    // mod1: idx:{6,7,25,8,9,10} -> pos:{5,6,nil,7,8,nil} => idx:{4,3,25,2,1,10}
    // mod2: idx:{4,3,25,2,1,10} -> pos:{nil,nil,nil,1,0,4}
    //                           => idx:{4,3,25,2,8,5}
    Xt_int ref_mpatch_idx[] = {4,3,25,2,8,5};

    Xt_idxlist mpatch_idxlist = xt_idxmod_new(patch_idxlist, mod, 2, mstate);

    check_idxlist(mpatch_idxlist, ref_mpatch_idx, patch_num);

    // check mstate:
    for(int i=0; i< patch_num; i++) {
      if (mstate[i] != ref_mstate[i])
        PUT_ERR("(mstate[i] != ref_mstate[i])\n");
    }

    xt_idxlist_delete(mpatch_idxlist);
    xt_idxlist_delete(patch_idxlist);
    xt_idxlist_delete(g2_dst_idxlist);
    xt_idxlist_delete(g2_src_idxlist);
    xt_idxlist_delete(g1_dst_idxlist);
    xt_idxlist_delete(g1_src_idxlist);
  }

  // finalize:
  xt_finalize();
  MPI_Finalize();
  return TEST_EXIT_CODE;
}
