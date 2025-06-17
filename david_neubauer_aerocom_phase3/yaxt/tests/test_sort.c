/**
 * @file test_sort.c
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

#include <assert.h>

//#define VERBOSE

#include <yaxt.h>

#include "tests.h"

static void
test1(void (*anysort_index)(Xt_int * a, int n, int * idx, int reset_index)) {
  int tsize = 8;
  Xt_int ivec[]        = {7,5,3,1,2,2,3,3};
  Xt_int sorted_ivec[] = {1,2,2,3,3,3,5,7};
  int pvec[]           = {1,2,3,4,5,6,7,8};
  int sorted_pvec[]    = {3,4,5,2,6,7,1,0};

  tsize = sizeof(ivec) /sizeof(ivec[0]);
  assert(tsize == sizeof(pvec) /sizeof(pvec[0]));

  anysort_index(ivec, tsize, pvec, 1);
  for(int i=0; i<tsize; i++) {
    if ( ivec[i] != sorted_ivec[i] )
      PUT_ERR("wrong sorting values\n");
    if ( pvec[i] != sorted_pvec[i] )
      PUT_ERR("wrong sorting positions\n");
  }
}

int main(void) {

  test1(xt_quicksort_index);
  test1(xt_mergesort_index);
  test1(xt_sort_index);

  return TEST_EXIT_CODE;
}
