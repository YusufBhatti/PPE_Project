/**
 * @file mergesort.c
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
#include <stdio.h>

#include "xt/quicksort.h"
#include "xt/mergesort.h"

static void insertionsort_idxpos(idxpos_type *w, int n);
static void insertionsort_index(Xt_int *val, int n, int *pos, int reset_pos);


static void insertionsort_idxpos(idxpos_type *w, int n) {

  int i, j;
  for(i = 1; i < n; i++) {

    idxpos_type tw = w[i];

    for(j = i;  j > 0 && w[j-1].idx > tw.idx; j--) {
      w[j] = w[j-1];
    }
    w[j] = tw;
  }
}

static void insertionsort_index(Xt_int *val, int n, int *pos, int reset_pos) {

  int i, j;

  if (pos != NULL) {

    if (reset_pos) {
      for(i = 0; i < n; i++) {
        pos[i]=i;
      }
    }

    for(i = 1; i < n; i++) {

      Xt_int tv = val[i];
      int tp = pos[i];

      for(j = i;  j > 0 && val[j-1] > tv; j--) {
        val[j] = val[j-1];
        pos[j] = pos[j-1];
      }

      val[j] = tv;
      pos[j] = tp;
    }

  } else {

    for(i = 1; i < n; i++) {

      Xt_int tv = val[i];

      for(j = i;  j > 0 && val[j-1] > tv; j--) {
        val[j] = val[j-1];
      }

      val[j] = tv;
    }

  }

}

static inline void
my_merge_idxpos(idxpos_type * restrict v, idxpos_type * restrict w, int lb1, int ub1, int lb2, int ub2) {

  int p  = lb1;
  int i1 = lb1;
  int i2 = lb2;

  while(i1<=ub1 && i2<=ub2) {
    if (v[i1].idx <= v[i2].idx) {
      w[p]=v[i1];
      i1++;
    } else {
      w[p]=v[i2];
      i2++;
    }
    p++;
  }
  while(i1<=ub1) {
    w[p]=v[i1];
    i1++;
    p++;
  }
  while(i2<=ub2) {
    w[p]=v[i2];
    i2++;
    p++;
  }
}

static void
my_mergesort_idxpos(idxpos_type * restrict v, idxpos_type * restrict w, int lb, int ub)
{

  int n = ub - lb + 1;
  if (n<2) return;

  if (n<9) {
    insertionsort_idxpos(&v[lb],n);
    return;
  }

  int m   = n/4;

  // 4-way subsort:
  int lb1 = lb;
  int ub1 = lb1 + m - 1;

  int lb2 = ub1+1;
  int ub2 = lb2 + m - 1;

  int lb3 = ub2+1;
  int ub3 = lb3 + m - 1;

  int lb4 = ub3+1;
  int ub4 = ub;

  my_mergesort_idxpos(v, w, lb1, ub1);
  my_mergesort_idxpos(v, w, lb2, ub2);
  my_mergesort_idxpos(v, w, lb3, ub3);
  my_mergesort_idxpos(v, w, lb4, ub4);

  // 2 x 2-way merge v -> w
  my_merge_idxpos(v, w, lb1, ub1, lb2, ub2);
  my_merge_idxpos(v, w, lb3, ub3, lb4, ub4);
  // final merge: w -> v
  my_merge_idxpos(w, v, lb1, ub2, lb3, ub4);

}

void xt_mergesort_idxpos(idxpos_type * restrict v, int n)
{

  idxpos_type *w;

  w = malloc((size_t)n * sizeof(w[0]));
  my_mergesort_idxpos(v, w, 0, n-1);
  free(w);

}

void xt_mergesort_index(Xt_int *val, int n, int *pos, int reset_pos) {

  // insertion sort for few elements:
  if (n<9) {
    insertionsort_index(val,n,pos,reset_pos);
    return;
  }

  idxpos_type *v;
  idxpos_type *w;

  v = malloc((size_t)n * sizeof(v[0]));
  w = malloc((size_t)n * sizeof(w[0]));

  if (reset_pos || pos == NULL) {
    for(int i=0; i<n; i++) {
      v[i].idx = val[i];
      v[i].pos = i;
    }
  } else {
    for(int i=0; i<n; i++) {
      v[i].idx = val[i];
      v[i].pos = pos[i];
    }
  }

  my_mergesort_idxpos(v,w, 0, n-1);
  if (pos) {
    for(int i=0; i<n; i++) {
      val[i] = v[i].idx;
      pos[i] = v[i].pos;
    }
  } else {
    for(int i=0; i<n; i++) {
      val[i] = v[i].idx;
    }
  }

  free(w);
  free(v);
}

