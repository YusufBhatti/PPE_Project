/**
 * @file quicksort.c
 * @brief Non-recursive stack version of Quicksort
 *
 *  based on N. Wirth's Pascal Book, 'Algorithms + Data Structures = Programms'.
 *  by Alan Miller ( 19 Jul 1995 )
 *
 *  based on:
 * - http://www.nag.com/nagware/examples.asp
 * - http://www.nag.com/nagware/Examples/nur.f90
 *
 *  see also:
 * - http://en.wikipedia.org/wiki/Quicksort
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

#include <stdlib.h>

#include "xt/quicksort.h"
#include "core/ppm_xfuncs.h"

void xt_quicksort_idxpos(idxpos_type  *v, int n) {

  enum { pointer_inc = 64 };
  unsigned pointer_size;

  int * stackl;
  int * stackr;

  int i, j, k, l, r;
  unsigned s;
  idxpos_type w, x;

  /*  Initialization */

  pointer_size = pointer_inc;

  stackl = xmalloc((size_t)pointer_size * sizeof(stackl[0]));
  stackr = xmalloc((size_t)pointer_size * sizeof(stackr[0]));

  s = 1;
  stackl[0] = 0;
  stackr[0] = n-1;

  /*  Start sorting

  ... keep taking the top request from the stack until s = 0.

  */

  while  ( s > 0 ) {

    --s;
    l = stackl[s];
    r = stackr[s];

    /* ... keep splitting a[l], ... ,a[r] until l>= r. */

    while ( l < r ) {

      i = l;
      j = r;
      k = (l+r) / 2;
      x = v[k];

      /* Search from lower end */

      while(1) {
        while(1)
          if ((v[i].idx < x.idx) || (v[i].idx == x.idx && v[i].pos < x.pos))
            i = i + 1;
          else
            break;

        /* Search from upper end */

        while(1)
          if ((x.idx < v[j].idx) || (x.idx == v[j].idx && x.pos < v[j].pos))
            j = j - 1;
          else
            break;

        /* Swap positions i & j */

        if ( i <= j ) {

          w      = v[i];
          v[i]   = v[j];
          v[j]   = w;

          i = i + 1;
          j = j - 1;

          if ( i > j ) break;
        } else
          break;
      }

      if ( j-l >= r-i ) {

        if ( l < j ) {

          if ( s >= pointer_size ) {
            pointer_size = pointer_size + pointer_inc;
            stackl = xrealloc(stackl, (size_t)pointer_size * sizeof(stackl[0]));
            stackr = xrealloc(stackr, (size_t)pointer_size * sizeof(stackr[0]));
          }

          stackl[s] = l;
          stackr[s] = j;
          ++s;
        }
        l = i;
      } else {

        if ( i < r ) {

          if ( s >= pointer_size ) {
            pointer_size = pointer_size + pointer_inc;
            stackl = xrealloc(stackl, (size_t)pointer_size * sizeof(stackl[0]));
            stackr = xrealloc(stackr, (size_t)pointer_size * sizeof(stackr[0]));
          }

          stackl[s] = i;
          stackr[s] = r;
          ++s;
        }
        r = j;
      }
    } /* ( l < r ) */
  } /* ( s /= 0 ) */

  free(stackl);
  free(stackr);
}

void xt_quicksort_index(Xt_int * v_idx, int n, int * v_pos, int reset_pos) {

  enum { pointer_inc = 64 };
  unsigned pointer_size;

  int * stackl;
  int * stackr;

  int i, j, k, l, r;
  unsigned s;
  Xt_int w_idx, x_idx;
  int w_pos, x_pos;

  /*  Initialization */
  int *v_pos_orig =  v_pos;

  if (!v_pos) v_pos = malloc((size_t)n  * sizeof(v_pos[0]));

  if (v_pos != v_pos_orig || reset_pos) {
    for(i=0; i<n; i++) {
      v_pos[i] = i;
    }
  }

  pointer_size = pointer_inc;

  stackl = xmalloc((size_t)pointer_size * sizeof(stackl[0]));
  stackr = xmalloc((size_t)pointer_size * sizeof(stackr[0]));

  s = 1;
  stackl[0] = 0;
  stackr[0] = n-1;

  /*  Start sorting

  ... keep taking the top request from the stack until s = 0.

  */

  while  ( s > 0 ) {

    --s;
    l = stackl[s];
    r = stackr[s];

    /* ... keep splitting a[l], ... ,a[r] until l>= r. */

    while ( l < r ) {

      i = l;
      j = r;
      k = (l+r) / 2;
      x_idx = v_idx[k];
      x_pos = v_pos[k];

      /* Search from lower end */

      while(1) {
        while(1)
          if ((v_idx[i] < x_idx) || (v_idx[i] == x_idx && v_pos[i] < x_pos))
            i = i + 1;
          else
            break;

        /* Search from upper end */

        while(1)
          if ((x_idx < v_idx[j]) || (x_idx == v_idx[j] && x_pos < v_pos[j]))
            j = j - 1;
          else
            break;

        /* Swap positions i & j */

        if ( i <= j ) {

          w_idx      = v_idx[i];
          v_idx[i]   = v_idx[j];
          v_idx[j]   = w_idx;

          w_pos      = v_pos[i];
          v_pos[i]   = v_pos[j];
          v_pos[j]   = w_pos;

          i = i + 1;
          j = j - 1;

          if ( i > j ) break;
        } else
          break;
      }

      if ( j-l >= r-i ) {

        if ( l < j ) {

          if ( s >= pointer_size ) {
            pointer_size = pointer_size + pointer_inc;
            stackl = xrealloc(stackl, (size_t)pointer_size * sizeof(stackl[0]));
            stackr = xrealloc(stackr, (size_t)pointer_size * sizeof(stackr[0]));
          }

          stackl[s] = l;
          stackr[s] = j;
          ++s;
        }
        l = i;
      } else {

        if ( i < r ) {

          if ( s >= pointer_size ) {
            pointer_size = pointer_size + pointer_inc;
            stackl = xrealloc(stackl, (size_t)pointer_size * sizeof(stackl[0]));
            stackr = xrealloc(stackr, (size_t)pointer_size * sizeof(stackr[0]));
          }

          stackl[s] = i;
          stackr[s] = r;
          ++s;
        }
        r = j;
      }
    } /* ( l < r ) */
  } /* ( s /= 0 ) */

  free(stackl);
  free(stackr);
  if (v_pos != v_pos_orig) free(v_pos);
}

