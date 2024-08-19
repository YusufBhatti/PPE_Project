/**
 * @file xt_handles.c
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
#include <stdlib.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt/xt_handles.h"


struct Xt_handle_set_type_ {
  unsigned num; // number of used handles
  unsigned cap; // capacity
  unsigned hnf; // hint for next free entry
  void **p; // vector of pointers to user data, p[handle] == pointer to user data
};

static const unsigned default_handle_set_cap = 16;

static void extend_handle_set(Xt_handle_set_type hset);

int xt_handle_is_valid(Xt_handle_set_type hset, int ih) {
  if(!hset) return 0;
  if (ih < 0 || (unsigned)ih >= hset->cap) return 0;
  if (hset->p[ih]) return 1;
  return 0;
}

void *xt_handle2pointer(Xt_handle_set_type hset, int ih) {
  if (!xt_handle_is_valid(hset, ih)) return NULL;
  return hset->p[ih];
}

Xt_handle_set_type xt_handle_set_new(int cap) {
  Xt_handle_set_type hset;
  hset = malloc(sizeof (*hset));
  assert(cap >= 0);
  hset->p = malloc((unsigned)cap * sizeof(void*));
  hset->num = 0;
  hset->cap = (unsigned)cap;
  if (hset->cap < 1) hset->cap = default_handle_set_cap;
  hset->hnf = 0;
  for (unsigned i=0; i<(unsigned)cap; i++) {
    hset->p[i] = NULL;
  }

  return hset;
}


static void extend_handle_set(Xt_handle_set_type hset) {
  if (!hset) {
    hset = xt_handle_set_new((int)default_handle_set_cap);
    return;
  }

  unsigned old_cap, new_cap;
  old_cap = hset->cap;
  if (old_cap < 2) {
    new_cap = default_handle_set_cap;
  } else {
    new_cap = old_cap*2;
  }

  hset->p = xrealloc(hset->p, (unsigned)new_cap * sizeof(void*));

  for (unsigned i=old_cap; i<new_cap; i++) {
    hset->p[i] = NULL;
  }
  hset->hnf = old_cap;
  hset->cap = new_cap;
}


void xt_handle_set_delete(Xt_handle_set_type hset) {
  free(hset->p);
  free(hset);
  hset=NULL;
}


int xt_handle_new(Xt_handle_set_type hset, void *p) {
  if (p == NULL) return -1;

  if (hset->num >= hset->cap) extend_handle_set(hset);
  for (unsigned j=0; j<hset->cap; j++) {
    unsigned i = (hset->hnf + j) % hset->cap;
    if (!hset->p[i]) {
      hset->p[i] = p;
      hset->hnf = (i+1) % hset->cap;
      hset->num++;
      return (int)i;
    }
  }

  die("internal error");
  /* GNU C realizes die does not return, other compilers lack this smart */
#ifndef __GNUC__
  return 0;
#endif
}

void xt_handle_delete(Xt_handle_set_type hset, int handle) {
  if (!xt_handle_is_valid(hset, handle)) die("delete_handle: invalid handle");
  hset->p[handle]=NULL;
  if (hset->num == 0) die("delete_handle: internal error - negative occupation count");
  hset->num--;
  if (hset->p[hset->hnf]) hset->hnf = (unsigned)handle;
}

