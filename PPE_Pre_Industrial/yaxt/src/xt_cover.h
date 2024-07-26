/**
 * @file xt_cover.h
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
#ifndef XT_COVERAGE_H
#define XT_COVERAGE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdbool.h>
#include <stdlib.h>

#include "xt/xt_core.h"

struct Xt_pos_ext_vec {
  size_t num_pos_ext, size_pos_ext;
  struct Xt_pos_ext *pos_ext;
};

struct Xt_pos_range {
  int start, end;
};

void
xt_cover_start(struct Xt_pos_ext_vec *restrict cover,
               size_t initial_size);

void
xt_cover_finish(struct Xt_pos_ext_vec *restrict cover);

/**
 * @param search_start_pos should be 0, if \a forward == true
 * or cover->num_pos_ext, if \a forward == false to search all of
 * cover, choose value to start search at if part of cover is known to
 * be non-matching
 * @return index of first position extent in cover that overlaps or is
 * adjacent to \a query
 */
size_t
xt_cover_search(struct Xt_pos_ext_vec *restrict cover,
                struct Xt_pos_range query, bool forward,
                size_t search_start_pos);

/**
 * append \range to cover
 *
 * @note user must ensure range actually is appendable, i.e. does not
 * overlap or precede an existing range in \a cover
 */
void
xt_cover_range_append(struct Xt_pos_ext_vec *restrict cover,
                      struct Xt_pos_ext range);

/**
 * @param search_start_pos should be 0, if \a forward == true
 * or cover->num_pos_ext, if \a forward == false to search all of
 * cover, choose value to start search at if part of cover is known to
 * be non-matching
 * @return SIZE_MAX if \a range could be fully integrated with \a
 * cover, position \a i of overlapping \a cover.pos_ext[i] otherwise
 */
size_t
xt_cover_insert_or_overlap(struct Xt_pos_ext_vec *restrict cover,
                           struct Xt_pos_range range, bool forward,
                           size_t search_start_pos);

/**
 * tests if sorted pos_ext in \a coverage do indeed fully cover \a idxlist
 */
bool
xt_idxlist_pos_ext_is_full_cover(Xt_idxlist idxlist,
                                 struct Xt_pos_ext_vec cover);


#endif
