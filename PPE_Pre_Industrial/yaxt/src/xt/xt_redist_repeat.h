/**
 * @file xt_redist_repeat.h
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

#ifndef XT_REDIST_REPEAT_H
#define XT_REDIST_REPEAT_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "xt/xt_redist.h"

/** \example test_redist_repeat.c
 */
/** \example test_redist_repeat_f.f90
 */
/** \example test_redist_repeat_parallel.c
 */
/** \example test_redist_repeat_parallel_f.f90
 */

/**
 * constructor for a redistribution that has a repetitive pattern, which
 * is described by the given redistribution
 * @param[in] redist          redistribution
 * @param[in] src_extent      extent that scales the given displacements
 *                            for the source data
 * @param[in] dst_extent      extent that scales the given displacements
 *                            for the destination data
 * @param[in] num_repetitions number of repetitions of the given
 *                            redistribution
 * @param[in] displacements   displacements for each repetition
 *
 * \remarks all redistributions need to be based on the same
 *          MPI communicator
 */
Xt_redist xt_redist_repeat_new(Xt_redist redist, MPI_Aint src_extent,
                               MPI_Aint dst_extent, int num_repetitions,
                               int displacements[num_repetitions]);

#endif // XT_REDIST_REPEAT_H