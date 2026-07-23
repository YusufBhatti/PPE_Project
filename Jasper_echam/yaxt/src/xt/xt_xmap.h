/**
 * @file xt_xmap.h
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

#ifndef XT_XMAP_H
#define XT_XMAP_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "xt/xt_core.h"

/** \example test_xmap_common.c
 */
/** \example test_xmap_common_f.f90
 */
/** \example test_xmap_common_parallel.c
 */
/** \example test_xmap_common_parallel_f.f90
 */

/**
 * \file xt_xmap.h
 * \brief exchange map declarations
 *
 * methods to be used with Xt_xmap
 */

typedef struct Xt_xmap_iter_ *Xt_xmap_iter;

/**
 * gets the MPI communicator of the exchange map
 * @param[in] xmap exchange map
 * @return         MPI communicator of the exchange map
 */
MPI_Comm xt_xmap_get_communicator(Xt_xmap xmap);

/**
 * gets the number of processes that will require data from
 * the local process in an exchange operation
 *
 * @param[in] xmap exchange map
 * @return         number of processes that will require data from
 *                 the local process in an exchange operation
 */
int xt_xmap_get_num_destinations(Xt_xmap xmap);

/**
 * gets the number of processes from which the local process will
 * require data in an exchange operation
 *
 * @param[in] xmap exchange map
 * @return number of prosses from which the local process will
 *         required data in an exchange operation
 */
int xt_xmap_get_num_sources(Xt_xmap xmap);

/**
 * gets the ranks of processes that will require data from
 * the local process in an exchange operation
 *
 * @param[in]  xmap  exchange map
 * @param[out] ranks ranks of processes that will require data
 *                   from the local process in an exchange
 *                   operation
 *
 * \remarks the user needs to provide an array for ranks that
 *          is sufficient in size to store all ranks
 */
void xt_xmap_get_destination_ranks(Xt_xmap xmap, int * ranks);

/**
 * gets the ranks of processes from which the local process will
 * require data in an exchange operation
 *
 * @param[in]  xmap  exchange map
 * @param[out] ranks ranks of processes from which the local process
 *                   will required data in an exchange operation
 *
 * \remarks the user needs to provide an array for ranks that
 *          is sufficient in size to store all ranks
 */
void xt_xmap_get_source_ranks(Xt_xmap xmap, int * ranks);

/**
 * gets an iterator for all local outgoing messages of an
 * exchange operation
 *
 * @param[in] xmap exchange map
 * @return         iterator for all local outgoing messages\n
 *                 is NULL if there are no local outgoing messages
 */
Xt_xmap_iter xt_xmap_get_out_iterator(Xt_xmap xmap);

/**
 * gets an iterator for all local incoming messages of an
 * exchange operation
 *
 * @param[in] xmap exchange map
 * @return         iterator for all local incoming messages\n
 *                 is NULL if there are no local incoming messages
 */
Xt_xmap_iter xt_xmap_get_in_iterator(Xt_xmap xmap);

/**
 * sets an exchange map iterator to the next item
 *
 * @param[in,out] iter exchange map iterator
 * @return             returns 0 if there was not further item
 */
int xt_xmap_iterator_next(Xt_xmap_iter iter);

/**
 * gets the rank of the current item of an exchange map iterator
 *
 * @param[in] iter exchange map iterator
 * @return         rank of the current item
 */
int xt_xmap_iterator_get_rank(Xt_xmap_iter iter);

/**
 * get the positions of all element of the current item of an exchange
 * map iterator
 *
 * @param[in] iter exchange map iterator
 * @return         exchange list of the current item
 */
int const * xt_xmap_iterator_get_transfer_pos(Xt_xmap_iter iter);

/**
 * get the number of positions returned by
 * \ref xt_xmap_iterator_get_transfer_pos of the current item of an exchange
 * map iterator
 * @param[in] iter exchange map iterator
 * @return         number of positions returned by
 *                 \ref xt_xmap_iterator_get_transfer_pos for the current item
 */
int xt_xmap_iterator_get_num_transfer_pos(Xt_xmap_iter iter);

/**
 * get the position ranges of all element of the current item of an exchange
 * map iterator as position extents
 *
 * @param[in] iter exchange map iterator
 * @return         exchange list of the current item
 */
const struct Xt_pos_ext* xt_xmap_iterator_get_transfer_pos_ext(Xt_xmap_iter iter);

/**
 * get the number of position extents returned by
 * \ref xt_xmap_iterator_get_transfer_pos_ext of the current item of an exchange
 * map iterator
 * @param[in] iter exchange map iterator
 * @return         number of positions returned by
 *                 \ref xt_xmap_iterator_get_transfer_pos_ext for the current item
 */
int xt_xmap_iterator_get_num_transfer_pos_ext(Xt_xmap_iter iter);

/**
 * destructor for an exchange map iterator
 *
 * @param[in] iter exchange map iterator
 */
void xt_xmap_iterator_delete(Xt_xmap_iter iter);

/**
 * destructor for an exchange map
 *
 * @param[in] xmap exchange map to be destructed
 */
void xt_xmap_delete(Xt_xmap xmap);

/**
 * gives an upper limit to the range of positions in the src part of the transfer_pos,
 *
 * @param[in] xmap exchange map
 * @return         upper limit of range of positions in the source list
 */
int xt_xmap_get_max_src_pos(Xt_xmap xmap);

/**
 * gives an upper limit to the range of positions in the dst part of the transfer_pos,
 *
 * @param[in] xmap exchange map
 * @return         upper limit of range of positions in the destination list
 */
int xt_xmap_get_max_dst_pos(Xt_xmap xmap);

/**
 * inquiry function to check if the inverse is defined \n
 * this call is collective among of processes that are
 * associated to the exchange map
 *
 * @param[in] xmap exchange map the is to be checked
 * @return         returns != 0 if the exchange map is bijective
 */
int xt_xmap_is_bijective(Xt_xmap xmap);

#endif // XT_XMAP_H
