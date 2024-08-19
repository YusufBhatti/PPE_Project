/**
 * @file xt_redist.h
 * @brief redistribution of data
 *
 * contains declaration the redistribution data structure, which
 * is derived from one or more xt_xmaps
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

#ifndef XT_REDIST_H
#define XT_REDIST_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>

#include "xt/xt_core.h"

/**
 * destructor
 *
 * @param[in,out] redist redistribution structure
 */
void xt_redist_delete(Xt_redist redist);

/**
 * synchronous redistribution of data
 *
 * @param[in]     redist     redistribution structure
 * @param[in]     num_arrays number of base addresses in src_data and dst_data
 * @param[in]     src_data   array containing the addresses of the first
 *                           elements of the input data
 * @param[in,out] dst_data   array containing the addresses of the first
 *                           elements of the output data
 *
 * @remark The above implies that NULL or any other invalid pointer
 * must not be used in either @a src_data or @a dst_data.
 */
void xt_redist_s_exchange(Xt_redist redist, int num_arrays,
                          const void **src_data, void **dst_data);

/**
 * synchronous redistribution of data - single array case
 *
 * @param[in]     redist   redistribution structure
 * @param[in]     src_data address of the first element of the input data
 * @param[in,out] dst_data address of the first element of the output data
 *
 * @remark The above implies that NULL or any other invalid pointer
 * must not be used in either @a src_data or @a dst_data.
 */
void xt_redist_s_exchange1(Xt_redist redist, const void *src_data, void *dst_data);

/**
 * gets a copy of the MPI_Datatype used for the data of the send operation with
 * the given rank
 *
 * @param[in] redist redistribution structure
 * @param[in] rank   MPI rank
 * \return MPI_Datatype for the data of the send operation with the given rank
 *
 * \remarks returns MPI_DATATYPE_NULL if there is no send operation with the given rank
 */
MPI_Datatype xt_redist_get_send_MPI_Datatype(Xt_redist redist, int rank);

/**
 * gets a copy of the MPI_Datatype used for the data of the recv operation with
 * the given rank
 *
 * @param[in] redist redistribution structure
 * @param[in] rank   MPI rank
 * \return MPI_Datatype for the data of the recv operation with the given rank
 *
 * \remarks returns MPI_DATATYPE_NULL if there is no recv operation with the given rank
 */
MPI_Datatype xt_redist_get_recv_MPI_Datatype(Xt_redist redist, int rank);

/**
 * returns a MPI communicator, which the redistribution is based on
 * @param[in] redist redistribution structure
 * \return MPI communicator, which the redistribution is based on
 */
MPI_Comm xt_redist_get_MPI_Comm(Xt_redist redist);

#endif // XT_REDIST_H
