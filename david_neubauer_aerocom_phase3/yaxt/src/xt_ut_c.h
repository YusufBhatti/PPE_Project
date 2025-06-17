/**
 * @file xt_ut_c.h
 * @brief supportes the unitrans interfaces/wrappers in xt_ut.f90 with yaxt functionality.
 *
 * @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
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

#ifndef XT_UT_C_H
#define XT_UT_C_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/**
 * @example test_ut.f90
 */

/**
 * Initialization of support data for unitrans (start sizes of allocated data).
 * We use the unitrans names. Unlike unitrans, the allocated data can grow.
 *
 * @param[in] decomp_size       number of decompositions
 * @param[in] comm_tmpl_size    number of communication templates
 * @param[in] comm_size         number of transpositions
 * @param[in] debug_lvl         debug level
 * @param[in] mode              default communication mode
 * @param[in] debug_unit        fortran file unit for debug output (not used)
 */
void xt_ut_init(int decomp_size, int comm_tmpl_size, int comm_size, int debug_lvl, int mode, int debug_unit);

/**
 * creates new decomposition and returns access handle
 *
 * @param[in] idx_vec           local index vector
 * @param[in] idx_vec_n         size of idx_vec
 * @return returns handle
 */
int xt_ut_init_decomposition_1d(Xt_int *idx_vec, int idx_vec_n);

/**
 * creates new unitrans compatible oneway transposition decomposition and returns access handle
 *
 * @param[in] decomp_handle_in   access handle for input decomposition
 * @param[in] decomp_handle_out  access handle for output decomposition
 * @param[in] mpi_world          mpi communicator
 * @param[in] icheck_unique      enable/disbale check of uniqueness
 * @return returns handle
 */
int xt_ut_init_oneway_transposition_template(int decomp_handle_in, int decomp_handle_out,
                                             int mpi_world, int icheck_unique);
/**
 * destroys transposition template
 *
 * @param[in] id                 transposition handle
 */
void xt_ut_destroy_transposition_template(int id);

#endif
