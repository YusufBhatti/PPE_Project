/**
 * @file xt_ut_c.c
 * @brief implementation of unitrans using yaxt
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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "xt/xt_mpi.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_idxvec.h"
#include "xt_ut_c.h"
#include "xt/xt_xmap.h"
#include "xt/xt_xmap_all2all.h"
#include "xt/xt_redist_p2p.h"

#include "xt/xt_handles.h"

typedef struct {
  int decomp_size;
  int comm_tmpl_size;
  int comm_trans_size;
  int debug_lvl;
  int mode;
  int debug_unit;
} xt_ut_init_type;

static xt_ut_init_type init_val = {0,0,0,0,0,0};


typedef struct {
  Xt_idxlist idxvec;
  int domain_size;
} xt_ut_deco_type;

static Xt_handle_set_type deco_set;

typedef struct {
  Xt_xmap xmap;
} xt_ut_template_type;

static Xt_handle_set_type template_set;


typedef struct {
  Xt_redist redist;
} xt_ut_trans_type;


static Xt_handle_set_type trans_set;

/* the external functions declared here are intended to be called from
 * Fortran, therefore there is no prototype declaration in a
 * corresponding .h file.
 */
#if (defined __GNUC__ && __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 5))\
  || (defined __clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-prototypes"
#endif

void xt_ut_abort(char *msg, char *source, int *line) __attribute__((noreturn));

void xt_ut_abort(char *msg, char *source, int *line) {
  Xt_abort(MPI_COMM_WORLD, msg, source, *line);
}

void xt_ut_init(int decomp_size, int comm_tmpl_size, int comm_trans_size, int debug_lvl, int mode, int debug_unit) {

  xt_initialize(MPI_COMM_WORLD);

  init_val.decomp_size     = decomp_size;
  init_val.comm_tmpl_size  = comm_tmpl_size;
  init_val.comm_trans_size = comm_trans_size;
  init_val.debug_lvl       = debug_lvl;
  init_val.mode            = mode;
  init_val.debug_unit      = debug_unit;

  deco_set = xt_handle_set_new(init_val.decomp_size);

  template_set = xt_handle_set_new(init_val.comm_tmpl_size);

  trans_set = xt_handle_set_new(init_val.comm_trans_size);
}

void xt_ut_finalize(void) {
  xt_finalize();
}

static int new_ideco(void) {
  xt_ut_deco_type *data = xmalloc(sizeof(*data));
  return xt_handle_new(deco_set, data);
}

static void delete_ideco(int id) {
  xt_ut_deco_type *deco = xt_handle2pointer(deco_set, id);
  free(deco);
  xt_handle_delete(deco_set, id);
}

void xt_ut_destroy_decomposition(int id) {
  xt_ut_deco_type *deco = xt_handle2pointer(deco_set, id);
  if (!deco) die("xt_ut_destroy_decomposition: invalid handle");
  xt_idxlist_delete(deco->idxvec);
  delete_ideco(id);
}

static int new_itemplate(void) {
  xt_ut_template_type *data = xmalloc(sizeof(*data));
  return xt_handle_new(template_set, data);
}

static void delete_itemplate(int id) {
  xt_ut_template_type *template = xt_handle2pointer(template_set, id);
  free(template);
  xt_handle_delete(template_set, id);
}

void xt_ut_destroy_transposition_template(int id) {
  xt_ut_template_type *template = xt_handle2pointer(template_set, id);
  if (!template) die("xt_ut_destroy_transposition_template: invalid handle");
  xt_xmap_delete(template->xmap);
  delete_itemplate(id);
}

static int new_itrans(void) {
  xt_ut_trans_type *data = xmalloc(sizeof(*data));
  return xt_handle_new(trans_set, data);
}

static void delete_itrans(int id) {
  xt_ut_trans_type *trans = xt_handle2pointer(trans_set, id);
  free(trans);
  xt_handle_delete(trans_set, id);
}

void xt_ut_destroy_transposition(int id) {
  xt_ut_trans_type *trans = xt_handle2pointer(trans_set, id);
  if (!trans) die("xt_ut_destroy_transposition: invalid handle");
  xt_redist_delete(trans->redist);
  delete_itrans(id);
}


MPI_Fint xt_ut_init_decomposition_1d(Xt_int *iv, int iv_n)
{
  MPI_Fint id = new_ideco();
  xt_ut_deco_type *deco = xt_handle2pointer(deco_set, id);
  if (!deco) die("xt_ut_init_decomposition_1d: internal error");
  deco->idxvec = xt_idxvec_new(iv, iv_n);
  return id;
}


MPI_Fint
xt_ut_init_oneway_transposition_template(int id_in, int id_out,
                                         int XT_UNUSED(mpi_world),
                                         int XT_UNUSED(icheck_unique)) {
  xt_ut_deco_type *deco_in = xt_handle2pointer(deco_set, id_in);
  if(!deco_in) die("xt_ut_init_oneway_transposition_template: invalid handle (first argument)");

  xt_ut_deco_type *deco_out = xt_handle2pointer(deco_set, id_out);
  if(!deco_out) die("xt_ut_init_oneway_transposition_template: invalid handle (second argument)");

  int id = new_itemplate();
  xt_ut_template_type *template = xt_handle2pointer(template_set, id);
  if (!template) die("xt_ut_init_oneway_transposition_template: internal error");
  template->xmap = xt_xmap_all2all_new(deco_in->idxvec, deco_out->idxvec, MPI_COMM_WORLD);

  return id;
}


MPI_Fint xt_ut_init_transposition_simple(MPI_Fint itemplate,
                                         MPI_Fint f_datatype) {
  xt_ut_template_type *template = xt_handle2pointer(template_set, itemplate);
  if (!template) die("xt_ut_init_transposition_simple: invalid handle");

  MPI_Datatype datatype = MPI_Type_f2c(f_datatype);
  MPI_Fint itrans = new_itrans();
  xt_ut_trans_type *trans = xt_handle2pointer(trans_set, itrans);
  if (!trans) die("xt_ut_init_transposition_simple: internal error");
  trans->redist = xt_redist_p2p_off_new(template->xmap, NULL, NULL, datatype);

  return itrans;
}

MPI_Fint xt_ut_init_transposition(MPI_Fint itemplate,
                                  MPI_Fint offset_in[],
                                  MPI_Fint XT_UNUSED(offset_in_size),
                                  MPI_Fint offset_out[],
                                  MPI_Fint XT_UNUSED(offset_out_size),
                                  MPI_Fint f_datatype) {

  xt_ut_template_type *template = xt_handle2pointer(template_set, itemplate);
  if (!template) die("xt_ut_init_transposition: invalid handle");

  MPI_Fint itrans = new_itrans();
  xt_ut_trans_type *trans = xt_handle2pointer(trans_set, itrans);
  if (!trans) die("xt_ut_init_transposition: internal error");
  MPI_Datatype datatype = MPI_Type_f2c(f_datatype);
  trans->redist = xt_redist_p2p_off_new(template->xmap, offset_in, offset_out, datatype);

  return itrans;
}


void xt_ut_transpose(const void **pt_in, int itrans, int XT_UNUSED(direction),
                     void **pt_out) {
  xt_ut_trans_type *trans = xt_handle2pointer(trans_set, itrans);
  if (!trans) die("xt_ut_transpose1: invalid handle");
  xt_redist_s_exchange(trans->redist, 1, pt_in, pt_out);
}
