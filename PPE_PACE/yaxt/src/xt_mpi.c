/**
 * @file xt_mpi.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *                                 Joerg Behrens <behrens@dkrz.de>
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 *         Joerg Behrens <behrens@dkrz.de>
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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <inttypes.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_core.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"

//! COMPACT_DT enables the anlysis of displacements in order to give a
//! more compact description to the datatype generators of MPI. For strong
//! enough MPI implementations this not required. Then you can undefine
//! COMPACT_DT and save some prcessing time within yaxt without losing communication
//! performance.
#define COMPACT_DT

#define zero_stripe ((struct Xt_offset_ext){ .start=0, .stride=0, .size=0 })

static MPI_Datatype
xt_mpi_generate_compact_datatype_block(const int *disp, const int *blocklengths,
                                       int count, MPI_Datatype old_type);
static MPI_Datatype
xt_mpi_generate_compact_datatype(int const *disp, int disp_len,
                                 MPI_Datatype old_type);


//taken from http://beige.ucs.indiana.edu/I590/node85.html
void xt_mpi_error(int error_code, MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  char error_string[MPI_MAX_ERROR_STRING];
  int length_of_error_string, error_class;

  MPI_Error_class(error_code, &error_class);
  MPI_Error_string(error_class, error_string, &length_of_error_string);
  fprintf(stderr, "%3d: %s\n", rank, error_string);
  MPI_Error_string(error_code, error_string, &length_of_error_string);
  fprintf(stderr, "%3d: %s\n", rank, error_string);
  MPI_Abort(comm, error_code);
}

#ifndef COMPACT_DT
static MPI_Datatype copy_mpi_datatype(MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_dup(old_type, &datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_simple(int displacement, MPI_Datatype old_type, MPI_Comm comm)
{
  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_create_indexed_block(1, 1, &displacement, old_type,
                                                &datatype), comm);
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_contiguous(int displacement, int blocklength,
                            MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  if (displacement == 0)
    xt_mpi_call(MPI_Type_contiguous(blocklength, old_type, &datatype),
                    comm);
  else
    xt_mpi_call(MPI_Type_create_indexed_block(1, blocklength,
                                                  &displacement, old_type,
                                                  &datatype), comm);

  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;

}

static MPI_Datatype
gen_mpi_datatype_vector(int count, int blocklength, int stride,
                        int offset, MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_vector(count, blocklength, stride, old_type,
                              &datatype), comm);
  if (offset != 0) {

    MPI_Datatype datatype_;
    int hindexed_blocklength = 1;
    MPI_Aint old_type_size, old_type_lb;

    xt_mpi_call(MPI_Type_get_extent(old_type, &old_type_lb,
                                    &old_type_size), comm);

    MPI_Aint displacement = offset * old_type_size;

    xt_mpi_call(MPI_Type_create_hindexed(1, &hindexed_blocklength,
                                         &displacement, datatype, &datatype_),
                comm);
    xt_mpi_call(MPI_Type_free(&datatype), comm);
    datatype = datatype_;
  }
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_indexed_block(int const * displacements, int blocklength,
                               int count, MPI_Datatype old_type, MPI_Comm comm)
{
  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_create_indexed_block(count, blocklength,
                                                (void *)displacements,
                                                old_type, &datatype), comm);
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_indexed(const int *displacements, const int *blocklengths,
                         int count, MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_indexed(count, (int*)blocklengths, (void*)displacements,
                                   old_type, &datatype), comm);
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static inline int
check_for_vector_type(const int *displacements, const int *blocklengths,
                      int count) {

  int blocklength = blocklengths[0];

  for (int i = 1; i < count; ++i)
    if (blocklengths[i] != blocklength)
      return 0;

  int stride = displacements[1] - displacements[0];

  for (int i = 1; i + 1 < count; ++i)
    if (displacements[i+1] - displacements[i] != stride)
      return 0;

  return 1;
}

static inline int check_for_indexed_block_type(const int *blocklengths,
                                               int count) {

  int blocklength = blocklengths[0];

  for (int i = 1; i < count; ++i)
    if (blocklengths[i] != blocklength)
      return 0;

  return 1;
}
#endif

MPI_Datatype
xt_mpi_generate_datatype_block(const int *displacements,
                               const int *blocklengths,
                               int count, MPI_Datatype old_type,
                               MPI_Comm comm) {
#ifdef COMPACT_DT
  (void)comm;
  return xt_mpi_generate_compact_datatype_block(displacements, blocklengths, count, old_type);
#else
  MPI_Datatype datatype;

  if (count == 0)
    datatype = MPI_DATATYPE_NULL;
  else if (count == 1 && blocklengths[0] == 1 && displacements[0] == 0)
    datatype = copy_mpi_datatype(old_type, comm);
  else if (count == 1 && blocklengths[0] == 1)
    datatype = gen_mpi_datatype_simple(displacements[0], old_type, comm);
  else if (count == 1)
    datatype = gen_mpi_datatype_contiguous(displacements[0], blocklengths[0],
                                           old_type, comm);
  else if (check_for_vector_type(displacements, blocklengths, count))
    datatype = gen_mpi_datatype_vector(count, blocklengths[0],
                                       displacements[1] - displacements[0],
                                       displacements[0], old_type, comm);
  else if (check_for_indexed_block_type(blocklengths, count))
    datatype = gen_mpi_datatype_indexed_block(displacements, blocklengths[0],
                                              count, old_type, comm);
  else
    datatype = gen_mpi_datatype_indexed(displacements, blocklengths, count,
                                        old_type, comm);

  return datatype;
#endif
}

MPI_Datatype xt_mpi_generate_datatype(int const * displacements, int count,
                                      MPI_Datatype old_type, MPI_Comm comm) {

  if (count <= 0)
    return MPI_DATATYPE_NULL;

#ifdef COMPACT_DT
  return xt_mpi_generate_compact_datatype(displacements, count, old_type);
#endif

  int * blocklengths = xmalloc((size_t)count * sizeof(*blocklengths));
  int new_count = 0;
  {
    int i = 0;
    do {
      int j = 1;
      while (i + j < count && displacements[i] + j == displacements[i + j])
        ++j;
      blocklengths[new_count++] = j;
      i += j;
    } while (i < count);
  }

  int * tmp_displ = NULL;
  const int *displ;

  if (new_count != count) {

    tmp_displ = xmalloc((size_t)new_count * sizeof(*tmp_displ));

    int offset = 0;

    for (int i = 0; i < new_count; ++i) {

      tmp_displ[i] = displacements[offset];
      offset += blocklengths[i];
    }

    displ = tmp_displ;
  } else
    displ = displacements;

  MPI_Datatype datatype;

  datatype = xt_mpi_generate_datatype_block(displ, blocklengths, new_count,
                                            old_type, comm);

  free(blocklengths);

  free(tmp_displ);

  return datatype;
}


static int
scan_stripe(const int *disp, int disp_len, struct Xt_offset_ext *restrict v,
            size_t vsize) {

  if (disp_len<1) return 0;

  struct Xt_offset_ext x = zero_stripe;
  size_t i = 0;
  int p = 0;
  while (p < disp_len) {

    if (!x.size) {
      x.start   = disp[p];
      x.stride = 1;
      x.size = 1;
      p++; continue;
    }

    if (x.size == 1) {
      x.stride = disp[p] - disp[p-1];
      x.size = 2;
      p++; continue;
    }

    // x.size >= 2:
    if (disp[p] - disp[p-1] == x.stride) {
      x.size++;
      p++; continue;
    }

    if (x.size > 2 || (x.size == 2 && x.stride == 1) ) {
      // we accept small contiguous vectors (nstrides==2, stride==1)
      if (i >= vsize) die("scan_stripe: vsize too small\n");
      v[i]= x;
      i++;
      x = zero_stripe;
      continue;
    }

    if (x.size == 2) {
      // break up trivial vec:
      if (i >= vsize) die("scan_stripe: vsize too small\n");
      v[i] = x;
      v[i].size = 1;
      v[i].stride = 1;
      i++;
      x.start += x.stride;
      x.size = 1;
      x.stride = 1;
      continue;
    }

  }

  // tail cases:
  if (x.size > 2 || (x.size == 2 && x.stride == 1) ) {
    if (i >= vsize) die("scan_stripe: vsize too small\n");
      v[i]= x;
      i++;
  } else if (x.size == 2) {
    if (i+1 >= vsize) die("scan_stripe: vsize too small\n");
      v[i] = x;
      v[i].size = 1;
      v[i].stride = 1;
      i++;
      v[i] = x;
      v[i].start += x.stride;
      v[i].size = 1;
      v[i].stride = 1;
      i++;
  } else if (x.size == 1) {
    if (i >= vsize) die("scan_stripe: vsize too small\n");
    v[i] = x;
    v[i].size = 1;
    v[i].stride = 1;
    i++;
  }

  size_t vn = i;

  // check:
  p = 0;
  for (size_t j = 0; j < vn; j++) {
    for (int k = 0; k < v[j].size; k++) {
      int d = v[j].start + k*v[j].stride;
      static const char die_msg[2][32] = { "scan_stripe: internal error (1)",
                                           "scan_stripe: internal error (2)" };
      int err_code = p > disp_len || (disp[p] != d) * 2;
      if (err_code) die(die_msg[err_code - 1]);
      p++;
    }
  }
  if (p != disp_len) die("scan_stripe: internal error (3)");

  return (int)vn;
}

static int
match_simple_vec(int *pstart, const struct Xt_offset_ext *v, int vlen,
                 MPI_Datatype old_type, int *disp, MPI_Datatype *dt) {
  // we only accept non-trivial matches (nsteps>2) with stride /= 1
  // using only one vector from v
  int p = *pstart;
  if (p >= vlen) return 0;
  int nstrides = v[p].size;
  int stride = v[p].stride;
  if (nstrides < 2 || stride == 1 ) return 0;

  (*pstart)++;

  *disp = vlen > 1 ? v[p].start : 0;

  MPI_Datatype dt1;
  xt_mpi_call(MPI_Type_vector(nstrides, 1, stride, old_type, &dt1),
              Xt_default_comm);

  int start = v[p].start - *disp;
  if (!start) {
    *dt = dt1;
    return nstrides;
  }

  // (start != 0) => add offset:
  MPI_Aint old_type_size, old_type_lb;
  xt_mpi_call(MPI_Type_get_extent(old_type, &old_type_lb,
                                  &old_type_size), Xt_default_comm);

  MPI_Aint displacement = start * old_type_size;
  int bl2 = 1;
  MPI_Datatype dt2;
  xt_mpi_call(MPI_Type_create_hindexed(1, &bl2, &displacement, dt1, &dt2),
              Xt_default_comm);

  xt_mpi_call(MPI_Type_free(&dt1), Xt_default_comm);

  *dt = dt2;
  return nstrides;
}

static int
match_block_vec(int *pstart, const struct Xt_offset_ext *v, int vlen,
                MPI_Datatype old_type, int *disp, MPI_Datatype *dt) {
  // using at least 3 vectors
  int p = *pstart;
  if (p+2 >= vlen) return 0;
  if (v[p].stride != 1 || v[p+1].stride != 1 ) return 0;
  int bl = v[p].size;
  if (bl < 1 || v[p+1].size != bl) return 0;

  int vstride = v[p+1].start - v[p].start;

  p += 2;
  while( p < vlen && v[p].stride == 1 && v[p].size == bl &&
         v[p].start - v[p-1].start == vstride ) {
    p++;
  }
  int n = p - *pstart;
  if (n<3) return 0;

  if (n == vlen) {
    *disp = 0;
  } else {
    *disp = v[*pstart].start;
  }

  MPI_Datatype dt1;
  xt_mpi_call(MPI_Type_vector(n, bl, vstride, old_type, &dt1),
              Xt_default_comm);

  int start = v[*pstart].start - *disp;

  *pstart = p;
  if (!start) {
    *dt = dt1;
    return n;
  }

  // (start != 0) => add offset:
  MPI_Aint old_type_size, old_type_lb;
  xt_mpi_call(MPI_Type_get_extent(old_type, &old_type_lb,
                                  &old_type_size), Xt_default_comm);

  MPI_Aint displacement = start * old_type_size;
  int bl2 = 1;
  MPI_Datatype dt2;
  xt_mpi_call(MPI_Type_create_hindexed(1, &bl2, &displacement, dt1, &dt2),
              Xt_default_comm);

  xt_mpi_call(MPI_Type_free(&dt1), Xt_default_comm);

  *dt = dt2;
  return n;
}

static int
match_contiguous(int *pstart, const struct Xt_offset_ext *v, int vlen,
                 MPI_Datatype old_type, int *disp, MPI_Datatype *dt) {
  int p = *pstart;
  if (p >= vlen) return 0;
  if (v[p].stride != 1) return 0;
  if (v[p].size < 2) return 0;

  if (vlen > 1) {
    *disp = v[*pstart].start;
  } else {
    *disp = 0;
  }

  int d = v[p].start - *disp;

  if (!d) {
    xt_mpi_call(MPI_Type_contiguous(v[p].size, old_type, dt), Xt_default_comm) ;
  } else {
    xt_mpi_call(MPI_Type_create_indexed_block(1, v[p].size, &d, old_type, dt),
                Xt_default_comm);
  }

  (*pstart)++;
  return v[p].size;
}

static int
match_indexed(int *pstart, const struct Xt_offset_ext *v, int vlen, MPI_Datatype old_type, int *disp, MPI_Datatype *dt) {
  // we only accept non-trivial matches
  int p = *pstart;
  if (p >= vlen) return 0;
  if (v[p].stride != 1) return 0;
  if (v[p].size < 2) return 0;

  while( p < vlen && v[p].stride == 1 ) {
    p++;
  }
  int n = p - *pstart;

  if (n < 2) return 0;

  int start;
  if (n == vlen) {
    start = 0;
  } else {
    start = v[*pstart].start;
  }

  *disp = start;
  int *restrict bl = xmalloc(2 * (size_t)n * sizeof (*bl)),
    *restrict d = bl + n;
  int hom_bl = 1;
  for (int i = 0; i < n; i++) {
    int iv = *pstart+i;
    d[i] = v[iv].start - start;
    bl[i] = v[iv].size;
    if (bl[i] != bl[0]) hom_bl = 0;
  }

  if (hom_bl) {
    xt_mpi_call(MPI_Type_create_indexed_block(n, bl[0], d, old_type, dt), Xt_default_comm);
  } else {
    xt_mpi_call(MPI_Type_indexed(n, bl, d, old_type, dt),
                Xt_default_comm);
  }

  *pstart = p;

  free(bl);
  return n;
}

static int
gen_fallback_type(int *istart, int *iend, const struct Xt_offset_ext *v,
                  int vlen, MPI_Datatype old_type, int *offset,
                  MPI_Datatype *dt) {
  if (*istart<0 || *iend < *istart || *iend >= vlen) return 0;

  int ia = *istart;
  int ib = *iend;

  int n = 0;
  for (int i=ia; i<=ib; i++) {
    n += v[i].size;
  }
  if (n<1) return 0;

  int start;
  if (ia == 0 && ib == vlen-1) {
    // generate absolute datatype
    start = 0;
  } else {
    // generate relative datatype that gets embedded by the caller
    start = v[ia].start;
  }

  *offset = start;

  int *d = xmalloc(sizeof (*d) * (size_t)n);
  int p=0;

  for (int i=ia; i<=ib; i++) {
    for (int k=0; k<v[i].size; k++) {
      d[p] = v[i].start + k * v[i].stride - start;
      p++;
    }
  }

  if (n==1 && d[0] == 0) {
    // At the moment we disable the embed aspect that was used to avoid MPI_Type_dup.
    // It turned out this creates more complexity when freeing the intermediate datatypes in the callin function.
    // maybe we just forget the whole idea - then we can also remove the embed logic here.
    //if (embed)
    //  *dt = old_type;
    //else
    xt_mpi_call(MPI_Type_dup(old_type, dt), Xt_default_comm);
  } else {
    xt_mpi_call(MPI_Type_create_indexed_block(n, 1, d, old_type, dt), Xt_default_comm);
  }
  free(d);

  *istart = *iend + 1;

  return n;
}

static MPI_Datatype
parse_stripe(const struct Xt_offset_ext *v, int vlen, MPI_Datatype old_type) {

  int set_start = -1;
  int set_end = -2;
  MPI_Datatype *restrict wdt = xmalloc(sizeof(*wdt) * (size_t)vlen);
  int *restrict wdisp = xmalloc(sizeof (*wdisp) * (size_t)vlen);
  MPI_Datatype match_dt;
  int match_disp;
  int p = 0;
  int m = 0;
  while (p<vlen) {
    if ( match_block_vec(&p, v, vlen, old_type, &match_disp, &match_dt) ||
         match_indexed(&p, v, vlen, old_type, &match_disp, &match_dt) ||
         match_simple_vec(&p, v, vlen, old_type, &match_disp, &match_dt) ||
         match_contiguous(&p, v, vlen, old_type, &match_disp, &match_dt) ) {
      if (set_start <= set_end) {
        gen_fallback_type(&set_start, &set_end, v, vlen, old_type, &wdisp[m], &wdt[m]);
        m++;
      }
      wdt[m] = match_dt;
      wdisp[m] = match_disp;
      m++;
    } else {
      if (set_start>set_end) set_start = p;
      set_end = p;
      p++;
    }

  }
  if (set_start <=  set_end) {
    gen_fallback_type(&set_start, &set_end, v, vlen, old_type, &wdisp[m], &wdt[m]);
    m++;
  }
  int wlen = m;
  MPI_Datatype result_dt;
  if (wlen == 1 ) {
    if (wdisp[0] == 0) {
      xt_mpi_call(MPI_Type_dup(wdt[0], &result_dt), Xt_default_comm);
    } else {
      die("parse_stripe: internal error; wlen == 1 && match_disp != 0\n");
    }
  } else {
    MPI_Aint old_type_lb, old_type_extent;
    MPI_Aint *restrict wbdisp = xmalloc((size_t)wlen * sizeof (*wbdisp));
    int *restrict wblocklength
      = xmalloc((size_t)wlen * sizeof (*wblocklength));;
    xt_mpi_call(MPI_Type_get_extent(old_type, &old_type_lb,
                                    &old_type_extent), Xt_default_comm);
    for(int i=0; i<wlen; i++) {
      wbdisp[i] = wdisp[i] * old_type_extent;
      wblocklength[i] = 1;
    }
    xt_mpi_call(MPI_Type_create_struct(wlen, wblocklength, wbdisp,
                                       wdt, &result_dt), Xt_default_comm);
    free(wblocklength);
    free(wbdisp);
  }
  xt_mpi_call(MPI_Type_commit(&result_dt), Xt_default_comm);
  for(m = 0; m<wlen; m++) {
    xt_mpi_call(MPI_Type_free(&wdt[m]), Xt_default_comm);
  }
  free(wdt);
  free(wdisp);
  return result_dt;
}

MPI_Datatype
xt_mpi_generate_datatype_stripe(const struct Xt_offset_ext *v,
                                int count, MPI_Datatype old_type,
                                MPI_Comm XT_UNUSED(comm))
{
  if (count < 1) return MPI_DATATYPE_NULL;
  return parse_stripe(v, count, old_type);
}


static MPI_Datatype
xt_mpi_generate_compact_datatype_block(const int *disp, const int *blocklengths,
                                       int count, MPI_Datatype old_type) {

  if (count < 1) return MPI_DATATYPE_NULL;
  struct Xt_offset_ext *restrict v = xmalloc(sizeof(*v) * (size_t)count);
  for (int i=0; i<count; i++) {
    v[i].start = disp[i];
    v[i].stride = 1;
    v[i].size = blocklengths[i];
  }
  MPI_Datatype dt = parse_stripe(v, count, old_type);
  free(v);
  return dt;
}

static MPI_Datatype
xt_mpi_generate_compact_datatype(int const *disp, int disp_len,
                                 MPI_Datatype old_type) {

  if (disp_len < 1) return MPI_DATATYPE_NULL;

  struct Xt_offset_ext *v = xmalloc(sizeof(*v) * (size_t)disp_len);
  int vlen = scan_stripe(disp, disp_len, v, (size_t)disp_len);
  MPI_Datatype dt = parse_stripe(v, vlen, old_type);
  free(v);
  return dt;
}

/* functions to handle optimizations on communicators */
int xt_mpi_comm_internal_keyval;

typedef unsigned long used_map_elem;

enum {
  used_map_elem_bits = sizeof (used_map_elem) * CHAR_BIT,
};

struct xt_mpi_comm_internal_attr {
  int refcount;
  unsigned used_map_size;
  used_map_elem used_map[];
};

static int
xt_mpi_comm_internal_keyval_copy(
  MPI_Comm XT_UNUSED(oldcomm), int XT_UNUSED(keyval),
  void *XT_UNUSED(extra_state), void *XT_UNUSED(attribute_val_in),
  void *attribute_val_out, int *flag)
{
  struct xt_mpi_comm_internal_attr *new_comm_attr
    = malloc(sizeof (struct xt_mpi_comm_internal_attr)
             + sizeof (used_map_elem));
  int retval;
  if (new_comm_attr)
  {
    new_comm_attr->refcount = 1;
    new_comm_attr->used_map_size = 1;
    new_comm_attr->used_map[0] = 1U;
    *(void **)attribute_val_out = new_comm_attr;
    *flag = 1;
    retval = MPI_SUCCESS;
  } else {
    *flag = 0;
    retval = MPI_ERR_NO_MEM;
  }
  return retval;
}

static int
xt_mpi_comm_internal_keyval_delete(
  MPI_Comm XT_UNUSED(comm), int XT_UNUSED(comm_keyval),
  void *attribute_val, void *XT_UNUSED(extra_state))
{
  free(attribute_val);
  return MPI_SUCCESS;
}

static int xt_mpi_tag_ub_val;

void
xt_mpi_init(void) {
  xt_mpi_call(MPI_Comm_create_keyval(xt_mpi_comm_internal_keyval_copy,
                                     xt_mpi_comm_internal_keyval_delete,
                                     &xt_mpi_comm_internal_keyval, NULL),
              Xt_default_comm);
  void *attr;
  int flag;
  xt_mpi_call(MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &attr, &flag),
              MPI_COMM_WORLD);
  assert(flag);
  xt_mpi_tag_ub_val = *(int *)attr;
}

void
xt_mpi_finalize(void) {
  xt_mpi_call(MPI_Comm_free_keyval(&xt_mpi_comm_internal_keyval),
              Xt_default_comm);
}

static struct xt_mpi_comm_internal_attr *
xt_mpi_comm_get_internal_attr(MPI_Comm comm)
{
  int attr_found;
  void *attr_val;
  xt_mpi_call(MPI_Comm_get_attr(comm, xt_mpi_comm_internal_keyval,
                                &attr_val, &attr_found),
              comm);
  return attr_found ? attr_val : NULL;
}

#if HAVE_DECL___BUILTIN_CTZL
#define ctzl(v) (__builtin_ctzl(v))
#elif HAVE_DECL___BUILTIN_CLZL
static inline int
ctzl(unsigned long v) {
  enum {
    ulong_bits = sizeof (unsigned long) * CHAR_BIT,
  };
  /* clear all but lowest 1 bit */
  v = v & ~(v - 1);
  int c = ulong_bits - __builtin_clzl(v) - 1;
  return c;
}
#else
static inline int
ctzl(unsigned long v) {
  enum {
    ulong_bits = sizeof (unsigned long) * CHAR_BIT,
  };
  // c will be the number of zero bits on the right
  unsigned int c = ulong_bits;
  v &= (unsigned long)-(long)v;
  if (v) c--;
#if SIZEOF_UNSIGNED_LONG * CHAR_BIT == 64
  if (v & UINT64_C(0x00000000ffffffff)) c -= 32;
  if (v & UINT64_C(0x0000ffff0000ffff)) c -= 16;
  if (v & UINT64_C(0x00ff00ff00ff00ff)) c -= 8;
  if (v & UINT64_C(0x0f0f0f0f0f0f0f0f)) c -= 4;
  if (v & UINT64_C(0x3333333333333333)) c -= 2;
  if (v & UINT64_C(0x5555555555555555)) c -= 1;
#elif SIZEOF_UNSIGNED_LONG * CHAR_BIT == 32
  if (v & 0x0000FFFFUL) c -= 16;
  if (v & 0x00FF00FFUL) c -= 8;
  if (v & 0x0F0F0F0FUL) c -= 4;
  if (v & 0x33333333UL) c -= 2;
  if (v & 0x55555555UL) c -= 1;
#else
  error "Unexpected size of long.\n"
#endif
  return (int)c;
}
#endif

MPI_Comm
xt_mpi_comm_smart_dup(MPI_Comm comm, int *tag_offset)
{
  MPI_Comm comm_dest;
  struct xt_mpi_comm_internal_attr *comm_xt_attr_val
    = xt_mpi_comm_get_internal_attr(comm);
  size_t position = 0;
  int refcount = comm_xt_attr_val ? comm_xt_attr_val->refcount : 0;
  if (comm_xt_attr_val
      && (refcount + 1) < xt_mpi_tag_ub_val / xt_mpi_num_tags) {
    comm_dest = comm;
    comm_xt_attr_val->refcount = ++refcount;
    size_t used_map_size = comm_xt_attr_val->used_map_size;
    while (position < used_map_size
           && comm_xt_attr_val->used_map[position] == ~(used_map_elem)0)
      ++position;
    if (position >= used_map_size) {
      /* sadly, we need to recreate the value to enlarge it */
      struct xt_mpi_comm_internal_attr *new_comm_xt_attr_val
        = xmalloc(sizeof (*new_comm_xt_attr_val)
                  + (used_map_size + 1) * sizeof (used_map_elem));
      new_comm_xt_attr_val->refcount = refcount;
      new_comm_xt_attr_val->used_map_size = (unsigned)(used_map_size + 1);
      for (size_t i = 0; i < used_map_size; ++i)
        new_comm_xt_attr_val->used_map[i] = comm_xt_attr_val->used_map[i];
      new_comm_xt_attr_val->used_map[used_map_size] = 1U;
      position *= used_map_elem_bits;
      xt_mpi_call(MPI_Comm_set_attr(comm_dest, xt_mpi_comm_internal_keyval,
                                    new_comm_xt_attr_val), comm_dest);
    } else {
      /* not all bits are set, find first unset position and insert */
      used_map_elem used_map_entry = comm_xt_attr_val->used_map[position],
        unset_lsb = ~used_map_entry & (used_map_entry + 1),
        bit_pos = (used_map_elem)ctzl(unset_lsb);
      comm_xt_attr_val->used_map[position] = used_map_entry | unset_lsb;
      position = position * used_map_elem_bits + (size_t)bit_pos;
    }
  } else {
    struct xt_mpi_comm_internal_attr *comm_attr
      = xmalloc(sizeof (*comm_attr) + sizeof (used_map_elem));
    comm_attr->refcount = 1;
    comm_attr->used_map_size = 1;
    comm_attr->used_map[0] = 1U;
    xt_mpi_call(MPI_Comm_dup(comm, &comm_dest), comm);
    xt_mpi_call(MPI_Comm_set_attr(comm_dest, xt_mpi_comm_internal_keyval,
                                  comm_attr), comm_dest);
  }
  *tag_offset = (int)(position * xt_mpi_num_tags);
  return comm_dest;
}

void
xt_mpi_comm_smart_dedup(MPI_Comm *comm, int tag_offset)
{
  struct xt_mpi_comm_internal_attr *comm_xt_attr_val
    = xt_mpi_comm_get_internal_attr(*comm);
  int refcount = comm_xt_attr_val ? --(comm_xt_attr_val->refcount) : 0;
  if (refcount < 1) {
    xt_mpi_call(MPI_Comm_free(comm), MPI_COMM_WORLD);
    *comm = MPI_COMM_NULL;
  } else {
    size_t position = (unsigned)tag_offset / xt_mpi_num_tags,
      map_elem = position / used_map_elem_bits,
      in_elem_bit = position % used_map_elem_bits;
    comm_xt_attr_val->used_map[map_elem] &= ~((used_map_elem)1 << in_elem_bit);
  }
}

void
xt_mpi_comm_mark_exclusive(MPI_Comm comm) {
  struct xt_mpi_comm_internal_attr *comm_attr
    = xmalloc(sizeof (*comm_attr) + sizeof (used_map_elem));
  comm_attr->refcount = 1;
  comm_attr->used_map_size = 1;
  comm_attr->used_map[0] = 1U;
  xt_mpi_call(MPI_Comm_set_attr(comm, xt_mpi_comm_internal_keyval,
                                comm_attr), comm);
}

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
