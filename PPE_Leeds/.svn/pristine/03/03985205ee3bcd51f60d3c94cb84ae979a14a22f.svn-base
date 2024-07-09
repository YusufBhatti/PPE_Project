/**
 * @file xt_arithmetic_util.h
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
#ifndef XT_ARITHMETIC_UTIL_H
#define XT_ARITHMETIC_UTIL_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <limits.h>

#include "xt/xt_core.h"

/* simple operations on Xt_int */
/**
 * @return -1 if x < 0, 1 otherwise
 */
static inline Xt_int
Xt_isign(Xt_int x)
{
#if (-1 >> 1) == -1
  return ((x >> (sizeof (Xt_int) * CHAR_BIT - 1)) |
          (Xt_int)((Xt_uint)(~x) >> (sizeof (Xt_int) * CHAR_BIT - 1)));
#else
  return (Xt_int)((x >= 0) - (x < 0));
#endif
}

/**
 * @return -1 if x < 0, 1 otherwise
 */
static inline int
isign(int x)
{
#if (-1 >> 1) == -1
  return ((x >> (sizeof (int) * CHAR_BIT - 1)) |
          (int)((unsigned int)(~x) >> (sizeof (int) * CHAR_BIT - 1)));
#else
  return (x >= 0) - (x < 0);
#endif
}

/**
 * @return ~0 if x < 0, 0 otherwise
 */
static inline int
isign_mask(int x)
{
#if (-1 >> 1) == -1
  return x >> (sizeof (int) * CHAR_BIT - 1);
#else
#warning Unusual behaviour of shift operator detected.
  return (x < 0) * ~0;
#endif
}


/**
 * @return ~(Xt_int)0 if x < 0, 0 otherwise
 */
static inline Xt_int
Xt_isign_mask(Xt_int x)
{
#if (-1 >> 1) == -1
  return x >> (sizeof (Xt_int) * CHAR_BIT - 1);
#else
#warning Unusual behaviour of shift operator detected.
  return (x < 0) * ~0;
#endif
}

/**
 * @return MIN(a, b)
 */
static inline int
imin(int a, int b)
{
  return a <= b ? a : b;
}

#endif
