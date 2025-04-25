/*
 * This file contains a set of routines to provide CRAY to IEEE binary and 
 * vice versa conversion for f90 on Unix machines. Unfortunatelly there is 
 * no implemented set of such functions in general available. To restrict
 * the f90 source code to f90 standard only this implementations are developed.
 * One basic convention is to name all routines with a leading util_ to make
 * this visible to the developer.
 *
 * Authors: L. Kornblueh      Max-Planck-Institute for Meteorology, Hamburg
 *          U. Schulzweida    Max-Planck-Institute for Meteorology, Hamburg
 *
 * Date:    7.5.1998
 *
 * $Id: util_convert.c,v 1.3 1999/05/20 08:19:00 m214003 Exp $
 *
 */

#include <string.h>
#include <limits.h>

#include "cfortran.h"

/****************************************************************************/
/* cfortran prototypes:                                                     */

/*
  PDOUBLE on old CRAY PVPs has been 128bit, but that is not true anymore
*/
/* #if  defined  (CRAY) */
/* #  define  PDOUBLE  PFLOAT */
/* #endif */

void cf_util_cray2ieee(double *dcrayf, double *dieeef, int nf);
FCALLSCSUB3(cf_util_cray2ieee, UTIL_CRAY2IEEE, util_cray2ieee,
	    PDOUBLE, PDOUBLE, INT)

void cf_util_ieee2cray(double *dieeef, double *dcrayf, int nf);
FCALLSCSUB3(cf_util_ieee2cray, UTIL_IEEE2CRAY, util_ieee2cray,
	    PDOUBLE, PDOUBLE, INT)

void cf_util_i8toi4(double *dint64, int *int32, int ni);
FCALLSCSUB3(cf_util_i8toi4, UTIL_I8TOI4, util_i8toi4,
	    PDOUBLE, PINT, INT)

void cf_util_i4toi8(int *int32, double *dint64, int ni);
FCALLSCSUB3(cf_util_i4toi8, UTIL_I4TOI8, util_i4toi8,
	    PINT, PDOUBLE, INT)



#if LONG_MAX >= 9223372036854775807L
#  define INT64  long
#  define UINT64 unsigned long
#else
#  define INT64  long long
#  define UINT64 unsigned long long
#endif


/****************************************************************************/

/*
 *  Convert Floating point, Cray to IEEE 64-bit
 *
 *
 *	input : cf	Cray Floating point numbers
 *		nf	Number of elements in cf
 *	output: ieeef	IEEE Floating point numbers
 *
 *  Format :
 *            sign  exponent  mantissa
 *  IEEE :     1      11        52
 *  Cray :     1      15        48
 */

void cf_util_cray2ieee(double *dcrayf, double *dieeef, int nf)
{
  UINT64 *crayf = (UINT64 *) dcrayf;
  UINT64 *ieeef = (UINT64 *) dieeef;

  const UINT64 cray_expbias = 0x4000;   /* Cray exponent bias */
  const UINT64 ieee_expbias = 0x3ff;    /* IEEE exponent bias */
#if LONG_MAX >= 9223372036854775807L
  const UINT64 mask1 = 0x8000000000000000;
#else
  const UINT64 mask1 = 0x8000000000000000LL;
#endif
  const UINT64 mask2 = 0x7fff;
#if LONG_MAX >= 9223372036854775807L
  const UINT64 mask3 = 0xfffffffffffff;
#else
  const UINT64 mask3 = 0xfffffffffffffLL;
#endif
  const UINT64 mask4 = 0x7ff;
#if LONG_MAX >= 9223372036854775807L
  const UINT64 indef = 0xfff8000000000000;   /* IEEE indefinite, 64 bit form */
#else
  const UINT64 indef = 0xfff8000000000000LL; /* IEEE indefinite, 64 bit form */
#endif
  int i;

  /* Set sign bit, exponent and mantissa in one vector loop */

  for (i = 0; i < nf; i++) {
    *(ieeef+i) = (*(crayf+i) & mask1)
       |
        ((((*(crayf+i) >> 48) & mask2)-cray_expbias+ieee_expbias-1) << 52)
       |
    	((*(crayf+i) << 5) & mask3);
  }
  /*
   * handle 0.0, underflow and overflow :
   *
   * if the exponent field goes negative then the Cray number was
   * either 0.0 or too small to represent in IEEE, in either case
   * set the IEEE result to all 0's which is the IEEE representation for 0.0.
   *
   * if the exponent field is too large for the IEEE field, set
   * the result to the indefinite value.
   */

  /* Cray's internal conversion seems to skip this underflow part and instead
     of indef sets 0 .... */

  for (i = 0; i < nf; i++) {
    if ((((*(crayf+i) >> 48) & mask2)-cray_expbias+ieee_expbias-1) < 0) {
      *(ieeef+i) = 0;
    } else if ((((*(crayf+i) >> 48) & mask2)-cray_expbias+ieee_expbias-1) 
	       > mask4) {
/*    *(ieeef+i) = indef;  */
      *(ieeef+i) = 0;
    }
  }


  return;
}

/*
 *  Convert Floating point, IEEE 64-bit, to Cray Floating point
 *
 *	input : ieeef	IEEE Floating point numbers (double precision)
 *		nf	Number of elements in ieeef
 *	output: crayf	Cray Floating point numbers
 *
 *  Format :
 *            sign  exponent  mantissa	unused
 *  IEEE :     1      11        52
 *  Cray :     1      15        48
 */

void cf_util_ieee2cray(double *dieeef, double *dcrayf, int nf)
{
  UINT64 *crayf = (UINT64 *) dcrayf;
  UINT64 *ieeef = (UINT64 *) dieeef;

  const UINT64 cray_expbias = 0x4000;      /* Cray exponent bias */
  const UINT64 ieee_expbias = 0x3ff;       /* IEEE exponent bias */
#if LONG_MAX >= 9223372036854775807L
  const UINT64 implied = 0x10000000000000; /* implied bit in IEEE mantissa */
  const UINT64 mask1 = 0x8000000000000000;
  const UINT64 mask3 = 0xfffffffffffff;
#else
  const UINT64 implied = 0x10000000000000LL; /* implied bit in IEEE mantissa */
  const UINT64 mask1 = 0x8000000000000000LL;
  const UINT64 mask3 = 0xfffffffffffffLL;
#endif
  const UINT64 mask4 = 0x7ff;

  int i;

  /* Set sign bit, exponent and mantissa in one vector loop */

  for (i = 0; i < nf; i++) {
    if (*(ieeef+i) == 0) {
      *(crayf+i) = 0;
    } else {
      *(crayf+i) = (*(ieeef+i) & mask1)
	 |
 	  ((((*(ieeef+i) >> 52) & mask4)-ieee_expbias+cray_expbias+1) << 48)
 	 |
 	  (((*(ieeef+i) & mask3)+implied) >> 5);
    }
  }

  return;
}

void cf_util_i8toi4(double *dint64, int *int32, int ni)
{
  INT64 *int64 = (INT64 *) dint64;
#if SIZEOF_INT == 8
  memcpy (int32, int64, (size_t) ni*8);
#else
  int i;

  for (i = 0; i < ni; i++) {
#ifdef DEBUG
    if ( (*(int64+i) > INT_MAX) ||  (*(int64+i) < INT_MIN)) {
      fprintf (stderr, 
	       "Integer out of range. Converted number contains nonsens\n");
    }
#endif
    /*
    *(int32+i) = (INT) *(int64+i);
    */
    *(int32+i) = (int) *(int64+i);
  }
#endif
  return;
}

void cf_util_i4toi8(int *int32, double *dint64, int ni)
{
  UINT64 *int64 = (UINT64 *) dint64;
#if SIZEOF_INT == 8
  memcpy (int64, int32, (size_t) ni*8);
#else

#if LONG_MAX >= 9223372036854775807L
  const UINT64 mask1 = 0xffffffff00000000;
#else
  const UINT64 mask1 = 0xffffffff00000000LL;
#endif
  const UINT64 mask2 = 0x0;
  int i;

  for (i = 0; i < ni; i++) {
    if ( *(int32+i) < 0 ) {
      *(int64+i) =  *(int32+i) | mask1;
    } else {
      *(int64+i) =  *(int32+i) | mask2;
    }
  }
#endif
  return;
}

