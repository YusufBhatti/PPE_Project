#ifndef PIO_CDF_INT_H
#define PIO_CDF_INT_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#ifdef HAVE_LIBNETCDF
#include "cdf_int.h"

void
cdiPioEnableNetCDFParAccess(void);

#endif
#endif

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
