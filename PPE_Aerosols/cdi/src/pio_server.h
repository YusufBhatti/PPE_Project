#ifndef PIO_SERVER_
#define PIO_SERVER_

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <mpi.h>

#include "pio_conf.h"

void cdiPioCollectorMessageLoop(const struct cdiPioConf *conf);

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
