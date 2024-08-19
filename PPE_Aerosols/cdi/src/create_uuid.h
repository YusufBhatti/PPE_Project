#ifndef CREATE_UUID_H
#define CREATE_UUID_H

#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include "cdi.h"

void
create_uuid(unsigned char uuid[CDI_UUID_SIZE]);

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
