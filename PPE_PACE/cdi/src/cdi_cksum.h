#ifndef CDI_CKSUM_H_
#define CDI_CKSUM_H_

#include <inttypes.h>

uint32_t cdiCheckSum(int type, int count, const void *data);

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
