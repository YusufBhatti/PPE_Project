#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdi_int.h"
#include "pio.h"
#include "cdipio.h"
#include "pio_impl.h"

const char *const cdiPioCmdStrTab[] = {
  "IO_Open_file",
  "IO_Close_file",
  "IO_Get_fp",
  "IO_Set_fp",
  "IO_Send_buffer",
  "IO_Finalize"
};

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
