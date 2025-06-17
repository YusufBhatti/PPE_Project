#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <math.h>
#include <stdbool.h>

#include "cdi.h"
#include "resource_handle.h"
/* FIXME: no longer needed when taxis updates are sent as separate data */
#include "taxis.h"

#include "pio_interface.h"
#include "pio_rpc.h"



struct clientBufSize
computeClientStreamBufSize(int streamID, const struct collSpec *collector)
{
  /* 1 record is filled in last to indicate number of records in total */
  struct clientBufSize rmaSizeSpec
    = { .bufSize = sizeof (struct winHeaderEntry), .numDataRecords = 1, .numRPCRecords = 0 };
  int vlistID = streamInqVlist(streamID);
  size_t nvars = (size_t)vlistNvars(vlistID);
  for (size_t varID = 0; varID < nvars; ++varID)
    {
      size_t chunkSize;
      {
        int varSize = vlistInqVarSize(vlistID, (int)varID);
        chunkSize = (size_t)ceilf(cdiPIOpartInflate_
                                    * (float)varSize
                                    / (float)collector->numClients);
      }
      rmaSizeSpec.numDataRecords += 2;
      rmaSizeSpec.bufSize += chunkSize * sizeof (double)
        /* re-align chunk to multiple of double size */
        + sizeof (double) - 1
        /* one header for data record, one for corresponding part
         * descriptor*/
        + 2 * sizeof (struct winHeaderEntry)
        /* FIXME: heuristic for size of packed Xt_idxlist */
        + sizeof (Xt_int) * chunkSize * 3;
    }

  // memory required for the function calls encoded
  // for remote execution
  // once per stream and timestep for each collector process
  // from one model process
  if (collector->sendRPCData)
    {
      rmaSizeSpec.numRPCRecords = numRPCFuncs;
      rmaSizeSpec.bufSize +=
        numRPCFuncs * sizeof (struct winHeaderEntry)
        /* data part of streamDefTimestep */
        + (2 * CDI_MAX_NAME + sizeof (taxis_t));
    }
  rmaSizeSpec.bufSize = roundUpToMultiple(rmaSizeSpec.bufSize, PIO_WIN_ALIGN);
  return rmaSizeSpec;
}
