#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <limits.h>

#include "dmemory.h"

#include "cdi.h"
#include "cdi_int.h"


static
void tstepsInitEntry(stream_t *streamptr, size_t tsID)
{
  streamptr->tsteps[tsID].curRecID     = CDI_UNDEFID;
  streamptr->tsteps[tsID].position     = 0;
  streamptr->tsteps[tsID].records      = NULL;
  streamptr->tsteps[tsID].recordSize   = 0;
  streamptr->tsteps[tsID].nallrecs     = 0;
  streamptr->tsteps[tsID].recIDs       = NULL;
  streamptr->tsteps[tsID].nrecs        = 0;
  streamptr->tsteps[tsID].next         = 0;

  ptaxisInit(&streamptr->tsteps[tsID].taxis);
}


int tstepsNewEntry(stream_t *streamptr)
{
  size_t tsID            = (size_t)streamptr->tstepsNextID++;
  size_t tstepsTableSize = (size_t)streamptr->tstepsTableSize;
  tsteps_t *tstepsTable  = streamptr->tsteps;

  /*
    If the table overflows, double its size.
  */
  if ( tsID == tstepsTableSize )
    {
      if ( tstepsTableSize == 0 ) tstepsTableSize = 1;
      if ( tstepsTableSize <= INT_MAX / 2)
        tstepsTableSize *= 2;
      else if ( tstepsTableSize < INT_MAX)
        tstepsTableSize = INT_MAX;
      else
        Error("Resizing of tstep table failed!");
      tstepsTable = (tsteps_t *) Realloc(tstepsTable,
                                         tstepsTableSize * sizeof (tsteps_t));
    }

  streamptr->tstepsTableSize = (int)tstepsTableSize;
  streamptr->tsteps          = tstepsTable;

  tstepsInitEntry(streamptr, tsID);

  streamptr->tsteps[tsID].taxis.used = TRUE;

  return (int)tsID;
}


void cdiCreateTimesteps(stream_t *streamptr)
{
  long ntsteps;
  long tsID;

  if ( streamptr->ntsteps < 0 || streamptr->tstepsTableSize > 0 )
    return;

  if ( streamptr->ntsteps == 0 ) ntsteps = 1;    /* <<<<<-------- */
  else ntsteps = streamptr->ntsteps;

  streamptr->tsteps = (tsteps_t *) Malloc((size_t)ntsteps*sizeof(tsteps_t));

  streamptr->tstepsTableSize = (int)ntsteps;
  streamptr->tstepsNextID    = (int)ntsteps;

  for ( tsID = 0; tsID < ntsteps; tsID++ )
    {
      tstepsInitEntry(streamptr, (size_t)tsID);
      streamptr->tsteps[tsID].taxis.used = TRUE;
    }
}
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
