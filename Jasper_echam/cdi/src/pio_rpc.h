#ifndef PIO_RPC_H
#define PIO_RPC_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdbool.h>

#include <mpi.h>
#include <yaxt.h>

enum collectorCommandTags {
  FINALIZE,
  RESOURCES,
  STREAMOPEN,
  STREAMCLOSE,
  STREAMDEFVLIST,
  WRITETS,
  BLOCK_XFER,
};

#define MAXWINBUFFERSIZE ((size_t)2048 * 1024 * 1024)

enum
{
  numRPCFuncs = 1,
  STREAMDEFTIMESTEP = -1,
  HEADERSIZEMARKER = -numRPCFuncs - 1,
  PARTDESCMARKER = -numRPCFuncs - 2,
};
enum {
  MINFUNCID = -numRPCFuncs,
  MAXFUNCID = -1,
  defVlistNInts = 2,
};
extern const char * const funcMap[numRPCFuncs];

struct headerSize
{
  int numDataEntries, numRPCEntries;
};

struct dataRecord
{
  int varID, nmiss;
};

union funcArgs
{
  struct
  {
    int streamID, tsID;
  } streamNewTimestep;
};

/* Describes offset and ID of serialized partition descriptor.
 * partDescMarker == PARTDESCMARKER, always. */
struct partDescRecord
{
  Xt_uid uid;
};

struct winHeaderEntry
{
  int id;
  union
  {
    struct headerSize headerSize;
    struct dataRecord dataRecord;
    union funcArgs funcArgs;
    struct partDescRecord partDesc;
  }  specific;
  int offset;
};

/* round size to next multiple of factor */
static inline size_t
roundUpToMultiple(size_t size, size_t factor)
{
  return (size + factor - 1)/factor * factor;
}

enum
{
  /* align window base addresses and sizes to this value */
  PIO_WIN_ALIGN = sizeof (double),
};

struct clientBufSize
{
  size_t bufSize;
  int numDataRecords, numRPCRecords;
};

struct collSpec
{
  int numClients;
  int numServers;
  bool sendRPCData;
};

struct clientBufSize
computeClientStreamBufSize(int streamID, const struct collSpec *collector);

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
