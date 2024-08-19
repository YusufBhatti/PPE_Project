#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <ctype.h>

#include <yaxt.h>

#include "cdi.h"
#include "cdi_int.h"
#include "dmemory.h"
#include "namespace.h"
#include "taxis.h"

#include "cdipio.h"
#include "pio.h"
#include "pio_client.h"
#include "pio_comm.h"
#include "pio_interface.h"
#include "pio_rpc.h"
#include "pio_util.h"
#include "pio_serialize.h"

static int
cdiPioClientStreamOpen(const char *filename, char filemode,
                       int filetype, stream_t *streamptr,
                       int recordBufIsToBeCreated)
{
  (void)streamptr; (void)recordBufIsToBeCreated;
  int fileID;
  if ( filemode == 'w' )
    {
      MPI_Comm comm = cdiPioInqInterComm();
      int clientRank = commInqRankModel(),
        numClients = cdiPioCommInqSizeClients(),
        numColl = commInqSizeColl(),
        collRank = cdiPioCollRank(clientRank, numClients, numColl);
      int streamID = streamptr->self;
      if (clientRank
          == cdiPioClientRangeStart(collRank, numClients, numColl))
        {
          reshSetStatus(streamID, &streamOps,
                        reshGetStatus(streamID, &streamOps)
                        & ~RESH_SYNC_BIT);
          char *msgBuffer;
          int msgSize = 0;
          int msgBufPos = reshPackBufferCreate(&msgBuffer, &msgSize, &comm);
          int size;
          size_t filename_len = strlen(filename);
          xassert(filename_len < (size_t)(INT_MAX - msgBufPos));
          int soHdr[3] = { streamptr->self, filetype, (int)filename_len };
          xmpi(MPI_Pack_size(3, MPI_INT, comm, &size));
          msgSize += size;
          xmpi(MPI_Pack_size(1, MPI_CHAR, comm, &size));
          msgSize += size;
          xmpi(MPI_Pack_size((int)filename_len, MPI_CHAR, comm, &size));
          msgSize += size;
          /* optimize to pos + size */
          msgBuffer = (char *)Realloc(msgBuffer, (size_t)msgSize);
          xmpi(MPI_Pack(soHdr, 3, MPI_INT,
                        msgBuffer, msgSize, &msgBufPos, comm));
          xmpi(MPI_Pack(&filemode, 1, MPI_CHAR,
                        msgBuffer, msgSize, &msgBufPos, comm));
          xmpi(MPI_Pack((void *)filename, (int)filename_len, MPI_CHAR,
                        msgBuffer, msgSize, &msgBufPos, comm));
          xmpi(MPI_Sendrecv(msgBuffer, msgSize, MPI_PACKED, collRank,
                            STREAMOPEN,
                            &fileID, 1, MPI_INT, collRank, STREAMOPEN,
                            comm, MPI_STATUS_IGNORE));
          Free(msgBuffer);
        }
      else
        xmpi(MPI_Recv(&fileID, 1, MPI_INT, collRank, STREAMOPEN,
                      comm, MPI_STATUS_IGNORE));
      if (fileID >= 0)
        {
          streamptr->filetype = filetype;
          cdiPioClientStreamWinInit(streamID);
        }
    }
  else
    Error("cdiPIO read support not implemented");
  return fileID;
}

static void
cdiPioClientStreamDefVlist_(int streamID, int vlistID)
{
  cdiStreamDefVlist_(streamID, vlistID);
  int clientRank = commInqRankModel(),
    numClients = cdiPioCommInqSizeClients(),
    numColl = commInqSizeColl(),
    collRank = cdiPioCollRank(clientRank, numClients, numColl);
  int sendRPCData
    = (clientRank
       == cdiPioClientRangeStart(collRank, numClients, numColl));
  if (sendRPCData)
    {
      MPI_Comm comm = cdiPioInqInterComm();
      reshSetStatus(streamID, &streamOps,
                    reshGetStatus(streamID, &streamOps) & ~RESH_SYNC_BIT);
      char *msgBuffer;
      int msgSize = 0;
      int msgBufPos = reshPackBufferCreate(&msgBuffer, &msgSize, &comm);
      {
        int size;
        xmpi(MPI_Pack_size(defVlistNInts, MPI_INT, comm, &size));
        msgSize += size;
      }
      /* optimize: pos + size */
      msgBuffer = Realloc(msgBuffer, (size_t)msgSize);
      int msgData[defVlistNInts] = { streamID, streamInqVlist(streamID) };
      xmpi(MPI_Pack(&msgData, defVlistNInts, MPI_INT,
                    msgBuffer, msgSize, &msgBufPos, comm));
      xmpi(MPI_Send(msgBuffer, msgBufPos, MPI_PACKED, collRank,
                    STREAMDEFVLIST, comm));
      Free(msgBuffer);
    }
  struct collSpec cspec = { .numClients = numClients,
                            .numServers = numColl,
                            .sendRPCData = sendRPCData };
  cdiPioClientStreamWinCreate(streamID, &cspec);
}

static void
cdiPioClientStreamWriteVar_(int streamID, int varID, int memtype,
                            const void *data, int nmiss)
{
  (void)streamID; (void)varID; (void)memtype; (void)data; (void)nmiss;
  xabort("parallel writing requires explicit partition information,"
         " use streamWriteVarPart!");
}

static void
cdiPioClientStreamWriteVarChunk_(int streamID, int varID, int memtype,
                                 const int rect[][2],
                                 const void *data, int nmiss)
{
  /* todo: handle transmission of float data */
  if (memtype != MEMTYPE_DOUBLE)
    Error("Writing of non-double type data not implemented!");
  int vlistID = streamInqVlist(streamID);
  int size = vlistInqVarSize(vlistID, varID),
    varShape[3];
  unsigned ndims = (unsigned)cdiPioQueryVarDims(varShape, vlistID, varID);
  Xt_int varShapeXt[3], origin[3] = { 0, 0, 0 };
  int chunkShape[3] = { 1, 1, 1 };
  /* FIXME: verify xt_int ranges are good enough */
  for (unsigned i = 0; i < 3; ++i)
    varShapeXt[i] = varShape[i];
  for (unsigned i = 0; i < ndims; ++i)
    chunkShape[i] = rect[i][1] - rect[i][0] + 1;
  int varSize = varShape[0] * varShape[1] * varShape[2];
  xassert(varSize == size);
  Xt_idxlist chunkDesc
    = xt_idxsection_new(0, (int)ndims, varShapeXt, chunkShape, origin);
  pioBufferPartData(streamID, varID, data, nmiss, chunkDesc);
  xt_idxlist_delete(chunkDesc);
}

static void
cdiPioClientStreamWriteVarPart(int streamID, int varID, const void *data,
                               int nmiss, Xt_idxlist partDesc)
{
  pioBufferPartData(streamID, varID, data, nmiss, partDesc);
}

static void
cdiPioClientStreamWriteScatteredVarPart(int streamID, int varID,
                                        const void *data,
                                        int numBlocks, const int blocklengths[],
                                        const int displacements[],
                                        int nmiss, Xt_idxlist partDesc)
{
  cdiPioBufferPartDataGather(streamID, varID, data, numBlocks,
                             blocklengths, displacements, nmiss, partDesc);
}

#if defined HAVE_LIBNETCDF
static void
cdiPioCdfDefTimestepNOP(stream_t *streamptr, int tsID)
{
  (void)streamptr; (void)tsID;
}
#endif

static void
cdiPioClientStreamNOP(stream_t *streamptr)
{
  (void)streamptr;
}


static void
cdiPioClientStreamClose(stream_t *streamptr, int recordBufIsToBeDeleted)
{
  (void)recordBufIsToBeDeleted;
  int streamID = streamptr->self;
  int clientRank = commInqRankModel(),
    numClients = cdiPioCommInqSizeClients(),
    numColl = commInqSizeColl(),
    collRank = cdiPioCollRank(clientRank, numClients, numColl);
  if (clientRank
      == cdiPioClientRangeStart(collRank, numClients, numColl))
    {
      MPI_Comm comm = cdiPioInqInterComm();
      reshSetStatus(streamID, &streamOps,
                    reshGetStatus(streamID, &streamOps) & ~RESH_SYNC_BIT);
      char *msgBuffer;
      int msgSize = 0;
      int msgBufPos = reshPackBufferCreate(&msgBuffer, &msgSize, &comm);
      {
        int size;
        xmpi(MPI_Pack_size(1, MPI_INT, comm, &size));
        msgSize += size;
      }
      /* optimize: pos + size */
      msgBuffer = Realloc(msgBuffer, (size_t)msgSize);
      xmpi(MPI_Pack(&streamptr->self, 1, MPI_INT,
                    msgBuffer, msgSize, &msgBufPos, comm));
      xmpi(MPI_Send(msgBuffer, msgBufPos, MPI_PACKED, collRank,
                    STREAMCLOSE, comm));
      Free(msgBuffer);
    }
  cdiPioClientStreamWinDestroy(streamID);
}

static void
cdiPioTaxisPackWrap(void *data, void *buf, int size, int *pos,
                    void *context)
{
  int taxisID = (int)(intptr_t)data;
  reshPackResource(taxisID, &taxisOps, buf, size, pos, context);
}

static int
cdiPioClientStreamDefTimestep_(stream_t *streamptr, int tsID)
{
  int taxisID = vlistInqTaxis(streamptr->vlistID);
  struct winHeaderEntry header = (struct winHeaderEntry){
    .id = STREAMDEFTIMESTEP,
    .specific.funcArgs.streamNewTimestep = { streamptr->self, tsID } };
  xassert(sizeof (void *) >= sizeof (int));
  pioBufferFuncCall(streamptr->self, header,
                    (void *)(intptr_t)taxisID, cdiPioTaxisPackWrap);
  return cdiStreamDefTimestep_(streamptr, tsID);
}

void
cdiPioClientSetup(int *pioNamespace_, int *pioNamespace)
{
  *pioNamespace_ = *pioNamespace = namespaceNew();
  int callerCDINamespace = namespaceGetActive();
  namespaceSetActive(*pioNamespace_);
  cdiPioSerializeSetMPI();
  namespaceSwitchSet(NSSWITCH_ABORT, NSSW_FUNC(cdiAbortC_MPI));
  namespaceSwitchSet(NSSWITCH_WARNING, NSSW_FUNC(cdiPioWarning));
  namespaceSwitchSet(NSSWITCH_STREAM_OPEN_BACKEND,
                     NSSW_FUNC(cdiPioClientStreamOpen));
  namespaceSwitchSet(NSSWITCH_STREAM_DEF_VLIST_,
                     NSSW_FUNC(cdiPioClientStreamDefVlist_));
  namespaceSwitchSet(NSSWITCH_STREAM_WRITE_VAR_,
                     NSSW_FUNC(cdiPioClientStreamWriteVar_));
  namespaceSwitchSet(NSSWITCH_STREAM_WRITE_VAR_CHUNK_,
                     NSSW_FUNC(cdiPioClientStreamWriteVarChunk_));
  namespaceSwitchSet(NSSWITCH_STREAM_WRITE_VAR_PART_,
                     NSSW_FUNC(cdiPioClientStreamWriteVarPart));
  namespaceSwitchSet(NSSWITCH_STREAM_WRITE_SCATTERED_VAR_PART_,
                     NSSW_FUNC(cdiPioClientStreamWriteScatteredVarPart));
  namespaceSwitchSet(NSSWITCH_STREAM_CLOSE_BACKEND,
                     NSSW_FUNC(cdiPioClientStreamClose));
  namespaceSwitchSet(NSSWITCH_STREAM_DEF_TIMESTEP_,
                     NSSW_FUNC(cdiPioClientStreamDefTimestep_));
  namespaceSwitchSet(NSSWITCH_STREAM_SYNC,
                     NSSW_FUNC(cdiPioClientStreamNOP));
#ifdef HAVE_LIBNETCDF
  namespaceSwitchSet(NSSWITCH_CDF_DEF_TIMESTEP,
                     NSSW_FUNC(cdiPioCdfDefTimestepNOP));
  namespaceSwitchSet(NSSWITCH_CDF_STREAM_SETUP,
                     NSSW_FUNC(cdiPioClientStreamNOP));
#endif
  namespaceSetActive(callerCDINamespace);
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
