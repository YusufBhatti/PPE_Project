#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <inttypes.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

#include "cdi.h"
#include "dmemory.h"
#include "namespace.h"
#include "pio.h"
#include "pio_comm.h"
#include "pio_impl.h"
#include "pio_util.h"

typedef struct
{
  MPI_File fh;
  int fileID;
  MPI_Offset pos;
  char *name;
  int collWriteSize[];
} aFiledataM;

static listSet *bibAFiledataM;

static int
fileIDTest(void *a, void *fileID)
{
  return ((aFiledataM *)a)->fileID == (int)(intptr_t)fileID;
}


/***************************************************************/

static aFiledataM *
initAFiledataFileWriteAtAll(const char *filename, size_t bs)
{
  MPI_Comm commPio = commInqCommPio();
  int sizePio = commInqSizePio();
  size_t nameSize = strlen(filename) + 1;
  aFiledataM *of = Malloc(sizeof (*of)
                           + sizeof (of->collWriteSize[0]) * (size_t)sizePio
                           + nameSize);
  of->name = (char *)((unsigned char *)of + sizeof (*of)
                      + sizeof (of->collWriteSize[0]) * (size_t)sizePio);
  memcpy(of->name, filename, nameSize);

  MPI_Info open_info = MPI_INFO_NULL;
  xmpi(MPI_Info_create(&open_info));
  xmpi(MPI_Info_set(open_info, "access_style", "write_once"));
  xmpi(MPI_Info_set(open_info, "collective_buffering", "true"));
  /* tell IBM PE to buffer just as much as one buffer holds */
  {
    char buf_size_str[3*sizeof(size_t)*CHAR_BIT/8+1];
    snprintf(buf_size_str, sizeof (buf_size_str), "%zu", bs);
    xmpi(MPI_Info_set(open_info, "IBM_io_buffer_size", buf_size_str));
    xmpi(MPI_Info_set(open_info, "IBM_largeblock_io", "false"));
  }
  xmpi(MPI_File_open(commPio, of->name,
                     MPI_MODE_CREATE|MPI_MODE_WRONLY|MPI_MODE_UNIQUE_OPEN,
                     open_info, &of->fh));
  xmpi(MPI_Info_free(&open_info));
  of->pos = 0;

  return of;
}

/***************************************************************/

static int
destroyAFiledataFileWriteAtAll(void *v)
{
  aFiledataM *of = v;

  /* close file */
  MPI_Offset endpos, fsize;
  endpos = of->pos;
  xmpi(MPI_File_get_size(of->fh, &fsize));
  /* does the file need to be truncated? */
  MPI_Comm commPio = commInqCommPio();
  int trailingOctets = fsize > endpos;
  xmpi(MPI_Allreduce(MPI_IN_PLACE, &trailingOctets, 1, MPI_INT, MPI_LOR,
                     commPio));
  if (trailingOctets)
    xmpi(MPI_File_set_size(of->fh, endpos));
  int iret = MPI_File_close(&of->fh);

  Free(of);

  return iret == MPI_SUCCESS ? 0 : -1;
}

/***************************************************************/

static bool
compareNamesFileWriteAtAll(void *v1, void *v2)
{
  aFiledataM *afm1 = v1, *afm2 = v2;
  return !strcmp(afm1->name, afm2->name);
}

/***************************************************************/

static size_t
fwFileWriteAtAll(int fileID, const void *buffer, size_t len)
{
  aFiledataM *of
    = listSetGet(bibAFiledataM, fileIDTest, (void *)(intptr_t)fileID);
  xassert(of && len <= INT_MAX);
  MPI_Comm commPio = commInqCommPio();
  int sizePio = commInqSizePio(),
    rankPio = commInqRankPio();
  /* find position to write to */
  of->collWriteSize[rankPio] = (int)len;
  xmpi(MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                     of->collWriteSize, 1, MPI_INT, commPio));
  MPI_Offset myPos = of->pos, nextWritePos;
  for (size_t i = 0; i < (size_t)rankPio; ++i)
    myPos += of->collWriteSize[i];
  nextWritePos = myPos;
  for (size_t i = (size_t)rankPio; i < (size_t)sizePio; ++i)
    nextWritePos += of->collWriteSize[i];
  /* write buffer */
  xassert(len <= INT_MAX);
  xmpi(MPI_File_write_at_all(of->fh, myPos, (void *)buffer, (int)len,
                             MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE));
  of->pos = nextWritePos;
  return len;
}

/***************************************************************/

static int fcFileWriteAtAll(int fileID)
{
  aFiledataM *of
    = listSetGet(bibAFiledataM, fileIDTest, (void *)(intptr_t)fileID);
  if (!of)
    xabort("listSet, fileID=%d not found", fileID);
  int iret = listSetRemove(bibAFiledataM, fileIDTest, (void *)(intptr_t)fileID);
  return iret;
}

/***************************************************************/
static void
elemCheck(void *q, void *nm)
{
  aFiledataM *afm = q;
  const char *name = nm;

  if (!strcmp(name, afm->name))
    xabort("Filename %s has already been added to set\n", name);
}


static int
fowFileWriteAtAll(const char *filename, const char *mode)
{
  static unsigned long buffersize = 0;
  int id;
  enum {
    bcastRoot = 0
  };
  MPI_Comm commPio = commInqCommPio ();
  int rankPio = commInqRankPio ();

  if ((mode[0] != 'w' && mode[0] != 'W') || mode[0] == 0 || mode[1] != 0)
    xabort("Unsupported mode \"%s\" in parallel file open.", mode);

  /* broadcast buffersize to collectors ( just once, for all files )*/

  if (!buffersize)
    {
      if (rankPio == bcastRoot)
        buffersize = findWriteAccumBufsize();
      xmpi(MPI_Bcast(&buffersize, 1, MPI_UNSIGNED_LONG, bcastRoot, commPio));
    }

  xdebug("buffersize=%lu", buffersize);

  listSetForeach(bibAFiledataM, elemCheck, (void *)filename);
  aFiledataM *of = initAFiledataFileWriteAtAll(filename, (size_t)buffersize);

  if ((id = listSetAdd(bibAFiledataM, of)) < 0 )
    xabort("filename %s not unique", of->name);

  of->fileID = id;
  return id;
}

/***************************************************************/

static void finalizeFileWriteAtAll(void)
{
  if (!listSetIsEmpty(bibAFiledataM))
    xabort("set bibAFiledataM not empty");
  else
    {
      xdebug("%s", "destroy set");
      listSetDelete(bibAFiledataM);
    }
}

/***************************************************************/

void
cdiPioFileWriteAtAllInit(void)
{
  bibAFiledataM = listSetNew( destroyAFiledataFileWriteAtAll, compareNamesFileWriteAtAll );

  namespaceSwitchSet(NSSWITCH_FILE_OPEN, NSSW_FUNC(fowFileWriteAtAll));
  namespaceSwitchSet(NSSWITCH_FILE_CLOSE, NSSW_FUNC(fcFileWriteAtAll));
  namespaceSwitchSet(NSSWITCH_FILE_WRITE, NSSW_FUNC(fwFileWriteAtAll));
  cdiPioFileWritingFinalize = finalizeFileWriteAtAll;

  if ( bibAFiledataM == NULL )
    xabort ( "listSetNew did not succeed" );
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
