#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>

#include "dmemory.h"
#include "service.h"
#include "error.h"
#include "file.h"
#include "binary.h"
#include "exse.h"
#include "swap.h"


enum {
  SRV_HEADER_LEN = 8,
};


static int initSrvLib      = 0;
static int srvDefaultHprec = 0;
static int srvDefaultDprec = 0;


/*
 * A version string.
 */
#undef  LIBVERSION
#define LIBVERSION      1.4.0
#define XSTRING(x)	#x
#define STRING(x)	XSTRING(x)
static const char srv_libvers[] = STRING(LIBVERSION) " of " __DATE__" " __TIME__;

const char *srvLibraryVersion(void)
{
  return srv_libvers;
}


static int SRV_Debug = 0;    /* If set to 1, debugging */


void srvDebug(int debug)
{
  SRV_Debug = debug;

  if ( SRV_Debug )
    Message("debug level %d", debug);
}

static
void srvLibInit()
{
  const char *envName = "SRV_PRECISION";

  char *envString = getenv(envName);
  if ( envString )
    {
      int nrun;
      if ( strlen(envString) == 2 ) nrun = 1;
      else                          nrun = 2;

      int pos = 0;
      while ( nrun-- )
	{
	  switch ( tolower((int) envString[pos]) )
	    {
	    case 'i':
	      {
		switch ( (int) envString[pos+1] )
		  {
		  case '4': srvDefaultHprec = EXSE_SINGLE_PRECISION; break;
		  case '8': srvDefaultHprec = EXSE_DOUBLE_PRECISION; break;
		  default:
		    Message("Invalid digit in %s: %s", envName, envString);
		  }
		break;
	      }
	    case 'r':
	      {
		switch ( (int) envString[pos+1] )
		  {
		  case '4': srvDefaultDprec = EXSE_SINGLE_PRECISION; break;
		  case '8': srvDefaultDprec = EXSE_DOUBLE_PRECISION; break;
		  default:
		    Message("Invalid digit in %s: %s", envName, envString);
		  }
		break;
	      }
	    default:
              {
                Message("Invalid character in %s: %s", envName, envString);
                break;
              }
            }
	  pos += 2;
	}
    }

  initSrvLib = 1;
}

static
void srvInit(srvrec_t *srvp)
{
  srvp->checked    = 0;
  srvp->byteswap   = 0;
  srvp->hprec      = 0;
  srvp->dprec      = 0;
  srvp->datasize   = 0;
  srvp->buffersize = 0;
  srvp->buffer     = NULL;
}


void *srvNew(void)
{
  if ( ! initSrvLib ) srvLibInit();

  srvrec_t *srvp = (srvrec_t *) Malloc(sizeof(srvrec_t));

  srvInit(srvp);

  return (void*)srvp;
}


void srvDelete(void *srv)
{
  srvrec_t *srvp = (srvrec_t *) srv;

  if ( srvp )
    {
      if ( srvp->buffer ) Free(srvp->buffer);
      Free(srvp);
    }
}


int srvCheckFiletype(int fileID, int *swap)
{
  size_t data = 0;
  size_t dimx = 0, dimy = 0;
  size_t fact = 0;
  int found = 0;
  unsigned char buffer[72], *pbuf;

  if ( fileRead(fileID, buffer, 4) != 4 ) return found;

  size_t blocklen  = (size_t) get_UINT32(buffer);
  size_t sblocklen = (size_t) get_SUINT32(buffer);

  if ( SRV_Debug )
    Message("blocklen = %d sblocklen = %d", blocklen, sblocklen);

  if ( blocklen == 32 )
    {
     *swap = 0;
      fact = blocklen>>3;
      if ( fileRead(fileID, buffer, blocklen+8) != blocklen+8 ) return found;
      pbuf = buffer+4*fact;      dimx = (size_t) get_UINT32(pbuf);
      pbuf = buffer+5*fact;      dimy = (size_t) get_UINT32(pbuf);
      pbuf = buffer+blocklen+4;  data = (size_t) get_UINT32(pbuf);
    }
  else if ( blocklen == 64 )
    {
     *swap = 0;
      fact = blocklen>>3;
      if ( fileRead(fileID, buffer, blocklen+8) != blocklen+8 ) return found;
      pbuf = buffer+4*fact;      dimx = (size_t) get_UINT64(pbuf);
      pbuf = buffer+5*fact;      dimy = (size_t) get_UINT64(pbuf);
      pbuf = buffer+blocklen+4;  data = (size_t) get_UINT32(pbuf);
    }
  else if ( sblocklen == 32 )
    {
     *swap = 1;
      fact = sblocklen>>3;
      if ( fileRead(fileID, buffer, sblocklen+8) != sblocklen+8 ) return found;
      pbuf = buffer+4*fact;       dimx = (size_t) get_SUINT32(pbuf);
      pbuf = buffer+5*fact;       dimy = (size_t) get_SUINT32(pbuf);
      pbuf = buffer+sblocklen+4;  data = (size_t) get_SUINT32(pbuf);
    }
  else if ( sblocklen == 64 )
    {
     *swap = 1;
      fact = sblocklen>>3;
      if ( fileRead(fileID, buffer, sblocklen+8) != sblocklen+8 ) return found;
      pbuf = buffer+4*fact;       dimx = (size_t) get_SUINT64(pbuf);
      pbuf = buffer+5*fact;       dimy = (size_t) get_SUINT64(pbuf);
      pbuf = buffer+sblocklen+4;  data = (size_t) get_SUINT32(pbuf);
    }

  fileRewind(fileID);

  if      ( data && dimx*dimy*fact == data ) found = 1;
  else if ( data && dimx*dimy*8    == data ) found = 1;

  if ( SRV_Debug )
    {
      Message("swap = %d fact = %d", *swap, fact);
      Message("dimx = %lu dimy = %lu data = %lu", dimx, dimy, data);
    }

  return found;
}


int srvInqHeader(void *srv, int *header)
{
  srvrec_t *srvp = (srvrec_t *) srv;

  for ( size_t i = 0; i < SRV_HEADER_LEN; i++ )
    header[i] = srvp->header[i];

  if ( SRV_Debug )
    Message("datasize = %lu", srvp->datasize);

  return 0;
}


int srvDefHeader(void *srv, const int *header)
{
  srvrec_t *srvp = (srvrec_t *) srv;

  for ( size_t i = 0; i < SRV_HEADER_LEN; i++ )
    srvp->header[i] = header[i];

  srvp->datasize = (size_t)(header[4] * header[5]);

  if ( SRV_Debug )
    Message("datasize = %lu", srvp->datasize);

  return 0;
}

static
int srvInqData(srvrec_t *srvp, int prec, void *data)
{
  size_t i;
  int ierr = 0;
  int byteswap = srvp->byteswap;
  size_t datasize = srvp->datasize;
  void *buffer = srvp->buffer;
  int dprec = srvp->dprec;

  switch ( dprec )
    {
    case EXSE_SINGLE_PRECISION:
      {
	if ( sizeof(FLT32) == 4 )
	  {
	    if ( byteswap ) swap4byte(buffer, datasize);

	    if ( dprec == prec )
	      memcpy(data, buffer, datasize*sizeof(FLT32));
	    else
	      for ( i = 0; i < datasize; i++ )
		((double *) data)[i] = (double) ((float *) buffer)[i];
	  }
	else
	  {
	    Error("not implemented for %d byte float", sizeof(FLT32));
	  }
	break;
      }
    case EXSE_DOUBLE_PRECISION:
	if ( sizeof(FLT64) == 8 )
	  {
	    if ( byteswap ) swap8byte(buffer, datasize);

	    if ( dprec == prec )
	      memcpy(data, buffer, datasize*sizeof(FLT64));
	    else
	      for ( i = 0; i < datasize; i++ )
		((float *) data)[i] = (float) ((double *) buffer)[i];
	  }
	else
	  {
	    Error("not implemented for %d byte float", sizeof(FLT64));
	  }
	break;
    default:
      {
	Error("unexpected data precision %d", dprec);
        break;
      }
    }

  return ierr;
}


int srvInqDataSP(void *srv, float *data)
{
  return srvInqData((srvrec_t *)srv, EXSE_SINGLE_PRECISION, (void *) data);
}


int srvInqDataDP(void *srv, double *data)
{
  return srvInqData((srvrec_t *)srv, EXSE_DOUBLE_PRECISION, (void *) data);
}


static int
srvDefData(void *srv, int prec, const void *data)
{
  srvrec_t *srvp = (srvrec_t *) srv;
  size_t i;
  int dprec, hprec;
  void *buffer;

  if ( srvDefaultDprec ) dprec = srvDefaultDprec;
  else                   dprec = srvp->dprec;

  if ( ! dprec ) dprec = prec;

  srvp->dprec = dprec;

  if ( srvDefaultHprec ) hprec = srvDefaultHprec;
  else                   hprec = srvp->hprec;

  if ( ! hprec ) hprec = dprec;

  srvp->hprec = hprec;

  int *header = srvp->header;

  size_t datasize = (size_t)(header[4] * header[5]);
  size_t blocklen = datasize * (size_t)dprec;

  srvp->datasize = datasize;

  size_t buffersize = srvp->buffersize;

  if ( buffersize != blocklen )
    {
      buffersize = blocklen;
      buffer = srvp->buffer;
      buffer = Realloc(buffer, buffersize);
      srvp->buffer = buffer;
      srvp->buffersize = buffersize;
    }
  else
    buffer = srvp->buffer;

  switch ( dprec )
    {
    case EXSE_SINGLE_PRECISION:
      {
	if ( dprec == prec )
	  memcpy(buffer, data, datasize*sizeof(FLT32));
	else
	  for ( i = 0; i < datasize; i++ )
	    ((float *) buffer)[i] = (float) ((double *) data)[i];

	break;
      }
    case EXSE_DOUBLE_PRECISION:
      {
	if ( dprec == prec )
	  memcpy(buffer, data, datasize*sizeof(FLT64));
	else
	  for ( i = 0; i < datasize; i++ )
	    ((double *) buffer)[i] = (double) ((float *) data)[i];

	break;
      }
    default:
      {
	Error("unexpected data precision %d", dprec);
        break;
      }
    }

  return 0;
}


int srvDefDataSP(void *srv, const float *data)
{
  return srvDefData(srv, EXSE_SINGLE_PRECISION, (void *) data);
}


int srvDefDataDP(void *srv, const double *data)
{
  return srvDefData(srv, EXSE_DOUBLE_PRECISION, (void *) data);
}


int srvRead(int fileID, void *srv)
{
  srvrec_t *srvp = (srvrec_t *) srv;
  size_t i;
  union {
    INT32 i32[SRV_HEADER_LEN];
    INT64 i64[SRV_HEADER_LEN];
  } tempheader;
  void *buffer;
  int status;

  if ( ! srvp->checked )
    {
      status = srvCheckFiletype(fileID, &srvp->byteswap);
      if ( status == 0 ) Error("Not a SERVICE file!");
      srvp->checked = 1;
    }

  int byteswap = srvp->byteswap;

  /* read header record */
  size_t blocklen = binReadF77Block(fileID, byteswap);

  if ( fileEOF(fileID) ) return -1;

  if ( SRV_Debug )
    Message("blocklen = %lu", blocklen);

  size_t hprec = blocklen / SRV_HEADER_LEN;

  srvp->hprec = (int)hprec;

  switch ( hprec )
    {
    case EXSE_SINGLE_PRECISION:
      {
	binReadInt32(fileID, byteswap, SRV_HEADER_LEN, tempheader.i32);

	for ( i = 0; i < SRV_HEADER_LEN; i++ )
          srvp->header[i] = (int)tempheader.i32[i];

	break;
      }
    case EXSE_DOUBLE_PRECISION:
      {
	binReadInt64(fileID, byteswap, SRV_HEADER_LEN, tempheader.i64);

	for ( i = 0; i < SRV_HEADER_LEN; i++ )
          srvp->header[i] = (int)tempheader.i64[i];

	break;
      }
    default:
      {
	Error("Unexpected header precision %d", hprec);
        break;
      }
    }

  size_t blocklen2 = binReadF77Block(fileID, byteswap);

  if ( blocklen2 != blocklen )
    {
      Warning("Header blocklen differ (blocklen1=%d; blocklen2=%d)!", blocklen, blocklen2);
      if ( blocklen2 != 0 ) return -1;
    }

  srvp->datasize = (size_t)(srvp->header[4] * srvp->header[5]);

  if ( SRV_Debug )
    Message("datasize = %lu", srvp->datasize);

  blocklen = binReadF77Block(fileID, byteswap);

  size_t buffersize = srvp->buffersize;

  if ( buffersize < blocklen )
    {
      buffersize = blocklen;
      buffer = srvp->buffer;
      buffer = Realloc(buffer, buffersize);
      srvp->buffer = buffer;
      srvp->buffersize = buffersize;
    }
  else
    buffer = srvp->buffer;

  size_t datasize = srvp->datasize;

  size_t dprec = blocklen / datasize;

  srvp->dprec = (int)dprec;

  if ( dprec != EXSE_SINGLE_PRECISION && dprec != EXSE_DOUBLE_PRECISION )
    {
      Warning("Unexpected data precision %d", dprec);
      return -1;
    }

  fileRead(fileID, buffer, blocklen);

  blocklen2 = binReadF77Block(fileID, byteswap);

  if ( blocklen2 != blocklen )
    {
      Warning("Data blocklen differ (blocklen1=%d; blocklen2=%d)!", blocklen, blocklen2);
      if ( blocklen2 != 0 ) return -1;
    }

  return 0;
}


void srvWrite(int fileID, void *srv)
{
  srvrec_t *srvp = (srvrec_t *) srv;
  size_t i;
  union
  {
    INT32 i32[SRV_HEADER_LEN];
    INT64 i64[SRV_HEADER_LEN];
  } tempheader;
  int byteswap = srvp->byteswap;
  int dprec  = srvp->dprec;
  int hprec  = srvp->hprec;
  int *restrict header = srvp->header;

  /* write header record */
  size_t blocklen = SRV_HEADER_LEN * (size_t)hprec;

  binWriteF77Block(fileID, byteswap, blocklen);

  switch ( hprec )
    {
    case EXSE_SINGLE_PRECISION:
      {
	for (i = 0; i < SRV_HEADER_LEN; i++)
          tempheader.i32[i] = (INT32) header[i];

	binWriteInt32(fileID, byteswap, SRV_HEADER_LEN, tempheader.i32);

	break;
      }
    case EXSE_DOUBLE_PRECISION:
      {
	for (i = 0; i < SRV_HEADER_LEN; i++)
          tempheader.i64[i] = (INT64) header[i];

	binWriteInt64(fileID, byteswap, SRV_HEADER_LEN, tempheader.i64);

	break;
      }
    default:
      {
	Error("unexpected header precision %d", hprec);
        break;
      }
    }

  binWriteF77Block(fileID, byteswap, blocklen);

  size_t datasize = (size_t)(header[4] * header[5]);
  blocklen = datasize * (size_t)dprec;

  binWriteF77Block(fileID, byteswap, blocklen);

  srvp->datasize = datasize;

  void *buffer = srvp->buffer;

  switch ( dprec )
    {
    case EXSE_SINGLE_PRECISION:
      {
	binWriteFlt32(fileID, byteswap, datasize, (FLT32 *) buffer);
	break;
      }
    case EXSE_DOUBLE_PRECISION:
      {
	binWriteFlt64(fileID, byteswap, datasize, (FLT64 *) buffer);
	break;
      }
    default:
      {
	Error("unexpected data precision %d", dprec);
        break;
      }
    }

  binWriteF77Block(fileID, byteswap, blocklen);
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
