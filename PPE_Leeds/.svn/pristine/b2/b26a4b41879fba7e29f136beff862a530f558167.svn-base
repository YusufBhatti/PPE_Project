#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>

#include "dmemory.h"
#include "extra.h"
#include "error.h"
#include "file.h"
#include "binary.h"
#include "exse.h"
#include "swap.h"


enum {
  EXT_HEADER_LEN = 4,
};


static int initExtLib       = 0;
static int extDefaultPrec   = 0;
static int extDefaultNumber = EXT_REAL;


/*
 * A version string.
 */
#undef  LIBVERSION
#define LIBVERSION      1.4.0
#define XSTRING(x)	#x
#define STRING(x)	XSTRING(x)
static const char ext_libvers[] = STRING(LIBVERSION) " of " __DATE__ " " __TIME__ ;

const char *extLibraryVersion(void)
{
  return ext_libvers;
}


static int EXT_Debug = 0;    /* If set to 1, debugging */


void extDebug(int debug)
{
  EXT_Debug = debug;

  if ( EXT_Debug )
    Message("debug level %d", debug);
}


static
void extLibInit()
{
  const char *envName = "EXT_PRECISION";

  char *envString = getenv(envName);
  if ( envString )
    {
      int pos = 0;

      if ( strlen(envString) == 2  )
	{
	  switch ( tolower((int) envString[pos]) )
	    {
	    case 'r':
	      {
		extDefaultNumber = EXT_REAL;
		switch ( (int) envString[pos+1] )
		  {
		  case '4': extDefaultPrec = EXSE_SINGLE_PRECISION; break;
		  case '8': extDefaultPrec = EXSE_DOUBLE_PRECISION; break;
		  default:
		    Message("Invalid digit in %s: %s", envName, envString);
		  }
		break;
	      }
	    case 'c':
	      {
		extDefaultNumber = EXT_COMP;
		switch ( (int) envString[pos+1] )
		  {
		  case '4': extDefaultPrec = EXSE_SINGLE_PRECISION; break;
		  case '8': extDefaultPrec = EXSE_DOUBLE_PRECISION; break;
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
	}
    }

  initExtLib = 1;
}

static
void extInit(extrec_t *extp)
{
  extp->checked    = 0;
  extp->byteswap   = 0;
  extp->prec       = 0;
  extp->number     = extDefaultNumber;
  extp->datasize   = 0;
  extp->buffersize = 0;
  extp->buffer     = NULL;
}


void *extNew(void)
{
  if ( ! initExtLib ) extLibInit();

  extrec_t *extp = (extrec_t *) Malloc(sizeof(extrec_t));

  extInit(extp);

  return (void*)extp;
}


void extDelete(void *ext)
{
  extrec_t *extp = (extrec_t *) ext;

  if ( extp )
    {
      if ( extp->buffer ) Free(extp->buffer);
      Free(extp);
    }
}


int extCheckFiletype(int fileID, int *swap)
{
  size_t fact = 0;
  size_t data =  0;
  size_t dimxy = 0;
  int found = 0;
  unsigned char buffer[40], *pbuf;

  if ( fileRead(fileID, buffer, 4) != 4 ) return found;

  size_t blocklen  = (size_t) get_UINT32(buffer);
  size_t sblocklen = (size_t) get_SUINT32(buffer);

  if ( EXT_Debug )
    Message("blocklen = %d sblocklen = %d", blocklen, sblocklen);

  if ( blocklen == 16 )
    {
     *swap = 0;
      fact = blocklen/4;
      if ( fileRead(fileID, buffer, blocklen+8) != blocklen+8 ) return found;
      pbuf = buffer+3*fact;      dimxy = (size_t) get_UINT32(pbuf);
      pbuf = buffer+blocklen+4;  data  = (size_t) get_UINT32(pbuf);
    }
  else if ( blocklen == 32 )
    {
     *swap = 0;
      fact = blocklen/4;
      if ( fileRead(fileID, buffer, blocklen+8) != blocklen+8 ) return found;
      pbuf = buffer+3*fact;      dimxy = (size_t) get_UINT64(pbuf);
      pbuf = buffer+blocklen+4;  data  = (size_t) get_UINT32(pbuf);
    }
  else if ( sblocklen == 16 )
    {
     *swap = 1;
      fact = sblocklen/4;
      if ( fileRead(fileID, buffer, sblocklen+8) != sblocklen+8 ) return found;
      pbuf = buffer+3*fact;       dimxy = (size_t) get_SUINT32(pbuf);
      pbuf = buffer+sblocklen+4;  data  = (size_t) get_SUINT32(pbuf);
    }
  else if ( sblocklen == 32 )
    {
     *swap = 1;
      fact = sblocklen/4;
      if ( fileRead(fileID, buffer, sblocklen+8) != sblocklen+8 ) return found;
      pbuf = buffer+3*fact;       dimxy = (size_t) get_SUINT64(pbuf);
      pbuf = buffer+sblocklen+4;  data  = (size_t) get_SUINT32(pbuf);
    }

  fileRewind(fileID);

  if      ( data && dimxy*fact   == data ) found = 1;
  else if ( data && dimxy*fact*2 == data ) found = 1;

  if ( EXT_Debug )
    {
      Message("swap = %d fact = %d", *swap, fact);
      Message("dimxy = %lu data = %lu", dimxy, data);
    }

  return found;
}


int extInqHeader(void *ext, int *header)
{
  extrec_t *extp = (extrec_t *) ext;

  for ( size_t i = 0; i < EXT_HEADER_LEN; i++ )
    header[i] = extp->header[i];

  if ( EXT_Debug ) Message("datasize = %lu", extp->datasize);

  return 0;
}


int extDefHeader(void *ext, const int *header)
{
  extrec_t *extp = (extrec_t *) ext;

  for ( size_t i = 0; i < EXT_HEADER_LEN; i++ )
    extp->header[i] = header[i];

  extp->datasize = (size_t)header[3];
  if ( extp->number == EXT_COMP ) extp->datasize *= 2;

  if ( EXT_Debug ) Message("datasize = %lu", extp->datasize);

  return 0;
}

static
int extInqData(extrec_t *extp, int prec, void *data)
{
  size_t i;
  int ierr = 0;
  int byteswap = extp->byteswap;
  size_t datasize = extp->datasize;
  void *buffer   = extp->buffer;
  int rprec    = extp->prec;

  switch ( rprec )
    {
    case EXSE_SINGLE_PRECISION:
      {
	if ( sizeof(FLT32) == 4 )
	  {
	    if ( byteswap ) swap4byte(buffer, datasize);

	    if ( rprec == prec )
	      memcpy(data, buffer, datasize*sizeof(FLT32));
	    else
	      for ( i = 0; i < datasize; ++i )
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

	    if ( rprec == prec )
	      memcpy(data, buffer, datasize*sizeof(FLT64));
	    else
	      for ( i = 0; i < datasize; ++i )
		((float *) data)[i] = (float) ((double *) buffer)[i];
	  }
	else
	  {
	    Error("not implemented for %d byte float", sizeof(FLT64));
	  }
	break;
    default:
      {
	Error("unexpected data precision %d", rprec);
	break;
      }
    }

  return ierr;
}


int extInqDataSP(void *ext, float *data)
{
  return extInqData((extrec_t *)ext, EXSE_SINGLE_PRECISION, (void *) data);
}


int extInqDataDP(void *ext, double *data)
{
  return extInqData((extrec_t *)ext, EXSE_DOUBLE_PRECISION, (void *) data);
}


static int extDefData(void *ext, int prec, const void *data)
{
  extrec_t *extp = (extrec_t *) ext;
  size_t i;
  int rprec;
  void *buffer;

  if ( extDefaultPrec ) rprec = extDefaultPrec;
  else                  rprec = extp->prec;

  if ( ! rprec ) rprec = prec;

  extp->prec = rprec;

  int *header = extp->header;

  size_t datasize = (size_t)header[3];
  if ( extp->number == EXT_COMP ) datasize *= 2;
  size_t blocklen = datasize * (size_t)rprec;

  extp->datasize = datasize;

  size_t buffersize = extp->buffersize;

  if ( buffersize != blocklen )
    {
      buffersize = blocklen;
      buffer = extp->buffer;
      buffer = Realloc(buffer, buffersize);
      extp->buffer = buffer;
      extp->buffersize = buffersize;
    }
  else
    buffer = extp->buffer;

  switch ( rprec )
    {
    case EXSE_SINGLE_PRECISION:
      {
	if ( rprec == prec )
	  memcpy(buffer, data, datasize*sizeof(FLT32));
	else
	  for ( i = 0; i < datasize; i++ )
	    ((float *) buffer)[i] = (float) ((double *) data)[i];

	break;
      }
    case EXSE_DOUBLE_PRECISION:
      {
	if ( rprec == prec )
	  memcpy(buffer, data, datasize*sizeof(FLT64));
	else
	  for ( i = 0; i < datasize; i++ )
	    ((double *) buffer)[i] = (double) ((float *) data)[i];

	break;
      }
    default:
      {
	Error("unexpected data precision %d", rprec);
        break;
      }
    }

  return 0;
}


int extDefDataSP(void *ext, const float *data)
{
  return extDefData(ext, EXSE_SINGLE_PRECISION, (void *) data);
}


int extDefDataDP(void *ext, const double *data)
{
  return extDefData(ext, EXSE_DOUBLE_PRECISION, (void *) data);
}


int extRead(int fileID, void *ext)
{
  extrec_t *extp = (extrec_t *) ext;
  size_t i;
  void *buffer;
  int status;

  if ( ! extp->checked )
    {
      status = extCheckFiletype(fileID, &extp->byteswap);
      if ( status == 0 ) Error("Not a EXTRA file!");
      extp->checked = 1;
    }

  int byteswap = extp->byteswap;

  /* read header record */
  size_t blocklen = binReadF77Block(fileID, byteswap);

  if ( fileEOF(fileID) ) return -1;

  if ( EXT_Debug )
    Message("blocklen = %lu", blocklen);

  size_t hprec = blocklen / EXT_HEADER_LEN;

  extp->prec = (int)hprec;

  switch ( hprec )
    {
    case EXSE_SINGLE_PRECISION:
      {
        INT32 tempheader[4];
	binReadInt32(fileID, byteswap, EXT_HEADER_LEN, tempheader);

	for ( i = 0; i < EXT_HEADER_LEN; i++ )
          extp->header[i] = (int)tempheader[i];

	break;
      }
    case EXSE_DOUBLE_PRECISION:
      {
        INT64 tempheader[4];
	binReadInt64(fileID, byteswap, EXT_HEADER_LEN, tempheader);

	for ( i = 0; i < EXT_HEADER_LEN; i++ )
          extp->header[i] = (int)tempheader[i];

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

  extp->datasize = (size_t)extp->header[3];

  if ( EXT_Debug ) Message("datasize = %lu", extp->datasize);

  blocklen = binReadF77Block(fileID, byteswap);

  size_t buffersize = (size_t)extp->buffersize;

  if ( buffersize < blocklen )
    {
      buffersize = blocklen;
      buffer = extp->buffer;
      buffer = Realloc(buffer, buffersize);
      extp->buffer = buffer;
      extp->buffersize = buffersize;
    }
  else
    buffer = extp->buffer;

  size_t dprec = blocklen / extp->datasize;

  if ( dprec == hprec )
    {
      extp->number = EXT_REAL;
    }
  else if ( dprec == 2*hprec )
    {
      dprec /= 2;
      extp->datasize *= 2;
      extp->number = EXT_COMP;
    }

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


int extWrite(int fileID, void *ext)
{
  extrec_t *extp = (extrec_t *) ext;
  size_t i;
  union { INT32 i32[EXT_HEADER_LEN]; INT64 i64[EXT_HEADER_LEN]; } tempheader;
  int byteswap = extp->byteswap;
  int rprec  = extp->prec;
  int number = extp->number;
  int *header = extp->header;

  /* write header record */
  size_t blocklen = EXT_HEADER_LEN * (size_t)rprec;

  binWriteF77Block(fileID, byteswap, blocklen);

  switch ( rprec )
    {
    case EXSE_SINGLE_PRECISION:
      {
	for ( i = 0; i < EXT_HEADER_LEN; i++ )
          tempheader.i32[i] = (INT32) header[i];

	binWriteInt32(fileID, byteswap, EXT_HEADER_LEN, tempheader.i32);

	break;
      }
    case EXSE_DOUBLE_PRECISION:
      {
	for ( i = 0; i < EXT_HEADER_LEN; i++ )
          tempheader.i64[i] = (INT64) header[i];

	binWriteInt64(fileID, byteswap, EXT_HEADER_LEN, tempheader.i64);

	break;
      }
    default:
      {
	Error("unexpected header precision %d", rprec);
        break;
      }
    }

  binWriteF77Block(fileID, byteswap, blocklen);

  size_t datasize = (size_t)header[3];
  if ( number == EXT_COMP ) datasize *= 2;
  blocklen = datasize * (size_t)rprec;

  binWriteF77Block(fileID, byteswap, blocklen);

  extp->datasize = datasize;

  void *buffer = extp->buffer;

  switch ( rprec )
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
	Error("unexpected data precision %d", rprec);
        break;
      }
    }

  binWriteF77Block(fileID, byteswap, blocklen);

  return 0;
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
