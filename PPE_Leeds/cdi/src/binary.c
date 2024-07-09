#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <inttypes.h>

#include "cdi.h"
#include "error.h"
#include "file.h"
#include "swap.h"
#include "binary.h"


UINT32 get_UINT32(unsigned char *x)
{
  switch (HOST_ENDIANNESS)
    {
    case CDI_BIGENDIAN:
      return((UINT32)(((UINT32)x[0]<<24)+((UINT32)x[1]<<16)+((UINT32)x[2]<< 8)+ (UINT32)x[3]));
    case CDI_LITTLEENDIAN:
      return ((UINT32)(((UINT32)x[3]<<24)+((UINT32)x[2]<<16)+((UINT32)x[1]<< 8)+ (UINT32)x[0]));
    default:
      Error("unhandled endianness %d", HOST_ENDIANNESS);
      return UINT32_C(0xFFFFFFFF);
    }
}


UINT32 get_SUINT32(unsigned char *x)
{
  switch (HOST_ENDIANNESS)
    {
    case CDI_BIGENDIAN:
      return((UINT32)(((UINT32)x[3]<<24)+((UINT32)x[2]<<16)+((UINT32)x[1]<< 8)+ (UINT32)x[0]));
    case CDI_LITTLEENDIAN:
      return((UINT32)(((UINT32)x[0]<<24)+((UINT32)x[1]<<16)+((UINT32)x[2]<< 8)+ (UINT32)x[3]));
    default:
      Error("unhandled endianness %d", HOST_ENDIANNESS);
      return UINT32_C(0xFFFFFFFF);
    }
}


UINT64 get_UINT64(unsigned char *x)
{
  switch (HOST_ENDIANNESS)
    {
    case CDI_BIGENDIAN:
      return((UINT64)(((UINT64)x[0]<<56)+((UINT64)x[1]<<48)+((UINT64)x[2]<<40)+((UINT64)x[3]<<32)+
                      ((UINT64)x[4]<<24)+((UINT64)x[5]<<16)+((UINT64)x[6]<< 8)+ (UINT64)x[7]));
    case CDI_LITTLEENDIAN:
      return((UINT64)(((UINT64)x[7]<<56)+((UINT64)x[6]<<48)+((UINT64)x[5]<<40)+((UINT64)x[4]<<32)+
                      ((UINT64)x[3]<<24)+((UINT64)x[2]<<16)+((UINT64)x[1]<< 8)+ (UINT64)x[0]));
    default:
      Error("unhandled endianness %d", HOST_ENDIANNESS);
      return UINT64_C(0xFFFFFFFFFFFFFFFF);
    }
}


UINT64 get_SUINT64(unsigned char *x)
{
  switch (HOST_ENDIANNESS)
    {
    case CDI_BIGENDIAN:
      return((UINT64)(((UINT64)x[7]<<56)+((UINT64)x[6]<<48)+((UINT64)x[5]<<40)+((UINT64)x[4]<<32)+
                      ((UINT64)x[3]<<24)+((UINT64)x[2]<<16)+((UINT64)x[1]<< 8)+ (UINT64)x[0]));
    case CDI_LITTLEENDIAN:
      return((UINT64)(((UINT64)x[0]<<56)+((UINT64)x[1]<<48)+((UINT64)x[2]<<40)+((UINT64)x[3]<<32)+
                      ((UINT64)x[4]<<24)+((UINT64)x[5]<<16)+((UINT64)x[6]<< 8)+ (UINT64)x[7]));
    default:
      Error("unhandled endianness %d", HOST_ENDIANNESS);
      return UINT64_C(0xFFFFFFFFFFFFFFFF);
    }
}


size_t binReadF77Block(int fileID, int byteswap)
{
  unsigned char f77block[4];
  size_t blocklen = 0;

  if ( fileRead(fileID, f77block, 4) == 4 )
    {
      if ( byteswap )
	blocklen = get_SUINT32(f77block);
      else
	blocklen =  get_UINT32(f77block);
    }

  return (blocklen);
}


void binWriteF77Block(int fileID, int byteswap, size_t blocksize)
{
  unsigned char f77block[4];

  switch (HOST_ENDIANNESS)
    {
    case CDI_BIGENDIAN:
      if ( byteswap )
	{
	  f77block[0] = (unsigned char) (blocksize);
	  f77block[1] = (unsigned char) (blocksize >>  8);
	  f77block[2] = (unsigned char) (blocksize >> 16);
	  f77block[3] = (unsigned char) (blocksize >> 24);
	}
      else
	{
	  f77block[3] = (unsigned char) (blocksize);
	  f77block[2] = (unsigned char) (blocksize >>  8);
	  f77block[1] = (unsigned char) (blocksize >> 16);
	  f77block[0] = (unsigned char) (blocksize >> 24);
	}
      break;
    case CDI_LITTLEENDIAN:
      if ( byteswap )
	{
	  f77block[3] = (unsigned char) (blocksize);
	  f77block[2] = (unsigned char) (blocksize >>  8);
	  f77block[1] = (unsigned char) (blocksize >> 16);
	  f77block[0] = (unsigned char) (blocksize >> 24);
	}
      else
	{
	  f77block[0] = (unsigned char) (blocksize);
	  f77block[1] = (unsigned char) (blocksize >>  8);
	  f77block[2] = (unsigned char) (blocksize >> 16);
	  f77block[3] = (unsigned char) (blocksize >> 24);
	}
      break;
    default:
      Error("unhandled endianness %d", HOST_ENDIANNESS);
    }

  if ( fileWrite(fileID, f77block, 4) != 4 )
    Error("write failed on %s", fileInqName(fileID));
}


int binReadInt32(int fileID, int byteswap, size_t size, INT32 *ptr)
{
  if ( sizeof(INT32) == 4 )
    {
      fileRead(fileID, (void *) ptr, 4*size);
      if ( byteswap ) swap4byte(ptr, size);
    }
  else
    {
      Error("not implemented for %d byte integer", sizeof(INT32));
    }

  return (0);
}


int binReadInt64(int fileID, int byteswap, size_t size, INT64 *ptr)
{
  if ( sizeof(INT64) == 8 )
    {
      fileRead(fileID, (void *) ptr, 8*size);
      if ( byteswap ) swap8byte(ptr, size);
    }
  else
    {
      Error("not implemented for %d byte integer", sizeof(INT64));
    }

  return (0);
}


int binReadFlt32(int fileID, int byteswap, size_t size, FLT32 *ptr)
{
  if ( sizeof(FLT32) == 4 )
    {
      fileRead(fileID, (void *) ptr, 4*size);
      if ( byteswap ) swap4byte(ptr, size);
    }
  else
    {
      Error("not implemented for %d byte float", sizeof(FLT32));
    }

  return (0);
}


int binReadFlt64(int fileID, int byteswap, size_t size, FLT64 *ptr)
{
  if ( sizeof(FLT64) == 8 )
    {
      fileRead(fileID, (void *) ptr, 8*size);
      if ( byteswap ) swap8byte(ptr, size);
    }
  else
    {
      Error("not implemented for %d byte float", sizeof(FLT64));
    }

  return (0);
}


int binWriteInt32(int fileID, int byteswap, size_t size, INT32 *ptr)
{
  if ( sizeof(INT32) == 4 )
    {
      if ( byteswap ) swap4byte(ptr, size);
      fileWrite(fileID, (void *) ptr, 4*size);
    }
  else
    {
      Error("not implemented for %d byte integer", sizeof(INT32));
    }

  return (0);
}


int binWriteInt64(int fileID, int byteswap, size_t size, INT64 *ptr)
{
  if ( sizeof(INT64) == 8 )
    {
      if ( byteswap ) swap8byte(ptr, size);
      fileWrite(fileID, (void *) ptr, 8*size);
    }
  else
    {
      Error("not implemented for %d byte integer", sizeof(INT64));
    }

  return (0);
}


int binWriteFlt32(int fileID, int byteswap, size_t size, FLT32 *ptr)
{
  if ( sizeof(FLT32) == 4 )
    {
      if ( byteswap ) swap4byte(ptr, size);
      fileWrite(fileID, (void *) ptr, 4*size);
    }
  else
    {
      Error("not implemented for %d byte float", sizeof(FLT32));
    }

  return (0);
}


int binWriteFlt64(int fileID, int byteswap, size_t size, FLT64 *ptr)
{
  if ( sizeof(FLT64) == 8 )
    {
      if ( byteswap ) swap8byte(ptr, size);
      fileWrite(fileID, (void *) ptr, 8*size);
    }
  else
    {
      Error("not implemented for %d byte float", sizeof(FLT64));
    }

  return (0);
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
