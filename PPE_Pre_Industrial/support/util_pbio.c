/*
 * This files contains the binary C-IO interface to read and write unblocked
 * files originally provided through the EMOS library from ECMWF. 
 *
 * Authors: L. Kornblueh      Max-Planck-Institute for Meteorology, Hamburg
 *
 * Date:   30.5.1998
 *
 * Changes:
 *
 * Uwe Schulzweida 24.1.2007
 *     Set default buffersize to 4*blksize
 *
 * $Id: util_pbio.c,v 1.4 1998/12/11 13:23:58 m214003 Exp $
 *
 */

/* LK, 20070525 does not work ... needs to be checked */
/* #ifdef CRAY */
/* #define FFIO */
/* #endif */

#ifndef EMOS


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef FFIO
#include <ffio.h>
#include <fcntl.h>
#endif

#ifndef HAVE_VALLOC
#define valloc malloc
#endif


#ifdef PTHREADS
#include <pthread.h>
#endif

#include "cfortran.h"

/****************************************************************************/
/* cfortran prototypes:                                                     */

void cf_pbopen(int *unit, char *name, char *mode, int *iret);
FCALLSCSUB4(cf_pbopen, PBOPEN, pbopen,
	    PINT, STRING, STRING, PINT)
void cf_pbseek(int *unit, int *offset, int *whence, int *iret);
FCALLSCSUB4(cf_pbseek, PBSEEK, pbseek,
	    PINT, PINT, PINT, PINT)
void cf_pbread(int *unit, void *buffer, int *nbytes, int *iret);
FCALLSCSUB4(cf_pbread, PBREAD, pbread,
	    PINT, PVOID, PINT, PINT)
void cf_pbwrite(int *unit, void *buffer, int *nbytes, int *iret);
FCALLSCSUB4(cf_pbwrite, PBWRITE, pbwrite,
	    PINT, PVOID, PINT, PINT)
void cf_pbclose(int *unit, int *iret);
FCALLSCSUB2(cf_pbclose, PBCLOSE, pbclose,
	    PINT, PINT)
void cf_pbflush(int *unit);
FCALLSCSUB1(cf_pbflush, PBFLUSH, pbflush,
	    PINT)

/****************************************************************************/

static FILE** fptable = NULL;
static int fptableSize = 0;
#ifdef PTHREADS
static pthread_mutex_t fpTableBusy = PTHREAD_MUTEX_INITIALIZER;
#endif
FILE * pbfp(long index) {
  if( (fptable == NULL) || ((int)index < 0) || ((int)index >= fptableSize) )
    return (FILE *) NULL;
  else
    return fptable[index];
}

/*
// Default buffer size for I/O operations (set via setvbuf)
*/
static long size = 0;
static int sizeSet = 0;
static char * envSize;
static char** fileBuffer = NULL;

/*
// Debug flags.
*/
#define DEBUGOFF 1
#define DEBUG (debugSet > DEBUGOFF )
static char * debugLevel;
static int debugSet = 0;
 
 
#define CURRENT_FILE (fptable[*unit])

void cf_pbopen(int *unit, char *name, char *mode, int *iret)
{
/*
// Purpose:
//  Opens file, returns the index of a UNIX FILE pointer held in 
//  an internal table (fptable).
//
// First time through, reads value in environment variable PBIO_BUFSIZE
// (if it is set) and uses it as the size to be used for internal file
// buffers; the value is passed to setvbuf. If PBIO_BUFSIZE is not set,
// a default value is used.
//
// Function  accepts:
//    name = filename
//    mode = r, w
//
//    Note: l1 and l2 are the lengths of the FORTRAN character strings
//          in name and mode.
//
// Function returns:
//   INTEGER iret:
//     -1 = Could not open file.
//     -2 = Invalid file name.
//     -3 = Invalid open mode specified
//      0 = OK.
*/
  int n;
  char *p;
  char flags[4];
#ifdef FFIO
  int fflags=0;
#endif

/*
// See if DEBUG switched on.
*/
    if( ! debugSet ) {
      debugLevel = getenv("PBIO_DEBUG");
      if( debugLevel == NULL )
        debugSet = DEBUGOFF;              /* off */
      else {
        int loop;
        for( loop = 0; loop < strlen(debugLevel) ; loop++ ) {
          if( ! isdigit(debugLevel[loop]) ) {
            printf("Invalid number string in PBIO_DEBUG: %s\n", debugLevel);
            printf("PBIO_DEBUG must comprise only digits [0-9].\n");
            debugSet = DEBUGOFF;
          }
        }
        debugSet = DEBUGOFF + atol( debugLevel );
      }
      if( DEBUG ) printf("PBIO_PBOPEN: debug switched on\n");
    }

  strcpy (flags, "");

  *unit = 0;
  *iret = 0;
                
  if( DEBUG ) printf("PBIO_PBOPEN: filename = %s\n", name);

/* build open flags from "modes" */

  p = mode;
  while (*p && (strlen(flags) < 3 )) {
    switch (*p) {
    case 'a':
    case 'A':
      strcat (flags, "a");
#ifdef FFIO
      fflags = O_RDWR | O_APPEND;
#endif
      break;
    case 'c':
    case 'C':
    case 'w':
    case 'W':
      strcat (flags, "w");
#ifdef FFIO
      fflags = O_WRONLY | O_CREAT | O_TRUNC;
#endif
      break;
    case 'r':
    case 'R':
      strcat (flags, "r");
#ifdef FFIO
      fflags = O_RDONLY;
#endif
      break;
    default:
      *iret = -3;
      return;
    }
    p++;
  }

/* if read/write change flags */

  if (!strcmp (flags, "wr") || !strcmp (flags, "rw")) {
    strcpy (flags, "r+w" );
#ifdef FFIO
    fflags = O_RDWR;
#endif
  }

  if( DEBUG ) printf("PBIO_PBOPEN: file open mode = %s\n", flags);

#ifdef FFIO
  *unit = ffopen (name, fflags, 0644);
  
  if (*unit < 0) {
    perror (name);
    *iret = -1;
  }
#else

/*
// Look for a free slot in fptable.
// (Create the table the first time through).
*/
#ifdef PTHREADS
/*
// Wait if another thread opening a file
*/
    pthread_mutex_lock(&fpTableBusy);
#endif

    n = 0;
    if( fptableSize == 0 ) {
      int i;
      fptableSize = 2;
      fptable = (FILE **) malloc(fptableSize*sizeof(FILE *));
      if( fptable == NULL ) {
        perror("Unable to allocate space for table of FILE pointers");
        exit(1);
      }

      fileBuffer = (char **) malloc(fptableSize*sizeof(char *));
      if( fileBuffer == NULL ) {
        perror("Unable to allocate space for FILE buffers");
        exit(1);
      }

      for( i = 0; i < fptableSize; i++ ) {
        fptable[i] = 0;
        fileBuffer[i] = NULL;
      }
    }
    else {
      while( n < fptableSize ) {
        if(fptable[n]==0) {
          *unit = n;
          break;
        }
        n++;
      }
    }
/*
// If the table overflows, double its size.
*/
    if( n == fptableSize) {
      int i;
      fptableSize = 2*fptableSize;
      fptable = (FILE **) realloc(fptable, fptableSize*sizeof(FILE *));
      if( fptable == NULL ) {
        perror("Unable to reallocate space for table of FILE pointers");
        exit(1);
      }
      n = fptableSize/2;

      fileBuffer = (char **) realloc(fileBuffer, fptableSize*sizeof(char *));
      if( fileBuffer == NULL ) {
        perror("Unable to allocate space for FILE buffers");
        exit(1);
      }

      n = fptableSize/2;
      for( i = n; i < fptableSize; i++ ) {
        fptable[i] = 0;
        fileBuffer[i] = NULL;
      }

      *unit = n;
    }

    if( DEBUG ) printf("PBIO_PBOPEN: fptable slot = %d\n", *unit);

    fptable[n] = fopen(name, flags );

    if(fptable[n] == NULL) {
      perror(name);
      *iret = -1;

#ifdef PTHREADS
      pthread_mutex_unlock(&fpTableBusy);
#endif
      return;
    }

/*
// Now allocate a buffer for the file, if necessary.
*/
    if( ! sizeSet ) {
      envSize = getenv("PBIO_BUFSIZE");
      if( envSize )
	{
	  int loop;
	  for( loop = 0; loop < strlen(envSize) ; loop++ ) {
	    if( ! isdigit(envSize[loop]) ) {
	      printf("Invalid number string in PBIO_BUFSIZE: %s\n", envSize);
	      printf("PBIO_BUFSIZE must comprise only digits [0-9].\n");
	      exit(1);
	    }
	  }
	  size = atol( envSize );
	}
      else
	{
	  struct stat filestat;
	  if ( stat(name, &filestat) != 0 )
	    perror(name);
	  size = 4*filestat.st_blksize;
	}
      if( size < 0 ) {
        printf("Invalid buffer size in PBIO_BUFSIZE: %s\n", envSize);
        printf("Buffer size defined by PBIO_BUFSIZE must be positive.\n");
        exit(1);
      }
      sizeSet = 1;
    }

    if( DEBUG ) printf("PBIO_PBOPEN: file buffer size = %ld\n", size);

    if( fileBuffer[n] == NULL ) {
      fileBuffer[n] = (char *) malloc(size);
      if( setvbuf(CURRENT_FILE, fileBuffer[*unit], _IOFBF, size) ) {
        perror("setvbuf failed");
        *iret = -1;
      }
    }

#ifdef PTHREADS
    pthread_mutex_unlock(&fpTableBusy);
#endif

#endif /* FFIO */

  return;
}

void cf_pbseek(int *unit, int *offset, int *whence, int *iret)
{
/*
//
// Purpose:
//   Seeks to a specified location in file.
//
//  Function  accepts:
//    unit = the index of a UNIX FILE pointer held in
//           an internal table (fptable).
//
//    offset = byte count
//
//    whence  = 0, from start of file
//            = 1, from current position
//            = 2, from end of file.  
//
//  Returns:
//    iret:
//      -2 = error in handling file,
//      -1 = end-of-file
//      otherwise,  = byte offset from start of file.
*/
  int my_offset = *offset;
  int my_whence = *whence;

/* must use negative offset if working from end-of-file */

  if ( *whence == 2 ) my_offset = - labs(my_offset);

#ifdef FFIO
  *iret = ffseek(*unit, my_offset, *whence);
  if (*iret == -1) *iret = -2;
#else

/*
// Must use negative offset if working from end-of-file
*/
  if( DEBUG ) { 
    printf("PBIO_PBSEEK: fptable slot = %d\n", *unit);
    printf("PBIO_PBSEEK: Offset = %d\n", my_offset);
    printf("PBIO_PBSEEK: Type of offset = %d\n", my_whence);
  }

  if( my_whence == 2 ) my_offset = - abs(my_offset);

  *iret = ftell(CURRENT_FILE);
  if( DEBUG ) printf("PBIO_PBSEEK: current position = %d\n", *iret);
  if( *iret == my_offset && my_whence == 0)
    *iret = 0;
  else
    *iret = fseek(CURRENT_FILE, my_offset, my_whence);

  if( DEBUG ) printf("PBIO_PBSEEK: fileSeek return code = %d\n",*iret);

  if( *iret != 0 ) {
    if( ! feof(CURRENT_FILE) ) {
      *iret = -2;             /* error in file-handling */
      perror("pbseek");
    }
    else
      *iret = -1;             /* end-of-file  */

    clearerr(CURRENT_FILE);
    return;
  }

/*
// Return the byte offset from start of file
*/
    *iret = ftell(CURRENT_FILE);

    if( DEBUG )
      printf("PBIO_PBSEEK: byte offset from start of file = %d\n",*iret);

#endif

  return;
}

void cf_pbread(int *unit, void *buffer, int *nbytes, int *iret)
{
/*
// Purpose:
//  Reads a block of bytes from a file..
//
//  Function  accepts:
//    unit = the index of a UNIX FILE pointer held in
//           an internal table (fptable).
//
//    nbytes = number of bytes to read.
//
//  Returns:
//    iret:
//      -2 = error in reading file,
//      -1 = end-of-file,
//      otherwise, = number of bytes read.
*/
 
#ifdef FFIO
  *iret = ffread(*unit,buffer,*nbytes);
  if(*iret != *nbytes) *iret = -1;
#else
  if( DEBUG ) {
    printf("PBIO_READ: fptable slot = %d. ", *unit);
    printf("Number of bytes to read = %d\n", *nbytes);
  }

  if( (*iret = fread(buffer, 1, *nbytes, CURRENT_FILE) ) != *nbytes) {
    if( ! feof(CURRENT_FILE) ) {
      *iret = -2;             /*  error in file-handling  */
      perror("pbread");
      clearerr(CURRENT_FILE);
      return;
    }
    else {
      *iret = -1;             /*  end-of-file */
      clearerr(CURRENT_FILE);
    }
  }

  if( DEBUG ) {
    printf("PBIO_READ: fptable slot = %d. ", *unit);
    printf("Number of bytes read = %d\n", *nbytes);
  }
#endif

  return;
}

void cf_pbwrite(int *unit, void *buffer, int *nbytes, int *iret)
{
/*
// Purpose:
//  Writes a block of bytes to a file.
//
//  Function  accepts:
//    unit = the index of a UNIX FILE pointer held in
//           an internal table (fptable).
//
//    nbytes = number of bytes to write.
//
//  Returns:
//    iret:
//      -1 = Could not write to file.
//     >=0 = Number of bytes written.
*/
#ifdef FFIO
  *iret = ffwrite(*unit,buffer,*nbytes);
  if(*iret != *nbytes) *iret = -1;
#else
  if( DEBUG ) {
    printf("PBIO_WRITE: fptable slot = %d. ", *unit);
    printf("Number of bytes to write = %d\n", *nbytes);
  }

  if( (*iret = fwrite(buffer, 1, *nbytes, CURRENT_FILE) ) != *nbytes) {
    perror("pbwrite");
    *iret = -1;
  }

  if( DEBUG ) {
    printf("PBIO_WRITE: fptable slot = %d. ", *unit);
    printf("PBIO_WRITE: number of bytes written = %d\n", *iret);
  }
#endif

  return;
}

void cf_pbclose(int *unit, int *iret)
{
/*
// Purpose:
//  Closes file.
//
//  Function  accepts:
//    unit = the index of a UNIX FILE pointer held in
//           an internal table (fptable).
////  Returns:
//    iret:
//      0 = OK.
//      otherwise = error in handling file.
*/
  if( DEBUG )
    printf("PBIO_CLOSE: fptable slot = %d\n", *unit);
#ifdef FFIO
  *iret = ffclose(*unit);
#else
  *iret = fclose(CURRENT_FILE);
#endif
  
  if (*iret != 0) perror("pbclose");
  CURRENT_FILE = 0;

  return;
}

/*
//------------------------------------------------------------------------
//  PBFLUSH - flush (from FORTRAN)
//------------------------------------------------------------------------
*/
void cf_pbflush(int *unit)
{
/*
// Purpose:	Flushes file.
*/
  if( DEBUG )
    printf("PBIO_FLUSH: fptable slot = %d\n", *unit);
#ifdef FFIO
  ffflush(*unit);
#else
  fflush(CURRENT_FILE);
#endif
}

#else
#define UNUSED
#endif
