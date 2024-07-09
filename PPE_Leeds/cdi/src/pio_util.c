#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <unistd.h>

#include "pio_util.h"
#include "cdi.h"
#include "dmemory.h"

static
char commands[][13] = { "FINALIZE",
                        "RESOURCES",
                        "WRITETS"};

/*****************************************************************************/
void
cdiAbortC_MPI(const char *caller, const char *filename,
              const char *functionname, int line,
              const char *errorString, va_list ap)
{
  {
    int rank = getMPICommWorldRank();
    fprintf(stderr, "ERROR, pe%d in %s, %s, line %d%s"
            "%s\nerrorString: \"", rank, functionname, filename, line,
            caller?", called from ":"", caller?caller:"");
  }
  vfprintf(stderr, errorString, ap);
  fputs("\"\n", stderr);
  if (callsToMPIAreAllowed())
    MPI_Abort(MPI_COMM_WORLD, 1);
  else
    abort();
  exit(EXIT_FAILURE);
  va_end(ap);
}

void cdiPioWarning(const char *caller, const char *fmt, va_list ap)
{
  int rank = getMPICommWorldRank();
  fprintf(stderr, "pe%d: Warning (%s) : ", rank, caller);
  vfprintf(stderr, fmt, ap);
  fputc('\n', stderr);
}

/*****************************************************************************/

/***************************************************************/

void pcdiXMPI(int iret, const char *filename, int line)
{
  char errorString[2][MPI_MAX_ERROR_STRING + 1];
  int len, errorClass, rank = getMPICommWorldRank();
  MPI_Error_class(iret, &errorClass);
  MPI_Error_string(errorClass, errorString[0], &len);
  errorString[0][len] = '\0';
  MPI_Error_string(iret, errorString[1], &len);
  errorString[1][len] = '\0';
  fprintf(stderr, "MPI ERROR, pe%d, %s, line %d,"
          "errorClass: \"%s\""
          "errorString: \"%s\"\n",
          rank, filename, line,
          errorString[0], errorString[1]);
  MPI_Abort(MPI_COMM_WORLD, 1);
}

/*****************************************************************************/

void pcdiXMPIStat ( int iret, const char *filename, int line, MPI_Status *status )
{
  char errorString[MPI_MAX_ERROR_STRING + 1];
  int len, rank = getMPICommWorldRank();

  if ( iret == MPI_ERR_IN_STATUS )
    {
      fprintf ( stderr, "------- checking error in request ----------\n" );
      switch ( status->MPI_ERROR )
        {
        case MPI_SUCCESS :
          fprintf ( stderr, "-------- mpi_success -----------\n" );
          break;
        case MPI_ERR_PENDING:
          fprintf ( stderr, "-------- mpi_err_pending ----------\n");
          break;
        default:
          MPI_Error_string ( status->MPI_ERROR, errorString, &len );
          errorString[len] = '\0';
          fprintf ( stderr,"MPI ERROR in request, pe%d, %s, line %d,"
                    "return value: %d, error_string: %s\n",
                    rank, filename, line, iret, errorString );
          MPI_Abort ( MPI_COMM_WORLD, iret );
        }
    }
  else
    xmpi ( iret );

  return;
}

/****************************************************/

void pcdiDebugMsg ( const char * cdiPioDebugString, const char *filename,
                    const char *functionname, int line, int tag, int source,
                    int nfinished )
{
  int rank = getMPICommWorldRank();

  fprintf ( stdout,
            "%s pe%d in %s, %s, line %d: command %s, source %d, finalized=%d\n",
            cdiPioDebugString, rank, functionname, filename, line,
            &commands[tag][0], source, nfinished );
}

/****************************************************/

void pcdiDebugMsg2 ( const char *filename, const char *functionname, int line,
                   int tag, int source, char * text )
{
  int rank = getMPICommWorldRank();

  fprintf ( stdout,
            "pe%d in %s, %s, line %d: command %s, source %d, %s\n",
            rank, functionname, filename, line,
            &commands[tag][0], source, text );
}

/****************************************************/

void printArray ( const char *cdiPioDebugString, const char *ps, const void *array, int n,
                  int datatype, const char * funname, const char * filename, int line )
{
  int i;
  int * iArray;
  double * dArray;

  {
    int rank = getMPICommWorldRank();
    fprintf ( stdout, "%s pe%d in %s, %s, line %d: %s = ",
              cdiPioDebugString, rank, funname, filename, line, ps );
  }

  switch ( datatype )
    {
    case DATATYPE_INT:
      iArray = ( int * ) array;
      for ( i = 0; i < n-1; i++ )
	fprintf ( stdout, "%d ", * ( iArray + i ));
      fprintf ( stdout, "%d\n", * ( iArray + n - 1 ));
      break;
    case DATATYPE_FLT:
      dArray = ( double * ) array;
      for ( i = 0; i < n-1; i++ )
	fprintf ( stdout, "%.2f ", * ( dArray + i ));
      fprintf ( stdout, "%.2f\n", * ( dArray + n-1 ));
      break;
    default:
      fprintf ( stdout, " ... no datatype defined\n" );
    }

  return;
}

int
cdiPioQueryVarDims(int varShape[3], int vlistID, int varID)
{
  int gridID = vlistInqVarGrid(vlistID, varID);
  int zaxisID = vlistInqVarZaxis(vlistID, varID);
  int gridType = gridInqType(gridID);
  switch (gridType)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
      varShape[0] = gridInqXsize(gridID);
      varShape[1] = gridInqYsize(gridID);
      break;
    case GRID_SPECTRAL:
      varShape[0] = 2;
      varShape[1] = gridInqSize(gridID) / 2;
      break;
    case GRID_GENERIC:
    case GRID_LCC:
    case GRID_GME:
    case GRID_CURVILINEAR:
    case GRID_UNSTRUCTURED:
      xabort("unimplemented grid type: %d", gridType);
      break;
    }
  varShape[2] = zaxisInqSize(zaxisID);
  /* FIXME: other grids have different dimensionality */
  return (varShape[2] == 1)?2:3;
}

size_t
cdiPioGetPageSize(bool largePageAlign)
{
  long pagesize = -1L;
  bool nameAssigned = false;
#if HAVE_DECL__SC_LARGE_PAGESIZE || HAVE_DECL__SC_PAGE_SIZE || HAVE_DECL__SC_PAGESIZE
  int name;
#if HAVE_DECL__SC_LARGE_PAGESIZE
  if (largePageAlign)
    {
      name = _SC_LARGE_PAGESIZE;
      nameAssigned = true;
    }
  else
#endif
    {
#if HAVE_DECL__SC_PAGESIZE || HAVE_DECL__SC_PAGE_SIZE
      name =
#if HAVE_DECL__SC_PAGESIZE
        _SC_PAGESIZE
#elif HAVE_DECL__SC_PAGE_SIZE
        _SC_PAGE_SIZE
#endif
        ;
      nameAssigned = true;
#endif
    }
  if (nameAssigned)
    pagesize = sysconf(name);
#endif
  if (pagesize == -1L)
    pagesize =
#if HAVE_DECL_PAGESIZE
      PAGESIZE
#elif HAVE_DECL_PAGE_SIZE
      PAGE_SIZE
#else
      commonPageSize
#endif
      ;
  return (size_t)pagesize;
}



void
deco1D_CCP(size_t nelems, const size_t *restrict weightPfxSums,
           size_t nparts, size_t *restrict separators)
{
  separators[0] = 0;
  separators[nparts] = nelems;
  size_t i = 0, k = 1;
  size_t weightTotal = weightPfxSums[nelems];
  while (k < nparts)
    {
      size_t target = k * weightTotal / nparts;
      do {
        ++i;
      } while (i < nelems && weightPfxSums[i] < target);
      separators[k] = i;
      ++k;
      --i;
    }
  /* todo: implement h2 and dp+ algorithms from
   * A. Pınar, C. Aykanat
   * Fast optimal load balancing algorithms for 1D partitioning
   */
}



/****************************************************/
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
