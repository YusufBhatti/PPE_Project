/*
 * Debian etch requires the following define
 * for posix_memalign.
 * It is not required for lenny and Suse 10.2 !
 */
#define _XOPEN_SOURCE 600

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <unistd.h>

/* #ifndef NOMPI */
/* #include "mpi.h" */
/* #endif */

typedef enum 
{ 
/*   ALLOC_mpi,  */
  ALLOC_aligned, 
  ALLOC_std 
} ALLOC_Type;

/* description of each single array entry */

struct buffer_entry_t {
  int size;
  int dims;
  double *start;
  double *buffer_position;
};

/* description of a single buffer pool per I/O server */

struct buffer_pool_t {
  int entries;
  int last_entry;
  int size;
  struct buffer_entry_t *entry;
  double *buffer;
};

/* global variables */

int handler = 0;

int pools = 0;
int max_pools = 32;

struct buffer_pool_t *pool = NULL;

/*****************************************************************************/

int init_buffer_pool(int entries) 
{
  int i;

  /* should combine malloc and realloc here to realloc ... */

  if (pool == NULL) 
    {
      pool = (struct  buffer_pool_t *) 
	malloc(max_pools*sizeof(struct  buffer_pool_t));
    }

  pools++;
  if (pools > max_pools) 
    {
      max_pools *= 2;
      pool = (struct  buffer_pool_t *) 
	realloc(pool, max_pools*sizeof(struct  buffer_pool_t)); 
    }

  pools++;
  handler++;

  pool[handler].entries = entries;
  pool[handler].last_entry = -1;
  pool[handler].size = 0;
  pool[handler].entry = (struct buffer_entry_t *) 
    malloc(entries*sizeof(struct buffer_entry_t));
  pool[handler].buffer = NULL; 

  for(i = 0; i < entries; i++) 
    {
      pool[handler].entry[i].size = 0;
      pool[handler].entry[i].dims = 0;
      pool[handler].entry[i].start = NULL;
    }
  
  return(handler);
}

/*****************************************************************************/

int setup_buffer_pool(int pd)
{
  int i;
  double *current;

  current = pool[pd].buffer;
  for (i = 0; i < pool[pd].entries; i++)
    {
      pool[pd].entry[i].buffer_position = current;
      current += pool[pd].entry[i].size;
    }

  return (0);
}

/*****************************************************************************/

int allocate_buffer_pool(int pd, int allocation_type)
{
  int status;
  size_t pagesize;
/* #ifndef NOMPI */
/*   char error_string[MPI_MAX_ERROR_STRING+1]; */
/*   int len; */
/* #endif */
  pagesize = (size_t) sysconf(_SC_PAGESIZE);

  switch (allocation_type)
    {
/*     case ALLOC_mpi: */
/* #ifndef NOMPI */
/*       status = MPI_Alloc_mem((MPI_Aint)(pool[pd].size*sizeof(int)),  */
/* 			     MPI_INFO_NULL, &(pool[pd].buffer)); */
/*       if (status != MPI_SUCCESS) */
/* 	{ */
/* 	  MPI_Error_string(status, error_string, &len); */
/* 	  error_string[len] = '\0'; */
/* 	  fprintf(stderr,"MPI error %4d: %s\n\n", status, error_string); */
/* 	} */
/* #endif */
/*       break; */
    case ALLOC_aligned:
#ifdef _SX
      pool[pd].buffer = malloc( (pool[pd].size*sizeof(double)));
      if ( pool[pd].buffer == NULL ) {
        status = -1;
      }
      else {
        status = 0;
      }
#else
      status = posix_memalign((void **) &(pool[pd].buffer), 
			      pagesize, (pool[pd].size*sizeof(double)));
#endif
      if (status != 0) 
	{
	  switch (status)
	    {
	    case EINVAL:
	      fprintf(stderr, "The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
	      break;
	    case ENOMEM:
	      fprintf(stderr, "There was insufficient memory to fulfill the allocation request.\n");
	      break;
	    }
	}
      break;
    case ALLOC_std:
      pool[pd].buffer = (double *) 
	malloc(pool[pd].size*sizeof(double));
      if (pool[pd].buffer == NULL) 
	{
	  fprintf(stderr, "There was insufficient memory to fulfill the allocation request.\n");
	}
      break;
    default:
      fprintf(stderr,"allocate_buffer_pool: allocation type not supported.\n");
      exit (1);
    }

  
  return (0);
}

/*****************************************************************************/

int copy_chunk_to_buffer_pool(int pd, int chunk_no, double *chunk)
{
  memcpy(pool[pd].entry[chunk_no].buffer_position,
	 pool[pd].entry[chunk_no].start,
	 pool[pd].entry[chunk_no].size*sizeof(double));

  return (0);
}

/*****************************************************************************/

int buffer_pool_entry(int pd, int size, int dims, double *chunk) 
{
  int n;
	  
  pool[pd].last_entry++;
  if (pool[pd].last_entry+1 > pool[pd].entries) 
    {
      fprintf(stderr,
	      "buffer_pool_entry: something went wrong adding a new entry.\n");
      exit(1);
    }

  pool[pd].size += size;

  n = pool[pd].last_entry;

  pool[pd].entry[n].size = size;
  pool[pd].entry[n].dims = dims;
  pool[pd].entry[n].start = chunk;

/*   fprintf(stderr,  */
/* 	  "register: handler: %3d entry: %4d start: 0x%p chunk: 0x%p\n",  */
/* 	  pd, n, pool[pd].entry[n].start, chunk); */

  return (pool[pd].last_entry);
}

/*****************************************************************************/

#ifdef _UTIL_BUFFER_POOL

void print_test(void)
{
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[handler].buffer[    0], pool[handler].buffer[ 1023]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[handler].buffer[ 1024], pool[handler].buffer[ 2047]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[handler].buffer[ 2048], pool[handler].buffer[ 3071]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[handler].buffer[ 3072], pool[handler].buffer[ 4095]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[handler].buffer[ 4096], pool[handler].buffer[ 5119]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[handler].buffer[ 5120], pool[handler].buffer[ 6143]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[handler].buffer[ 6144], pool[handler].buffer[ 7167]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[handler].buffer[ 7168], pool[handler].buffer[ 8191]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[handler].buffer[ 8192], pool[handler].buffer[ 9215]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[handler].buffer[ 9216], pool[handler].buffer[10239]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[handler].buffer[10240], pool[handler].buffer[11263]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[handler].buffer[11264], pool[handler].buffer[12287]);
}

int main (int argc, char *argv[])
{
  int pd;
  int i, n;

  double **chunk;

  chunk = (double **) malloc(12*sizeof(double *));
  for (i = 0; i < 12; i++)
    {
      chunk[i] = (double *) malloc(1024*sizeof(double));
      for (n = 0; n < 1024; n++)
	{
	  chunk[i][n] = (double) i + 0.0001* (double) (n + 1);
	}
    }

  /* initialize basic buffer with 12 entries */

  pd = init_buffer_pool(12);
  
  /* add entries into buffer pool */

  for (i = 0; i < 12; i++) 
    {
      buffer_pool_entry(pd, chunk[i], 1024, 2);
      fprintf(stderr,"entry: %4d %5d %4d 0x%p\n", 
	      i, 
	      pool[pd].size, 
	      pool[pd].entry[i].size, 
	      pool[pd].entry[i].start); 
    }

  fprintf(stderr,"Buffer size for all chunks: %d elements.\n", 
	  pool[pd].size);

  /* allocate buffer */

  allocate_buffer_pool(pd, alloc_std);

  /* set pointer entries into buffer pool */

  setup_buffer_pool(pd);

  /* copy data to buffer pool */

  for (i = 0; i < 12; i++) 
    {
      copy_chunk_to_buffer_pool(pd, i, chunk[i]);
    }

  /* check buffer */

  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[pd].buffer[    0], pool[pd].buffer[ 1023]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[pd].buffer[ 1024], pool[pd].buffer[ 2047]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[pd].buffer[ 2048], pool[pd].buffer[ 3071]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[pd].buffer[ 3072], pool[pd].buffer[ 4095]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[pd].buffer[ 4096], pool[pd].buffer[ 5119]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[pd].buffer[ 5120], pool[pd].buffer[ 6143]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[pd].buffer[ 6144], pool[pd].buffer[ 7167]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[pd].buffer[ 7168], pool[pd].buffer[ 8191]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[pd].buffer[ 8192], pool[pd].buffer[ 9215]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[pd].buffer[ 9216], pool[pd].buffer[10239]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[pd].buffer[10240], pool[pd].buffer[11263]);
  fprintf(stderr, "%8.4f %8.4f\n", 
	  pool[pd].buffer[11264], pool[pd].buffer[12287]);

  /* no cleanup yet */

  return (0);
}

#endif

