#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined (__APPLE__) && defined (__MACH__)
#include <sys/types.h>
#include <sys/sysctl.h>
#endif
#ifdef __linux
#include <malloc.h>
#include <unistd.h>
#endif
#ifdef _AIX
#include <malloc.h>
#include <sys/systemcfg.h>
#endif

#include "cfortran.h"

void cf_util_set_malloc_alignment(void);
FCALLSCSUB0(cf_util_set_malloc_alignment, UTIL_SET_MALLOC_ALIGNMENT, util_set_malloc_alignment)

void cf_util_mem_zero_double(void *ptr, int count); 
FCALLSCSUB2(cf_util_mem_zero_double, UTIL_MEM_ZERO_DOUBLE, util_mem_zero_double, PVOID, INT)

void cf_util_mem_zero_float(void *ptr, int count); 
FCALLSCSUB2(cf_util_mem_zero_float, UTIL_MEM_ZERO_FLOAT, util_mem_zero_float, PVOID, INT)

void cf_util_mem_zero_long(void *ptr, int count); 
FCALLSCSUB2(cf_util_mem_zero_long, UTIL_MEM_ZERO_LONG, util_mem_zero_long, PVOID, INT)

void cf_util_mem_zero_int(void *ptr, int count); 
FCALLSCSUB2(cf_util_mem_zero_int, UTIL_MEM_ZERO_INT, util_mem_zero_int, PVOID, INT)

const int ldebug = 0; /* false by default, set to 1 to enable */

int get_cachelinesize(int *cachelinesize)
{
  int ret;
#if defined (__APPLE__) && defined (__MACH__)
  int mib[2] = { CTL_HW, HW_CACHELINE }; 
  unsigned int miblen = 2;
  size_t len = 4;

  /* returns cache linesize in bytes */
  if ((ret=sysctl(mib, miblen, cachelinesize, &len, NULL, 0)) != 0)
    {
      perror("get_cachelinesize");
    }
#endif  
#ifdef __linux
  *cachelinesize = (int) sysconf (_SC_LEVEL1_DCACHE_LINESIZE);
  ret = 0;
#endif
#ifdef _AIX
  /* returns cache linesize in bytes */
  *cachelinesize = _system_configuration.dcache_line; 
  ret = 0;
#endif
  return ret;
}

void cf_util_set_malloc_alignment(void)
{
#ifdef _AIX
  static int cacheLineSize = 0;
  int ret;
  if (! cacheLineSize)
    {
      (void) get_cachelinesize(&cacheLineSize);
    }
  if ((ret=mallopt(M_MALIGN, cacheLineSize)) != 0) 
    {
      perror("util_set_malloc_alignment");
    }
  else
    {
      if (ldebug)
	fprintf(stderr,"util_set_malloc_alignment: changed system default memory allocation alignment to L1 data cache line bounds.\n");
    }
#endif
  return;
}

void cf_util_mem_zero_double(void *ptr, int count) 
{
#ifdef _AIX
  static long cacheLineSize = 0;
  char *addr = ptr;
  char *limit;
  long length = (long) (count * sizeof(double));
  
  if (! cacheLineSize)
    {
      int tmp;
      (void) get_cachelinesize(&tmp);
      cacheLineSize = (long) tmp;
    }
  if (ldebug)
    fprintf(stderr, "L1 data cache line size: %d bytes\n", cacheLineSize);
  
  /* zero any initial portion to first cache line boundary */
  if ((long)addr & (cacheLineSize-1)) 
    {
      limit = (char*) (length < cacheLineSize ?
		       ((long)addr+length) :
		       (((long)addr + cacheLineSize) & ~(cacheLineSize-1)));
      if (ldebug)
	fprintf(stderr, "initial portion: %ld ...\n", limit-addr); 
     /* code assumes ptr will have long alignment! */
#pragma unroll(4)
      for (; addr < limit; addr += sizeof(long))
	*((long *)addr) = 0;
    }
  
  limit = (char *) (((long)ptr + length) & ~(cacheLineSize-1));
  if (ldebug)
    fprintf(stderr, "dcbz portion: %ld ...\n", limit-addr);
  for (; addr < limit; addr += cacheLineSize)
    {
      __dcbz(addr);
    }
  /* zero final portion smaller than a cache line */
  limit = (char *) ((long)ptr + length);
  if (addr < limit) {
    if (ldebug)
      fprintf(stderr, "final portion: %ld ...\n", limit-addr);
    /* code assumes length will be a multiple of sizeof(long)! */
#pragma execution_frequency(very_high)
#pragma unroll(4)
    for (; addr < limit; addr += sizeof(long)) {
      *((long *)addr) = 0;
    }
  }

#else

  /* fallback solution, memset used instead */

  static long cacheLineSize = 0;
  char *addr = ptr;
  
  if (! cacheLineSize)
    {
      int tmp;
      (void) get_cachelinesize(&tmp);
      cacheLineSize = (long) tmp;
    }
  if (ldebug)
    fprintf(stderr, "L1 data cache line size: %ld bytes\n", cacheLineSize);
  
  /* zero any initial portion to first cache line boundary */
  if (ldebug)
    if ((long)addr & (cacheLineSize-1))
      { 
	fprintf(stderr, "Array not L1 data cache line size aligned\n");
      }
  
  if (ldebug)
    fprintf(stderr, "memset fallback ...\n");
   
  if (ldebug)
    fprintf(stderr, "memset count: %d\n", count);
  
  memset(ptr, '\0', (size_t) count);
  
#endif
}

void cf_util_mem_zero_float(void *ptr, int count) 
{
#ifdef _AIX
  static long cacheLineSize = 0;
  char *addr = ptr;
  char *limit;
  long length = (long) (count * sizeof(float));
  
  if (! cacheLineSize)
    {
      int tmp;
      (void) get_cachelinesize(&tmp);
      cacheLineSize = (long) tmp;
    }
  if (ldebug)
    fprintf(stderr, "L1 data cache line size: %d bytes\n", cacheLineSize);
  
  /* zero any initial portion to first cache line boundary */
  if ((long)addr & (cacheLineSize-1)) 
    {
      limit = (char*) (length < cacheLineSize ?
		       ((long)addr+length) :
		       (((long)addr + cacheLineSize) & ~(cacheLineSize-1)));
      if (ldebug)
	fprintf(stderr, "initial portion: %ld ...\n", limit-addr); 
     /* code assumes ptr will have long alignment! */
#pragma unroll(4)
      for (; addr < limit; addr += sizeof(long))
	*((long *)addr) = 0;
    }
  
  limit = (char *) (((long)ptr + length) & ~(cacheLineSize-1));
  if (ldebug)
    fprintf(stderr, "dcbz portion: %ld ...\n", limit-addr);
  for (; addr < limit; addr += cacheLineSize)
    {
      __dcbz(addr);
    }
  /* zero final portion smaller than a cache line */
  limit = (char *) ((long)ptr + length);
  if (addr < limit) {
    if (ldebug)
      fprintf(stderr, "final portion: %ld ...\n", limit-addr);
    /* code assumes length will be a multiple of sizeof(long)! */
#pragma execution_frequency(very_high)
#pragma unroll(4)
    for (; addr < limit; addr += sizeof(long)) {
      *((long *)addr) = 0;
    }
  }

#else

  /* fallback solution, memset used instead */

  static long cacheLineSize = 0;
  char *addr = ptr;
  
  if (! cacheLineSize)
    {
      int tmp;
      (void) get_cachelinesize(&tmp);
      cacheLineSize = (long) tmp;
    }
  if (ldebug)
    fprintf(stderr, "L1 data cache line size: %ld bytes\n", cacheLineSize);
  
  /* zero any initial portion to first cache line boundary */
  if (ldebug)
    if ((long)addr & (cacheLineSize-1))
      { 
	fprintf(stderr, "Array not L1 data cache line size aligned\n");
      }
  
  if (ldebug)
    fprintf(stderr, "memset fallback ...\n");
   
  if (ldebug)
    fprintf(stderr, "memset count: %d\n", count);
  
  memset(ptr, '\0', (size_t) count);
  
#endif
}

void cf_util_mem_zero_long(void *ptr, int count) 
{
#ifdef _AIX
  static long cacheLineSize = 0;
  char *addr = ptr;
  char *limit;
  long length = (long) (count * sizeof(int64_t));
  
  if (! cacheLineSize)
    {
      int tmp;
      (void) get_cachelinesize(&tmp);
      cacheLineSize = (long) tmp;
    }
  if (ldebug)
    fprintf(stderr, "L1 data cache line size: %d bytes\n", cacheLineSize);
  
  /* zero any initial portion to first cache line boundary */
  if ((long)addr & (cacheLineSize-1)) 
    {
      limit = (char*) (length < cacheLineSize ?
		       ((long)addr+length) :
		       (((long)addr + cacheLineSize) & ~(cacheLineSize-1)));
      if (ldebug)
	fprintf(stderr, "initial portion: %ld ...\n", limit-addr); 
     /* code assumes ptr will have long alignment! */
#pragma unroll(4)
      for (; addr < limit; addr += sizeof(long))
	*((long *)addr) = 0;
    }
  
  limit = (char *) (((long)ptr + length) & ~(cacheLineSize-1));
  if (ldebug)
    fprintf(stderr, "dcbz portion: %ld ...\n", limit-addr);
  for (; addr < limit; addr += cacheLineSize)
    {
      __dcbz(addr);
    }
  /* zero final portion smaller than a cache line */
  limit = (char *) ((long)ptr + length);
  if (addr < limit) {
    if (ldebug)
      fprintf(stderr, "final portion: %ld ...\n", limit-addr);
    /* code assumes length will be a multiple of sizeof(long)! */
#pragma execution_frequency(very_high)
#pragma unroll(4)
    for (; addr < limit; addr += sizeof(long)) {
      *((long *)addr) = 0;
    }
  }

#else

  /* fallback solution, memset used instead */

  static long cacheLineSize = 0;
  char *addr = ptr;
  
  if (! cacheLineSize)
    {
      int tmp;
      (void) get_cachelinesize(&tmp);
      cacheLineSize = (long) tmp;
    }
  if (ldebug)
    fprintf(stderr, "L1 data cache line size: %ld bytes\n", cacheLineSize);
  
  /* zero any initial portion to first cache line boundary */
  if (ldebug)
    if ((long)addr & (cacheLineSize-1))
      { 
	fprintf(stderr, "Array not L1 data cache line size aligned\n");
      }
  
  if (ldebug)
    fprintf(stderr, "memset fallback ...\n");
   
  if (ldebug)
    fprintf(stderr, "memset count: %d\n", count);
  
  memset(ptr, '\0', (size_t) count);
  
#endif
}

void cf_util_mem_zero_int(void *ptr, int count) 
{
#ifdef _AIX
  static long cacheLineSize = 0;
  char *addr = ptr;
  char *limit;
  long length = (long) (count * sizeof(int32_t));
  
  if (! cacheLineSize)
    {
      int tmp;
      (void) get_cachelinesize(&tmp);
      cacheLineSize = (long) tmp;
    }
  if (ldebug)
    fprintf(stderr, "L1 data cache line size: %d bytes\n", cacheLineSize);
  
  /* zero any initial portion to first cache line boundary */
  if ((long)addr & (cacheLineSize-1)) 
    {
      limit = (char*) (length < cacheLineSize ?
		       ((long)addr+length) :
		       (((long)addr + cacheLineSize) & ~(cacheLineSize-1)));
      if (ldebug)
	fprintf(stderr, "initial portion: %ld ...\n", limit-addr); 
     /* code assumes ptr will have long alignment! */
#pragma unroll(4)
      for (; addr < limit; addr += sizeof(long))
	*((long *)addr) = 0;
    }
  
  limit = (char *) (((long)ptr + length) & ~(cacheLineSize-1));
  if (ldebug)
    fprintf(stderr, "dcbz portion: %ld ...\n", limit-addr);
  for (; addr < limit; addr += cacheLineSize)
    {
      __dcbz(addr);
    }
  /* zero final portion smaller than a cache line */
  limit = (char *) ((long)ptr + length);
  if (addr < limit) {
    if (ldebug)
      fprintf(stderr, "final portion: %ld ...\n", limit-addr);
    /* code assumes length will be a multiple of sizeof(long)! */
#pragma execution_frequency(very_high)
#pragma unroll(4)
    for (; addr < limit; addr += sizeof(long)) {
      *((long *)addr) = 0;
    }
  }

#else

  /* fallback solution, memset used instead */

  static long cacheLineSize = 0;
  char *addr = ptr;
  
  if (! cacheLineSize)
    {
      int tmp;
      (void) get_cachelinesize(&tmp);
      cacheLineSize = (long) tmp;
    }
  if (ldebug)
    fprintf(stderr, "L1 data cache line size: %ld bytes\n", cacheLineSize);
  
  /* zero any initial portion to first cache line boundary */
  if (ldebug)
    if ((long)addr & (cacheLineSize-1))
      { 
	fprintf(stderr, "Array not L1 data cache line size aligned\n");
      }
  
  if (ldebug)
    fprintf(stderr, "memset fallback ...\n");
   
  if (ldebug)
    fprintf(stderr, "memset count: %d\n", count);
  
  memset(ptr, '\0', (size_t) count);
  
#endif
}

#ifdef __MEMZERO_TEST__
/* test program */
int main(int argc, char *argv[])
{
  double *d;
  int i;

  cf_util_set_malloc_alignment();

  d = (double *) malloc(sizeof(double)*1048576);

  for (i = 0; i < 1048576; i++)
    {
      d[i] = 123.4;
    }
  
  cf_util_mem_zero_double(d, 1048576);

  for (i = 0; i < 1048576; i++)
    {
      if (d[i] != 0.0)
	{
	  fprintf(stderr, "%d %f\n", i, d[i]);
	}
    }
  
  free(d);

  return 0;
}
#endif
