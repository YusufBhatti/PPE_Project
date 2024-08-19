#ifndef _DMEMORY_H
#define _DMEMORY_H

#include <stdio.h>

// if DEBUG_MEMORY is defined setenv MEMORY_DEBUG to debug memory
#define  DEBUG_MEMORY

#ifndef  WITH_FUNCTION_NAME
#define  WITH_FUNCTION_NAME
#endif

extern size_t  memTotal(void);
extern void    memDebug(int debug);
extern void    memExitOnError(void);

#if  defined  DEBUG_MEMORY

extern void   *memRealloc(void *ptr, size_t size, const char *file, const char *functionname, int line);
extern void   *memCalloc (size_t nmemb, size_t size, const char *file, const char *functionname, int line);
extern void   *memMalloc (size_t size, const char *file, const char *functionname, int line);
extern void    memFree   (void *ptr, const char *file, const char *functionname, int line);

#if  defined  WITH_FUNCTION_NAME
#  define  Realloc(p, s)  memRealloc((p), (s), __FILE__, __func__, __LINE__)
#  define   Calloc(n, s)   memCalloc((n), (s), __FILE__, __func__, __LINE__)
#  define   Malloc(s)      memMalloc((s), __FILE__, __func__, __LINE__)
#  define     Free(p)        memFree((p), __FILE__, __func__, __LINE__)
#else
#  define  Realloc(p, s)  memRealloc((p), (s), __FILE__, (void *) NULL, __LINE__)
#  define   Calloc(n, s)   memCalloc((n), (s), __FILE__, (void *) NULL, __LINE__)
#  define   Malloc(s)      memMalloc((s), __FILE__, (void *) NULL, __LINE__)
#  define     Free(p)        memFree((p), __FILE__, (void *) NULL, __LINE__)
#endif

#endif /* DEBUG_MEMORY */

#endif /* _DMEMORY_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
