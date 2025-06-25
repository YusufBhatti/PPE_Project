/*
 *
 * Generate backtrace, need to get more working solutions. The current once
 * working are:
 *
 *     Linux/glibc  gcc/gfortran
 *     Linux/glibc  gcc/NAG f95
 *     MacOSX       gcc/gfortran
 *
 * Luis Kornblueh, MPIM, January 2008
 *
 * Luis Kornblueh, MPIM, April 2010, added MacOSX
 *
 */

#include "config.h"

#include "cfortran.h"

/*
 * cfortran prototypes:
 */

void cf_util_backtrace(void);
FCALLSCSUB0(cf_util_backtrace, UTIL_BACKTRACE, util_backtrace) 

#if defined (HAVE_EXECINFO_H) && defined (__GNUC__)

#include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>
#if defined (__APPLE__) && defined (__MACH__)
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <mach-o/dyld.h>
#endif

#if defined (__linux)

void cf_util_backtrace(void)
{
  void *callstack[32];
  int frames;
  char **symbols;
  int i;

  frames = backtrace(callstack, 32);
  symbols = backtrace_symbols(callstack, frames);

  for (i = 0; i < frames; i++) {
    fprintf(stderr,"%s\n", symbols[i]);
  }

  fprintf(stderr, "\nUse addr2line for address to line number conversion.\n\n");

  free(symbols);

  return;
}

#elif defined (__APPLE__) && defined (__MACH__)

static char *pipeGets(char *buf, size_t  buflen, FILE  *f)  
{
  char *p = buf;
  size_t len = buflen;
  size_t datalen = 0;
  
  *buf = '\0';
  while (datalen < (buflen - 1)) 
    {
      fgets(p, len, f);
      if (feof(f)) break;
      if (ferror(f) && (errno != EINTR)) break;
      if (strchr(buf, '\n')) break;
      datalen = strlen(buf);
      p = buf + datalen;
      len -= datalen;
    }
  
  return (buf[0] ? buf : NULL);
}

static void GetNameOfAndPathToThisProcess(char *executable)
{
  char path[PATH_MAX];
  uint32_t size = sizeof(path);

  if (_NSGetExecutablePath(path, &size) == 0)
    {
      realpath(path, executable);
      /* fprintf(stderr, "programs name: %s\n", executable);   */
    }
  else
    {
      printf("buffer too small; need size %u\n", size);
    }

  return;
}

void cf_util_backtrace(void)
{
  void *callstack[128];
  int frames;
  char **symbols;
  int i;

  char saved_env[128], *env;

  char pipeBuf[PATH_MAX];
  char pathToThisProcess[PATH_MAX];

#ifdef __LP64__
  char arch[] = "x86_64";
#else
  char arch[] = "i386";
#endif

  FILE *f;

  /* fprintf(stderr,"arch: %s\n", arch); */

  env = getenv("NSUnbufferedIO");
  if (env) {
    strlcpy(saved_env, env, sizeof(saved_env));
  }
  setenv("NSUnbufferedIO", "YES", 1);

  GetNameOfAndPathToThisProcess(pathToThisProcess);

  /* fprintf(stderr,"executable: %s\n", pathToThisProcess); */

  frames = backtrace (callstack, 128);
  symbols = backtrace_symbols (callstack, frames);

  snprintf(pipeBuf, sizeof(pipeBuf), "/usr/bin/atos -o \"%s\" -arch \"%s\"", pathToThisProcess, arch);

  fprintf(stderr, "\n");
  f = popen(pipeBuf, "r+");
  if (f) {
    setbuf(f, 0);
    for (i = 0; i < frames; i++) {
      fprintf(f, "%#lx\n", (long)callstack[i]);
      /* fprintf(stderr, "%#lx\n", (long)callstack[i]); */
      pipeGets(pipeBuf, sizeof(pipeBuf), f);
#ifdef __LP64__
      fprintf(stderr, "%3d  0x%016llx  %s", i, (unsigned long long)callstack[i], pipeBuf);
#else
      fprintf(stderr, "%3d  0x%08lx  %s", i, (unsigned long)callstack[i], pipeBuf);
#endif
    }
    pclose(f);
  }
  fprintf(stderr, "\n");
  
  if (env) {
    setenv("NSUnbufferedIO", saved_env, 1);
  } else {
    unsetenv("NSUnbufferedIO");
  }

  free (symbols);

  return;
}
     
#else

void cf_util_backtrace(void)
{
  fprintf (stderr, "Stack traceback not available.\n");

  return;
}

#endif

#elif  defined (HAVE_UCONTEXT_H)

#include <ucontext.h>

#if (defined(sun) || defined(__sun)) && (defined(__SVR4) || defined(__svr4__))

void cf_util_backtrace(void)
{
  (void) printstack(2); /* fd 2 is always stderr */

  return;
}

#else

void cf_util_backtrace(void)
{
  fprintf (stderr, "Stack traceback not available.\n");

  return;
}

#endif
 
#else

void cf_util_backtrace(void)
{
  fprintf (stderr, "Stack traceback not available.\n");

  return;
}

#endif
