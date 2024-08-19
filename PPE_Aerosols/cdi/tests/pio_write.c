#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#ifdef USE_MPI
#include <mpi.h>
#include <yaxt.h>
#else
typedef int MPI_Comm;
#endif

#include "cdi.h"
#include "dmemory.h"
#include "pio_write.h"
#ifdef USE_MPI
#include "cdipio.h"
#include "pio_util.h"
#endif

static const struct model_config default_setup
#ifdef __cplusplus
= { 12, 6, 3, 5, 5, FILETYPE_GRB, DATATYPE_PACK24,1,"grb"};
#else
  = { .nlon = 12, .nts = 3, .nlat = 6, .nvars = 5,
      .filetype = FILETYPE_GRB, .datatype = DATATYPE_PACK24,
      .compute_checksum = 1,
      .suffix = "grb",
      .max_nlev = 5,
};
#endif

static const struct {
  char suffix[4];
  int type, defaultDT, defaultGrid;
} suffix2type[] = {
  { "nc", FILETYPE_NC, DATATYPE_FLT64, GRID_LONLAT },
  { "grb",  FILETYPE_GRB, DATATYPE_PACK24, GRID_LONLAT },
  { "nc2", FILETYPE_NC2, DATATYPE_FLT64, GRID_LONLAT },
  { "nc4", FILETYPE_NC4, DATATYPE_FLT64, GRID_LONLAT },
  { "ext", FILETYPE_EXT, DATATYPE_FLT64, GRID_GENERIC, },
  { "svc", FILETYPE_SRV, DATATYPE_FLT64, GRID_GENERIC, },
  { "ieg", FILETYPE_IEG, DATATYPE_FLT64, GRID_LONLAT },
};

static void
invalidOptionDie(const char *format, ...)
{
  va_list ap;
  va_start(ap, format);
  vfprintf(stderr, format, ap);
  exit(EXIT_FAILURE);
}

static int
parse_intarg(const char msg[])
{
  char *end;
  long temp = strtol(optarg, &end, 0);
  if ((errno == ERANGE && (temp == LONG_MAX || temp == LONG_MIN))
      || (errno != 0 && temp == 0)) {
    perror(msg);
    exit(EXIT_FAILURE);
  }
  if (temp > INT_MAX || temp < INT_MIN)
    invalidOptionDie("range error: %ld\n", temp);
  return (int)temp;
}

static unsigned
parse_unsignedarg(const char msg[])
{
  char *end;
  unsigned long temp = strtoul(optarg, &end, 0);
  if ((errno == ERANGE && (temp == ULONG_MAX))
      || (errno != 0 && temp == 0))
    {
      perror(msg);
      exit(EXIT_FAILURE);
    }
  if (temp > UINT_MAX)
    invalidOptionDie("range error: %ld\n", temp);
  return (unsigned)temp;
}

typedef int (*pioRoleFunc)(MPI_Comm commSuper, int IOMode, int nProcsIO);

static void
parse_long_option(int pioConfHandle, pioRoleFunc *pioRoleAssign,
                  const char *str)
{
  static const char cacheRedistStr[] = "cache-redists",
    pioRoleSchemeOptionStr[] = "pio-role-scheme";
  if (!strncmp(str, cacheRedistStr, sizeof (cacheRedistStr) - 1))
    {
#ifdef USE_MPI
      if (str[sizeof (cacheRedistStr) - 1] == '\0'
          || !strcmp(str + sizeof (cacheRedistStr) - 1, "=true"))
        cdiPioConfSetRedistCache(pioConfHandle, 1);
      else if (!strcmp(str + sizeof (cacheRedistStr) - 1, "=false"))
        cdiPioConfSetRedistCache(pioConfHandle, 0);
      else
        invalidOptionDie("invalid option argument to -qcache-redists: %s\n",
                         optarg + sizeof (cacheRedistStr) - 1);
#else
      invalidOptionDie("CDI-PIO option -q%s ignored in non-MPI mode\n",
                       cacheRedistStr);
#endif
    }
  else if (!strncmp(str, pioRoleSchemeOptionStr,
                    sizeof (pioRoleSchemeOptionStr) - 1))
    {
#ifdef USE_MPI
      static const char pioRoleSchemeLastN[]="last",
        pioRoleSchemeFirstN[]="first", pioRoleSchemeBalanced[]="balanced";
      if (str[sizeof (pioRoleSchemeOptionStr) - 1] == '=')
        {
          const char *optargarg = str + sizeof (pioRoleSchemeOptionStr);
          if (!strcmp(optargarg, pioRoleSchemeLastN))
            *pioRoleAssign = cdiPioCSRLastN;
          else if (!strcmp(optargarg, pioRoleSchemeFirstN))
            *pioRoleAssign = cdiPioCSRFirstN;
          else if (!strcmp(optargarg, pioRoleSchemeBalanced))
            *pioRoleAssign = cdiPioCSRBalanced;
          else
            invalidOptionDie("unknown scheme argument \"%s\" to -q%s",
                             optargarg, pioRoleSchemeOptionStr);
        }
      else
        invalidOptionDie("long option %s needs argument\n",
                         pioRoleSchemeOptionStr);
#else
      invalidOptionDie("CDI-PIO option -q%s ignored in non-MPI mode\n",
                       pioRoleSchemeOptionStr);
#endif
    }
  else
    invalidOptionDie("unknown long option: %s\n", str);
}


int main(int argc, char *argv[])
{
  struct model_config setup = default_setup;

  MPI_Comm commModel;
  int pioConfHandle = 0;
  pioRoleFunc pioRoleAssign = 0;
#ifdef USE_MPI
  MPI_Comm commGlob;
  int sizeGlob;
  int rankGlob;
  int IOMode = PIO_MPI;
  int nProcsIO = 2;

  xmpi ( MPI_Init ( &argc, &argv));
  commGlob = MPI_COMM_WORLD;
  xt_initialize(commGlob);
  xmpi ( MPI_Comm_set_errhandler ( commGlob, MPI_ERRORS_RETURN ));
  xmpi ( MPI_Comm_size ( commGlob, &sizeGlob ));
  xmpi ( MPI_Comm_rank ( commGlob, &rankGlob ));

  pioConfHandle = cdiPioConfCreate();
  pioRoleAssign = cdiPioCSRLastN;
#endif

  /* seed random generator */
  unsigned randSeed;
  {
#ifdef USE_MPI
    if (rankGlob == 0)
#endif
      {
        struct timeval tv;
        int status = gettimeofday(&tv, NULL);
        if (status != 0)
          {
            perror("failed seed generation!");
            exit(1);
          }
        randSeed = (unsigned)(tv.tv_sec ^ tv.tv_usec);
      }
  }

  {
    int opt;
    while ((opt = getopt(argc, argv, "f:m:n:z:t:y:cs:q:"
#ifdef USE_MPI
                         "p:w:"
#endif
                         )) != -1)
      switch (opt) {
#ifdef USE_MPI
      case 'p':
        IOMode = cdiPioStr2IOMode(optarg);
        if (IOMode < 0)
          {
            fprintf(stderr, "Unsupported PIO mode requested: %s\n", optarg);
            exit(EXIT_FAILURE);
          }
        break;
      case 'w':
        {
          long temp = strtol(optarg, NULL, 0);
          if (temp < 0 || temp > INT_MAX/2)
            {
              fprintf(stderr, "Unsupported number of I/O servers: %ld\n", temp);
              exit(EXIT_FAILURE);
            }
          nProcsIO = (int)temp;
        }
        break;
#endif
      case 'q':
        parse_long_option(pioConfHandle, &pioRoleAssign, optarg);
        break;
      case 'f':
        {
          int found = 0;
          for (size_t i = 0;
               i < sizeof (suffix2type) / sizeof (suffix2type[0]);
               ++i)
            if (!strcmp(optarg, suffix2type[i].suffix))
              {
                found = 1;
                setup.filetype = suffix2type[i].type;
                setup.suffix = suffix2type[i].suffix;
                setup.datatype = suffix2type[i].defaultDT;
                break;
              }
          if (!found)
            {
              fprintf(stderr, "Unsupported format requested: %s\n", optarg);
              exit(EXIT_FAILURE);
            }
        }
        break;
      case 'm':
        setup.nlon = parse_intarg("error parsing number of longitudes");
        break;
      case 'n':
        setup.nlat = parse_intarg("error parsing number of latitudes");
        break;
      case 'y':
        setup.nvars = parse_intarg("error parsing number of variables");
        if (setup.nvars < 1)
          {
            fputs("number of levels must be greater than zero!\n",
                  stderr);
            exit(EXIT_FAILURE);
          }
        if (setup.nvars > 127)
          {
            fputs("number of variables must not exceed 127!\n",
                  stderr);
            exit(EXIT_FAILURE);
          }
        break;
      case 'z':
        setup.max_nlev = parse_intarg("error parsing number of levels");
        if (setup.max_nlev < 1)
          {
            fputs("number of levels must be greater than zero!\n",
                  stderr);
            exit(EXIT_FAILURE);
          }
        break;
      case 't':
        setup.nts = parse_intarg("error parsing number of timesteps");
        break;
      case 'c':
        setup.compute_checksum = 0;
        break;
      case 's':
        randSeed = parse_unsignedarg("error parsing random seed");
        break;
      default: /* '?' */
        fprintf(stderr, "Usage: %s "
                "[-m nlon] [-n nlat] [-z nlev] [-t nts] [-y num_vars]"
#ifdef USE_MPI
                " [-p PIO_MODE] [-w NIOSERVERS] [-c]"
#endif
                "\n", argv[0]);
        exit(EXIT_FAILURE);
      }

  }

#ifdef USE_MPI
  MPI_Bcast(&randSeed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
#endif
  srandom(randSeed);
#ifdef USE_MPI
  if (rankGlob == 0)
#endif
    fprintf(stderr, "random seed=%u\n", randSeed);

#ifdef USE_MPI
  int pioNamespace;
  cdiPioConfSetIOMode(pioConfHandle, IOMode);
  cdiPioConfSetCSRole(pioConfHandle, pioRoleAssign(commGlob, IOMode,
                                                   nProcsIO));
  cdiPioConfSetPartInflate(pioConfHandle, 1.0f);
  int initNamespace = namespaceGetActive();
  commModel = cdiPioInit(commGlob, pioConfHandle, &pioNamespace);
  if (commModel != MPI_COMM_NULL)
    {
      namespaceSetActive(pioNamespace);
#else
      commModel = -1;
#endif

      modelRun (setup, commModel);

#ifdef USE_MPI
    }
  namespaceSetActive(initNamespace);
  cdiPioConfDestroy(pioConfHandle);
  pioFinalize ();
  xt_finalize();
  MPI_Finalize ();
#endif
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
