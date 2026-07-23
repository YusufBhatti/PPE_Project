#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include <yaxt.h>

#include "dmemory.h"
#include "error.h"
#include "resource_handle.h"

#include "cdipio.h"
#include "pio_conf.h"

static int cdiPioConfCompareP(void *conf1, void *conf2);
static void cdiPioConfDestroyP(void *conf);
static void cdiPioConfPrintP(void *cdiPioConfPtr, FILE * fp);

struct intCodeStrMap
{
  const char text[28];
  int code;
};

static size_t
mapSearchStr(const struct intCodeStrMap map[], size_t mapSize,
             const char *str)
{
  size_t retval = SIZE_MAX;
  for (size_t i = 0; i < mapSize; ++i)
    if (!strcmp(str, map[i].text))
      {
        retval = i;
        break;
      }
  return retval;
}

static size_t
mapSearchCode(const struct intCodeStrMap map[], size_t mapSize,
              int code)
{
  size_t retval = SIZE_MAX;
  for (size_t i = 0; i < mapSize; ++i)
    if (code == map[i].code)
      {
        retval = i;
        break;
      }
  return retval;
}

static const struct intCodeStrMap modeMap[] = {
  { "PIO_NONE", PIO_NONE },
  { "PIO_MPI", PIO_MPI },
  { "PIO_FPGUARD", PIO_FPGUARD },
  { "PIO_ASYNCH", PIO_ASYNCH },
  { "PIO_WRITER", PIO_WRITER },
  { "PIO_MPI_FW_ORDERED", PIO_MPI_FW_ORDERED },
  { "PIO_MPI_FW_AT_ALL", PIO_MPI_FW_AT_ALL },
  { "PIO_MPI_FW_AT_REBLOCK", PIO_MPI_FW_AT_REBLOCK },
};

int
cdiPioStr2IOMode(const char *modeStr)
{
  size_t idx = mapSearchStr(modeMap, sizeof (modeMap) / sizeof (modeMap[0]),
                            modeStr);
  int mode = (idx != SIZE_MAX) ? modeMap[idx].code : -1;
  return mode;
}

const char *
cdiPioIOMode2Str(int IOMode)
{
  if (IOMode < PIO_MINIOMODE || IOMode > PIO_MAXIOMODE)
    return "";
  else
    return modeMap[IOMode].text;
}

static const struct intCodeStrMap roleMap[] = {
  { "PIO_ROLE_CLIENT", PIO_ROLE_CLIENT },
  { "PIO_ROLE_COLLECTOR", PIO_ROLE_COLLECTOR },
  { "PIO_ROLE_WRITER", PIO_ROLE_WRITER },
  { "PIO_ROLE_WRITER_COLLECTOR", PIO_ROLE_WRITER_COLLECTOR },
  { "PIO_ROLE_FPGUARD", PIO_ROLE_FPGUARD },
};

#if 0
static int
cdiPioStr2CSRole(const char *roleStr)
{
  size_t idx = mapSearchStr(roleMap, sizeof (roleMap) / sizeof (roleMap[0]),
                            roleStr);
  int role = (idx != SIZE_MAX) ? roleMap[idx].code : -1;
  return role;
}
#endif

static const char *
cdiPioCSRole2Str(int role)
{
  size_t pos = mapSearchCode(roleMap, sizeof (roleMap) / sizeof (roleMap[0]),
                             role);
  const char *roleStr = (pos != SIZE_MAX) ? roleMap[pos].text : "";
  return roleStr;
}

const resOps cdiPioConfOps = {
  cdiPioConfCompareP,
  cdiPioConfDestroyP,
  cdiPioConfPrintP,
  /* serialization of configuration is not supported */
  0, 0, 0
};

static int
cdiPioConfCompareP(void *p1, void *p2)
{
  struct cdiPioConf *a = p1, *b = p2;
  return (a->IOMode != b->IOMode)
    | (a->clientServerRole != b->clientServerRole)
    | (a->partInflate != b->partInflate)
    | (a->postCommSetupActions != b->postCommSetupActions);
}

static void
cdiPioConfDestroyP(void *conf)
{
  Free(conf);
}

static void
cdiPioConfPrintP(void *cdiPioConfPtr, FILE * fp)
{
  struct cdiPioConf *conf = cdiPioConfPtr;
  const char *iomodeStr = cdiPioIOMode2Str(conf->IOMode),
    *CSRoleStr = cdiPioCSRole2Str(conf->clientServerRole);
  if (!iomodeStr[0])
    iomodeStr = "(invalid!)";
  if (!CSRoleStr[0])
    CSRoleStr = "(invalid!)";
  fprintf(fp, "configuration object %p\n"
          "IOMode = %s\n"
          "client/server = %s\n"
          "part data imbalance = %f\n"
          "aligning of block buffers to large pages is %sabled\n"
          "record aggregation buffer size %zu\n"
          "callback after setup of communication = %p\n",
          cdiPioConfPtr, iomodeStr, CSRoleStr, conf->partInflate,
          conf->largePageAlign ? "en" : "dis",
          (size_t)conf->recordAggBufLimMB * 1024 * 1024,
          (void *)conf->postCommSetupActions);
}


int cdiPioConfCreate(void)
{
  struct cdiPioConf *conf = Malloc(sizeof (*conf));
  conf->IOMode = PIO_NONE;
  conf->clientServerRole = PIO_ROLE_CLIENT;
  conf->partInflate = 1.1f;
  conf->postCommSetupActions = cdiPioNoPostCommSetup;
  conf->largePageAlign = false;
  conf->cacheRedists = true;
  conf->recordAggBufLimMB = 128;
  conf->xmap_new = xt_xmap_dist_dir_new;
  conf->stripify = true;
  int resH = reshPut(conf, &cdiPioConfOps);
  /* configuration objects are never forwarded */
  reshSetStatus(resH, &cdiPioConfOps,
                reshGetStatus(resH, &cdiPioConfOps) & ~RESH_SYNC_BIT);
  return resH;
}

void cdiPioConfDestroy(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  cdiPioConfDestroyP(conf);
  reshRemove(confResH, &cdiPioConfOps);
}

void cdiPioConfSetPartInflate(int confResH, float partInflate)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  xassert(partInflate >= 1.0f);
  conf->partInflate = partInflate;
}

float cdiPioConfGetPartInflate(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return conf->partInflate;
}

void cdiPioConfSetIOMode(int confResH, int IOMode)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  conf->IOMode = IOMode;
}

int cdiPioConfGetIOMode(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return conf->IOMode;
}

void cdiPioConfSetCSRole(int confResH, int CSRole)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  conf->clientServerRole = CSRole;
}

int cdiPioConfGetCSRole(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return conf->clientServerRole;
}

void
cdiPioConfSetPostCommSetupActions(int confResH,
                                  void (*postCommSetupActions)(void))
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  conf->postCommSetupActions = postCommSetupActions;
}

void (*cdiPioConfGetPostCommSetupActions(int confResH))(void)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return conf->postCommSetupActions;
}

void cdiPioConfSetLargePageAlign(int confResH, int largePageAlign)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  conf->largePageAlign = largePageAlign != 0;
}

int cdiPioConfGetLargePageAlign(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return conf->largePageAlign;
}

void cdiPioConfSetRedistCache(int confResH, int doCache)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  conf->cacheRedists = doCache;
}

int cdiPioConfGetRedistCache(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return conf->cacheRedists;
}

void cdiPioConfSetRecordAggBufLim(int confResH, int lim_mb)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  if (lim_mb > 0)
    conf->recordAggBufLimMB = (unsigned)lim_mb;
  else
    Error("unexpected negative buffer size value %d requested", lim_mb);
}

int cdiPioConfGetRecordAggBufLim(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return (int)conf->recordAggBufLimMB;
}

void cdiPioConfSetXmapNew(int confResH, xmap_new_func_ptr xmap_new)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  conf->xmap_new = xmap_new;
}

xmap_new_func_ptr cdiPioConfGetXmapNew(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return conf->xmap_new;
}

void cdiPioConfSetStripeConversion(int confResH, int doConversion)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  conf->stripify = doConversion;
}

int cdiPioConfGetStripeConversion(int confResH)
{
  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);
  return conf->stripify;
}

