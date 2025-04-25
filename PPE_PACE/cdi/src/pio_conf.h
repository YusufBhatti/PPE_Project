#ifndef PIO_CONF_H
#define PIO_CONF_H


#include <stdbool.h>

#include <yaxt.h>

/*
 * declare data structures and functions to manipulate/query CDI-PIO
 * configuration object.
 */

#include "resource_handle.h"

typedef Xt_xmap (*xmap_new_func_ptr)(Xt_idxlist src_idxlist,
                                     Xt_idxlist dst_idxlist,
                                     MPI_Comm comm);

/*
 * cdiPioConf is meant to be internal to the library and not to be
 * used directly if possible, CDI-PIO users should rely on the
 * corresponding functions instead.
 */
struct cdiPioConf {
  int IOMode;
  int clientServerRole;
  float partInflate;
  unsigned recordAggBufLimMB;
  void (*postCommSetupActions)(void);
  xmap_new_func_ptr xmap_new;
  bool largePageAlign;
  bool cacheRedists;
  bool stripify;
};

extern const resOps cdiPioConfOps;

int cdiPioConfCreate(void);

void cdiPioConfDestroy(int confResH);

void cdiPioConfSetPartInflate(int confResH, float partInflate);

float cdiPioConfGetPartInflate(int confResH);

void cdiPioConfSetIOMode(int confResH, int IOMode);

int cdiPioConfGetIOMode(int confResH);

void cdiPioConfSetCSRole(int confResH, int CSRole);

int cdiPioConfGetCSRole(int confResH);

void
cdiPioConfSetpostCommSetupActions(int confResH,
                                  void (*postCommSetupActions)(void));

void (*cdiPioConfGetpostCommSetupActions(int confResH))(void);

int
cdiPioIOModeStr2Code(const char *modeStr);

const char *
cdiPioIOMode2Str(int IOMode);

int cdiPioConfGetLargePageAlign(int confResH);

void cdiPioConfSetLargePageAlign(int confResH, int largePageAlign);

void cdiPioConfSetRedistCache(int confResH, int doCache);

int cdiPioConfGetRedistCache(int confResH);

void cdiPioConfSetRecordAggBufLim(int confResH, int lim_mb);

int cdiPioConfGetRecordAggBufLim(int confResH);

void cdiPioConfSetXmapNew(int confResH, xmap_new_func_ptr xmap_new);

xmap_new_func_ptr cdiPioConfGetXmapNew(int confResH);

void cdiPioConfSetStripeConversion(int confResH, int doConversion);

int cdiPioConfGetStripeConversion(int confResH);

#endif
