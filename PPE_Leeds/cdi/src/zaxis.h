#ifndef _ZAXIS_H
#define _ZAXIS_H

void zaxisGetTypeDescription(int zaxisType, int* outPositive, const char** outName, const char** outLongName, const char** outStdName, const char** outUnit);  //The returned const char* point to static storage. Don't free or modify them.

unsigned cdiZaxisCount(void);

void cdiZaxisGetIndexList(unsigned numIDs, int *IDs);

void
zaxisUnpack(char * unpackBuffer, int unpackBufferSize,
            int * unpackBufferPos, int originNamespace, void *context,
            int force_id);

void zaxisDefLtype2(int zaxisID, int ltype2);

const resOps *getZaxisOps(void);

#endif

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
