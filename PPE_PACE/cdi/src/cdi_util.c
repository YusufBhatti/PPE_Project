#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include "cdi.h"
#include "cdi_int.h"

void cdiDecodeParam(int param, int *pnum, int *pcat, int *pdis)
{
  unsigned uparam = (unsigned)param;
  unsigned upnum;

  *pdis = 0xff   & uparam;
  *pcat = 0xff   & uparam >> 8;
  upnum = 0xffff & uparam >> 16;
  if ( upnum > 0x7fffU ) upnum = 0x8000U - upnum;
  *pnum = (int)upnum;
}


int cdiEncodeParam(int pnum, int pcat, int pdis)
{
  unsigned uparam, upnum;

  if ( pcat < 0 || pcat > 255 ) pcat = 255;
  if ( pdis < 0 || pdis > 255 ) pdis = 255;

  upnum = (unsigned)pnum;
  if ( pnum < 0 ) upnum = (unsigned)(0x8000 - pnum);

  uparam = upnum << 16 | (unsigned)(pcat << 8) | (unsigned)pdis;

  return ((int)uparam);
}


void cdiDecodeDate(int date, int *year, int *month, int *day)
{

  int iyear = date / 10000;
  *year = iyear;
  int idate = abs(date - iyear * 10000),
    imonth = idate / 100;
  *month = imonth;
  *day   = idate - imonth * 100;
}


int cdiEncodeDate(int year, int month, int day)
{
  int iyear = abs(year),
    date = iyear * 10000 + month * 100 + day;
  if ( year < 0 ) date = -date;
  return (date);
}


void cdiDecodeTime(int time, int *hour, int *minute, int *second)
{
  int ihour = time / 10000,
    itime = time - ihour * 10000,
    iminute = itime / 100;
  *hour   = ihour;
  *minute = iminute;
  *second = itime - iminute * 100;
}


int cdiEncodeTime(int hour, int minute, int second)
{
  int time = hour*10000 + minute*100 + second;

  return time;
}


void cdiParamToString(int param, char *paramstr, int maxlen)
{
  int dis, cat, num;
  int len;

  cdiDecodeParam(param, &num, &cat, &dis);

  size_t umaxlen = maxlen >= 0 ? (unsigned)maxlen : 0U;
  if ( dis == 255 && (cat == 255 || cat == 0 ) )
    len = snprintf(paramstr, umaxlen, "%d", num);
  else  if ( dis == 255 )
    len = snprintf(paramstr, umaxlen, "%d.%d", num, cat);
  else
    len = snprintf(paramstr, umaxlen, "%d.%d.%d", num, cat, dis);

  if ( len >= maxlen || len < 0)
    fprintf(stderr, "Internal problem (%s): size of input string is too small!\n", __func__);
}


const char *cdiUnitNamePtr(int cdi_unit)
{
  const char *cdiUnits[] = {
    /*  0 */  "undefined",
    /*  1 */  "Pa",
    /*  2 */  "hPa",
    /*  3 */  "mm",
    /*  4 */  "cm",
    /*  5 */  "dm",
    /*  6 */  "m",
  };
  enum { numUnits = (int) (sizeof(cdiUnits)/sizeof(char *)) };
  const char *name = ( cdi_unit > 0 && cdi_unit < numUnits ) ?
    cdiUnits[cdi_unit] : NULL;
  return name;
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
