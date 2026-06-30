#ifndef CDI_DATETIME_H
#define CDI_DATETIME_H

typedef struct
{
  long date;
  long time;
}
DateTime;

static inline int
datetimeCmp(DateTime dt1, DateTime dt2)
{
  return 2 * ((dt1.date > dt2.date) - (dt1.date < dt2.date))
    + (dt1.time > dt2.time) - (dt1.time < dt2.time);
}


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
