#include <stdlib.h>

#include "cfortran.h"

/****************************************************************************/
/* cfortran prototypes:                                                     */

void cf_util_exit(int exit_no);
FCALLSCSUB1(cf_util_exit, UTIL_EXIT, util_exit, 
	    INT)

void cf_util_abort(void);
FCALLSCSUB0(cf_util_abort, UTIL_ABORT, util_abort)

int cf_util_system(char *s);
FCALLSCFUN1(INT, cf_util_system, UTIL_SYSTEM, util_system,
	    PSTRING)

#ifdef __XT3__ 
void cf_util_base_iobuf(void);
FCALLSCSUB0(cf_util_base_iobuf, UTIL_BASE_IOBUF, util_base_iobuf) 
#endif

/****************************************************************************/

#ifdef __XT3__ 
static int bsize=65536;
static char *obuf, *ebuf;

void cf_util_base_iobuf(void)
{
  int stat1, stat2;

  obuf = malloc(bsize);
  ebuf = malloc(bsize);

  stat1=setvbuf( stdout, obuf, _IOFBF, bsize );
  stat2=setvbuf( stderr, ebuf, _IOFBF, bsize );

  return;
} 
#endif

void cf_util_exit(int exit_no)
{
  exit(exit_no);
}

void cf_util_abort(void)
{
  exit(1);
}

int cf_util_system(char *s)
{
#ifdef __XT3__ 
  return(-1);  
#else
  return(system(s));
#endif
}
