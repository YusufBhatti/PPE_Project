#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <pwd.h>

#include <sys/types.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#elif HAVE_SYS_UNISTD_H
#include <sys/unistd.h>
#endif

#include <sys/types.h>
#include <sys/utsname.h>

#ifdef _AIX
#include <sys/systemcfg.h>
#endif

#include "cfortran.h"

/****************************************************************************/
/* cfortran prototypes:                                                     */

void cf_util_os_system(char *name, int *actual_len);
FCALLSCSUB2(cf_util_os_system, UTIL_OS_SYSTEM, util_os_system,
            PSTRING, PINT)
void cf_util_user_name(char *name, int *actual_len);
FCALLSCSUB2(cf_util_user_name, UTIL_USER_NAME, util_user_name,
            PSTRING, PINT)
void cf_util_node_name(char *name, int *actual_len);
FCALLSCSUB2(cf_util_node_name, UTIL_NODE_NAME, util_node_name,
            PSTRING, PINT)

/****************************************************************************/

/* funcion implemetations */ 


void cf_util_user_name(char *name, int *actual_len)
{
  struct passwd *current;
 
#if defined (__XT3__) || defined (__CRAYXC)
  current = NULL;
#else
  current = getpwuid(getuid());
#endif
  if (current == NULL) {
#if defined (__XT3__) || defined (__CRAYXC)
    strcpy(name, "user name not available");
#else
    strcpy(name, "unknown user name");
#endif
  } else {
    if (strlen(current->pw_name) == 0) {
      strcpy(name, "unknown user name");
    } else {
      strcpy(name, current->pw_gecos);
      strcat(name, " (");
      strcat(name, current->pw_name);
      strcat(name, ")");
    }
  }    

  *actual_len = strlen(name);
}

void cf_util_os_system(char *name, int *actual_len)
{
  struct utsname system_info;
  int iret;
#ifdef _AIX
  int imp;
#endif

  iret = uname(&system_info);

  /* Cray doesn't handle operating system information proper */

#ifdef _UNICOSMP
  strcpy(name, "UNICOS");
#else
  strcpy(name, system_info.sysname);
#endif
  strcat(name, " ");
#ifdef _AIX
  strcat(name, system_info.version);
  strcat(name, ".");
  strcat(name, system_info.release);
#else
  strcat(name, system_info.release);
#endif
  strcat(name, " ");
#ifdef _AIX
  imp = _system_configuration.implementation;
  if (imp == POWER_4)
    strcat(name, "Power4");
  else if (imp == POWER_5)
    strcat(name, "Power5");
  else if (imp == POWER_6)
    strcat(name, "Power6");
  else if (imp == POWER_7)
    strcat(name, "Power7");
  else {
    strcat(name, "unknown");
  }
#else
  strcat(name, system_info.machine);
#endif

  *actual_len = strlen(name);
}

void cf_util_node_name(char *name, int *actual_len)
{
  struct utsname host_info;
  int iret;
  char *fqdn;

  iret = uname(&host_info);
  if (iret != EFAULT) 
    {
      strcpy (name, host_info.nodename);
      
    } 
  else 
    {
      strcpy(name, "unknown");
    }

  /* if FQDN is returned remove domainname */
  if ((fqdn = index(name,'.')) != NULL) 
    {
      *fqdn = '\0';
    }
  
  *actual_len = strlen(name);

  return;
}

#ifdef _UNIT_TEST

int main(int argc, char *argv[])
{
  char hostname[72];
  char osname[72];
  char username[72];
  int len;

  cf_util_node_name(hostname, &len);

  fprintf(stderr,"hostname: <%s> len: %d\n", hostname, len);

  cf_util_os_system(osname, &len);

  fprintf(stderr,"osname: <%s> len: %d\n", osname, len);

  cf_util_user_name(username, &len);

  fprintf(stderr,"username: <%s> len: %d\n", username, len);

  return 0;
}

#endif
