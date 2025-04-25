#define _POSIX_C_SOURCE 200112L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <unistd.h>

#include <sys/types.h>
#include <sys/stat.h>

#define PID_LEN 8

int util_islink(char *path)
{
  char *buf;
  ssize_t len;
  long max_buf_len;
  int iret;

  iret = -1;

  max_buf_len = pathconf("/", _PC_NAME_MAX);

  buf = (char *) malloc(max_buf_len*sizeof(char));
  
  if ((len = readlink(path, buf, max_buf_len-1)) != -1)
    {
      buf[len] = '\0';
      /* file is a link ..., return C true */
      iret = 1;
    }
  else
    {
      if (errno == EINVAL)
	{
	  /* file is not a link ... */
	  iret = 0;
	}
    }

  free(buf);
  
  /* something else went wrong ... */
  return iret;
}

int util_tmpnam_len(void)
{
  return (L_tmpnam + 1 + PID_LEN);
}

int util_tmpnam(char *filename)
{
  char *ptr;
  char pid_string[PID_LEN + 1]; 
  pid_t pid;

  pid = getpid();

  sprintf(pid_string, "%8.8ld", (long) pid);

  ptr = tmpnam(filename);
  
  strcat(filename, pid_string); 	  

  return strlen(filename);
}

int util_filesize(char *filename)
{
  struct stat statbuf;

  if (stat(filename, &statbuf) == -1)
    {
      return 0;
    }

  return ((int) statbuf.st_size);
}
