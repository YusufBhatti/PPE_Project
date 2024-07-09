#ifdef __xlc__
#pragma comment(user,",,")
#endif
#include <string.h>

const char remote_url[] = "";
const char branch[] = "";
const char revision[] = "";

void repository_url(char *name, int *actual_len)
{
  if (strlen(remote_url) > *actual_len)
    {
      *actual_len = 0;
    }
  else
    {
      strcpy(name, remote_url);
      *actual_len = strlen(name);
    }

  return;
}

void branch_name(char *name, int *actual_len)
{
  if (strlen(branch) > *actual_len)
    {
      *actual_len = 0;
    }
  else
    {
      strcpy(name, branch);
      *actual_len = strlen(name);
    }

  return;
}

void revision_key(char *name, int *actual_len)
{
  if (strlen(revision) > *actual_len)
    {
      *actual_len = 0;
    }
  else
    {
      strcpy(name, revision);
      *actual_len = strlen(name);
    }

  return;
}

