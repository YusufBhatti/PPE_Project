#include <assert.h>
#include <stdlib.h>

#include "cdi.h"

int main()
{
  assert(tableInqParName(-1, -1, NULL) != 0);
  return EXIT_SUCCESS;
}
