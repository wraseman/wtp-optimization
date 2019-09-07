/* fdoublef.c -- September 16, 1993 */
#include "wtp.h"

int fdoublef(FILE *fp, const char *fmt, double x)
/*
* Purpose: Print a double value if greater than or equal to zero else
*          print  spaces.
*
* Inputs:  Simplar to fprintf.
*    fp  = File pointer to disk file or output device.
*    fmt = Format string. The format string should start with %w.pf where
*          'w' is the width and 'p' is the precision.
*    x   = value to be printed.
*
* Returned:
*    Error code as defined by fprintf().
*/
{
  register int i;
  register int e = 0;

  if (x >= 0.0)
  {
    e |= fprintf(fp, fmt, x);
  }
  else
  {
    for (i = atoi(&fmt[1]); i > 0; i--)
      e |= fprintf(fp, " ");
  }
  return (e);
}
