/* ftab.c */
#include "wtp.h"

int ftab(FILE *fp, register short tab)
/*
*  Purpose: print 'tab' spaces on output device.
*/
{
  register int e = 0;
  for (; tab > 0; tab--)
    e |= fprintf(fp, " ");
  return (e);
}
