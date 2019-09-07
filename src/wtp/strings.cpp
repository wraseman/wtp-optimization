/* strings.c */

/* Functions which alter strings from uppercase to lowercase and vice versa. 
 * source: http://stackoverflow.com/questions/23618316/undefined-reference-to-strlwr 
 */

#include "wtp.h"

char *strupr(char *string)
{
  unsigned char *p = (unsigned char *)string;

  while (*p)
  {
    *p = toupper(*p);
    p++;
  }

  return string;
}
char *strlwr(char *string)
{
  unsigned char *p = (unsigned char *)string;

  while (*p)
  {
    *p = tolower(*p);
    p++;
  }

  return string;
}