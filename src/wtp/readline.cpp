/* readLine.c */
#include "wtp.h"

int read_line(FILE *fp, char **id, char **value)
/*
*  Purpose: Support function to reading WTP data file.  This function will
*           read the next line from fp and return a string pointer to the 
*           data element id and a string pointer to the value.
*/
{
  register char *id1, *value1;
  static char buffer[80];

  if (fp == NULL)
    return (FALSE);

  memset(buffer, 0, sizeof(buffer)); /* clear memory on stack */
  fgets(buffer, sizeof(buffer), fp); /* read in first line of file to buffer */

  /* Remove newline character from buffer. */
  id1 = buffer;
  if (id1 != '\0')
  {
    while (*id1 != '\0')
      id1++;            /* Move to end of buffer */
    while (*id1 <= ' ') /* compares ASCII values [0-255] of *id1 and a space */
    {                   /* source: https://stackoverflow.com/questions/22736348/char-comparison-in-c */
      /* source: http://www.asciitable.com/ */
      *id1 = '\0';
      id1--;
    }
  }

  /* Position 'id1' to point to first character in buffer. */
  id1 = buffer;
  while (((*id1 == ' ') || (*id1 == '\t')) && (*id1 != '\0'))
    id1++; // WJR: added tab functionality

  value1 = id1;
  if ((*value1 != '\0'))
  {
    while ((*value1 != ' ') && (*value1 != '\t') && (*value1 != '=') && (*value1 != '\0'))
      value1++; /* Move past letters */

    if (*value1 == '=')
    {
      *value1 = '\0'; /* Mark end of id */
      value1++;
    }
    else if ((*value1 == ' ') || (*value1 == '\t')) /* If there is whitespace... */
    {
      *value1 = '\0';
      value1++;
      while ((*value1 != '=') && (*value1 != '\0'))
        value1++;
      if (*value1 == '=')
        value1++;
    }
    else
    {
      /* *value = \0 */
    }
    while (*value1 <= ' ' && *value1 != '\0')
      value1++;
  }

  *id = id1;
  *value = value1;

  if (strcmp(id1, "End") == 0)
    return (TRUE);
  else if (*id1 != '\0' && *value1 != '\0')
    return (TRUE);
  else
    return (FALSE);
}
