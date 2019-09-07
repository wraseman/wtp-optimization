/* save_wtp.c */
#include "wtp.h"

int save_wtp(char *file_name, struct ProcessTrain *train, FILE *ferr)
/*
*  Purpose: Save process train data to a data file..
*
*  Input:
*    file_name : The calling routine is responsibal for insuring that the
*                file name is valid... such as having a .wtp extension.
*
*  Return:
*    TRUE : file was successfuly saved.
*    FALSE: an error occoured.
*
*  Note:
*   1. This routine was orginally developed to support cli_wtp(),
*      specifically the -s<filename> switch.
*
*  Michael D. Cummins
*     July 1993
*/
{
  FILE *fp;
  int success = FALSE;

  if (file_name != NULL && strlen(file_name) > 0)
  {
    fp = fopen(file_name, "w");
    if (fp == NULL)
    {
      if (ferr != NULL)
        fprintf(ferr, "Could not open %s\n", file_name);
      success = FALSE;
    }
    else
    {
      success = writewtp(fp, train, ferr);
      if (success == TRUE)
      {
        strncpy(train->file_name, file_name, sizeof(train->file_name) - 1);
      }
      fflush(fp);
      fclose(fp);
    }
  }

  return (success);
}
