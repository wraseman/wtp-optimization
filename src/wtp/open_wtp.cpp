/* open_wtp.c */
#include "wtp.h"

int open_wtp(
    const char *file_name,
    struct ProcessTrain *train,
    FILE *fout
    // void    (*wait)(char *msg),
    // long max_line_count
)
/*
*  Purpose: Read process train data from a file into 'train'.
*
*  Input:
*    file_name : Must be a valid existing file.
*
*  Return:
*    TRUE : file was successfully loaded into train.
*    FALSE: an error occurred.
*
*    train:
*      1. contains process train data.
*      2. file_name is copyed into train->file_name.
*
*  Note:
*   1. This routine was orginally developed to support cli_wtp(),
*      specifically the -o<filename> switch.
*
*  Michael D. Cummins
*     April 1993
*/
{
  FILE *fp;
  int success = FALSE;

  if (file_name != NULL && strlen(file_name) > 0)
  {
    fp = fopen(file_name, "r");
    if (fp == NULL)
    {
      if (fout != NULL)
        fprintf(fout, "Could not open %s\n", file_name);
      success = FALSE;
    }
    else
    {
      strncpy(train->file_name, file_name, sizeof(train->file_name) - 1);
      // success = read_wtp( fp, train, fout, wait, max_line_count );
      success = read_wtp(fp, train, fout);
      fclose(fp);
    }
  }

  return (success);
}
