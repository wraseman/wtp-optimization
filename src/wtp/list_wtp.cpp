/* List_WTP.c */
#include "wtp.h"

int list_wtp(
    struct ProcessTrain *train, /* Process train control structure    */
    FILE *fp                    /* stdprn, stdout, ...                */
    // void      (*wait_gui)(char *msg),/* User interface function            */
    // long max_line_count
)
/*
*  Purpose: Print Process Train and Unit Process data to *fp.
*
*  Inputs:
*   *train              Process train control structure.
*   *fp                 Opened output device or file.
*   (*wait)(char *msg)  User interface call back function for listing to CRT.
*   max_line_count      Number of lines the CRT can display.
*
*  Return:
*    TRUE/FALSE for success/fail. Very little can go wrong in this function.
*
*  Notes:
*   1. list_wtp() is ANSI.
*   2. list_wtp() is called from list_gui() to print the process train
*      and unit process data to an output device.
*   3. The \r is needed on Laser printers.
*
*  Michael D. Cummins
*      July 1993
*/
{
  int e = 0;     /* fprintf() error flag */
  int count = 0; /* Printed line counter */
  int type;      /* unit process index   */
  int i;         /* data element index   */
  struct DataID *data;
  struct UnitProcess *unit;
  char buffer[120];
  char line[120];

  if (train == NULL || fp == NULL)
    return (FALSE);

  /* Print title at top of page */
  e |= fprintf(fp, "Process train data for %s\r\n", train->file_name);
  count++;

  /* Main loop over process train */
  for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
  {
    /* Pull type into register. */
    type = unit->type;

    // if( (wait_gui != NULL) && (count >= max_line_count - 5) )
    //   {
    //     count = 0;
    //     wait_gui( "Press Continue" );
    //   }

    /* Print Name of unit process. */
    e |= fprintf(fp, "%s\r\n", UnitProcessTable[type].name);
    count++;
    count++;

    /* Loop over unit process data. */
    for (data = UnitProcessTable[type].data, i = 0; data->id; data++, i++)
    {
      /* Name of data element. */
      sprintf(line, "        %s ", data->name);

      /* Add '.'s to column 40 */
      while (strlen(line) < 50)
        strcat(line, ".");

      /* Add value of data element. */
      UnitProcessTable[type].write(buffer, i, unit);
      strcat(line, " ");
      strcat(line, buffer);

      /* Add ' 's to column 48 */
      while (strlen(line) < 58)
        strcat(line, " ");

      /* Add units of measurment. */
      strcat(line, data->units);

      /* Output 'line' to 'fp' */
      e |= fprintf(fp, "%s\r\n", line);
      count++;
    }
  }

  // if( wait_gui!=NULL ) wait_gui("Listing complete" );

  if (e < 0)
    return (FALSE);
  else
    return (TRUE);
}
