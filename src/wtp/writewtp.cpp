/* WriteWTP.c */
#include "wtp.h"

int writewtp(FILE *fp, struct ProcessTrain *train, FILE *ferr)
/*
*  Purpoes: Write process train and unit treatment process data to a file.
*
*  Inputs:
*    fp    = Calling routine must fopen(xxx,"w") file.
*    train = Process train data to be saved.   
*    ferr  = stderr or file to report errors to.
*
*  Return: TRUE/FALSE for success/fail.
*
*  Notes:
*   1. writewtp() is ANSI.
*   2. This function is intended to support save_wtp().
*   3. The return value from fprintf() and the write_xxx() support functions
*      is the number of characters transmited or a negative value if an 
*      error occours.  The return values are or(ed) together into an error
*      flag called 'e'.  If any return value is negative then the sign bit
*      of 'e' will be set.  The sign bit is checked at the conclusion of 
*      writting each unit process data packet and an error message printed
*      to 'ferr' if an error occours.
*
*  Michael D. Cummins
*     July 1993
*/
{
  struct UnitProcess *unit;
  int e; /* error flag from fprintf() */
  int success = TRUE;
  int i;
  short type;
  short tab = 0;
  struct DataID *data;
  int max_unit_types;
  char buffer[80];

  if (fp == NULL || train == NULL)
    return (FALSE);

  max_unit_types = 0;
  while (UnitProcessTable[max_unit_types].name != NULL)
    max_unit_types++;

  e = fprintf(fp, "U.S. EPA WTP Model:%s\n", VERSION);
  if (e < 0)
  {
    if (ferr != NULL)
      fprintf(ferr, "Could not write to data file.\n");
    return (FALSE);
  }

  for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
  {
    e = 0;
    type = unit->type;
    if (type >= 0 && type < max_unit_types)
    {
      e |= fprintf(fp, "%s\n", UnitProcessTable[unit->type].name);
      tab += 2;

      if (UnitProcessTable[type].write == NULL && unit->type != WTP_EFFLUENT)
      /* WTP_EFFLUENT caveat needed because it now has no inputs */
      {
        if (ferr != NULL)
        {
          fprintf(ferr,
                  "Output function not installed in UnitProcessTable[] for %s\n",
                  UnitProcessTable[type].name);
        }
        success = FALSE;
      }
      else
      {
        for (data = UnitProcessTable[unit->type].data, i = 0;
             data->id != NULL;
             data++, i++)
        {
          e |= UnitProcessTable[type].write(buffer, i, unit);
          e |= ftab(fp, tab);
          e |= fprintf(fp, "%-8s= %s\n", data->id, buffer);
        }
      }
      e |= ftab(fp, tab);
      e |= fprintf(fp, "End\n");
      tab -= 2;
    }
    if (e < 0)
    {
      if (ferr != NULL)
      {
        fprintf(ferr,
                "Error in writting %s\n", UnitProcessTable[type].name);
      }
      success = FALSE; /* Keep on going?  */
    }
  } /* end for(unit...) */

  return (success);
}
