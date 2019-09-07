/* Read_WTP.c  -- September 20, 1993 */

#include "wtp.h"

int strip_buffer(char buffer[], FILE *fout)
{

    /* Purpose: remove comments, whitespace, and return characters from buffer.
     *
     * Inputs:
     *  buffer - character array which contains name of unit process or simulation mode.
     *  fout - where to print.
     *
     * Return:
     *  continue - when exiting strip_buffer() issue a 'continue' statement to the for-loop.
     *
     */

    int cont = FALSE;
    char *ptr;

    /* Eliminate whitespace */
    while ((*buffer == ' ') || (*buffer == '\t'))
        buffer++;

    /* Skip any lines with comments (signified by '#') */
    if (buffer[0] == '#')
    {
        return cont = TRUE;
    }

    /* Skip any blank lines*/
    if (buffer[0] == '\n')
    { // blank lines
        return cont = TRUE;
    }

    /* Strip out any '\n' character(s) */
    ptr = strchr(buffer, '\n');
    if (ptr != NULL)
    { /* If '\n' is found... */
        *ptr = '\0';
    }

    /* Strip out any comments */
    if (strchr(buffer, '#') != NULL)
    { /* If '#' is found... */
        ptr = strchr(buffer, '#');
        *ptr = '\0';
    }

    /* Remove any trailing symbols from right to left until it reaches a letter or a period (for abbreviations) */
    while (isalpha(*ptr) == 0 && (*ptr != '.'))
    {
        *ptr = '\0';
        ptr--;
    }

    /* Echo modified buffer */
    if (fout != NULL)
    {
        fprintf(fout, "%s\n", buffer);
    }

    return cont = FALSE;
}

int read_wtp(
    FILE *fin,                  /* Process Train data file         */
    struct ProcessTrain *train, /* Process train control structure */
    FILE *fout                  /* stdout or NULL                  */
                                //   void             (*wait)(char *msg),/* User intrface function          */
    // long max_line_count         /* Lines the crt will display      */
)
/*
*  Purpose: Read process train and unit treatment process data from a file
*           into a process train control structure.  Any existing data in
*           the process train control structure are lost.
*
*  Inputs:
*    fin  = The calling routine must fopen(..."r") the file containg
*           process train data.
*    fout = where to print loading status & error messages. Can be NULL.
*    wait = User interface call back function that is called when
*           'max_line_count' lines have been sent to fout.  read_wtp()
*           supplies a 'msg' to wait() that states if the listing is
*           complete or there is more to come.
*    max_line_count
*           Number of line that can be displayed on the CRT.
*  Return:
*    TRUE  = success, *train contains new data
*    FALSE = read error... *train will contain as much as possible.
*
*  Note:
*    1. This function is intended to support open_wtp().
*    2. This function will identify previously used unit process names
*       in a data file and convert to the current name.
*    3. Any error messages are sent to 'fout'.
*
*  Michael D. Cummins
*     July 1993
*/
{
    int type;    /* Unit process type */
    int i;       /* Data element index */
                 //  int                 mm=0;      /* mm: Printed line counter */
    int success; /* returned value */
    int ok;
    int done;
    int found;

    char *id;
    char *val;
    struct UnitProcess *unit;
    FILE *ferr; /* send error messages to *fout or *stderr */
    struct DataID *data;
    int max_unit_types;
    char buffer[256];
    register char *unit_name;

    /* Determine how many entries in UnitProcessTable */
    max_unit_types = 0;
    while (UnitProcessTable[max_unit_types].name != NULL)
        max_unit_types++;

    /* clear memory on stack */
    memset(buffer, 0, sizeof(buffer));

    /* This may change. */
    ferr = fout;

    /* self protection */
    if (fin == NULL || train == NULL)
    {
        fprintf(ferr, "Error: failed self-protection in read_wtp()\n");
        exit(EXIT_FAILURE);
    }

    /* Initialize process train to an empty list. */
    while (FirstUnitProcess(train))
        RemoveUnitProcess(FirstUnitProcess(train));

    /* The following is the code to the current version of the reader. */
    if (fout != NULL)
    {
        fprintf(fout, "%s", buffer); /* Echo Title line from data file. */
                                     //      mm++;
    }

    success = TRUE;

    /* Log progress */
    if (fout != NULL)
    {
        fprintf(fout, "Treatment train processes:\n");
    }

    while (!feof(fin))
    {
        /* Read name of unit process. */
        if (fgets(buffer, 120, fin) != NULL)
        {
            unit_name = buffer;
            if (strip_buffer(unit_name, fout) == TRUE)
            { /* Remove comments, whitespace, and returns from buffer */
                continue;
            }

            /* Identify unit process and add to end of process train. */
            type = 0;
            unit = NULL;
            while (type < max_unit_types && unit == NULL)
            {
                if (strcmp(unit_name, UnitProcessTable[type].name) == 0)
                {
                    unit = AddUnitProcess(train, type);
                }
                else
                {
                    type++;
                }
            }

            type = 0; /* use 'unit->type' */

            if (unit == NULL)
            {
                fprintf(ferr, "Error: unknown unit process name '%s' found in %s\n", unit_name, train->file_name);
                exit(EXIT_FAILURE); /* error handling added by WJR, 8/17 */
                                    //                success = FALSE; /* Unknown unit process. */
            }
            else
            {

                /* Unit Process is valid... Read the unit process data. */
                type = unit->type;

                if (UnitProcessTable[type].read == NULL && unit->type != WTP_EFFLUENT)
                { /* WTP_EFFLUENT caveat needed because it now has no inputs */
                    /* Oops, a program development error. */
                    if (ferr != NULL)
                    {
                        fprintf(ferr,
                                "read function not installed in UnitProcessTable[] for %s.\n",
                                UnitProcessTable[type].name);
                    }

                    success = FALSE;
                }
                else
                { /* read function exists... read data packet. */
                    done = FALSE;
                    ok = TRUE;

                    while (!done) /* Iterate through unit process data entries */
                    {
                        if (read_line(fin, &id, &val) != TRUE)
                        {
                            /* This should not happen. */
                            if (ferr != NULL)
                            {
                                fprintf(ferr,
                                        "Unexpected end of file or blank line.\n");
                            }
                            done = TRUE; /* End of file or read error. */
                            ok = FALSE;
                        }
                        else
                        {
                            /* We have a valid line of data... now process it */
                            if (strcmp(id, "End") == 0)
                            {
                                done = TRUE; /* Normal end of data packet. */
                            }
                            else
                            {
                                /* Determine which data element. */
                                found = FALSE;
                                i = 0;
                                for (data = UnitProcessTable[type].data; data->id && !found; data++)
                                {
                                    if (strcmp(data->id, id) == 0)
                                        found = TRUE;
                                    else
                                        i++;
                                }

                                /* At this point 'i' is an index to the data element */

                                if (found)
                                {
                                    UnitProcessTable[type].read(val, i, unit);
                                }
                                else
                                {
                                    /* Incorrect ID in data file. */
                                    if (ferr != NULL)
                                    {
                                        ok = FALSE;
                                        fprintf(ferr, "%s not member of %s\n", id, UnitProcessTable[type].name);
                                    }
                                }
                            } /* End of else ( strcmp( id, "End")==0 ) */
                        }     /* End of else (read_line(fin, &id, &val)!=TRUE) */
                    }         /* End of while(!done) */
                }             /* End of else (UnitProcessTable[type].read==NULL && unit->type!=WTP_EFFLUENT) */

                if (ok == FALSE)
                {
                    fprintf(ferr, "%s has corrupt %s data.\n", train->file_name, UnitProcessTable[unit->type].name);
                    success = FALSE;
                }

            } /* End of else ( unit == NULL ) */

            //          /* Check line count */
            //          if( (fout!=NULL) && (wait!=NULL) && (mm>=max_line_count) )
            //            {
            //              mm = 0;
            //              wait( "Press Continue" );
            //            }
        }
    } /* End of while(!eof(fin)) */

    //   if( wait != NULL ) wait("Listing complete");

    return (success);
}