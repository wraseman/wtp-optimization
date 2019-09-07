/* THM_WTP.c   -- September 16, 1993
*    thm_wtp() is called from thm_gui() and prints "THM & disinfection"
*    to an opened output device.
*/
#include "wtp.h"

short thm_labels(const char *name, FILE *fp);

int thm_wtp(
    struct ProcessTrain *train, /* Process train control structure */
    FILE *fout                  /* stdout,stdprn, or fopen(.."w")  */
    // void         (*wait_gui)(char *msg),/* CRT pause function              */
    // long max_line_count
)
/*
*  Purpose: Print a one screen summary of Disinfection and DBP parameters.
*
*  Notes:
*   1. thm_wtp() is ANSI.
*   2. Called from thm_gui() to support "THM & Disinfection" on main menu.
*   3. The \r is needed on Laser printers.
*
*  Michael D. Cummins
*     April 1993
*/
{
  int numprocs;
  double *CT_saved;
  double log_required_init;
  double init_ratio_g;
  double init_ratio_v;
  double init_ratio_c;
  int ec_step_1_comply = 2;
  int ec_exempt = 2;
  int CT_index;
  struct UnitProcess *unit;
  struct Effluent *eff;

  /* Self protection */
  if (train == NULL || fout == NULL)
    return (FALSE);

  /* Count number of unit processes */
  numprocs = 0;
  for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
  {
    numprocs++;
  }

  /* Allocate memory to save CT_Ratio */
  CT_saved = (double *)calloc(3 * numprocs + 10, sizeof(double));
  CT_index = 0;

  coldflag = TRUE;

  runmodel(train);
  for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
  {
    if (unit->type == INFLUENT)
    {
      log_required_init = unit->eff.log_required_g;
      init_ratio_g = unit->eff.ct_ratio;
      init_ratio_v = unit->eff.ct_ratio_v;
      init_ratio_c = unit->eff.ct_ratio_c;
    }
    CT_saved[CT_index] = unit->eff.ct_ratio;
    CT_index++;
    CT_saved[CT_index] = unit->eff.ct_ratio_v;
    CT_index++;
    CT_saved[CT_index] = unit->eff.ct_ratio_c;
    CT_index++;
  }

  coldflag = FALSE;
  if (runmodel(train) == FALSE)
    return (FALSE);

  runmodel(train);
  CT_index = 0;

  thm_labels(train->file_name, fout);
  for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
  {
    eff = &unit->eff;
    fprintf(fout, "%-17s", UnitProcessTable[unit->type].name);
    /*      fdoublef(fout,"% 5.1f", eff->pH                 );
      fdoublef(fout,"% 7.1f", eff->FreeCl2 * MW_Cl2   );
      fdoublef(fout,"% 8.1f", eff->NH2Cl   * MW_Cl2   );  */
    if (log_required_init != 999999.0)
    {
      fdoublef(fout, "%6.1f", CT_saved[CT_index]);
      CT_index++;
      fdoublef(fout, "%7.1f", CT_saved[CT_index]);
      CT_index++;
      fdoublef(fout, "%7.1f", CT_saved[CT_index]);
      CT_index++;
    }
    else
    {
      fprintf(fout, "    na");
      fprintf(fout, "     na");
      fprintf(fout, "     na");
    }
    fdoublef(fout, "%6.1f", eff->FreeCl2 * MW_Cl2);
    fdoublef(fout, "%6.1f", eff->NH2Cl * MW_Cl2);
    fdoublef(fout, "%7.1f", eff->chlorite);
    fdoublef(fout, "%6.0f", eff->BrO3);
    fdoublef(fout, "%6.1f", eff->TTHM);
    fdoublef(fout, "%6.1f", eff->HAA5);

    fprintf(fout, "\r\n");

    if (unit->type == WTP_EFFLUENT)
    {
      ec_step_1_comply = eff->ec_meeting_step1;
      ec_exempt = eff->ec_exempt;
    }

  } //End "for(unit=FirstUnitProcess(train); unit; unit=NextUnitProcess(unit) )"

  fprintf(fout, "\r\n");

  //E.C. outputs
  if (log_required_init == 999999.0)
  { // GW is source
    fprintf(fout, "E.C. requirements do not apply - Groundwater Is Source\r\n");
    //count += 1;
  }
  else if (ec_exempt == 2)
  {
    fprintf(fout, "E.C. compliance not determined - no WTP_EFFLUENT\r\n");
    //count += 1;
  }
  else //EC applies and compliance can be determined
  {
    if (ec_exempt == TRUE)
      fprintf(fout, "E.C. not required - raw TOC, raw SUVA, and/or finished TOC <= 2\r\n");
    else
      fprintf(fout, "E.C. raw TOC, raw SUVA, and finished TOC <= 2 exemptions do not apply\r\n");
    if (ec_step_1_comply == TRUE)
      fprintf(fout, "E.C. Step 1 TOC removal requirement ACHIEVED\r\n");
    else
      fprintf(fout, "E.C. Step 1 TOC removal requirement NOT ACHIEVED\r\n");
    //count +=2;
  }

  //CT Ratio Notes
  if (log_required_init == 999999.0)
  {
    fprintf(fout, "na = not applicable - no CT requirements exist for groundwater\r\n");
  }

  if (init_ratio_v == 1.0)
  {
    fprintf(fout, "Virus CT Ratio = 1 because physical removal provided all disinfection\r\n");
  }
  if (init_ratio_g == 1.0)
  {
    fprintf(fout, "Giardia CT Ratio = 1 because physical removal provided all disinfection\r\n");
  }
  if (init_ratio_c == 1.0)
  {
    fprintf(fout, "Crypto. CT Ratio = 1 because physical removal provided all disinfection\r\n");
  }

  // if( wait_gui!=NULL ) wait_gui("Listing complete");
  free(CT_saved);

  // if (max_line_count > 0)
  // {
  //   /* Do Nothing */
  // }
  return (TRUE);
}

short thm_labels(const char *name, FILE *fp)
{
  static const char *header[] = {
      "                                 Table 11                                 ",
      "          Predicted Inactivation at Minimum Temperature and Peak Flow     ",
      "               and DBPs at Plant Flow and Influent Temperature            ",
      "--------------------------------------------------------------------------",
      "                 _____ CT Ratios ____  Cl2   NH2Cl  ClO2- BrO3- TTHM  HAA5",
      "Location         Giardia Virus Crypto (mg/L as Cl2)(mg/L) <----(ug/L)---->",
      "--------------------------------------------------------------------------",
      NULL};
  register int count = 0;
  register const char **ptr;

  fprintf(fp, "Process Train: %s\r\n", name);
  fprintf(fp, "\r\n");

  for (ptr = header; *ptr != NULL; ptr++, count++)
    fprintf(fp, "%s\r\n", *ptr);

  return (count + 2); /* Number of lines printed */
}
