/* Mfuf_Rmv.c */
#include "wtp.h"

void mfuf_rmv(struct UnitProcess *unit)
/*
*   This subroutine calculates changes in flow, etc. through an MFUF_UP
*
*
* Note: Developed for SWAT
*/
{
  //FILE *fptr;

  /* Inputs: */
  double recovery; /* Effluent/Influent flow           (%)  */

  struct Mfuf *mfuf;
  struct Effluent *eff;

  if (unit == NULL || unit->data.ptr == NULL || unit->type != MFUF_UP)
    return;

  /* Get inputs from UnitProcess structure */
  mfuf = unit->data.mfuf;
  recovery = mfuf->recover;
  eff = &unit->eff;

  /* Some Protection/Range-Checking */
  if (recovery < 1.0)
    recovery = 1.0;
  if (recovery > 100)
    recovery = 100;

  /* Set "eff->Flow to account for recovery" */
  eff->Flow *= recovery / 100.0;

  /* Set this CLO2 flag */
  //  eff->clo2rawflag = FALSE;

  /* Eliminate any ozone residual */
  eff->o3_res = 0.0;

} /* End subroutine "mfuf_rmv()" */