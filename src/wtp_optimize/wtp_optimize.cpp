/* wtp_optimize.cpp */

#include "wtp_optimize.h"
#include "wtp.h"
#include "auto_dose.h"
#include "borg.h"
#include <iostream>

// wtp_optimize.cpp contains functions which implement the various run modes in the wtp-optimize project

/* These run modes include:
* - single simulation (with automatic chemical dosing)
* - multi-simulation (with automatic chemical dosing)
* - simulation-optimization (with automatic chemical dosing)
* - validation
*/

/* Single simulation mode (with automatic chemical dosing) */
int single_sim_mode(struct ProcessTrain *train,               /* Process train control structure    */
                    struct OperationalParameters *operations, // operational parameters
                    struct InfluentParameters *influent_wq,
                    FILE *fout /* Output file */
)
/*
* Purpose:
*   Run complete WTP model and send output to *fout.  Multiple outputtables are
*   produced by this function.  
*
* Inputs:
*   train     = Process train control structure.
*   fout      = stdout, stdprn, or fopen(...,"w") file
*
* Return: TRUE/FALSE for success/fail.
*
* Notes:
*  1. This function is intended to support the "Run model" operation on
*     the user interface.
*  2. run_wtp() calls runmodel().
*  3. The \r is needed on Laser printers.
*
* Michael D. Cummins
*    July 1993
*/
{
    int e = 0;              /* fprintf() error flag */
    int bankclosed = FALSE; /*Ensures that output for bankfilter credit shows only once*/
    int presed_cntr = 0;    /* Ensures that output for presed credit shows only once */
    int soft2_cntr = 0;     /* Ensures that output for 2nd stage softening credit shows only once */
    short eff_cntr = 0;     /*counts number of WTP Effluents in train */
    short tap_cntr = 0;     /*counts number of AVG_TAP points in train*/
    short eos_cntr = 0;     /*counts number of END_OF_SYSTEM points in train*/
                            //  short             co2_cntr=0; /*counts number of CARBON_DIOXIDE points in train*/
    short count = 0;        /* Printed line counter */
    double toc_rem;         /* toc removal % */
    double lt2_wscp_credit = 0.0;
    //  double		    giardia_credit=0.0;
    //  double   		    virus_credit=0.0;
    //  double           	    crypto_credit=0.0;
    double ife_lr_credit = 1.0;
    double cfe_lr_credit = 0.5;
    double wscp_credit = 0.5;
    double soft2_credit = 0.5;
    double zero_lr = 0.0;
    const char *optfilter = "Turbidity-Optimized Filters";
    const char *wscp = "Watershed Control Program";
    const char *soft2basin = "2-Stage Softening";
    register struct UnitProcess *unit;
    register struct Effluent *eff;
    struct Effluent *influent;
    struct Effluent *effluent;
    struct Effluent *tap;
    struct Effluent *eos;

    char buffer[120];
    //int                         loc;  // WJR, location (0 or 1) of pH and alkalinity adjustment, WJR (6/10/17)

    /* Self protection: */
    if (train == NULL || fout == NULL)
        return (FALSE);

    coldflag = FALSE;

    if (runmodel(train) == FALSE)
        return (FALSE);

    e |= fprintf(fout, "\nProcess Train: %s\n", train->file_name);
    e |= fprintf(fout, "\n");
    count += 2;

    //  unit = LastUnitProcess(train)->eff.influent;
    //  e |= fprintf( fout,"influent pH:%9.1f\r\n",unit->data.influent->pH);

    /*********************************************************************************************************************************************/

    /* Influent, effluent, and avg_tap water quality parameters first */

    for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
    {
        switch (unit->type)
        {
        case INFLUENT:

            // change influent water quality parameters based on input file values
            unit->data.influent->alkalinity = influent_wq->alkalinity;
            unit->data.influent->nh3 = influent_wq->nh3;
            unit->data.influent->bromide = influent_wq->bromide;
            unit->data.influent->calcium = influent_wq->calcium;
            unit->data.influent->hardness = influent_wq->hardness;
            unit->data.influent->pH = influent_wq->pH;
            unit->data.influent->temp = influent_wq->temp;
            unit->data.influent->toc = influent_wq->toc;
            unit->data.influent->ntu = influent_wq->ntu;
            unit->data.influent->uv254 = influent_wq->uv254;

            influent = &unit->eff;
            break;

            //            case CARBON_DIOXIDE:
            //                if(co2_cntr==0) co2_1 = &unit->eff;
            //                co2_cntr++; break;

            //            case RAPID_MIX:
            //                rapid_mix = &unit->eff;break;

        case WTP_EFFLUENT:
            eff_cntr++;
            effluent = &unit->eff;
            break;

        case AVG_TAP:
            tap_cntr++;
            tap = &unit->eff;
            break;

        case END_OF_SYSTEM:
            eos_cntr++;
            eos = &unit->eff;
            break;

        default:
            break;
        }
    }

    // read in decision variables
    double alk_setpt_1 = operations->alk_setpt;
    double pH_setpt_1 = operations->pH_setpt;
    double cl2_setpt = 0.2;
    //double alk_setpt_2 = 20.0;
    double pH_setpt_2 = 8.0;
    double DBP_safety_factor = operations->DBPsf;

    auto_dose(train, pH_setpt_1, alk_setpt_1, cl2_setpt, pH_setpt_2, DBP_safety_factor, fout);

    int debug = TRUE;

    if (debug == TRUE)
    {
        //Calculate toc removal through plant here
        if (eff_cntr > 0)
        {
            if (influent->TOC > 0.0)
                toc_rem = (influent->TOC - effluent->TOC) / influent->TOC * 100.0;
            else
                toc_rem = 0.0;
        }

        e |= fprintf(fout, "                                     Table 1                                 \r\n");
        e |= fprintf(fout, "         Water Quality Summary for Raw, Finished, and Distributed Water      \r\n");
        e |= fprintf(fout, "          At Plant Flow (%4.1f MGD) and Influent Temperature (%3.1f C)       \r\n", influent->Flow, influent->DegK - 273.15);
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        e |= fprintf(fout, "Parameter             Units      Raw Water   Effluent   Avg. Tap   End of Sys\r\n");
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        count += 6;

        if (eff_cntr > 0 && tap_cntr == 0 && eos_cntr == 0)
        { /* Print-out influent and effluent information only */

            e |= fprintf(fout, "pH                    (-)       %8.1f%11.1f\r\n", influent->pH, effluent->pH);
            e |= fprintf(fout, "Alkalinity      (mg/L as CaCO3) %8.0f%11.0f\r\n", influent->Alk * MW_CaCO3 / 2, effluent->Alk * MW_CaCO3 / 2);
            e |= fprintf(fout, "TOC                  (mg/L)     %8.1f%11.1f\r\n", influent->TOC, effluent->TOC);
            e |= fprintf(fout, "UVA                  (1/cm)     %8.3f%11.3f\r\n", influent->UV, effluent->UV_out);
            e |= fprintf(fout, "(T)SUVA              (1/cm)     %8.1f%11.1f\r\n", influent->SUVA, effluent->SUVA);
            e |= fprintf(fout, "Ca Hardness     (mg/L as CaCO3) %8.0f%11.0f\r\n", influent->Ca_aq * MW_CaCO3, effluent->Ca_aq * MW_CaCO3);
            e |= fprintf(fout, "Mg Hardness     (mg/L as CaCO3) %8.0f%11.0f\r\n", influent->Mg_aq * MW_CaCO3, effluent->Mg_aq * MW_CaCO3);
            e |= fprintf(fout, "Ammonia-N            (mg/L)     %8.2f%11.2f\r\n", influent->NH3 * MW_NH3, effluent->NH3 * MW_NH3);
            e |= fprintf(fout, "Bromide              (ug/L)     %8.0f%11.0f\r\n", influent->Br * MW_Br * 1000, effluent->Br * MW_Br * 1000);
            e |= fprintf(fout, "Free Cl2 Res.     (mg/L as Cl2) %8.1f%11.1f\r\n", influent->FreeCl2 * MW_Cl2, effluent->FreeCl2 * MW_Cl2);
            e |= fprintf(fout, "Chloramine Res.   (mg/L as Cl2) %8.1f%11.1f\r\n", influent->NH2Cl * MW_Cl2, effluent->NH2Cl * MW_Cl2);
            e |= fprintf(fout, "TTHMs                (ug/L)     %8.0f%11.0f\r\n", influent->TTHM, effluent->TTHM);
            e |= fprintf(fout, "HAA5                 (ug/L)     %8.0f%11.0f\r\n", influent->HAA5, effluent->HAA5);
            e |= fprintf(fout, "HAA6                 (ug/L)     %8.0f%11.0f\r\n", influent->HAA6, effluent->HAA6);
            e |= fprintf(fout, "HAA9                 (ug/L)     %8.0f%11.0f\r\n", influent->HAA9, effluent->HAA9);
            e |= fprintf(fout, "TOX                  (ug/L)     %8.0f%11.0f\r\n", influent->tox, effluent->tox);
            e |= fprintf(fout, "Bromate              (ug/L)     %8.0f%11.0f\r\n", influent->BrO3, effluent->BrO3);
            e |= fprintf(fout, "Chlorite             (mg/L)     %8.1f%11.1f\r\n", influent->chlorite, effluent->chlorite);
            e |= fprintf(fout, "TOC Removal         (percent)   %19.0f\r\n", toc_rem);
            //E.C. outputs
            if (influent->log_required_g == 999999.0)
            { // GW is source
                e |= fprintf(fout, "E.C. requirements do not apply - Groundwater Is Source\r\n");
                count += 1;
            }
            else
            {
                if (effluent->ec_exempt == TRUE)
                    e |= fprintf(fout, "E.C. not required - raw TOC, raw SUVA, and/or finished TOC <= 2\r\n");
                else
                    e |= fprintf(fout, "E.C. raw TOC, raw SUVA, and finished TOC <= 2 exemptions do not apply\r\n");
                if (effluent->ec_meeting_step1 == TRUE)
                    e |= fprintf(fout, "E.C. Step 1 TOC removal requirement ACHIEVED\r\n");
                else
                    e |= fprintf(fout, "E.C. Step 1 TOC removal requirement NOT ACHIEVED\r\n");
                count += 2;
            }
            //CT Ratio Outputs
            if (swflag == FALSE /*influent->log_required_g == 999999.0*/ && gw_virus_flag == FALSE)
            {
                e |= fprintf(fout, "CT Ratios Not Applicable - Groundwater Is Source\r\n");
                count += 1;
            }
            else if (swflag == TRUE)
            {
                //   if(influent->log_required
                e |= fprintf(fout, "CT Ratios                        \r\n");
                e |= fprintf(fout, "  Virus                (-)      %8.1f%11.1f\r\n", influent->ct_ratio_v, effluent->ct_ratio_v);
                e |= fprintf(fout, "  Giardia              (-)      %8.1f%11.1f\r\n", influent->ct_ratio, effluent->ct_ratio);
                e |= fprintf(fout, "  Cryptosporidium*     (-)      %8.1f%11.1f\r\n", influent->ct_ratio_c, effluent->ct_ratio_c);
                e |= fprintf(fout, "*Crypto. Disinfection Calcs. Based on Proposed LT2 Rule\r\n");
                count += 5;
            }
            else // (swflag==FALSE && gw_virus_flag==TRUE) **/
            {
                e |= fprintf(fout, "CT Ratios                        \r\n");
                e |= fprintf(fout, "  Virus*               (-)      %8.1f%11.1f\r\n", influent->ct_ratio_v, effluent->ct_ratio_v);
                e |= fprintf(fout, "*CT Ratios N/A for Giardia & Crypto. in a GW source\r\n");
                count += 3;
            }
        } /* End of "if (eff_cntr > 0 && tap_cntr == 0 && eos_cntr == 0)" */

        else if (eff_cntr == 0)

        { /* There's no WTP_EFFLUENT (and we'll assume there's also no AVG_TAP)*/

            e |= fprintf(fout, "pH                     (-)      %8.1f\r\n", influent->pH);
            e |= fprintf(fout, "Alkalinity      (mg/L as CaCO3) %8.0f\r\n", influent->Alk * MW_CaCO3 / 2);
            e |= fprintf(fout, "TOC                  (mg/L)     %8.1f\r\n", influent->TOC);
            e |= fprintf(fout, "UV                   (1/cm)     %8.3f\r\n", influent->UV);
            e |= fprintf(fout, "(T)SUVA              (1/cm)     %8.1f\r\n", influent->SUVA);
            e |= fprintf(fout, "Ca Hardness     (mg/L as CaCO3) %8.0f\r\n", influent->Ca_aq * MW_CaCO3);
            e |= fprintf(fout, "Mg Hardness     (mg/L as CaCO3) %8.0f\r\n", influent->Mg_aq * MW_CaCO3);
            e |= fprintf(fout, "Ammonia-N            (mg/L)     %8.2f\r\n", influent->NH3 * MW_NH3);
            e |= fprintf(fout, "Bromide              (ug/L)     %8.0f\r\n", influent->Br * MW_Br * 1000);
            e |= fprintf(fout, "Free Cl2 Res.     (mg/L as Cl2) %8.1f\r\n", influent->FreeCl2 * MW_Cl2);
            e |= fprintf(fout, "Chloramine Res.   (mg/L as Cl2) %8.1f\r\n", influent->NH2Cl * MW_Cl2);
            e |= fprintf(fout, "TTHMs                (ug/L)     %8.0f\r\n", influent->TTHM);
            e |= fprintf(fout, "HAA5                 (ug/L)     %8.0f\r\n", influent->HAA5);
            e |= fprintf(fout, "HAA6                 (ug/L)     %8.0f\r\n", influent->HAA6);
            e |= fprintf(fout, "HAA9                 (ug/L)     %8.0f\r\n", influent->HAA9);
            e |= fprintf(fout, "TOX                  (ug/L)     %8.0f\r\n", influent->tox);
            e |= fprintf(fout, "Bromate              (ug/L)     %8.0f\r\n", influent->BrO3);
            e |= fprintf(fout, "Chlorite             (mg/L)     %8.1f\r\n", influent->chlorite);
            e |= fprintf(fout, "TOC Removal         (percent)              not determined\r\n");
            //EC Compliance
            if (influent->log_required_g == 999999.0)
            { // GW is source
                e |= fprintf(fout, "E.C. requirements do not apply - Groundwater Is Source\r\n");
                count += 1;
            }
            else
            {
                e |= fprintf(fout, "No WTP EFFLUENT - E.C. compliance cannot be determined\r\n");
                count += 1;
            }
            //CT Ratio Outputs
            if (swflag == FALSE /*influent->log_required_g == 999999.0*/ && gw_virus_flag == FALSE)
            {
                e |= fprintf(fout, "CT Ratios Not Applicable - Groundwater Is Source\r\n");
                count += 1;
            }
            else if (swflag == TRUE)
            {
                //   if(influent->log_required
                e |= fprintf(fout, "CT Ratios                        \r\n");
                e |= fprintf(fout, "  Virus                (-)      %8.1f\r\n", influent->ct_ratio_v);
                e |= fprintf(fout, "  Giardia              (-)      %8.1f\r\n", influent->ct_ratio);
                e |= fprintf(fout, "  Cryptosporidium*     (-)      %8.1f\r\n", influent->ct_ratio_c);
                e |= fprintf(fout, "*Crypto. Disinfection Calcs. Based on Proposed LT2 Rule\r\n");
                count += 5;
            }
            else // (swflag==FALSE && gw_virus_flag==TRUE) **/
            {
                e |= fprintf(fout, "CT Ratios                        \r\n");
                e |= fprintf(fout, "  Virus                (-)      %8.1f\r\n", influent->ct_ratio_v);
                e |= fprintf(fout, "CT Ratios N/A for for Giardia & Crypto.\r\n");
                count += 3;
            }

        } /* End of else if (eff_cntr == 0) */

        else if (eff_cntr > 0 && tap_cntr > 0 && eos_cntr == 0)

        { /* Print-out influent, effluent, and avg_tap information */

            e |= fprintf(fout, "pH                     (-)      %8.1f%11.1f%11.1f\r\n", influent->pH, effluent->pH, tap->pH);
            e |= fprintf(fout, "Alkalinity      (mg/L as CaCO3) %8.0f%11.0f%11.0f\r\n", influent->Alk * MW_CaCO3 / 2, effluent->Alk * MW_CaCO3 / 2, tap->Alk * MW_CaCO3 / 2);
            e |= fprintf(fout, "TOC                  (mg/L)     %8.1f%11.1f%11.1f\r\n", influent->TOC, effluent->TOC, tap->TOC);
            e |= fprintf(fout, "UV                   (1/cm)     %8.3f%11.3f%11.3f\r\n", influent->UV, effluent->UV_out, tap->UV_out);
            e |= fprintf(fout, "(T)SUVA              (1/cm)     %8.1f%11.1f%11.1f\r\n", influent->SUVA, effluent->SUVA, tap->SUVA);
            e |= fprintf(fout, "Ca Hardness     (mg/L as CaCO3) %8.0f%11.0f%11.0f\r\n", influent->Ca_aq * MW_CaCO3, effluent->Ca_aq * MW_CaCO3, tap->Ca_aq * MW_CaCO3);
            e |= fprintf(fout, "Mg Hardness     (mg/L as CaCO3) %8.0f%11.0f%11.0f\r\n", influent->Mg_aq * MW_CaCO3, effluent->Mg_aq * MW_CaCO3, tap->Mg_aq * MW_CaCO3);
            e |= fprintf(fout, "Ammonia-N            (mg/L)     %8.2f%11.2f%11.2f\r\n", influent->NH3 * MW_NH3, effluent->NH3 * MW_NH3, tap->NH3 * MW_NH3);
            e |= fprintf(fout, "Bromide              (ug/L)     %8.0f%11.0f%11.0f\r\n", influent->Br * MW_Br * 1000, effluent->Br * MW_Br * 1000, tap->Br * MW_Br * 1000);
            e |= fprintf(fout, "Free Cl2 Res.     (mg/L as Cl2) %8.1f%11.1f%11.1f\r\n", influent->FreeCl2 * MW_Cl2, effluent->FreeCl2 * MW_Cl2, tap->FreeCl2 * MW_Cl2);
            e |= fprintf(fout, "Chloramine Res.   (mg/L as Cl2) %8.1f%11.1f%11.1f\r\n", influent->NH2Cl * MW_Cl2, effluent->NH2Cl * MW_Cl2, tap->NH2Cl * MW_Cl2);
            e |= fprintf(fout, "TTHMs                (ug/L)     %8.0f%11.0f%11.0f\r\n", influent->TTHM, effluent->TTHM, tap->TTHM);
            e |= fprintf(fout, "HAA5                 (ug/L)     %8.0f%11.0f%11.0f\r\n", influent->HAA5, effluent->HAA5, tap->HAA5);
            e |= fprintf(fout, "HAA6                 (ug/L)     %8.0f%11.0f%11.0f\r\n", influent->HAA6, effluent->HAA6, tap->HAA6);
            e |= fprintf(fout, "HAA9                 (ug/L)     %8.0f%11.0f%11.0f\r\n", influent->HAA9, effluent->HAA9, tap->HAA9);
            e |= fprintf(fout, "TOX                  (ug/L)     %8.0f%11.0f%11.0f\r\n", influent->tox, effluent->tox, tap->tox);
            e |= fprintf(fout, "Bromate              (ug/L)     %8.0f%11.0f%11.0f\r\n", influent->BrO3, effluent->BrO3, tap->BrO3);
            e |= fprintf(fout, "Chlorite             (mg/L)     %8.1f%11.1f%11.1f\r\n", influent->chlorite, effluent->chlorite, tap->chlorite);
            e |= fprintf(fout, "TOC Removal         (percent)   %19.0f\r\n", toc_rem);
            //E.C. outputs
            if (influent->log_required_g == 999999.0)
            { // GW is source
                e |= fprintf(fout, "E.C. requirements do not apply - Groundwater Is Source\r\n");
                count += 1;
            }
            else
            {
                if (effluent->ec_exempt == TRUE)
                    e |= fprintf(fout, "E.C. not required - raw TOC, raw SUVA, and/or finished TOC <= 2\r\n");
                else
                    e |= fprintf(fout, "E.C. raw TOC, raw SUVA, and finished TOC <= 2 exemptions do not apply\r\n");
                if (effluent->ec_meeting_step1 == TRUE)
                    e |= fprintf(fout, "E.C. Step 1 TOC removal requirement ACHIEVED\r\n");
                else
                    e |= fprintf(fout, "E.C. Step 1 TOC removal requirement NOT ACHIEVED\r\n");
                count += 2;
            }
            //CT Ratio Outputs
            if (swflag == FALSE /*influent->log_required_g == 999999.0*/ && gw_virus_flag == FALSE)
            {
                e |= fprintf(fout, "CT Ratios Not Applicable - Groundwater Is Source\r\n");
                count += 1;
            }
            else if (swflag == TRUE)
            {
                //   if(influent->log_required
                e |= fprintf(fout, "CT Ratios                        \r\n");
                e |= fprintf(fout, "  Virus                (-)      %8.1f%11.1f%11.1f\r\n", influent->ct_ratio_v, effluent->ct_ratio_v, tap->ct_ratio_v);
                e |= fprintf(fout, "  Giardia              (-)      %8.1f%11.1f%11.1f\r\n", influent->ct_ratio, effluent->ct_ratio, tap->ct_ratio);
                e |= fprintf(fout, "  Cryptosporidium*     (-)      %8.1f%11.1f%11.1f\r\n", influent->ct_ratio_c, effluent->ct_ratio_c, tap->ct_ratio_c);
                e |= fprintf(fout, "*Crypto. Disinfection Calcs. Based on Proposed LT2 Rule\r\n");
                count += 5;
            }
            else // (swflag==FALSE && gw_virus_flag==TRUE) **/
            {
                e |= fprintf(fout, "CT Ratios                        \r\n");
                e |= fprintf(fout, "  Virus                (-)      %8.1f%11.1f%11.1f\r\n", influent->ct_ratio_v, effluent->ct_ratio_v, tap->ct_ratio_v);
                e |= fprintf(fout, "CT Ratios N/A for for Giardia & Crypto.\r\n");
                count += 3;
            }

        } /* End of "else if (eff_cntr > 0 && tap_cntr > 0 && eos_cntr ==0)") */

        else if (eff_cntr > 0 && tap_cntr > 0 && eos_cntr > 0)

        { /* Print-out influent, effluent, avg_tap, and end-of-system information */

            e |= fprintf(fout, "pH                     (-)      %8.1f%11.1f%11.1f%11.1f\r\n", influent->pH, effluent->pH, tap->pH, eos->pH);
            e |= fprintf(fout, "Alkalinity      (mg/L as CaCO3) %8.0f%11.0f%11.0f%11.0f\r\n", influent->Alk * MW_CaCO3 / 2, effluent->Alk * MW_CaCO3 / 2, tap->Alk * MW_CaCO3 / 2, eos->Alk * MW_CaCO3 / 2);
            e |= fprintf(fout, "TOC                  (mg/L)     %8.1f%11.1f%11.1f%11.1f\r\n", influent->TOC, effluent->TOC, tap->TOC, eos->TOC);
            e |= fprintf(fout, "UV                   (1/cm)     %8.3f%11.3f%11.3f%11.3f\r\n", influent->UV, effluent->UV_out, tap->UV_out, eos->UV_out);
            e |= fprintf(fout, "(T)SUVA              (1/cm)     %8.1f%11.1f%11.1f%11.1f\r\n", influent->SUVA, effluent->SUVA, tap->SUVA, eos->SUVA);
            e |= fprintf(fout, "Ca Hardness     (mg/L as CaCO3) %8.0f%11.0f%11.0f%11.0f\r\n", influent->Ca_aq * MW_CaCO3, effluent->Ca_aq * MW_CaCO3, tap->Ca_aq * MW_CaCO3, eos->Ca_aq * MW_CaCO3);
            e |= fprintf(fout, "Mg Hardness     (mg/L as CaCO3) %8.0f%11.0f%11.0f%11.0f\r\n", influent->Mg_aq * MW_CaCO3, effluent->Mg_aq * MW_CaCO3, tap->Mg_aq * MW_CaCO3, eos->Mg_aq * MW_CaCO3);
            e |= fprintf(fout, "Ammonia-N            (mg/L)     %8.2f%11.2f%11.2f%11.2f\r\n", influent->NH3 * MW_NH3, effluent->NH3 * MW_NH3, tap->NH3 * MW_NH3, eos->NH3 * MW_NH3);
            e |= fprintf(fout, "Bromide              (ug/L)     %8.0f%11.0f%11.0f%11.0f\r\n", influent->Br * MW_Br * 1000, effluent->Br * MW_Br * 1000, tap->Br * MW_Br * 1000, eos->Br * MW_Br * 1000);
            e |= fprintf(fout, "Free Cl2 Res.     (mg/L as Cl2) %8.1f%11.1f%11.1f%11.1f\r\n", influent->FreeCl2 * MW_Cl2, effluent->FreeCl2 * MW_Cl2, tap->FreeCl2 * MW_Cl2, eos->FreeCl2 * MW_Cl2);
            e |= fprintf(fout, "Chloramine Res.   (mg/L as Cl2) %8.1f%11.1f%11.1f%11.1f\r\n", influent->NH2Cl * MW_Cl2, effluent->NH2Cl * MW_Cl2, tap->NH2Cl * MW_Cl2, eos->NH2Cl * MW_Cl2);
            e |= fprintf(fout, "TTHMs                (ug/L)     %8.0f%11.0f%11.0f%11.0f\r\n", influent->TTHM, effluent->TTHM, tap->TTHM, eos->TTHM);
            e |= fprintf(fout, "HAA5                 (ug/L)     %8.0f%11.0f%11.0f%11.0f\r\n", influent->HAA5, effluent->HAA5, tap->HAA5, eos->HAA5);
            e |= fprintf(fout, "HAA6                 (ug/L)     %8.0f%11.0f%11.0f%11.0f\r\n", influent->HAA6, effluent->HAA6, tap->HAA6, eos->HAA6);
            e |= fprintf(fout, "HAA9                 (ug/L)     %8.0f%11.0f%11.0f%11.0f\r\n", influent->HAA9, effluent->HAA9, tap->HAA9, eos->HAA9);
            e |= fprintf(fout, "TOX                  (ug/L)     %8.0f%11.0f%11.0f%11.0f\r\n", influent->tox, effluent->tox, tap->tox, eos->tox);
            e |= fprintf(fout, "Bromate              (ug/L)     %8.0f%11.0f%11.0f%11.0f\r\n", influent->BrO3, effluent->BrO3, tap->BrO3, eos->BrO3);
            e |= fprintf(fout, "Chlorite             (mg/L)     %8.1f%11.1f%11.1f%11.1f\r\n", influent->chlorite, effluent->chlorite, tap->chlorite, eos->chlorite);
            e |= fprintf(fout, "TOC Removal         (percent)   %19.0f\r\n", toc_rem);
            //E.C. outputs
            if (influent->log_required_g == 999999.0)
            { // GW is source
                e |= fprintf(fout, "E.C. requirements do not apply - Groundwater Is Source\r\n");
                count += 1;
            }
            else
            {
                if (effluent->ec_exempt == TRUE)
                    e |= fprintf(fout, "E.C. not required - raw TOC, raw SUVA, and/or finished TOC <= 2\r\n");
                else
                    e |= fprintf(fout, "E.C. raw TOC, raw SUVA, and finished TOC <= 2 exemptions do not apply\r\n");
                if (effluent->ec_meeting_step1 == TRUE)
                    e |= fprintf(fout, "E.C. Step 1 TOC removal requirement ACHIEVED\r\n");
                else
                    e |= fprintf(fout, "E.C. Step 1 TOC removal requirement NOT ACHIEVED\r\n");
                count += 2;
            }
            //CT Ratio Outputs
            if (swflag == FALSE /*influent->log_required_g == 999999.0*/ && gw_virus_flag == FALSE)
            {
                e |= fprintf(fout, "CT Ratios Not Applicable - Groundwater Is Source\r\n");
                count += 1;
            }
            else if (swflag == TRUE)
            {
                //   if(influent->log_required
                e |= fprintf(fout, "CT Ratios                        \r\n");
                e |= fprintf(fout, "  Virus                (-)      %8.1f%11.1f%11.1f%11.1f\r\n", influent->ct_ratio_v, effluent->ct_ratio_v, tap->ct_ratio_v, eos->ct_ratio_v);
                e |= fprintf(fout, "  Giardia              (-)      %8.1f%11.1f%11.1f%11.1f\r\n", influent->ct_ratio, effluent->ct_ratio, tap->ct_ratio, eos->ct_ratio);
                e |= fprintf(fout, "  Cryptosporidium*     (-)      %8.1f%11.1f%11.1f%11.1f\r\n", influent->ct_ratio_c, effluent->ct_ratio_c, tap->ct_ratio_c, eos->ct_ratio_c);
                e |= fprintf(fout, "*Crypto. Disinfection Calcs. Based on Proposed LT2 Rule\r\n");
                count += 5;
            }
            else // (swflag==FALSE && gw_virus_flag==TRUE) **/
            {
                e |= fprintf(fout, "CT Ratios                        \r\n");
                e |= fprintf(fout, "  Virus               (-)      %8.1f%11.1f%11.1f%11.1f\r\n", influent->ct_ratio_v, effluent->ct_ratio_v, tap->ct_ratio_v, eos->ct_ratio_v);
                e |= fprintf(fout, "CT Ratios N/A for for Giardia & Crypto.\r\n");
                count += 3;
            }

        } /* End of "else if (eff_cntr > 0 && tap_cntr > 0 && eos_cntr > 0)") */

        else if (eff_cntr > 0 && tap_cntr == 0 && eos_cntr > 0)

        { /* Print-out influent, effluent, and end-of-system information */

            e |= fprintf(fout, "pH                     (-)      %8.1f%11.1f%22.1f\r\n", influent->pH, effluent->pH, eos->pH);
            e |= fprintf(fout, "Alkalinity      (mg/L as CaCO3) %8.0f%11.0f%22.0f\r\n", influent->Alk * MW_CaCO3 / 2, effluent->Alk * MW_CaCO3 / 2, eos->Alk * MW_CaCO3 / 2);
            e |= fprintf(fout, "TOC                  (mg/L)     %8.1f%11.1f%22.1f\r\n", influent->TOC, effluent->TOC, eos->TOC);
            e |= fprintf(fout, "UV                   (1/cm)     %8.3f%11.3f%22.3f\r\n", influent->UV, effluent->UV_out, eos->UV_out);
            e |= fprintf(fout, "(T)SUVA              (1/cm)     %8.1f%11.1f%22.1f\r\n", influent->SUVA, effluent->SUVA, eos->SUVA);
            e |= fprintf(fout, "Ca Hardness     (mg/L as CaCO3) %8.0f%11.0f%22.0f\r\n", influent->Ca_aq * MW_CaCO3, effluent->Ca_aq * MW_CaCO3, eos->Ca_aq * MW_CaCO3);
            e |= fprintf(fout, "Mg Hardness     (mg/L as CaCO3) %8.0f%11.0f%22.0f\r\n", influent->Mg_aq * MW_CaCO3, effluent->Mg_aq * MW_CaCO3, eos->Mg_aq * MW_CaCO3);
            e |= fprintf(fout, "Ammonia-N            (mg/L)     %8.2f%11.2f%22.2f\r\n", influent->NH3 * MW_NH3, effluent->NH3 * MW_NH3, eos->NH3 * MW_NH3);
            e |= fprintf(fout, "Bromide              (ug/L)     %8.0f%11.0f%22.0f\r\n", influent->Br * MW_Br * 1000, effluent->Br * MW_Br * 1000, eos->Br * MW_Br * 1000);
            e |= fprintf(fout, "Free Cl2 Res.     (mg/L as Cl2) %8.1f%11.1f%22.1f\r\n", influent->FreeCl2 * MW_Cl2, effluent->FreeCl2 * MW_Cl2, eos->FreeCl2 * MW_Cl2);
            e |= fprintf(fout, "Chloramine Res.   (mg/L as Cl2) %8.1f%11.1f%22.1f\r\n", influent->NH2Cl * MW_Cl2, effluent->NH2Cl * MW_Cl2, eos->NH2Cl * MW_Cl2);
            e |= fprintf(fout, "TTHMs                (ug/L)     %8.0f%11.0f%22.0f\r\n", influent->TTHM, effluent->TTHM, eos->TTHM);
            e |= fprintf(fout, "HAA5                 (ug/L)     %8.0f%11.0f%22.0f\r\n", influent->HAA5, effluent->HAA5, eos->HAA5);
            e |= fprintf(fout, "HAA6                 (ug/L)     %8.0f%11.0f%22.0f\r\n", influent->HAA6, effluent->HAA6, eos->HAA6);
            e |= fprintf(fout, "HAA9                 (ug/L)     %8.0f%11.0f%22.0f\r\n", influent->HAA9, effluent->HAA9, eos->HAA9);
            e |= fprintf(fout, "TOX                  (ug/L)     %8.0f%11.0f%22.0f\r\n", influent->tox, effluent->tox, eos->tox);
            e |= fprintf(fout, "Bromate              (ug/L)     %8.0f%11.0f%22.0f\r\n", influent->BrO3, effluent->BrO3, eos->BrO3);
            e |= fprintf(fout, "Chlorite             (mg/L)     %8.1f%11.1f%22.1f\r\n", influent->chlorite, effluent->chlorite, eos->chlorite);
            e |= fprintf(fout, "TOC Removal         (percent)   %19.0f\r\n", toc_rem);
            //E.C. outputs
            if (influent->log_required_g == 999999.0)
            { // GW is source
                e |= fprintf(fout, "E.C. requirements do not apply - Groundwater Is Source\r\n");
                count += 1;
            }
            else
            {
                if (effluent->ec_exempt == TRUE)
                    e |= fprintf(fout, "E.C. not required - raw TOC, raw SUVA, and/or finished TOC <= 2\r\n");
                else
                    e |= fprintf(fout, "E.C. raw TOC, raw SUVA, and finished TOC <= 2 exemptions do not apply\r\n");
                if (effluent->ec_meeting_step1 == TRUE)
                    e |= fprintf(fout, "E.C. Step 1 TOC removal requirement ACHIEVED\r\n");
                else
                    e |= fprintf(fout, "E.C. Step 1 TOC removal requirement NOT ACHIEVED\r\n");
                count += 2;
            }
            //CT Ratio Outputs
            if (swflag == FALSE /*influent->log_required_g == 999999.0*/ && gw_virus_flag == FALSE)
            {
                e |= fprintf(fout, "CT Ratios Not Applicable - Groundwater Is Source\r\n");
                count += 1;
            }
            else if (swflag == TRUE)
            {
                //   if(influent->log_required
                e |= fprintf(fout, "CT Ratios                        \r\n");
                e |= fprintf(fout, "  Virus                (-)      %8.1f%11.1f%22.1f\r\n", influent->ct_ratio_v, effluent->ct_ratio_v, eos->ct_ratio_v);
                e |= fprintf(fout, "  Giardia              (-)      %8.1f%11.1f%22.1f\r\n", influent->ct_ratio, effluent->ct_ratio, eos->ct_ratio);
                e |= fprintf(fout, "  Cryptosporidium*     (-)      %8.1f%11.1f%22.1f\r\n", influent->ct_ratio_c, effluent->ct_ratio_c, eos->ct_ratio_c);
                e |= fprintf(fout, "*Crypto. Disinfection Calcs. Based on Proposed LT2 Rule\r\n");
                count += 5;
            }
            else // (swflag==FALSE && gw_virus_flag==TRUE) **/
            {
                e |= fprintf(fout, "CT Ratios                        \r\n");
                e |= fprintf(fout, "  Virus               (-)      %8.1f%11.1f%22.1f\r\n", influent->ct_ratio_v, effluent->ct_ratio_v, eos->ct_ratio_v);
                e |= fprintf(fout, "CT Ratios N/A for for Giardia & Crypto.\r\n");
                count += 3;
            }

        } /* End of "else if (eff_cntr > 0 && tap_cntr == 0 && eos_cntr ==0)") */

        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        count += 21;

        /*This code describes if CT_Ratio =1 due to physical removal*/
        if (influent->ct_ratio_v == 1.0)
        {
            e |= fprintf(fout, "Virus CT Ratio = 1 because other credits met full disinfection requirements\r\n");
            count += 1;
        }
        if (influent->ct_ratio == 1.0)
        {
            e |= fprintf(fout, "Giardia CT Ratio = 1 because other credits met full disinfection requirements\r\n");
            count += 1;
        }
        if (influent->ct_ratio_c == 1.0)
        {
            e |= fprintf(fout, "Crypto. CT Ratio = 1 because other credits met full disinfection requirements\r\n");
            count += 1;
        }
        if (influent->ct_ratio == 1.0 || influent->ct_ratio_v == 1.0 || influent->ct_ratio_c == 1.0)
        {
            e |= fprintf(fout, "\r\n");
            count += 1;
        }

        // if( (wait_gui!=NULL) && (count > (max_line_count - 10)) )
        //   {
        //     wait_gui("Press Continue");
        //     count = 0;
        //   }

        /*********************************************************************************************************************************************/

        e |= fprintf(fout, "\r\n\n");
        e |= fprintf(fout, "                                      Table 2                                \r\n");
        e |= fprintf(fout, "                             Selected Input Parameters                       \r\n");
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        e |= fprintf(fout, "Parameter                                     Value  Units                   \r\n");
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        if (FirstUnitProcess(train)->type == INFLUENT)
        {
            e |= fprintf(fout, "TEMPERATURES\r\n");
            e |= fprintf(fout, "     Average%38.1f  (deg. C)\r\n", FirstUnitProcess(train)->data.influent->temp);
            e |= fprintf(fout, "     Minimum%38.1f  (deg. C)\r\n", FirstUnitProcess(train)->data.influent->low_temp);
            e |= fprintf(fout, "PLANT FLOW RATES\r\n");
            e |= fprintf(fout, "     Average%38.3f  (mgd)\r\n", FirstUnitProcess(train)->data.influent->avg_flow);
            e |= fprintf(fout, "     Peak Hourly%34.3f  (mgd)\r\n", FirstUnitProcess(train)->data.influent->peak_flow);
            e |= fprintf(fout, "\r\n");
            e |= fprintf(fout, "DISINFECTION INPUTS/CALCULATED VALUES\r\n");
            if (swflag == TRUE)
                e |= fprintf(fout, "Surface Water Plant?                           TRUE\r\n");
            else
                e |= fprintf(fout, "Surface Water Plant?                          FALSE\r\n");

            if (swflag == FALSE && gw_virus_flag == TRUE)
                e |= fprintf(fout, "Groundwater Disinfection Required?             TRUE\r\n");
            else if (swflag == TRUE) /*do nothing*/
                ;
            else
                e |= fprintf(fout, "Groundwater Disinfection Required?            FALSE\r\n");

            if (swflag == TRUE || (swflag == FALSE && gw_virus_flag == TRUE))
            {
                e |= fprintf(fout, "Giardia: Total Disinfection Credit Required%7.1f  (logs)\r\n", tot_dis_req_g);
                e |= fprintf(fout, "Giardia: Credit Achieved (other than by CT)%7.1f  (logs)\r\n", tot_giardia_lr);
                if (tot_dis_req_g - tot_giardia_lr >= 0.0)
                    e |= fprintf(fout, "Giardia: Inactivation Credit by CT Required%7.1f  (logs)\r\n", tot_dis_req_g - tot_giardia_lr);
                else
                    e |= fprintf(fout, "Giardia: Inactivation Credit by CT Required%7.1f  (logs)\r\n", 0.0);
                e |= fprintf(fout, "\r\n");

                e |= fprintf(fout, "Virus: Total Disinfection Credit Required%9.1f  (logs)\r\n", tot_dis_req_v);
                e |= fprintf(fout, "Virus: Credit Achieved (other than by CT)%9.1f  (logs)\r\n", tot_virus_lr);
                if (tot_dis_req_v - tot_virus_lr >= 0.0)
                    e |= fprintf(fout, "Virus: Inactivation Credit by CT Required%9.1f  (logs)\r\n", tot_dis_req_v - tot_virus_lr);
                else
                    e |= fprintf(fout, "Virus: Inactivation Credit by CT Required%9.1f  (logs)\r\n", 0.0);
                e |= fprintf(fout, "\r\n");

                if (FirstUnitProcess(train)->data.influent->lt2_wscp_flag == TRUE)
                    lt2_wscp_credit = 0.5;
                e |= fprintf(fout, "Crypto.: Total Disinfection Credit Required%7.1f  (logs)\r\n", tot_dis_req_c);
                e |= fprintf(fout, "Crypto.: Credit Achieved (other than by CT)%7.1f  (logs)\r\n", tot_crypto_lr + lt2_wscp_credit);

                if ((tot_dis_req_c - tot_crypto_lr - lt2_wscp_credit) >= 0.0 && bin34_inactreqd == 0.0)
                    e |= fprintf(fout, "Crypto.: Inactivation Credit by CT Required%7.1f  (logs)\r\n", tot_dis_req_c - tot_crypto_lr - lt2_wscp_credit);
                else if ((tot_dis_req_c - tot_crypto_lr - lt2_wscp_credit) < 0.0 && bin34_inactreqd == 0.0)
                    e |= fprintf(fout, "Crypto.: Inactivation Credit by CT Required%7.1f  (logs)\r\n", 0.0);
                else //(bin34_inactreqd > 0.0)
                {
                    e |= fprintf(fout, "Crypto.: Inactivation Credit by CT Required%7.1f  (logs)\r\n", bin34_inactreqd);
                    e |= fprintf(fout, "Crypto. CT Credit Set by 1.0-log Requirement of Bins 3,4\r\n");
                    count++;
                }
                e |= fprintf(fout, "\r\n");
                count += 16;
            }
            else
            {
                e |= fprintf(fout, "Source is a Groundwater Not Requiring Disinfection CT\r\n");
                e |= fprintf(fout, "\r\n");
                count += 2;
            }
        }
        else
        {
            e |= fprintf(fout, "Input Error: 'Influent' must be entered as first unit process\r\n");
            e |= fprintf(fout, "\r\n");
            count += 2;
        }

        /* Echo Log Removals (or, in the case of UV_DIS, Inactivations) Inputted/Creditted Here*/
        e |= fprintf(fout, "\nDISINFECT. CREDITS (not incl. CT): Giardia    Virus      Crypto\r\n");
        e |= fprintf(fout, "(in order of appearance)\r\n");
        for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
        {
            switch (unit->type)
            {
            case FILTER: //////////////////////////////////////////////////////////////FILTER////////////////////
                if (unit->data.filter->filt_stage == 1 && conv_filtflag == TRUE)
                {
                    e |= fprintf(fout, "%-16s%23.1f  %9.1f    %8.1f\r\n", UnitProcessTable[unit->type].name,
                                 unit->data.filter->giardia_lr_conv, unit->data.filter->virus_lr_conv, unit->data.filter->crypto_lr_conv);
                    count++;
                }
                else if (unit->data.filter->filt_stage == 1 && dir_filtflag == TRUE)
                {
                    e |= fprintf(fout, "%-16s%23.1f  %9.1f    %8.1f\r\n", UnitProcessTable[unit->type].name,
                                 unit->data.filter->giardia_lr_df, unit->data.filter->virus_lr_df, unit->data.filter->crypto_lr_df);
                    count++;
                }
                else if (unit->data.filter->filt_stage == 2)
                {
                    e |= fprintf(fout, "%-16s -2nd stage filt. %5.1f  %9.1f    %8.1f\r\n", UnitProcessTable[unit->type].name,
                                 zero_lr, zero_lr, unit->data.filter->crypto_lr_2);
                    count++;
                }
                else /*filt_stage==0*/
                    ;

                if (unit->data.filter->filt_stage == 1 && unit->data.filter->ife_turb_flag == TRUE)
                {
                    e |= fprintf(fout, "%-27s%12.1f  %9.1f    %8.1f\r\n", optfilter,
                                 zero_lr, zero_lr, ife_lr_credit);
                    count++;
                }
                else if (unit->data.filter->filt_stage == 1 && unit->data.filter->cfe_turb_flag == TRUE)
                {
                    e |= fprintf(fout, "%-27s%12.1f  %9.1f    %8.1f\r\n", optfilter,
                                 zero_lr, zero_lr, cfe_lr_credit);
                    count++;
                }
                else /*filt_stage==0 or 2 || xfe_turb_flag==FALSE*/
                    ;
                break;

            case DE_FILTER: /////////////////////////////////////////////////////DE_FILTER////////////////////
                if (unit->data.def->filt_stage == 1)
                {
                    e |= fprintf(fout, "%-16s%23.1f  %9.1f    %8.1f\r\n", UnitProcessTable[unit->type].name,
                                 unit->data.def->giardia_lr, unit->data.def->virus_lr, unit->data.def->crypto_lr_1);
                    count++;
                }
                else if (unit->data.def->filt_stage == 2)
                {
                    e |= fprintf(fout, "%-16s -2nd stage filt. %5.1f  %9.1f    %8.1f\r\n", UnitProcessTable[unit->type].name,
                                 zero_lr, zero_lr, unit->data.def->crypto_lr_2);
                    count++;
                }
                else /*filt_stage==0*/
                    ;
                break;

            case SLOW_FILTER: /////////////////////////////////////////////////////SLOW_FILTER////////////////////
                if (unit->data.ssf->filt_stage == 1)
                {
                    e |= fprintf(fout, "%-16s%23.1f  %9.1f    %8.1f\r\n", UnitProcessTable[unit->type].name,
                                 unit->data.ssf->giardia_lr, unit->data.ssf->virus_lr, unit->data.ssf->crypto_lr_1);
                    count++;
                }
                else if (unit->data.ssf->filt_stage == 2)
                {
                    e |= fprintf(fout, "%-16s -2nd stage filt. %5.1f  %9.1f    %8.1f\r\n", UnitProcessTable[unit->type].name,
                                 zero_lr, zero_lr, unit->data.ssf->crypto_lr_2);
                    count++;
                }
                else /*filt_stage==0*/
                    ;
                break;

            case CART_FILTER:
            case BAG_FILTER: //////////////////////////////////////////////////////BAG_FILTER OR CART_FILTER////////////////////
                if (unit->data.altf->filt_stage == 1)
                {
                    e |= fprintf(fout, "%-16s%23.1f  %9.1f    %8.1f\r\n", UnitProcessTable[unit->type].name,
                                 unit->data.altf->giardia_lr, unit->data.altf->virus_lr, unit->data.altf->crypto_lr_1);
                    count++;
                }
                else if (unit->data.altf->filt_stage == 2)
                {
                    e |= fprintf(fout, "%-16s -2nd stage filt. %5.1f  %9.1f    %8.1f\r\n", UnitProcessTable[unit->type].name,
                                 zero_lr, zero_lr, unit->data.altf->crypto_lr_2);
                    count++;
                }
                else /*filt_stage==0*/
                    ;
                break;

            case MFUF_UP: /////////////////////////////////////////////////////////////MFUF_UP////////////////////
                if (unit->data.mfuf->filt_stage == 1)
                {
                    e |= fprintf(fout, "%-16s%23.1f  %9.1f    %8.1f\r\n", UnitProcessTable[unit->type].name,
                                 unit->data.mfuf->giardia_lr, unit->data.mfuf->virus_lr, unit->data.mfuf->crypto_lr_1);
                    count++;
                }
                else if (unit->data.mfuf->filt_stage == 2)
                {
                    e |= fprintf(fout, "%-16s -2nd stage filt. %5.1f  %9.1f    %8.1f\r\n", UnitProcessTable[unit->type].name,
                                 zero_lr, zero_lr, unit->data.mfuf->crypto_lr_2);
                    count++;
                }
                else /*filt_stage==0*/
                    ;
                break;

            case NF_UP: /////////////////////////////////////////////////////////////////////NF////////////////////
                if (unit->data.nf->filt_stage == 1 && unit->data.nf->treat_fraction >= 1.0)
                {
                    e |= fprintf(fout, "%-16s%23.1f  %9.1f    %8.1f\r\n", UnitProcessTable[unit->type].name,
                                 unit->data.nf->giardia_lr, unit->data.nf->virus_lr, unit->data.nf->crypto_lr);
                    count++;
                }
                else if (unit->data.nf->filt_stage == 2 && unit->data.nf->treat_fraction >= 1.0)
                {
                    e |= fprintf(fout, "%-16s -2nd stage filt. %5.1f  %9.1f    %8.1f\r\n", UnitProcessTable[unit->type].name,
                                 zero_lr, zero_lr, unit->data.nf->crypto_lr);
                    count++;
                }
                else /*filt_stage==0 || treat_fraction < 1.0*/
                    ;
                break;

            case BANK_FILTER: ///////////////////////////////////////////////////////////BANK_FILTER//////////////
                if (bankfflag == 1 && bankclosed == FALSE)
                {
                    e |= fprintf(fout, "%-16s%23.1f  %9.1f    %8.1f\r\n", UnitProcessTable[unit->type].name,
                                 zero_lr, zero_lr, unit->data.bankf->crypto_lr_close);
                    count++;
                }
                else if (bankfflag == 2 && bankclosed == FALSE)
                {
                    e |= fprintf(fout, "%-16s%23.1f  %9.1f    %8.1f\r\n", UnitProcessTable[unit->type].name,
                                 zero_lr, zero_lr, unit->data.bankf->crypto_lr_far);
                    count++;
                }
                else /*bankfflag == 0 || bankclosed == TRUE */
                    ;
                bankclosed = TRUE;
                break;

            case GAC: /////////////////////////////////////////////////////////////GAC////////////////////
                if (unit->data.gac->filt_stage == 2)
                {
                    e |= fprintf(fout, "%-16s -2nd stage filt. %5.1f  %9.1f    %8.1f\r\n", UnitProcessTable[unit->type].name,
                                 zero_lr, zero_lr, unit->data.gac->crypto_lr_2);
                    count++;
                }
                break;

            case PRESED_BASIN: /////////////////////////////////////////////////////////PRESED_BASIN//////////////
                if (lt2presedflag == TRUE && presed_cntr < 1)
                {
                    e |= fprintf(fout, "%-16s%23.1f  %9.1f    %8.1f\r\n", UnitProcessTable[unit->type].name,
                                 zero_lr, zero_lr, unit->data.presed->crypto_lr);
                    count++;
                }
                presed_cntr++;
                break;

            case SETTLING_BASIN: ///////////////////////////////////////////////2ND STAGE SOFTENING//////////////
                if (soft2flag == TRUE && soft2_cntr < 1)
                {
                    e |= fprintf(fout, "%-17s%22.1f  %9.1f    %8.1f\r\n", soft2basin,
                                 zero_lr, zero_lr, soft2_credit);
                    count++;
                    soft2_cntr++;
                }
                break;

            case UV_DIS: /////////////////////////////////////////////////////////UV_DIS////////////////////
                e |= fprintf(fout, "%-16s%23.1f  %9.1f    %8.1f\r\n", UnitProcessTable[unit->type].name,
                             unit->data.uvdis->giardia_li, unit->data.uvdis->virus_li, unit->data.uvdis->crypto_li);
                count++;
                break;

            case INFLUENT: /////////////////////////////////////////////////INFLUENT-WSCP//////////////
                if (unit->data.influent->lt2_wscp_flag == TRUE)
                {
                    e |= fprintf(fout, "%-25s%14.1f  %9.1f    %8.1f\r\n", wscp,
                                 zero_lr, zero_lr, wscp_credit);
                    count++;
                }
                break;

            default:
                break;
            } /*End "switch" for unit processes*/
        }     /*End "for loop" over process train for echoing UP log removals and other lr credits*/

        e |= fprintf(fout, "CHEMICAL DOSES\r\n"); //********************************CHEM DOSE INPUTS
        e |= fprintf(fout, "(in order of appearance)\r\n");
        for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
        {
            switch (unit->type)
            {
            case ALUM:
                e |= fprintf(fout, "%-16s%34.1f  %10s\r\n", UnitProcessTable[unit->type].name,
                             unit->data.alum->dose, UnitProcessTable[unit->type].data->units);
                count++;
                break;
            case IRON:
                e |= fprintf(fout, "%-16s%34.1f  %10s\r\n", UnitProcessTable[unit->type].name,
                             unit->data.iron->dose, UnitProcessTable[unit->type].data->units);
                count++;
                break;
            case CHLORINE_DIOXIDE:
                e |= fprintf(fout, "%-16s%34.1f  %10s\r\n", UnitProcessTable[unit->type].name,
                             unit->data.clo2->dose, UnitProcessTable[unit->type].data->units);
                count++;
                break;
            case LIME:
                e |= fprintf(fout, "%-16s%34.1f  %10s\r\n", UnitProcessTable[unit->type].name,
                             unit->data.lime->dose, UnitProcessTable[unit->type].data->units);
                count++;
                break;
            case HYPOCHLORITE:
                e |= fprintf(fout, "%-16s%34.1f  %10s\r\n", UnitProcessTable[unit->type].name,
                             unit->data.chemical->naocl, UnitProcessTable[unit->type].data->units);
                count++;
                break;
            case CHLORINE:
                e |= fprintf(fout, "%-16s%34.1f  %10s\r\n", UnitProcessTable[unit->type].name,
                             unit->data.chemical->chlor, UnitProcessTable[unit->type].data->units);
                count++;
                break;
            case SULFURIC_ACID:
                e |= fprintf(fout, "%-16s%34.1f  %10s\r\n", UnitProcessTable[unit->type].name,
                             unit->data.chemical->h2so4, UnitProcessTable[unit->type].data->units);
                count++;
                break;
            case SODA_ASH:
                e |= fprintf(fout, "%-16s%34.1f  %10s\r\n", UnitProcessTable[unit->type].name,
                             unit->data.chemical->soda, UnitProcessTable[unit->type].data->units);
                count++;
                break;
            case AMMONIA:
                e |= fprintf(fout, "%-16s%34.1f  %10s\r\n", UnitProcessTable[unit->type].name,
                             unit->data.chemical->nh3, UnitProcessTable[unit->type].data->units);
                count++;
                break;
            case PERMANGANATE:
                e |= fprintf(fout, "%-16s%34.1f  %10s\r\n", UnitProcessTable[unit->type].name,
                             unit->data.chemical->kmno4, UnitProcessTable[unit->type].data->units);
                count++;
                break;
            case CARBON_DIOXIDE:
                e |= fprintf(fout, "%-16s%34.1f  %10s\r\n", UnitProcessTable[unit->type].name,
                             unit->data.chemical->co2, UnitProcessTable[unit->type].data->units);
                count++;
                break;
            case OZONE:
                e |= fprintf(fout, "%-16s%34.1f  %10s\r\n", UnitProcessTable[unit->type].name,
                             unit->data.chemical->o3, UnitProcessTable[unit->type].data->units);
                count++;
                break;
            case SODIUM_HYDROXIDE:
                e |= fprintf(fout, "%-16s%34.1f  %10s\r\n", UnitProcessTable[unit->type].name,
                             unit->data.chemical->naoh, UnitProcessTable[unit->type].data->units);
                count++;
                break;
            case SULFUR_DIOXIDE:
                e |= fprintf(fout, "%-16s%34.1f  %10s\r\n", UnitProcessTable[unit->type].name,
                             unit->data.chemical->so2, UnitProcessTable[unit->type].data->units);
                count++;
                break;
            case AMMONIUM_SULFATE:
                e |= fprintf(fout, "%-16s%34.1f  %10s\r\n", UnitProcessTable[unit->type].name,
                             unit->data.chemical->nh3, UnitProcessTable[unit->type].data->units);
                count++;
                break;
            default:
                break;
            }
        } //end of echoing chemical dose inputs

        /*Echo inputted T10 and T50 ratios and Volumes here*/
        e |= fprintf(fout, "\nPROCESS HYDRAULIC PARAMETERS:    T10/Tth    T50/Tth    VOL. (MG)\r\n");
        e |= fprintf(fout, "(in order of appearance)\r\n");
        for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
        {
            switch (unit->type)
            {
            case RAPID_MIX:
            case SLOW_MIX:
            case SETTLING_BASIN:
            case BASIN:
            case CONTACT_TANK:
            case CLEARWELL:
            case O3_CONTACTOR:
                e |= fprintf(fout, "%-16s%23.1f  %9.1f    %10.4f\r\n", UnitProcessTable[unit->type].name,
                             unit->data.basin->sb_t_ten, unit->data.basin->sb_mean, unit->data.basin->volume);
                count++;
                break;
            case FILTER:
                e |= fprintf(fout, "%-16s%23.1f  %9.1f    %10.4f\r\n", UnitProcessTable[unit->type].name,
                             unit->data.filter->fi_t_ten, unit->data.filter->fi_mean, unit->data.filter->volume);
                count++;
                break;
            case SLOW_FILTER:
                e |= fprintf(fout, "%-16s%23.1f  %9.1f    %10.4f\r\n", UnitProcessTable[unit->type].name,
                             unit->data.ssf->fi_t_ten, unit->data.ssf->fi_mean, unit->data.ssf->volume);
                count++;
                break;
            case DE_FILTER:
                e |= fprintf(fout, "%-16s%23.1f  %9.1f    %10.4f\r\n", UnitProcessTable[unit->type].name,
                             unit->data.def->fi_t_ten, unit->data.def->fi_mean, unit->data.def->volume);
                count++;
                break;
            case PRESED_BASIN:
                e |= fprintf(fout, "%-16s%23.1f  %9.1f    %10.4f\r\n", UnitProcessTable[unit->type].name,
                             unit->data.presed->sb_t_ten, unit->data.presed->sb_mean, unit->data.presed->volume);
                count++;
                break;
            default:
                break;
            } /*End "switch"*/
        }     /*End "for loop" over process train for echoing unit process hydraulic data*/

        /*Echo inputted GAC, MFUF_UP, and NF_UP operating/design parameters here*/
        for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
        {
            if (unit->type == GAC)
            {
                e |= fprintf(fout, "\nGAC OPERATION INPUTS\r\n");
                e |= fprintf(fout, "Empty Bed Contact Time%28.1f  (minutes)\r\n", unit->data.gac->ebct);
                e |= fprintf(fout, "Reactivation Frequency%28.1f  (days)\r\n", unit->data.gac->regen);
                e |= fprintf(fout, "Sys. Config. ('S'= Single, 'B'= Blended)%10c\r\n", unit->data.gac->config);
                count += 5;
                if (unit->data.gac->config == 'S')
                {
                    e |= fprintf(fout, "Eff. TOC Calc. ('M'= Max., 'A'= Avg.)%13c\r\n", unit->data.gac->toc_calc);
                    count += 1;
                } /*End "if config = S"*/
            }     /*End "if type == GAC"*/

            if (unit->type == MFUF_UP)
            {
                e |= fprintf(fout, "\nMF/UF OPERATION INPUTS\r\n");
                e |= fprintf(fout, "Recovery%42.1f  (%%)\r\n", unit->data.mfuf->recover);
                count += 3;
            } /*End "if type ==  MFUF_UP"*/

            if (unit->type == NF_UP)
            {
                e |= fprintf(fout, "\nNANOFILTRATION OPERATION INPUTS\r\n");
                e |= fprintf(fout, "Molecular Weight Cut-off%26.1f  (Daltons)\r\n", unit->data.nf->mwc);
                e |= fprintf(fout, "Recovery%42.1f  (%%)\r\n", unit->data.nf->recover);
                e |= fprintf(fout, "TOC Removal%39.1f  (%%)\r\n", unit->data.nf->toc_rem);
                e |= fprintf(fout, "UVA Removal%39.1f  (%%)\r\n", unit->data.nf->uva_rem);
                e |= fprintf(fout, "Bromide Removal%35.1f  (%%)\r\n", unit->data.nf->br_rem);
                e |= fprintf(fout, "Flow Fraction Treated%29.1f  (%%)\r\n", unit->data.nf->treat_fraction);
                count += 8;
            } /*End "if type ==  NF_UP"*/

            if (unit->type == SLOW_FILTER)
            {
                e |= fprintf(fout, "\nSLOW SAND FILTER OPERATION INPUTS\r\n");
                e |= fprintf(fout, "TOC Removal at Temp. < 10 C%23.1f  (%%)\r\n", unit->data.ssf->toc_rem_lowt);
                e |= fprintf(fout, "TOC Removal at 10 C  <= Temp. < 20 C%14.1f  (%%)\r\n", unit->data.ssf->toc_rem_midt);
                e |= fprintf(fout, "TOC Removal at 20 C <= Temp.%22.1f  (%%)\r\n", unit->data.ssf->toc_rem_hight);
                count += 5;
            } /*End "if type ==  NF_UP"*/

        } /*End for loop over process train*/

        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n\n");
        count += 12;

        /*************************************************************************************************************/

        /* Print pH,TOC,UV,Alk, Temp, Cl2, NH2Cl, NH3 */
        e |= fprintf(fout, "\r\n");
        e |= fprintf(fout, "                                  Table 3                                    \r\n");
        e |= fprintf(fout, "                       Predicted Water Quality Profile                       \r\n");
        e |= fprintf(fout, "         At Plant Flow (%4.1f MGD) and Influent Temperature (%3.1f C)        \r\n", influent->Flow, influent->DegK - 273.15);
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        e |= fprintf(fout, "                                                           | Residence Time |\r\n");
        e |= fprintf(fout, "                 pH   TOC     UVA   (T)SUVA   Cl2    NH2Cl | Process|  Cum. |\r\n");
        e |= fprintf(fout, "Location        (-)  (mg/L)  (1/cm) (L/mg-m) (mg/L) (mg/L) |  (hrs) | (hrs) |\r\n");
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        count += 9;

        /* what are the results at the end of each process unit each process has a
   eff pointer that points to a unit holding the effluent results for each
   unit as specified in the EFFLUENT STRUCTURE if we use the standard UNIT
   Structure. We establish a new variable here to use for chlorine residual
   for all subsequent process following WTP_Effluent Process called
   clres_aft_wtp  */

        for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
        {
            eff = &unit->eff;

            /*    if(unit->type == WTP_EFFLUENT)       new chlorine residual 
	   clres_aft_wtp = eff->FreeCl2 * MW_Cl2; */

            e |= fprintf(fout, "%-16s", UnitProcessTable[unit->type].name);
            e |= fprintf(fout, "%4.1f", eff->pH);
            e |= fprintf(fout, "%6.1f", eff->TOC);
            e |= fprintf(fout, "%9.3f", eff->UV_out);
            e |= fprintf(fout, "%7.1f", eff->SUVA);
            e |= fprintf(fout, "%8.1f", eff->FreeCl2 * MW_Cl2);
            e |= fprintf(fout, "%7.1f", eff->NH2Cl * MW_Cl2);
            e |= fprintf(fout, "%9.2f", eff->processtime);
            e |= fprintf(fout, "%9.2f", eff->traintime);
            e |= fprintf(fout, "\r\n");
            count++;
        }
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        count++;
        e |= fprintf(fout, "TOC Removal (percent):   %4.0f\r\n", toc_rem);
        count++;
        //E.C. outputs
        if (influent->log_required_g == 999999.0)
        { // GW is source
            e |= fprintf(fout, "E.C. requirements do not apply - Groundwater Is Source\r\n");
            count += 1;
        }
        else
        {
            if (effluent->ec_exempt == TRUE)
                e |= fprintf(fout, "E.C. not required - raw TOC, raw SUVA, and/or finished TOC <= 2\r\n");
            else
                e |= fprintf(fout, "E.C. raw TOC, raw SUVA, and finished TOC <= 2 exemptions do not apply\r\n");
            if (effluent->ec_meeting_step1 == TRUE)
                e |= fprintf(fout, "E.C. Step 1 TOC removal requirement ACHIEVED\r\n");
            else
                e |= fprintf(fout, "E.C. Step 1 TOC removal requirement NOT ACHIEVED\r\n");
            count += 2;
        }
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        count++;

        // if( (wait_gui!=NULL) && (count > (max_line_count - 10)) )
        //   {
        //     wait_gui("Press Continue");
        //     count = 0;
        //   }

        /* Print Br, Ca, Mg, Inactivation, Solids */
        e |= fprintf(fout, "\r\n\n");
        e |= fprintf(fout, "                                     Table 4                                 \r\n");
        e |= fprintf(fout, "                         Predicted Water Quality Profile                     \r\n");
        e |= fprintf(fout, "         At Plant Flow (%4.1f MGD) and Influent Temperature (%3.1f C)        \r\n", influent->Flow, influent->DegK - 273.15);
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        e |= fprintf(fout, "                               Calcium    Magnesium                          \r\n");
        e |= fprintf(fout, "                   pH   Alk    Hardness   Hardness   Solids  NH3-N  Bromide  \r\n");
        e |= fprintf(fout, "Location          (-)  (mg/L)   (mg/L)     (mg/L)    (mg/L)  (mg/L)  (ug/L)  \r\n");
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        count += 8;

        for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
        {
            eff = &unit->eff;
            e |= fprintf(fout, "%-16s", UnitProcessTable[unit->type].name);
            e |= fprintf(fout, "%5.1f", eff->pH);
            e |= fprintf(fout, "%6.1f", eff->Alk * MW_CaCO3 / 2);
            e |= fprintf(fout, "%10.0f", eff->Ca_aq * MW_CaCO3);
            e |= fprintf(fout, "%10.0f", eff->Mg_aq * MW_CaCO3);
            e |= fprintf(fout, "%11.1f", eff->solids);
            e |= fprintf(fout, "%8.1f", eff->NH3 * MW_NH3);
            e |= fprintf(fout, "%8.0f", eff->Br * MW_Br * 1000.0);
            e |= fprintf(fout, "\r\n");
            count++;
        }
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        count++;
        //  if( (wait_gui!=NULL) && (count > (max_line_count - 10)) )
        //    {

        //      wait_gui("Press Continue");
        //      count = 0;
        //    }

        /* Print THM's */
        e |= fprintf(fout, "\r\n");
        e |= fprintf(fout, "                                   Table 5                                   \r\n");
        e |= fprintf(fout, "                  Predicted Trihalomethanes and other DBPs                   \r\n");
        e |= fprintf(fout, "             At Average Flow (%4.1f MGD) and Temperature (%3.1f C)           \r\n", influent->Flow, influent->DegK - 273.15);
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        e |= fprintf(fout, "                  BrO3- ClO2-   TOX  |CHCl3   BDCM     DBCM    CHBr3   TTHMs \r\n");
        e |= fprintf(fout, "Location         (ug/L) (mg/L) (ug/L)|(ug/L) (ug/L)   (ug/L)   (ug/L) (ug/L) \r\n");
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        count += 8;

        for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
        {
            eff = &unit->eff;
            e |= fprintf(fout, "%-16s", UnitProcessTable[unit->type].name);
            e |= fprintf(fout, "%6.0f", eff->BrO3);
            e |= fprintf(fout, "%7.1f", eff->chlorite);
            e |= fprintf(fout, "%7.0f", eff->tox);
            e |= fprintf(fout, "%7.0f", eff->CHCl3);
            e |= fprintf(fout, "%7.0f", eff->CHBrCl2);
            e |= fprintf(fout, "%9.0f", eff->CHBr2Cl);
            e |= fprintf(fout, "%8.0f", eff->CHBr3);
            e |= fprintf(fout, "%8.0f", eff->TTHM);
            e |= fprintf(fout, "\r\n");
            count++;
        }
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        count++;

        //  if( (wait_gui!=NULL) && (count > (max_line_count - 10)) )
        //    {

        //      wait_gui("Press Continue");
        //      count = 0;
        //    }

        /* Print HaloAcetic Acids */
        e |= fprintf(fout, "\r\n");
        e |= fprintf(fout, "                                   Table 6                                   \r\n");
        e |= fprintf(fout, "                   Predicted Haloacetic Acids - through HAA5                 \r\n");
        e |= fprintf(fout, "             At Average Flow (%4.1f MGD) and Temperature (%3.1f C)           \r\n", influent->Flow, influent->DegK - 273.15);
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        e |= fprintf(fout, "                  MCAA   DCAA   TCAA   MBAA   DBAA   HAA5                    \r\n");
        e |= fprintf(fout, "Location         (ug/L) (ug/L) (ug/L) (ug/L) (ug/L) (ug/L)                   \r\n");
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        count += 8;

        for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
        {
            eff = &unit->eff;
            e |= fprintf(fout, "%-16s", UnitProcessTable[unit->type].name);
            e |= fprintf(fout, "%6.0f", eff->MCAA);
            e |= fprintf(fout, "%7.0f", eff->DCAA);
            e |= fprintf(fout, "%7.0f", eff->TCAA);
            e |= fprintf(fout, "%7.0f", eff->MBAA);
            e |= fprintf(fout, "%7.0f", eff->DBAA);
            e |= fprintf(fout, "%7.0f", eff->HAA5);
            e |= fprintf(fout, "\r\n");
            count++;
        }
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        count++;

        //  if( (wait_gui!=NULL) && (count > (max_line_count - 10)) )
        //    {

        //      wait_gui("Press Continue");
        //      count = 0;
        //    }

        /* Print HaloAcetic Acids (continued) */
        e |= fprintf(fout, "\r\n");
        e |= fprintf(fout, "                                     Table 7                                 \r\n");
        e |= fprintf(fout, "                  Predicted Haloacetic Acids (HAA6 through HAA9)             \r\n");
        e |= fprintf(fout, "            At Average Flow (%4.1f MGD) and Influent Temperature (%3.1f C)   \r\n", influent->Flow, influent->DegK - 273.15);
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        e |= fprintf(fout, "                  BCAA  BDCAA  DBCAA   TBAA   HAA6   HAA9                    \r\n");
        e |= fprintf(fout, "Location         (ug/L) (ug/L) (ug/L) (ug/L) (ug/L) (ug/L)                   \r\n");
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        count += 8;

        for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
        {
            eff = &unit->eff;
            e |= fprintf(fout, "%-16s", UnitProcessTable[unit->type].name);
            e |= fprintf(fout, "%6.0f", eff->BCAA);
            e |= fprintf(fout, "%7.0f", eff->BDCAA);
            e |= fprintf(fout, "%7.0f", eff->DBCAA);
            e |= fprintf(fout, "%7.0f", eff->TBAA);
            e |= fprintf(fout, "%7.0f", eff->HAA6);
            e |= fprintf(fout, "%7.0f", eff->HAA9);
            e |= fprintf(fout, "\r\n");
            count++;
        }
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        count++;

        //  if( (wait_gui!=NULL) && (count > (max_line_count - 10)) )
        //    {

        //      wait_gui("Press Continue");
        //      count = 0;
        //    }

        /* Print Disinfection-related parameters */
        e |= fprintf(fout, "\r\n");
        e |= fprintf(fout, "                                   Table 8                                   \n");
        e |= fprintf(fout, "         Predicted Disinfection Parameters - Residuals and CT Ratios         \n");
        e |= fprintf(fout, "         At Plant Flow (%4.1f MGD) and Influent Temperature (%3.1f C)        \r\n", influent->Flow, influent->DegK - 273.15);
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        e |= fprintf(fout, "                                                               CT Ratios     \r\n");
        e |= fprintf(fout, "                 Temp  pH    Cl2   NH2Cl  Ozone   ClO2  -------------------- \r\n");
        e |= fprintf(fout, "Location          (C)  (-)  (mg/L) (mg/L) (mg/L) (mg/L) Giardia Virus Crypto \r\n");
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        count += 9;

        for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
        {
            eff = &unit->eff;
            e |= fprintf(fout, "%-16s", UnitProcessTable[unit->type].name);
            e |= fprintf(fout, "%5.1f", eff->DegK - 273.15);
            e |= fprintf(fout, "%5.1f", eff->pH);
            e |= fprintf(fout, "%6.1f", eff->FreeCl2 * MW_Cl2);
            e |= fprintf(fout, "%7.1f", eff->NH2Cl * MW_Cl2);
            e |= fprintf(fout, "%7.2f", eff->o3_res);
            e |= fprintf(fout, "%7.2f", eff->clo2_res);
            /* e |= fprintf(fout,"%7.2f", eff->processtime        );  */
            if (influent->log_required_g == 999999.0 && gw_virus_flag == FALSE) //groundwater w/o virus dis. req'd.
            {
                e |= fprintf(fout, "    na");
                e |= fprintf(fout, "     na");
                e |= fprintf(fout, "    na");
            }
            else if (influent->log_required_g == 999999.0 && gw_virus_flag == TRUE) //groundwater w/virus dis. req'd.
            {
                e |= fprintf(fout, "    na");
                e |= fprintf(fout, "%8.1f", eff->ct_ratio_v);
                e |= fprintf(fout, "    na");
            }
            else //surface water
            {
                e |= fprintf(fout, "%7.1f", eff->ct_ratio);
                e |= fprintf(fout, "%8.1f", eff->ct_ratio_v);
                e |= fprintf(fout, "%7.1f", eff->ct_ratio_c);
            }
            e |= fprintf(fout, "\r\n");
            count++;
        }
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        count += 1;
        if (influent->log_required_g == 999999.0)
        {
            e |= fprintf(fout, "na = not applicable - no CT requirements exist for groundwater\r\n");
            count += 1;
        }

        /*This code describes if CT_Ratio = 1 due to physical removal*/
        if (influent->ct_ratio_v == 1.0)
        {
            e |= fprintf(fout, "Virus CT Ratio = 1 because other credits met full disinfection requirements\r\n");
            count += 1;
        }
        if (influent->ct_ratio == 1.0)
        {
            e |= fprintf(fout, "Giardia CT Ratio = 1 because other credits met full disinfection requirements\r\n");
            count += 1;
        }
        if (influent->ct_ratio_c == 1.0)
        {
            e |= fprintf(fout, "Crypto. CT Ratio = 1 because other credits met full disinfection requirements\r\n");
            count += 1;
        }
        if (influent->ct_ratio == 1.0 || influent->ct_ratio_v == 1.0 || influent->ct_ratio_c == 1.0)
        {
            e |= fprintf(fout, "\r\n");
            count += 1;
        }

        //  if( (wait_gui!=NULL) && (count > (max_line_count - 10)) )
        //    {

        //      wait_gui("Press Continue");
        //      count = 0;
        //    }

        /*********************************************************************************************************************************************/
        /* Print Disinfection-related parameters - part (b) */
        e |= fprintf(fout, "\r\n");
        e |= fprintf(fout, "                                    Table 9                                  \n");
        e |= fprintf(fout, "                Predicted Disinfection Parameters - CT Values                \n");
        e |= fprintf(fout, "         At Plant Flow (%4.1f MGD) and Influent Temperature (%3.1f C)        \r\n", influent->Flow, influent->DegK - 273.15);
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        e |= fprintf(fout, "                                                                             \r\n");
        e |= fprintf(fout, "                     Cl2     NH2Cl   Ozone   ClO2                            \r\n");
        e |= fprintf(fout, "Location             <-----(mg/L * minutes)----->                             \r\n");
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        count += 9;

        for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
        {
            eff = &unit->eff;
            e |= fprintf(fout, "%-16s", UnitProcessTable[unit->type].name);
            e |= fprintf(fout, "%8.1f", eff->ct_cl2);
            e |= fprintf(fout, "%9.1f", eff->ct_nh2cl);
            e |= fprintf(fout, "%8.1f", eff->ct_o3);
            e |= fprintf(fout, "%8.1f", eff->ct_clo2);
            e |= fprintf(fout, "\r\n");
            count++;
        }
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        count += 1;

        // if( (wait_gui!=NULL) && (count > (max_line_count - 10)) )
        //    {

        //      wait_gui("Press Continue");
        //      count = 0;
        //    }

        /*********************************************************************************************************************************************/

        /*
  e |= fprintf(fout,"\r\n");
  e |= fprintf(fout,"                                 Debugging Table                            \r\n");
  e |= fprintf(fout,"      Model Used to Predict DBPs in Each Unit Process - Avg. Conditions     \r\n");
  e |= fprintf(fout,"----------------------------------------------------------------------------\r\n");
  e |= fprintf(fout,"                 DBP Model                                                  \r\n");
  e |= fprintf(fout,"Location           Used                                                     \r\n");
  e |= fprintf(fout,"----------------------------------------------------------------------------\r\n");
  count += 7;

  for( unit=FirstUnitProcess(train); unit; unit=NextUnitProcess(unit) )
    {
      eff = &unit->eff;
      e |= fprintf( fout,"%-16s"  , UnitProcessTable[unit->type].name );
      e |= fprintf( fout,"%-50s"  , eff->dbpmodel                     );
      e |= fprintf( fout, "\r\n");
      count++;
    }

 if( (wait_gui!=NULL) && (count > (max_line_count - 10)) )
   {
    wait_gui("Press Continue");
    count = 0;
   }
  */
        /*
  e |= fprintf(fout,"\r\n");
  e |= fprintf(fout,"                                 Debugging Table                            \r\n");
  e |= fprintf(fout,"        Output of Process Train Description Flags - Avg. Conditions         \r\n");
  e |= fprintf(fout,"----------------------------------------------------------------------------\r\n");
  e |= fprintf(fout,"                                                                            \r\n");
  e |= fprintf(fout,"Process Description           YES/NO                                        \r\n");
  e |= fprintf(fout,"----------------------------------------------------------------------------\r\n");
  count += 7;

      e |= fprintf( fout,"Pre-Sedimentation");
      if   (pre_sedflag == TRUE) e |= fprintf(fout,"               YES\n");
      else                       e |= fprintf(fout,"                NO\n");

      e |= fprintf( fout,"Coagulation");
      if   (coagflag == TRUE) e |= fprintf(fout,"                     YES\n");
      else                    e |= fprintf(fout,"                      NO\n");

      e |= fprintf( fout,"Sedimentation");
      if   (sedflag == TRUE) e |= fprintf(fout,"                   YES\n");
      else                   e |= fprintf(fout,"                    NO\n");

      e |= fprintf( fout,"Filtration");
      if   (filtflag == TRUE) e |= fprintf(fout,"                      YES\n");
      else                    e |= fprintf(fout,"                       NO\n");

      e |= fprintf( fout,"  Direct Filtration");
      if   (dir_filtflag == TRUE) e |= fprintf(fout,"             YES\n");
      else                        e |= fprintf(fout,"              NO\n");

      e |= fprintf( fout,"  Biofiltration");
      if   (bio_filtflag == TRUE) e |= fprintf(fout,"                 YES\n");
      else                        e |= fprintf(fout,"                  NO\n");

      e |= fprintf( fout,"Ozonation");
      if   (o3flag == TRUE)   e |= fprintf(fout,"                       YES\n");
      else                    e |= fprintf(fout,"                        NO\n");

      e |= fprintf( fout,"  Pre-Ozonation");
      if   (pre_o3flag == TRUE) e |= fprintf(fout,"                 YES\n");
      else                      e |= fprintf(fout,"                  NO\n");

      e |= fprintf( fout,"  Intermediate-Ozonation");
      if   (int_o3flag == TRUE)     e |= fprintf(fout,"        YES\n");
      else                          e |= fprintf(fout,"         NO\n");

      e |= fprintf( fout,"  Post-Ozonation");
      if   (post_o3flag == TRUE) e |= fprintf(fout,"                YES\n");
      else                       e |= fprintf(fout,"                 NO\n");

      e |= fprintf( fout,"MF/UF");
      if   (mfufflag == TRUE) e |= fprintf(fout,"          YES\n");
      else                   e |= fprintf(fout,"           NO\n");

      e |= fprintf( fout,"NF");
      if   (nfflag == TRUE)  e |= fprintf(fout,"          YES\n");
      else                   e |= fprintf(fout,"           NO\n");

      e |= fprintf( fout,"GAC");
      if   (gacflag == TRUE) e |= fprintf(fout,"                             YES\n");
      else                   e |= fprintf(fout,"                              NO\n");

      e |= fprintf( fout,"  Biofiltration with GAC");
      if   (bio_gacflag == TRUE) e |= fprintf(fout,"        YES\n");
      else                       e |= fprintf(fout,"         NO\n");

      e |= fprintf( fout,"Chlorine Dioxide Addition");
      if   (clo2flag == TRUE)    e |= fprintf(fout,"       YES\n");
      else                       e |= fprintf(fout,"        NO\n");

      e |= fprintf( fout, "\r\n");

*/

        /* RUN MODEL AGAIN and Print Inactivation at Min Temp & Max Flow */
        coldflag = TRUE;
        if (runmodel(train) == FALSE)
            return (FALSE);
        coldflag = FALSE;

        if (swflag == TRUE)
            strcpy(buffer, "for Surface Water Plant ");
        else
            strcpy(buffer, "for Groundwater Plant ");
        if (coagflag)
        {
            if (filtflag)
                strcat(buffer, "with Coagulation and Filtration");
            else
                strcat(buffer, "with Coagulation");
        }
        else
        {
            if (filtflag)
                strcat(buffer, "with Filtration");
            else
                strcat(buffer, ""); /* without coagulation or filtration. */
            ;
        }

        e |= fprintf(fout, "\r\n");
        e |= fprintf(fout, "                                  Table 10                                   \r\n");
        e |= fprintf(fout, "                      Predicted Disinfection Parameters                      \r\n");
        e |= fprintf(fout, "          At Peak Flow (%4.1f MGD) and Minimum Temperature (%3.1f C)         \r\n", influent->Flow, influent->DegK - 273.15);
        for (size_t i = 0; i < (76 - strlen(buffer)) / 2; i++)
            e |= fprintf(fout, " ");
        e |= fprintf(fout, "%s\r\n", buffer);
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        e |= fprintf(fout, "                                                              CT Ratios      \r\n");
        e |= fprintf(fout, "                 Temp  pH    Cl2   NH2Cl  Ozone   ClO2  ---------------------\r\n");
        e |= fprintf(fout, "Location          (C)  (-)  (mg/L) (mg/L) (mg/L) (mg/L) Giardia Virus Crypto \r\n");
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        count += 10;

        for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
        {
            eff = &unit->eff;
            e |= fprintf(fout, "%-16s", UnitProcessTable[unit->type].name);
            e |= fprintf(fout, "%5.1f", eff->DegK - 273.15);
            e |= fprintf(fout, "%5.1f", eff->pH);
            e |= fprintf(fout, "%6.1f", eff->FreeCl2 * MW_Cl2);
            e |= fprintf(fout, "%7.1f", eff->NH2Cl * MW_Cl2);
            e |= fprintf(fout, "%7.2f", eff->o3_res);
            e |= fprintf(fout, "%7.2f", eff->clo2_res);
            /* e |= fprintf(fout,"%7.2f", eff->processtime        );    */
            if (influent->log_required_g == 999999.0 && gw_virus_flag == FALSE) //groundwater w/virus dis. req'd.
            {
                e |= fprintf(fout, "    na");
                e |= fprintf(fout, "     na");
                e |= fprintf(fout, "    na");
            }
            else if (influent->log_required_g == 999999.0 && gw_virus_flag == TRUE) //groundwater w/virus dis. req'd.
            {
                e |= fprintf(fout, "    na");
                e |= fprintf(fout, "%8.1f", eff->ct_ratio_v);
                e |= fprintf(fout, "    na");
            }
            else //surface water
            {
                e |= fprintf(fout, "%7.1f", eff->ct_ratio);
                e |= fprintf(fout, "%8.1f", eff->ct_ratio_v);
                e |= fprintf(fout, "%7.1f", eff->ct_ratio_c);
            }
            e |= fprintf(fout, "\r\n");
            count++;
        }
        e |= fprintf(fout, "-----------------------------------------------------------------------------\r\n");
        count += 1;

        if (influent->log_required_g == 999999.0)
        {
            e |= fprintf(fout, "na = not applicable - no CT requirements for groundwater\r\n");
            count += 1;
        }

        /*This code describes if CT_Ratio = 1 due to physical removal*/
        if (influent->ct_ratio_v == 1.0)
        {
            e |= fprintf(fout, "Virus CT Ratio = 1 because other credits met full disinfection requirements\r\n");
            count += 1;
        }
        if (influent->ct_ratio == 1.0)
        {
            e |= fprintf(fout, "Giardia CT Ratio = 1 because other credits met full disinfection requirements\r\n");
            count += 1;
        }
        if (influent->ct_ratio_c == 1.0)
        {
            e |= fprintf(fout, "Crypto. CT Ratio = 1 because other credits met full disinfection requirements\r\n");
            count += 1;
        }
        if (influent->ct_ratio == 1.0 || influent->ct_ratio_v == 1.0 || influent->ct_ratio_c == 1.0)
        {
            e |= fprintf(fout, "\r\n");
            count += 1;
        }

        //  if( (wait_gui!=NULL) && (count > (max_line_count - 10)) )
        //    {

        //      wait_gui("Press Continue");
        //      count = 0;
        //    }

        /*
  e |= fprintf(fout,"\r\n");
  e |= fprintf(fout,"                                 Utility Table                              \r\n");
  e |= fprintf(fout," Predicted Parameters at Minimum Temperature and Peak Hour Flow             \r\n");
  for( i=0; i< (76-strlen(buffer))/2; i++ ) e |= fprintf(fout," ");
  e |= fprintf(fout,"%s\r\n", buffer);
  e |= fprintf(fout,"----------------------------------------------------------------------------\r\n");
  e |= fprintf(fout,"                 Temp  pH    Alk    UV254                                   \r\n");
  e |= fprintf(fout,"Location          (C)  (-)  (mg/L)  (1/cm)                                  \r\n");
  e |= fprintf(fout,"----------------------------------------------------------------------------\r\n");
  count += 7;

  for( unit=FirstUnitProcess(train); unit; unit=NextUnitProcess(unit) )
    {
      eff = &unit->eff;
      e |= fprintf( fout,"%-16s"  , UnitProcessTable[unit->type].name );
      e |= fprintf(fout,"%5.1f", eff->DegK - 273.15      );
      e |= fprintf(fout,"%5.1f", eff->pH                 );
      e |= fprintf(fout,"%7.0f", eff->Alk * MW_CaCO3/2   );
      e |= fprintf(fout,"%9.3f", eff->UV                );
      e |= fprintf( fout, "\r\n");
      count++;
    } */
        /*  e |= fprintf(fout,"\n\n\nNote:  If Avg. Water Temp < 10 C then the model will adjust the \r\n");
      e |= fprintf(fout,"       Avg. Water Temp to 10 C for Trihalomethane Formation.  The\n");
      e |= fprintf(fout,"       input temperature is used for all other appropriate equations.\n");   */
        /*
  if( (wait_gui!=NULL) && (count > (max_line_count - 10)) )
    {
      wait_gui("Press Continue");
      count = 0;
    }

*/
    }
    //  if( wait_gui!=NULL)
    //	     wait_gui("Listing Complete");

    if (e < 0)
        return (FALSE);
    else
        return (TRUE);
}

/* Multi-simulation mode */
// void multi_sim_mode(ProcessTrain *train, MultiSimParameters params);

/* Multi-simulation mode (with automatic chemical dosing) */
// void multi_sim_mode(ProcessTrain *train, MultiSimAutoDose params);

/* Simulation-optimization mode */
void sim_opt_mode(ProcessTrain *train)
{
    /* Initialize simulation-optimization variables */
    // int n_timesteps = params.ts.n_years * params.ts.timesteps_per_year; // calcualte number of time steps
    // double ***samples;                                                  // 3D array of Monte Carlo time series data
    std::vector<std::vector<double>> mc_samples; // 2D vector of Monte Carlo samples
    // int n_simopt_params = 10;                                           /* Number of different water quality parameters generated by Monte Carlo sampling (listed below) */
    /* 1) Alkalinity (mg/L as CaCO3) */
    /* 2) Ammonia (mg/L as N) */
    /* 3) Bromide (ug/L) */
    /* 4) Calcium hardness (mg/L as CaCO3) */
    /* 5) Hardness (total) (mg/L as CaCO3) */
    /* 6) pH (-) */
    /* 7) Temperature (deg C) */
    /* 8) Total organic carbon (mg/L) */
    /* 9) Turbidity (NTU) */
    /* 10) UV absorbance at 254 nm (1/cm) */

    /* Initialize Borg functions and variables*/
    BORG_Archive result;
    BORG_Problem problem;

    const char *borg_header = "alk_Q1 pH_Q1 DBPsf_Q1 alk_Q2 pH_Q2 DBPsf_Q2 alk_Q3 pH_Q3 DBPsf_Q3 alk_Q4 pH_Q4 DBPsf_Q4 wc_freq_TTHM wc_freq_HAA5 expect_solids expect_lime expect_co2 lraa_const toc_const ct_const \n"; /* Header for Borg results file */

    /* Initialize problem formulation variables */
    // int nvars;   // number of decision variables
    // int nobjs;   // number of objectives
    // int nconsts; // number of constraints
    // int n_var_types = 3; // number of variable types (this multiplied by the number of timesteps/year to get total number of variables)

    double *vars;
    double *objs;
    double *consts;

    /* Read in Monte Carlo samples */
    // const char *file_montecarlo = "./in/monte_carlo/influent-wq-data.csv";
    // FILE *fin = fopen(file_montecarlo, "r"); // open input file
    // if (!fin)
    // { /* File Error handling */
    //     printf("Error opening %s \n, check that the proper directory and file name has been specified.\n", file_montecarlo);
    //     exit(EXIT_FAILURE);
    // }

    /* Allocate memory to a 3D array of doubles */
    /* source: http://www.geeksforgeeks.org/dynamically-allocate-2d-array-c/ (Using pointer to a pointer) */
    // samples = (double ***)malloc(n_timesteps * sizeof(double **));
    // for (int i = 0; i < n_timesteps; i++)
    // {
    //     samples[i] = (double **)malloc(n_simopt_params * sizeof(double *));
    //     for (int j = 0; j < n_simopt_params; j++)
    //     {
    //         samples[i][j] = (double *)malloc(params.mc.n_samples * sizeof(double));
    //     }
    // }
    // read_montecarlo(samples, params.mc.n_samples, n_timesteps, n_simopt_params, fin);

    // /* Select problem formulation */
    // if (params.problem == 1)
    // {
    //     nvars = params.ts.timesteps_per_year * n_var_types;
    //     // nvars = n_var_types;
    //     // nconsts = 1;
    //     // nconsts = 3;
    //     // nobjs = 6;

    //     vars = (double *)malloc(nvars * sizeof(double));
    //     objs = (double *)malloc(nobjs * sizeof(double));
    //     consts = (double *)malloc(nconsts * sizeof(double));
    // }

    /* Initialize problem formulation */
    // printf("number of variables: %d\n", N_VARS);
    // printf("number of objectives: %d\n", N_OBJS);
    // printf("number of constraints: %d\n", N_CONSTS);
    vars = (double *)malloc(N_VARS * sizeof(double));
    objs = (double *)malloc(N_OBJS * sizeof(double));
    consts = (double *)malloc(N_CONSTS * sizeof(double));

    /* Create the WTP problem, defining the number of decision variables, objectives and constraints.
                                 * The last argument, wtp, references the function that evaluates the WTP problem. */
    problem = BORG_Problem_create(N_VARS, N_OBJS, N_CONSTS, wtp);

    /* Set the lower and upper bounds for each decision variable. */
    // WJR: include problem logic for other problems and use timesteps_per_year to set how many quarters there are
    for (int i = 0; i < N_VAR_TYPES; i++)
    {
        for (int j = 0; j < N_QUARTERS_PER_YEAR; j++)
        {
            if (i == 0)
                BORG_Problem_set_bounds(problem, i + 3 * j, 0.0, 150.0); // set bounds of alkalinity decision variable
            if (i == 1)
                BORG_Problem_set_bounds(problem, i + 3 * j, 6.0, 9.5); // set bounds of pH decision variable
            if (i == 2)
                BORG_Problem_set_bounds(problem, i + 3 * j, 0.02, 0.75); // set bounds of disinfection byproduct safety factor
        }
    }

    /* Set the epsilon values used by the Borg MOEA.  Epsilons define the problem resolution, which controls
                                 * how many Pareto optimal solutions are generated and how far apart they are spaced. */
    for (int i = 0; i < N_OBJS; i++)
    {
        // WJR: include problem logic for other problems
        BORG_Problem_set_epsilon(problem, i, 0.01); // WJR: epsilon should change for each objective
    }

    /* Run the Borg MOEA on the WTP problem for <input #> function evaluations.*/
    // int nfe = params.num_func_evals;
    result = BORG_Algorithm_run(problem, simopt_params.num_func_evals);

    // /* Record current date and time */
    // time_t rawtime;
    // struct tm *info;
    // char buffer[256];

    // time(&rawtime);
    // info = localtime(&rawtime);
    // strftime(buffer, sizeof(buffer), "%m-%d-%y_%H-%M-%S", info);
    // std::string current_datetime(buffer); /* String which contains the date and time at the start of running the program */

    /* Print the Pareto optimal solutions.*/
    std::string result_dir; // result directory
    std::string result_ext; // result extension
    result_dir = "./out/sim_opt/result/";
    result_ext = ".result";
    std::string filename_result = result_dir + current_datetime + result_ext;

    FILE *fres = fopen(filename_result.c_str(), "w"); // open file stream

    if (!fres)
    { /* Error handling */
        printf("Error opening %s \n, check that the proper directory and file name has been specified.\n", filename_result.c_str());
        exit(EXIT_FAILURE);
    }

    fprintf(fres, borg_header); /* Print problem formulation header */

    BORG_Archive_print(result, fres);

    /* Free any allocated memory. */
    BORG_Archive_destroy(result);
    BORG_Problem_destroy(problem);
    free(vars);
    free(objs);
    free(consts);

    std::cout << "Complete! Check output file " << filename_result << " for results." << std::endl;

    /* Close file stream(s) */
    // fclose(fin);
    fclose(fres);
}

/* Validation mode */
// void validation_mode(ProcessTrain *train, ValidationParameters params);
