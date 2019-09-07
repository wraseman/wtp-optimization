/* auto_dose.c */

/* Purpose: automatically adjust chemical doses based on water quality set points.  */

#include "wtp.h"
#include "wtp_optimize.h"
#include "auto_dose.h"

void auto_dose(struct ProcessTrain *train,
               double pH_setpt_1, double alk_setpt_1, double cl2_setpt, double pH_setpt_2, double DBP_safety_factor, 
               FILE *fres)
{
    /* Purpose: automatically adjust chemical doses based on user defined set points. 
     * 
     * Inputs: 
     *  train       = Process train control structure.
     *  alk_setpt_1 = Alkalinity set point for raw water quality adjustment (aka location #1).
     *  pH_setpt_1  = pH set point for raw water quality adjustment (aka location #1).
     *  cl2_setpt   = Chlorine residual set point for the end of distribution system. 
     *  pH_setpt_2  = pH set point for effluent corrosion control (aka location #2).
     *  DBP_safety_factor = safety factor (as a fraction) of DBP group compared to maximum contaminant level.
     * 
     * Outputs: 
     *  None
     * */

    register struct UnitProcess *unit;

    struct Effluent *influent;
    struct Effluent *effluent;
    struct Effluent *eos;
    struct Effluent *co2_1; // effluent of CO2 addition, location #1
    struct Effluent *rapid_mix;

    int debug = TRUE; // (TRUE/FALSE): TRUE to print debugging statements to output file, FALSE for no output

    /* Put date/timestamp on treatment train/root-finding file */
    std::string log_dir; // log directory
    std::string log_ext; // log extension
    log_dir = "./out/single_sim/log/";
    log_ext = ".log";
    std::string log_filepath = log_dir + current_datetime + log_ext;

    /* Open file stream to log file */
    FILE *flog = fopen(log_filepath.c_str(), "w"); // print out the treatment train and rooting finding progress

    if (!flog)
    { /* Error handling */
        printf("Error opening %s \n, check that the proper directory and file name has been specified.\n", log_filepath.c_str());
        exit(EXIT_FAILURE);
    }

    /* Print out treatment train/root-finding file description */
    if (debug == TRUE)
    {
        fprintf(flog, "############################### TREATMENT TRAIN AND ROOT FINDING FILE ###############################\n");
        fprintf(flog, "# The following file has information about the treatment train parameters and reports the progress of \n");
        fprintf(flog, "# the root-finding algorithms (e.g. Secant and Bisection Methods) during the automatic chemical dose \n");
        fprintf(flog, "# routine, auto_dose(). This file is overwritten each time auto_dose() is called, so only the  \n");
        fprintf(flog, "# information about the last call of auto_dose() will be shown. \n");
        fprintf(flog, "# This file is particularly useful for checking the intial and final configurations of the treatment \n");
        fprintf(flog, "# train and for debugging the root-finding algorithms. \n");
        // fprintf(flog, "# The name train&rootfinding_MM-DD-YY_HH-MM-SS.out is used to uniquely identify each run by the date and time.\n");
        fprintf(flog, "\n");
    }

    RF_FUNC_PTR mod_dose_check_target_ptr = mod_dose_check_target; //set function pointer to pass mod_dose_check_target() function into root_find()

    int iter = 0;
    int max_iter = 50;
    static int run_number = 0; // number of model runs (maximum run_number = (number of function evaluations * Monte Carlo evaluations [or 1 if not in MC mode]))

    double lime_lo = 0.0;    // set lower bound for lime dosing search
    double lime_up = 1000.0; // set upper bound for lime dosing search
    double co2_lo = 0.0;     // set lower bound for CO2 dosing search
    double co2_up = 1000.0;  // set upper bound for CO2 dosing search
    double naocl_lo = 0.0;   // set lower bound for CO2 dosing search
    double naocl_up = 300.0; // set upper bound for CO2 dosing search
    double alum_lo = 18.0;   // set lower bound for CO2 dosing search
    double alum_up = 1000.0; // set upper bound for CO2 dosing search
    int co2_cntr = 0;        /* counts number of CARBON_DIOXIDE points in train, WJR */

    double lime_dose_pH;  // lime dose needed to meet pH setpoint
    double lime_dose_alk; // lime does needed to meet alkalinity setpoint

    int setptflag_pH_1 = FALSE;
    int setptflag_alk_1 = FALSE;
    int setptflag_DBP = FALSE;
    int setptflag_cl2 = FALSE;
    int setptflag_pH_2 = FALSE;

    double TTHM_MCL = 80.0; // maximum contaminant level for TTHM
    double HAA5_MCL = 60.0; // maximum contaminant level for HAA5
    double toc_rem_req;

    double TTHM_target = TTHM_MCL * (1 - DBP_safety_factor); // maximum contaminat level * a safety factor
    double HAA5_target = HAA5_MCL * (1 - DBP_safety_factor); // maximum contaminat level * a safety factor
    double toc_rem_target;
    double toc_safety_factor = 0.05;
    double pH_setpt_alum = 5.5;

    /* Initialize inputs to root-finding function */
    int dose_unit, dose_location, target_unit, target_location;
    char target_param;
    double target, x_lo, x_up;

    /* Run WTP model and define water quality sampling locations (influent, effluent, end of system, etc.) */
    if (debug == TRUE)
    {
        fprintf(flog, "============ START of auto_dose() ============\n");
        fprintf(flog, "model run number: %d \n", run_number);
        fprintf(flog, "alk_setpt_1: %f \n", alk_setpt_1);
        fprintf(flog, "pH_setpt_1: %f \n", pH_setpt_1);
        fprintf(flog, "End of system Cl2 setpoint: %f \n", cl2_setpt);
        fprintf(flog, "pH_setpt_2: %f \n", pH_setpt_2);
        fprintf(flog, "DBP_safety_factor: %f \n", DBP_safety_factor);
        fprintf(flog, "toc_safety_factor: %f \n", toc_safety_factor);
        fprintf(flog, "\n=========== Process Train =========== \n");
        writewtp(flog, train, flog); // write out process train to file
        fprintf(flog, "\n");
    }

    /* Run model */
    runmodel(train);

    /* Get unit process effluent values */
    for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
    {
        switch (unit->type)
        {
        case INFLUENT:
            influent = &unit->eff;
            break;

        case CARBON_DIOXIDE:
            if (co2_cntr == 0)
            {
                co2_1 = &unit->eff;
            }
            co2_cntr++; // keep track of number of CO2 additions in case there are multiple
            break;

        case RAPID_MIX:
            rapid_mix = &unit->eff;
            break;

        case WTP_EFFLUENT:
            effluent = &unit->eff;
            break;

        case END_OF_SYSTEM:
            eos = &unit->eff;
            break;

        default:
            break;
        }
    }

    while (setptflag_pH_1 == FALSE ||
           setptflag_alk_1 == FALSE ||
           setptflag_DBP == FALSE ||
           setptflag_cl2 == FALSE ||
           setptflag_pH_2 == FALSE)
    {

        // Set all setpoint flags to false at the beginning of each loop
        setptflag_pH_1 = FALSE;
        setptflag_alk_1 = FALSE;
        setptflag_DBP = FALSE;
        setptflag_cl2 = FALSE;
        setptflag_pH_2 = FALSE;

        // Automate raw water pH and alkalinity adjustment
        /*===================== START RAW WATER pH AND ALKALINITY ADJUSTMENT ==========================*/

        /* If alkalinity is too low, add lime */
        if ((co2_1->Alk * MW_CaCO3 / 2.0) < alk_setpt_1)
        {
            // adjust alkalinity by adding lime
            dose_unit = LIME;
            dose_location = 0;
            target_param = 'A';
            target_unit = CARBON_DIOXIDE;
            target_location = 0;
            target = alk_setpt_1;
            x_lo = lime_lo;
            x_up = lime_up;

            if (debug == TRUE)
            {
                fprintf(flog, "auto_dose: add lime to meet raw water alkalinity setpoint. Entering rootfind_and_mod_dose()! \n\n");
                write_dose_target(dose_unit, dose_location, target_param, target_unit, target_location, target, flog);
            }

            lime_dose_alk = rootfind_and_mod_dose(mod_dose_check_target_ptr, train, dose_unit, dose_location, target_param, target_unit, target_location, target, x_lo, x_up, flog);
            //bisect_method(mod_dose_check_target_ptr, train, LIME, 0, 'A', CARBON_DIOXIDE, 0, alk_setpt_1, lime_lo, lime_up, flog);
        }

        /* If pH is too low, add lime; if pH is too high, add carbon dioxide */
        if (co2_1->pH < pH_setpt_1)
        { // adjust pH by adding lime

            dose_unit = LIME;
            dose_location = 0;
            target_param = 'P';
            target_unit = CARBON_DIOXIDE;
            target_location = 0;
            target = pH_setpt_1;
            x_lo = lime_lo;
            x_up = lime_up;

            if (debug == TRUE)
            {
                fprintf(flog, "auto_dose: add lime to meet raw water pH setpoint. Entering rootfind_and_mod_dose()! \n\n");
                write_dose_target(dose_unit, dose_location, target_param, target_unit, target_location, target, flog);
            }
            lime_dose_pH = rootfind_and_mod_dose(mod_dose_check_target_ptr, train, dose_unit, dose_location, target_param, target_unit, target_location, target, x_lo, x_up, flog);

            /* Select whichever lime dose is higher (to acheive either pH or alkalinity setpoint) */
            if (lime_dose_alk > lime_dose_pH)
            {
                // adjust alkalinity by adding lime
                dose_unit = LIME;
                dose_location = 0;
                target_param = 'A';
                target_unit = CARBON_DIOXIDE;
                target_location = 0;
                target = alk_setpt_1;
                x_lo = lime_lo;
                x_up = lime_up;

                if (debug == TRUE)
                {
                    fprintf(flog, "auto_dose: add lime to meet raw water alkalinity setpoint. Entering rootfind_and_mod_dose()! \n\n");
                    write_dose_target(dose_unit, dose_location, target_param, target_unit, target_location, target, flog);
                }

                lime_dose_alk = rootfind_and_mod_dose(mod_dose_check_target_ptr, train, dose_unit, dose_location, target_param, target_unit, target_location, target, x_lo, x_up, flog);
            }
        }
        else
        { // adjust pH by adding carbon dioxide

            dose_unit = CARBON_DIOXIDE;
            dose_location = 0;
            target_param = 'P';
            target_unit = CARBON_DIOXIDE;
            target_location = 0;
            target = pH_setpt_1;
            x_lo = co2_lo;
            x_up = co2_up;

            if (debug == TRUE)
            {
                fprintf(flog, "auto_dose: add CO2 to meet raw water pH setpoint. Entering rootfind_and_mod_dose()! \n\n");
                write_dose_target(dose_unit, dose_location, target_param, target_unit, target_location, target, flog);
            }

            rootfind_and_mod_dose(mod_dose_check_target_ptr, train, dose_unit, dose_location, target_param, target_unit, target_location, target, x_lo, x_up, flog);
        }

        /*======================== END RAW WATER pH AND ALKALINITY ADJUSTMENT ==========================*/

        /*===================== START ALUM ADJUSTMENT TO ACHIEVE DBP AND TOC REGS ==========================*/
        /* For first iteration, set alum to minimum dose (18 mg/L) to control turbidity */
        if (debug == TRUE)
            fprintf(flog, "auto_dose: set alum to minimum dose (18 mg/L). \n");

        for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
        {
            switch (unit->type)
            {
            case ALUM:
                unit->data.alum->dose = 18.0;
                break;

            default:
                break;
            }
        }

        /* Adjust alum dose to achieve TOC regulations */
        if ((influent->Alk * MW_CaCO3 / 2.0) < 60.0)
        {
            if ((influent->TOC < 2.0))
            {
                // Do nothing if TOC < 2.0
                if (debug == TRUE)
                    fprintf(flog, "auto_dose: influent TOC =< 2.0, do not change alum dose. \n");
            }
            else
            {
                if ((influent->TOC < 4.0))
                {
                    toc_rem_req = 35.0;
                    toc_rem_target = toc_rem_req * (1 + toc_safety_factor);
                    if (debug == TRUE)
                        fprintf(flog, "auto_dose: influent TOC < 4.0, TOC removal target set to 35%%. \n");
                }
                else
                {
                    if ((influent->TOC < 8.0))
                    {
                        toc_rem_req = 45.0;
                        toc_rem_target = toc_rem_req * (1 + toc_safety_factor);
                        if (debug == TRUE)
                            fprintf(flog, "auto_dose: influent 4.0 <= TOC < 8.0, TOC removal target set to 45%%. \n");
                    }
                    else
                    {
                        toc_rem_req = 45.0;
                        toc_rem_target = toc_rem_req * (1 + toc_safety_factor);
                        if (debug == TRUE)
                            fprintf(flog, "auto_dose: influent TOC >= 8.0, TOC removal target set to 50%%. \n");
                    }
                }
                if (((influent->TOC - effluent->TOC) / influent->TOC * 100.0) < toc_rem_target)
                {

                    dose_unit = ALUM;
                    dose_location = 0;
                    target_param = 'O';
                    target_unit = END_OF_SYSTEM;
                    target_location = 0;
                    target = toc_rem_target;
                    x_lo = alum_lo;
                    x_up = alum_up;

                    if (debug == TRUE)
                    {
                        fprintf(flog, "auto_dose: add alum to meet TOC target. Entering rootfind_and_mod_dose()! \n\n");
                        write_dose_target(dose_unit, dose_location, target_param, target_unit, target_location, target, flog);
                    }

                    rootfind_and_mod_dose(mod_dose_check_target_ptr, train, dose_unit, dose_location, target_param, target_unit, target_location, target, x_lo, x_up, flog);
                }
            }
        }
        else
        {
            if (debug == TRUE)
                fprintf(stderr, "auto_dose: have not yet included logic for raw water alkalinity >= 60 mg/L! \n");
            exit(EXIT_FAILURE);
        }

        /* Adjust alum dose to achieve DBP regulations */
        if (eos->TTHM > TTHM_target)
        {

            dose_unit = ALUM;
            dose_location = 0;
            target_param = 'T';
            target_unit = END_OF_SYSTEM;
            target_location = 0;
            target = TTHM_target;
            x_lo = alum_lo;
            x_up = alum_up;

            if (debug == TRUE)
            {
                fprintf(flog, "auto_dose: add alum to meet TTHM. Entering rootfind_and_mod_dose()! \n\n");
                write_dose_target(dose_unit, dose_location, target_param, target_unit, target_location, target, flog);
            }
            rootfind_and_mod_dose(mod_dose_check_target_ptr, train, dose_unit, dose_location, target_param, target_unit, target_location, target, x_lo, x_up, flog);
        }

        if (eos->HAA5 > HAA5_target)
        {

            dose_unit = ALUM;
            dose_location = 0;
            target_param = 'H';
            target_unit = END_OF_SYSTEM;
            target_location = 0;
            target = HAA5_target;
            x_lo = alum_lo;
            x_up = alum_up;

            if (debug == TRUE)
            {
                fprintf(flog, "auto_dose: add alum to meet HAA5. Entering rootfind_and_mod_dose()! \n\n");
                write_dose_target(dose_unit, dose_location, target_param, target_unit, target_location, target, flog);
            }
            rootfind_and_mod_dose(mod_dose_check_target_ptr, train, dose_unit, dose_location, target_param, target_unit, target_location, target, x_lo, x_up, flog);
        }

        /* If achieving TOC and DBP regs reduces rapid mix pH below 5.5, adjust alum dose to achieve minimum pH. 
           Note: if this case occurs, either the TOC removal or DBP MCL regs will be violated. */
        if (rapid_mix->pH < pH_setpt_alum)
        { // ensure that pH does not drop too low in plant

            double alum_min = 0.0; // alum dose may need to be less than 18.0 mg/L, set lower guess to 0.0 mg/L

            dose_unit = ALUM;
            dose_location = 0;
            target_param = 'P';
            target_unit = RAPID_MIX;
            target_location = 0;
            target = pH_setpt_alum;
            x_lo = alum_min;
            x_up = alum_up;

            if (debug == TRUE)
            {
                fprintf(flog, "auto_dose: adjust alum to meet pH target. TOC removal or DBP MCL regs will be violated. \n");
                fprintf(flog, "auto_dose: alum dose minimum guess set to 0.0 mg/L to attempt to achieve pH >= 5.5. Entering rootfind_and_mod_dose()! \n\n");
                write_dose_target(dose_unit, dose_location, target_param, target_unit, target_location, target, flog);
            }
            rootfind_and_mod_dose(mod_dose_check_target_ptr, train, dose_unit, dose_location, target_param, target_unit, target_location, target, x_lo, x_up, flog);
        }

        /*===================== END ALUM ADJUSTMENT TO ACHIEVE DBP AND TOC REGS ==========================*/

        //Using root-finding algorithms, automate sodium hypochlorite dosing to maintain end of system chlorine residual
        /*===================== START EOS CL2 ADJUSTMENT TO MAINTAIN RESIDUAL ==========================*/

        if ((eos->FreeCl2 * MW_Cl2) != cl2_setpt)
        {

            dose_unit = HYPOCHLORITE;
            dose_location = 0;
            target_param = 'C';
            target_unit = END_OF_SYSTEM;
            target_location = 0;
            target = cl2_setpt;
            x_lo = naocl_lo;
            x_up = naocl_up;

            if (debug == TRUE)
            {
                fprintf(flog, "auto_dose: Add sodium hypochlorite to meet cl2 residual target.  Entering rootfind_and_mod_dose()! \n\n");
                write_dose_target(dose_unit, dose_location, target_param, target_unit, target_location, target, flog);
            }
            rootfind_and_mod_dose(mod_dose_check_target_ptr, train, dose_unit, dose_location, target_param, target_unit, target_location, target, x_lo, x_up, flog);
        }
        /*===================== END EOS CL2 ADJUSTMENT TO MAINTAIN RESIDUAL ==========================*/

        //Using root-finding algorithms, automate lime and/or CO2 dosing for corrosion control
        /*===================== START EFFLUENT pH AND ALKALINITY ADJUSTMENT FOR CORROSION CONTROL ==========================*/

        if (effluent->pH < pH_setpt_2)
        { // adjust pH by adding lime

            dose_unit = LIME;
            dose_location = 1;
            target_param = 'P';
            target_unit = WTP_EFFLUENT;
            target_location = 0;
            target = pH_setpt_2;
            x_lo = lime_lo;
            x_up = lime_up;

            if (debug == TRUE)
            {
                fprintf(flog, "auto_dose: add lime to meet corrosion control pH criteria. Entering rootfind_and_mod_dose()! \n\n");
                write_dose_target(dose_unit, dose_location, target_param, target_unit, target_location, target, flog);
            }

            rootfind_and_mod_dose(mod_dose_check_target_ptr, train, dose_unit, dose_location, target_param, target_unit, target_location, target, x_lo, x_up, flog);
        }
        else
        { // adjust pH by adding carbon dioxide

            dose_unit = CARBON_DIOXIDE;
            dose_location = 1;
            target_param = 'P';
            target_unit = WTP_EFFLUENT;
            target_location = 0;
            target = pH_setpt_2;
            x_lo = co2_lo;
            x_up = co2_up;

            if (debug == TRUE)
            {
                fprintf(flog, "auto_dose: add CO2 to meet corrosion control pH criteria. Entering rootfind_and_mod_dose()! \n\n");
            }

            rootfind_and_mod_dose(mod_dose_check_target_ptr, train, dose_unit, dose_location, target_param, target_unit, target_location, target, x_lo, x_up, flog);
        }
        /*======================== END EFFLUENT pH AND ALKALINITY ADJUSTMENT FOR CORROSION CONTROL ==========================*/

        // Check whether all setpoints have been met.
        if ((co2_1->Alk * MW_CaCO3 / 2.0) >= alk_setpt_1 - ERROR_TOL)
        {
            setptflag_alk_1 = TRUE;
            if (debug == TRUE)
                fprintf(flog, "auto_dose: alkalinity setpoint conditions met! Setting setptflag_alk_1 to TRUE. \n");
        }
        if (fabs(co2_1->pH - pH_setpt_1) < ERROR_TOL)
        {
            setptflag_pH_1 = TRUE;
            if (debug == TRUE)
                fprintf(flog, "auto_dose: pH setpoint conditions met! Setting setptflag_pH_1 to TRUE. \n");
        }
        if (fabs((eos->FreeCl2 * MW_Cl2) - cl2_setpt) < ERROR_TOL)
        {
            setptflag_cl2 = TRUE;
            if (debug == TRUE)
                fprintf(flog, "auto_dose: chlorine setpoint conditions met! Setting setptflag_cl2 to TRUE. \n");
        }

        if (fabs(effluent->pH - pH_setpt_2) < ERROR_TOL)  
        {
            setptflag_pH_2 = TRUE;
        }
        if ((eos->TTHM < (TTHM_MCL + ERROR_TOL) && eos->HAA5 < (HAA5_MCL + ERROR_TOL) &&
             ((influent->TOC - effluent->TOC) / influent->TOC * 100.0) >= (toc_rem_req - ERROR_TOL)))
        {
            setptflag_DBP = TRUE;
            if (debug == TRUE)
                fprintf(flog, "auto_dose: influent TOC >= 2.0 and TTHM, HAA5, and TOC removal conditions met! Setting setptflag_DBP to TRUE. \n");
        }
        else if ((eos->TTHM < (TTHM_MCL + ERROR_TOL) && eos->HAA5 < (HAA5_MCL + ERROR_TOL) &&
                  influent->TOC < 2.0))
        {
            if (debug == TRUE)
                fprintf(flog, "auto_dose: influent TOC < 2.0 (TOC removal exempt) and TTHM and HAA5 conditions met! Setting setptflag_DBP to TRUE. \n");
            setptflag_DBP = TRUE;
        }
        else if (fabs(rapid_mix->pH - pH_setpt_alum) < ERROR_TOL)
        {
            if (debug == TRUE)
                fprintf(flog, "auto_dose: TTHM, HAA5, or TOC could not be achieved! Alum dose set to get pH of 5.5 (the minimum). Setting setptflag_DBP to TRUE.  \n");
            setptflag_DBP = TRUE;
        }

        iter++;

        if (debug == TRUE)
        {
            fprintf(flog, "\n ERROR_TOL = %f\n", ERROR_TOL);
            fprintf(flog, "\n======== CHECK alkalinity setpt #1 =========\n");
            fprintf(flog, "Alkalinity: %f \n", co2_1->Alk * MW_CaCO3 / 2.0);
            fprintf(flog, "alk_setpt_1 - ERROR_TOL: %f \n", alk_setpt_1 - ERROR_TOL);
            fprintf(flog, "\n======== CHECK pH setpt #1 ========\n");
            fprintf(flog, "pH: %f \n", co2_1->pH);
            fprintf(flog, "(fabs(co2_1->pH - pH_setpt_1) < ERROR_TOL)? Answer: %d\n", (fabs(co2_1->pH - pH_setpt_1) < ERROR_TOL));
            fprintf(flog, "\n=========== CHECK DBP variables ==========\n");
            fprintf(flog, "End of system TTHM: %f \n", eos->TTHM);
            fprintf(flog, "End of system TTHM + ERROR_TOL: %f \n", TTHM_MCL + ERROR_TOL);
            fprintf(flog, "End of system HAA5: %f \n", eos->HAA5);
            fprintf(flog, "End of system HAA5: %f \n", HAA5_MCL + ERROR_TOL);
            fprintf(flog, "Influent TOC: %f \n", influent->TOC);
            fprintf(flog, "Effluent TOC: %f \n", effluent->TOC);
            fprintf(flog, "Percent TOC removal: %f\n", (influent->TOC - effluent->TOC) / influent->TOC * 100.0);
            fprintf(flog, "TOC removal requirement - ERROR_TOL: %f\n", toc_rem_req - ERROR_TOL);
            fprintf(flog, "Difference between rapid mix pH and pH setpoint: %f\n", fabs(rapid_mix->pH - pH_setpt_alum));
            fprintf(flog, "\n");
            fprintf(flog, "========== CHECK FLAGS for auto_dose iteration %d =========\n", iter);
            fprintf(flog, "setptflag_alk_1: %d \n", setptflag_alk_1);
            fprintf(flog, "setptflag_pH_1: %d \n", setptflag_pH_1);
            fprintf(flog, "setptflag_DBP: %d \n", setptflag_DBP);
            fprintf(flog, "setptflag_cl2: %d \n", setptflag_cl2);
            fprintf(flog, "setptflag_pH_2: %d \n", setptflag_pH_2);
            fprintf(flog, "\n");
        }

        if (iter > max_iter)
        {
            fprintf(stderr, "Maximum iterations exceeded in while loop for auto_dose()!\n");
            exit(EXIT_FAILURE);
        }

    } /* end while-loop */

    if (debug == TRUE)
    {
        fprintf(flog, "=========== Process Train =========== \n");
        writewtp(flog, train, flog); // write out process train to file
        fprintf(flog, "\n");
    }

    run_number++;
    fclose(flog);  // close log file

    /* Write water treatment plant process train to result file */
    if (fres) 
    {
        fprintf(fres, "=========== Process Train =========== \n");
        writewtp(fres, train, fres); // write out process train to file
        fprintf(fres, "\n");
    }

    // // Print out end of system pH
    // printf("Effluent pH: %f\n", effluent->pH);  // debugging
    // printf("End of system pH: %f\n", eos->pH);  // debugging

    return;
}

void write_dose_target(int dose_unit, int dose_location, char target_param, int target_unit, int target_location, double target, FILE *flog)
{

    /* Purpose: write information about chemical dose adjustments and water quality targets to output file */

    fprintf(flog, "dose_unit: %d \n", dose_unit);
    fprintf(flog, "dose_location: %d \n", dose_location);
    fprintf(flog, "target_param: %c \n", target_param);
    fprintf(flog, "target_unit: %d \n", target_unit);
    fprintf(flog, "target: %f \n", target);
    fprintf(flog, "target location: %d \n", target_location);
    fprintf(flog, "\n");

    return;
}