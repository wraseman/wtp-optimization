#include "wtp_optimize.h"
#include "wtp.h"
#include "auto_dose.h"
#include <vector>
#include <iostream>

/* Metadata for Monte Carlo water quality influent data */
#define HEADER true                              // header in Monte Carlo CSV (comma separated value) file
#define SIM_COL 0                                // simulation number column
#define YEAR_COL 1                               // year column
#define QUARTER_COL 2                            // quarter column
#define ALK_COL 3                                // alkalinity column
#define AMMONIA_COL 4                            // ammonia column
#define BROMIDE_COL 5                            // bromide column
#define CALCIUM_COL 6                            // calcium column
#define HARD_COL 7                               // hardness column
#define PH_COL 8                                 // pH column
#define TEMP_COL 9                               // temperature column
#define TOC_COL 10                               // total organic carbon column
#define TURB_COL 11                              // turbidity column
#define UV254_COL 12                             // UV254 absorbance column

void wtp(double *vars, double *objs, double *consts)
{
/*
* Purpose:
*   Define WTP problem formulation and alter treatment train to reflect Borg
*   decision variables.
*/
    /* Iteration variable */
    int i, j, k;

    /* Unit Process and Effluent structures */
    register struct UnitProcess *unit;
    // struct Effluent *influent;
    //    struct Effluent             *lime_1;  // first lime dosing location
    //    struct Effluent             *co2_1;  // first co2 dosing location
    //    struct Effluent             *rapid_mix;
    //    struct Effluent             *effluent;
    struct Effluent *eos;

    /* Initialize treatment train data structure */
    ProcessTrain *train;
    if ((train = AllocProcessTrain()) == NULL) // Allocate memory for process train pointer for WTP Model
    {
        fprintf(stderr, "Error: cannot allocate memory for process train.\n");
    }

    // int success = open_wtp("./in/wtp_train/conv.wtp", train, stdout);  // open and read in water treatment plant data. print treatment train to terminal.
    int success = open_wtp("./in/wtp_train/conv.wtp", train, NULL); // open and read in water treatment plant data. do not print treatment train to terminal.

    if (!success)
    {
        fprintf(stderr, "Error: treatment type either not specified or not of type string in config.json");
        exit(EXIT_FAILURE);
    }

    // /* Check that the treatment train was read correctly */
    // writewtp(stdout, train, stderr); // write out process train to file

    /* Read in Monte Carlo influent water quality data */
    static std::vector<std::vector<double>> monte_carlo; // 2D vector holding Monte Carlo data
    static int count = 0;                                // keep track of how many times wtp() has been called
    std::string filename = "./in/monte_carlo/influent-wq-data.csv";
    int num_wq_scenarios = simopt_params.num_wq_scenarios;  // number of water quality scenarios

    count++; // increment problem call count
    std::cout << "Starting WTP problem call number " << count << std::endl;

    // only read in the data the first time the wtp problem is called
    if (count == 1)
    {
        monte_carlo = read_montecarlo(filename.c_str(), HEADER, num_wq_scenarios);
    }

    /* Initialize accounting variables*/
    int lime_cntr;            /* Record how many times lime is added to water */
    int co2_cntr;             /* Record how many times co2 is added to water */
    int TTHM_exceed_cntr = 0; /* Count exceedances of TTHM safety threshold */
    int HAA5_exceed_cntr = 0; /* Count exceedances of HAA5 safety threshold */
    int toc_viol_cntr = 0;    /* Count violations of TOC removal threshold */
    int ct_viol_cntr = 0;     /* Count violations of contact time ratio */
    static int fe_count = 0;  /* Record number of function evaluations */
    fe_count++;

    /* Disinfection Byproduct Maximum Contaminant Level */
    double MCL_TTHM = 80.0; /* Maximum contaminant level for TTHMs */
    double MCL_HAA5 = 60.0; /* Maximum contaminant level for HAA5 */

    // /* Define unit costs for chemicals */
    // double unit_cost[] = {0.15, 0.2, 0.175, 1.5}; // unit cost for chemical doses [$/lb]: lime, carbon dioxide, alum, and sodium hypochlorite
    // double cost_conv_const = 8.34;                // cost conversion constant [lb*L/(mg*gallon)]
    // double year_to_days = 365.25;                 // number of days in a year

    /* Initialize automatic chemical dosing variables */
    double alk_setpt_1;       // alkalinity setpoint #1
    double pH_setpt_1;        // pH setpoint #1
    double DBP_safety_factor; // disinfection byproduct safety factor
    double cl2_setpt = 0.2;   // chlorine residual at end of distribution system (non-decision variable setpoint)
    double pH_setpt_2 = 8.00; // pH at end of the treatment plant (non-decision variable setpoint)

    double LRAA_TTHM = -DBL_MAX;    // locational running annual average for TTHMs
    double LRAA_HAA5 = -DBL_MAX;    // locational running annual average for HAA5s
    int LRAA_TTHM_exceed_count = 0; // keep track of how many times the LRAA regulation for TTHMs is exceeded, if any
    int LRAA_HAA5_exceed_count = 0; // keep track of how many times the LRAA regulation for HAA5s is exceeded, if any

    /* Initialize arrays for water quality parameters/chemical doses to be montiored for each timestep */
    double *eos_TTHM = (double *)malloc(N_TIMESTEPS * sizeof(double));   // maximum level of TTHM in the distribution system (assumed to be at the end of system)
    double *eos_HAA5 = (double *)malloc(N_TIMESTEPS * sizeof(double));   // maximum level of HAA5 in the distribution system (assumed to be at the end of system)
    double *eos_solids = (double *)malloc(N_TIMESTEPS * sizeof(double)); // amount of solids produced at the end of the treatment plant
    double *lime_dose = (double *)malloc(N_TIMESTEPS * sizeof(double));  // total lime dose used in the treatment plant
    double *co2_dose = (double *)malloc(N_TIMESTEPS * sizeof(double));   // total carbon dioxide dose used in the treatment plant
    double *alum_dose = (double *)malloc(N_TIMESTEPS * sizeof(double));  // alum dose used in the treatment plant

    /* Initialize values which are used to calculate the objective functions */
    double *freq_TTHM_exceed = (double *)malloc(num_wq_scenarios * sizeof(double)); // frequency of TTHM concentrations greater than MCL
    double *freq_HAA5_exceed = (double *)malloc(num_wq_scenarios * sizeof(double)); // frequency of HAA5 concentrations greater than MCL
    double *avg_solids = (double *)malloc(num_wq_scenarios * sizeof(double));       // average solids production
    double *avg_lime_dose = (double *)malloc(num_wq_scenarios * sizeof(double));    // average lime dose
    double *avg_co2_dose = (double *)malloc(num_wq_scenarios * sizeof(double));     // average carbon dioxide dose

    /* Reset arrays for each function evaluation */
    for (i = 0; i < N_TIMESTEPS; i++)
    {
        eos_TTHM[i] = -DBL_MAX;
        eos_HAA5[i] = -DBL_MAX;
        eos_solids[i] = -DBL_MAX;
        lime_dose[i] = -DBL_MAX;
        co2_dose[i] = -DBL_MAX;
        alum_dose[i] = -DBL_MAX;
    }

    for (i = 0; i < num_wq_scenarios; i++)
    {
        avg_solids[i] = -DBL_MAX;
        freq_TTHM_exceed[i] = -DBL_MAX;
        freq_HAA5_exceed[i] = -DBL_MAX;
        avg_lime_dose[i] = -DBL_MAX;
        avg_co2_dose[i] = -DBL_MAX;
    }

    /* Initialize summation variables */
    double sum_eos_solids = 0;
    double sum_eos_TTHM = 0;
    double sum_eos_HAA5 = 0;
    double sum_lime_dose = 0;
    double sum_co2_dose = 0;

    double sum_avg_solids = 0;
    double sum_avg_lime_dose = 0;
    double sum_avg_co2_dose = 0;

    /* Run for each Monte Carlo sample */
    for (k = 0; k < num_wq_scenarios; k++)
    {

        /* Reset accounting and summation variables for each Monte Carlo sample */
        TTHM_exceed_cntr = 0; /* Count exceedances of TTHM safety threshold */
        HAA5_exceed_cntr = 0; /* Count exceedances of HAA5 safety threshold */

        toc_viol_cntr = 0;  /* Count violations of TOC removal threshold */
        ct_viol_cntr = 0;     /* Count violations of contact time ratio */

        sum_eos_solids = 0;
        sum_eos_TTHM = 0;
        sum_eos_HAA5 = 0;
        sum_lime_dose = 0;
        sum_co2_dose = 0;

        /* Run time series of water qualities (timestep is a quarter of a year) */
        for (i = 0; i < N_YEARS; i++)
        {
            for (j = 0; j < N_QUARTERS_PER_YEAR; j++)
            {
                int row = j + i * N_QUARTERS_PER_YEAR + k * N_TIMESTEPS;

                /* Verify that data is read in correctly */
                int expected_sim = k + 1;
                int expected_year = i + 1;
                int expected_quarter = j + 1;

                int actual_sim = monte_carlo[row][SIM_COL];
                int actual_year = monte_carlo[row][YEAR_COL];
                int actual_quarter = monte_carlo[row][QUARTER_COL];

                validate_sim_year_quarter(expected_sim, actual_sim, "simulation number", SIM_COL, filename.c_str());
                validate_sim_year_quarter(expected_year, actual_year, "year", YEAR_COL, filename.c_str());
                validate_sim_year_quarter(expected_quarter, actual_quarter, "quarter", QUARTER_COL, filename.c_str());

                /* Define decision variables */
                alk_setpt_1 = vars[N_VAR_TYPES * j]; 
                pH_setpt_1 = vars[N_VAR_TYPES * j + 1];
                DBP_safety_factor = vars[N_VAR_TYPES * j + 2];

                /* Check that decision variables are being read in correctly */
                validate_vars_bounds(alk_setpt_1, "alkalinity setpoint #1", MIN_ALK, MAX_ALK);
                validate_vars_bounds(pH_setpt_1, "pH setpoint #1 ", MIN_PH, MAX_PH);
                validate_vars_bounds(DBP_safety_factor, "DBP safety factor", MIN_DBP_SF, MAX_DBP_SF);

                /* Verify that water quality data is within reasonable bounds */
                double pH = monte_carlo[row][PH_COL];
                validate_wq_bounds(pH, "pH", 0.0, 14.0); // note: pH can theoretically be below 0 and above 14, but it is highly unlikely for these conditions to arise
                double temp = monte_carlo[row][TEMP_COL];
                validate_wq_bounds(temp, "temperature", 0.0, 100.0);
                double toc = monte_carlo[row][TOC_COL];
                validate_wq_nonnegative(toc, "total organic carbon");
                double uv254 = monte_carlo[row][UV254_COL];
                validate_wq_nonnegative(uv254, "uv254 absorbance");
                double bromide = monte_carlo[row][BROMIDE_COL];
                validate_wq_nonnegative(bromide, "bromide");
                double alkalinity = monte_carlo[row][ALK_COL];
                validate_wq_nonnegative(alkalinity, "alkalinity");
                double calcium = monte_carlo[row][CALCIUM_COL];
                validate_wq_nonnegative(calcium, "calcium");
                double hardness = monte_carlo[row][HARD_COL];
                validate_wq_nonnegative(hardness, "total hardness");
                double nh3 = monte_carlo[row][AMMONIA_COL];
                validate_wq_nonnegative(nh3, "ammonia");
                double ntu = monte_carlo[row][TURB_COL];
                validate_wq_nonnegative(ntu, "turbidity");

                /* Reset chemical doses back to zero for process train */
                chem_reset(train);

                /* Set influent water quality of treatment plant based on Monte Carlo simulations */
                for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
                {
                    switch (unit->type)
                    {
                    case INFLUENT:
                        unit->data.influent->pH = pH;
                        unit->data.influent->temp = temp;
                        unit->data.influent->toc = toc;
                        unit->data.influent->uv254 = uv254;
                        unit->data.influent->bromide = bromide;
                        unit->data.influent->alkalinity = alkalinity;
                        unit->data.influent->calcium = calcium;
                        unit->data.influent->hardness = hardness;
                        unit->data.influent->nh3 = nh3;
                        unit->data.influent->ntu = ntu;
                        break;

                    default:
                        break;
                    }
                }

                /* Run model with automated chemical dosing based on water quality setpoints */
                auto_dose(train, pH_setpt_1, alk_setpt_1, cl2_setpt, pH_setpt_2, DBP_safety_factor, NULL);

                //					meet_TOC_req(train, cl2_setpt, pH_setpt_2);
                /* Gather info about water quality and chemical doses from WTP and distribution system for each model run */
                lime_cntr = 0; /* Reset lime dose counter */
                co2_cntr = 0;  /* Reset co2 dose counter */

                for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
                {
                    switch (unit->type)
                    {
                    // case INFLUENT:
                    //     influent = &unit->eff;
                    //     break;

                    case LIME:
                        if (lime_cntr == 0)
                        {
                            lime_dose[N_QUARTERS_PER_YEAR * i + j] = unit->data.lime->dose;
                            //								lime_1 = &unit->eff;
                        }
                        else if (lime_cntr == 1)
                        {
                            lime_dose[N_QUARTERS_PER_YEAR * i + j] += unit->data.lime->dose;
                        }
                        else
                        {
                            fprintf(stderr, "main: incorrect lime dosage accounting! \n");
                            exit(EXIT_FAILURE);
                        }
                        lime_cntr++;
                        break;

                    case CARBON_DIOXIDE:
                        if (co2_cntr == 0)
                        {
                            co2_dose[N_QUARTERS_PER_YEAR * i + j] = unit->data.chemical->co2;
                            //								co2_1 = &unit->eff;
                        }
                        else if (co2_cntr == 1)
                        {
                            co2_dose[N_QUARTERS_PER_YEAR * i + j] += unit->data.chemical->co2;
                        }
                        else
                        {
                            fprintf(stderr, "main: incorrect co2 dosage accounting! \n");
                            exit(EXIT_FAILURE);
                        }
                        co2_cntr++;
                        break;

                    case ALUM:
                        alum_dose[N_QUARTERS_PER_YEAR * i + j] = unit->data.alum->dose;
                        break;
                        //
                        //							case RAPID_MIX:
                        //								rapid_mix = &unit->eff;
                        //								break;

                        //						case HYPOCHLORITE:
                        //							naocl_dose[i] = unit->data.chemical->naocl;

                        //						case WTP_EFFLUENT:
                        //							effluent = &unit->eff;break;

                    case END_OF_SYSTEM:
                        eos = &unit->eff;
                        break;

                    default:
                        break;
                    }
                }

                /* Validation: record important parameters from each model run */
                // WJR include problem logic (1,2,3)
                //				fprintf(f_valid, "\n%d %d ", fe_count, i);  // function evaluations and Monte Carlo simulations
                //				fprintf(f_valid, "%f %f %f ", (influent->Alk * MW_CaCO3 / 2.0), influent->pH, influent->TOC);  // influent water quality
                //				fprintf(f_valid, "%f %f %f ", vars[0], vars[1], vars[2]);  // decision variables (chemical setpoints and safety factors)
                //				fprintf(f_valid, "%f %f %f %f ", lime_dose[i], co2_dose[i], alum_dose[i], naocl_dose[i]);  // chemical doses
                //				fprintf(f_valid, "%f %f %f %f ", lime_1->pH, co2_1->pH, rapid_mix->pH, eos->pH);  // // pH at various unit processes within plant
                //				fprintf(f_valid, "%f %f %f %f ", (lime_1->Alk * MW_CaCO3 / 2.0), (co2_1->Alk * MW_CaCO3 / 2.0), (rapid_mix->Alk * MW_CaCO3 / 2.0),
                //						(effluent->Alk * MW_CaCO3 / 2.0));  // alkalinity at varioius unit processes within plant
                //				fprintf(f_valid, "%f %f %f ", eos->solids, eos->TTHM, eos->HAA5);  // solids and DBP formation at end of system
                //				fprintf(f_valid, "%f %d", ((influent->TOC - effluent->TOC) / influent->TOC * 100.0), effluent->ec_meeting_step1);  // TOC removal/enhanced coagulation requirements

                /* Track DBP values and solids production by the end of the system */
                eos_TTHM[N_QUARTERS_PER_YEAR * i + j] = eos->TTHM;
                eos_HAA5[N_QUARTERS_PER_YEAR * i + j] = eos->HAA5;
                eos_solids[N_QUARTERS_PER_YEAR * i + j] = eos->solids;

                /* Calculate total solids produced, total chemical doses used, and cumulative DBP concentrations */
                sum_eos_TTHM += eos_TTHM[N_QUARTERS_PER_YEAR * i + j];
                sum_eos_HAA5 += eos_HAA5[N_QUARTERS_PER_YEAR * i + j];
                sum_eos_solids += eos_solids[N_QUARTERS_PER_YEAR * i + j];
                sum_lime_dose += lime_dose[N_QUARTERS_PER_YEAR * i + j];
                sum_co2_dose += co2_dose[N_QUARTERS_PER_YEAR * i + j];

                /* Record DBP safety threshold exceedances to calculate reliability objectives */
                if (eos_TTHM[N_QUARTERS_PER_YEAR * i + j] > MCL_TTHM - ERROR_TOL)
                {
                    TTHM_exceed_cntr++;
                }

                if (eos_HAA5[N_QUARTERS_PER_YEAR * i + j] > MCL_HAA5 - ERROR_TOL)
                {
                    HAA5_exceed_cntr++;
                }

                // Check total organic carbon constraint
                // printf("%f %d", ((influent->TOC - eos->TOC) / influent->TOC * 100.0), eos->ec_meeting_step1);  // debugging
                if ( (eos->ec_exempt == FALSE)  && (eos->ec_meeting_step1 == FALSE) )
                // if enhanced coagulation (EC) exemptions do not apply, and step 1 of EC is not met (i.e., TOC removal requirements)
                // then, there is a violation of the TOC removal contraint. 
                {
                    toc_viol_cntr += 1; 
                }

                // Check contact time ratio constraint
                // printf("ct_ratio: %f\n", eos->ct_ratio);  // debugging
                // printf("ct_ratio_c: %f\n", eos->ct_ratio_c);  // debugging
                // printf("ct_ratio_v: %f\n", eos->ct_ratio_v);  // debugging
                if ( (eos->ct_ratio < 1.0) || (eos->ct_ratio_c < 1.0) || (eos->ct_ratio_v < 1.0) )  
                // CT ratio must be >= 1.0 for giardia, cryptosporidium, and virus for compliance
                {
                    ct_viol_cntr += 1; 
                }
            }
        }

        // Locational running annual average DBPs constraint
        int n_LRAA_calcs = N_YEARS * N_QUARTERS_PER_YEAR - 3; // since LRAA_i=(Q_(i-3)+Q_(i-2)+Q_(i-1)+Q_i)/4,
                                                              // the LRAA can only be calculated the total number of timesteps - 3
        for (i = 0; i < (n_LRAA_calcs); i++)
        {

            LRAA_TTHM = (eos_TTHM[i] + eos_TTHM[i + 1] + eos_TTHM[i + 2] + eos_TTHM[i + 3]) / 4;

            if (LRAA_TTHM > MCL_TTHM)
            {
                LRAA_TTHM_exceed_count++;
            }

            LRAA_HAA5 = (eos_HAA5[i] + eos_HAA5[i + 1] + eos_HAA5[i + 2] + eos_HAA5[i + 3]) / 4;

            if (LRAA_HAA5 > MCL_HAA5)
            {
                LRAA_HAA5_exceed_count++;
            }
        }

        // Record objective values
        freq_TTHM_exceed[k] = (double)(TTHM_exceed_cntr) / (double)(N_TIMESTEPS); // Minimize TTHM exceedances of MCL
        freq_HAA5_exceed[k] = (double)(HAA5_exceed_cntr) / (double)(N_TIMESTEPS); // Minimize HAA5 exceedances of MCL
        avg_solids[k] = sum_eos_solids / (double)(N_TIMESTEPS);                   // Minimize average concentration of solids produced
        avg_lime_dose[k] = sum_lime_dose / (double)(N_TIMESTEPS);                 // Minimize average cumulative lime dose
        avg_co2_dose[k] = sum_co2_dose / (double)(N_TIMESTEPS);                   // Minimize average cumulative co2 dose

        sum_avg_solids += avg_solids[k];
        sum_avg_lime_dose += avg_lime_dose[k];
        sum_avg_co2_dose += avg_co2_dose[k];
    }

    /* Calculate objective function values */
    double worst_freq_TTHM_exceed = max(freq_TTHM_exceed, num_wq_scenarios);
    double worst_freq_HAA5_exceed = max(freq_HAA5_exceed, num_wq_scenarios);
    double avg_avg_eos_solids = sum_avg_solids / (double)(num_wq_scenarios);
    double avg_avg_lime_dose = sum_avg_lime_dose / (double)(num_wq_scenarios);
    double avg_avg_co2_dose = sum_avg_co2_dose / (double)(num_wq_scenarios);

    objs[0] = worst_freq_TTHM_exceed; // Minimize worst-case frequency of TTHM concentrations greater than MCL
    objs[1] = worst_freq_HAA5_exceed; // Minimize worst-case frequency of HAA5 concentrations greater than MCL
    objs[2] = avg_avg_eos_solids;     // Minimize expected solids production
    objs[3] = avg_avg_lime_dose;      // Minimize expected lime dose
    objs[4] = avg_avg_co2_dose;       // Minimize expected carbon dioxide dose

    /* Calculate constraint violations */
    // Note: all satisfied constraints must have a value of 0. Any non-zero value is considered a constraint violation. 
    consts[0] = LRAA_TTHM_exceed_count + LRAA_HAA5_exceed_count;  // locational running annual average DBP constraint
    consts[1] = toc_viol_cntr;  // total organic carbon removal constraint
    consts[2] = ct_viol_cntr;  // contact time ratio constraint

    /* Free dynamically allocated memory */
    free(eos_TTHM);
    free(eos_HAA5);
    free(eos_solids);
    free(lime_dose);
    free(co2_dose);
    //				free(alum_dose    );
    //				free(naocl_dose   );

    std::cout << "Finishing WTP problem call number " << count << std::endl;

    return;
}

//             /* Gather info about water quality and chemical doses from WTP and distribution system for each model run */
//             lime_cntr = 0; /* Reset lime dose counter */
//             co2_cntr = 0;  /* Reset co2 dose counter */

//             for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
//             {
//                 switch (unit->type)
//                 {
//                 case INFLUENT:
//                     influent = &unit->eff;
//                     break;

//                 case LIME:
//                     if (lime_cntr == 0)
//                     {
//                         lime_dose[N_QUARTERS_PER_YEAR * i + j] = unit->data.lime->dose;
//                         //								lime_1 = &unit->eff;
//                     }
//                     else if (lime_cntr == 1)
//                     {
//                         lime_dose[N_QUARTERS_PER_YEAR * i + j] += unit->data.lime->dose;
//                     }
//                     else
//                     {
//                         fprintf(stderr, "main: incorrect lime dosage accounting! \n");
//                         exit(EXIT_FAILURE);
//                     }
//                     lime_cntr++;
//                     break;

//                 case CARBON_DIOXIDE:
//                     if (co2_cntr == 0)
//                     {
//                         co2_dose[N_QUARTERS_PER_YEAR * i + j] = unit->data.chemical->co2;
//                         //								co2_1 = &unit->eff;
//                     }
//                     else if (co2_cntr == 1)
//                     {
//                         co2_dose[N_QUARTERS_PER_YEAR * i + j] += unit->data.chemical->co2;
//                     }
//                     else
//                     {
//                         fprintf(stderr, "main: incorrect co2 dosage accounting! \n");
//                         exit(EXIT_FAILURE);
//                     }
//                     co2_cntr++;
//                     break;

//                 case ALUM:
//                     alum_dose[N_QUARTERS_PER_YEAR * i + j] = unit->data.alum->dose;
//                     break;
//                     //
//                     //							case RAPID_MIX:
//                     //								rapid_mix = &unit->eff;
//                     //								break;

//                     //						case HYPOCHLORITE:
//                     //							naocl_dose[i] = unit->data.chemical->naocl;

//                     //						case WTP_EFFLUENT:
//                     //							effluent = &unit->eff;break;

//                 case END_OF_SYSTEM:
//                     eos = &unit->eff;
//                     break;

//                 default:
//                     break;
//                 }
//             }

//             /* Validation: record important parameters from each model run */
//             // WJR include problem logic (1,2,3)
//             //				fprintf(f_valid, "\n%d %d ", fe_count, i);  // function evaluations and Monte Carlo simulations
//             //				fprintf(f_valid, "%f %f %f ", (influent->Alk * MW_CaCO3 / 2.0), influent->pH, influent->TOC);  // influent water quality
//             //				fprintf(f_valid, "%f %f %f ", vars[0], vars[1], vars[2]);  // decision variables (chemical setpoints and safety factors)
//             //				fprintf(f_valid, "%f %f %f %f ", lime_dose[i], co2_dose[i], alum_dose[i], naocl_dose[i]);  // chemical doses
//             //				fprintf(f_valid, "%f %f %f %f ", lime_1->pH, co2_1->pH, rapid_mix->pH, eos->pH);  // // pH at various unit processes within plant
//             //				fprintf(f_valid, "%f %f %f %f ", (lime_1->Alk * MW_CaCO3 / 2.0), (co2_1->Alk * MW_CaCO3 / 2.0), (rapid_mix->Alk * MW_CaCO3 / 2.0),
//             //						(effluent->Alk * MW_CaCO3 / 2.0));  // alkalinity at varioius unit processes within plant
//             //				fprintf(f_valid, "%f %f %f ", eos->solids, eos->TTHM, eos->HAA5);  // solids and DBP formation at end of system
//             //				fprintf(f_valid, "%f %d", ((influent->TOC - effluent->TOC) / influent->TOC * 100.0), effluent->ec_meeting_step1);  // TOC removal/enhanced coagulation requirements

//             /* Track DBP values and solids production by the end of the system */
//             eos_TTHM[N_QUARTERS_PER_YEAR * i + j] = eos->TTHM;
//             eos_HAA5[N_QUARTERS_PER_YEAR * i + j] = eos->HAA5;
//             eos_solids[N_QUARTERS_PER_YEAR * i + j] = eos->solids;

//             /* Calculate total solids produced, total chemical doses used, and cumulative DBP concentrations */
//             sum_eos_TTHM += eos_TTHM[N_QUARTERS_PER_YEAR * i + j];
//             sum_eos_HAA5 += eos_HAA5[N_QUARTERS_PER_YEAR * i + j];
//             sum_eos_solids += eos_solids[N_QUARTERS_PER_YEAR * i + j];
//             sum_lime_dose += lime_dose[N_QUARTERS_PER_YEAR * i + j];
//             sum_co2_dose += co2_dose[N_QUARTERS_PER_YEAR * i + j];

//             /* Record DBP safety threshold exceedances to calculate reliability objectives */
//             if (eos_TTHM[N_QUARTERS_PER_YEAR * i + j] > MCL_TTHM - ERROR_TOL)
//             {
//                 TTHM_exceed_cntr++;
//             }

//             if (eos_HAA5[N_QUARTERS_PER_YEAR * i + j] > MCL_HAA5 - ERROR_TOL)
//             {
//                 HAA5_exceed_cntr++;
//             }
//         }
//     }

//     // Locational running annual average DBPs constraint
//     int n_LRAA_calcs = N_YEARS * N_QUARTERS_PER_YEAR - 3; // since LRAA_i=(Q_(i-3)+Q_(i-2)+Q_(i-1)+Q_i)/4,
//                                                           // the LRAA can only be calculated the total number of timesteps - 3
//     for (i = 0; i < (n_LRAA_calcs); i++)
//     {

//         LRAA_TTHM = (eos_TTHM[i] + eos_TTHM[i + 1] + eos_TTHM[i + 2] + eos_TTHM[i + 3]) / 4;

//         if (LRAA_TTHM > MCL_TTHM)
//         {
//             LRAA_TTHM_exceed_count++;
//         }

//         LRAA_HAA5 = (eos_HAA5[i] + eos_HAA5[i + 1] + eos_HAA5[i + 2] + eos_HAA5[i + 3]) / 4;

//         if (LRAA_HAA5 > MCL_HAA5)
//         {
//             LRAA_HAA5_exceed_count++;
//         }
//     }

//     // Record objective values
//     freq_TTHM_exceed[k] = (double)(TTHM_exceed_cntr) / (double)(N_TIMESTEPS); // Minimize TTHM exceedances of MCL
//     freq_HAA5_exceed[k] = (double)(HAA5_exceed_cntr) / (double)(N_TIMESTEPS); // Minimize HAA5 exceedances of MCL
//     avg_solids[k] = sum_eos_solids / (double)(N_TIMESTEPS);                   // Minimize average concentration of solids produced
//     avg_lime_dose[k] = sum_lime_dose / (double)(N_TIMESTEPS);                 // Minimize average cumulative lime dose
//     avg_co2_dose[k] = sum_co2_dose / (double)(N_TIMESTEPS);                   // Minimize average cumulative co2 dose

//     sum_avg_solids += avg_solids[k];
//     sum_avg_lime_dose += avg_lime_dose[k];
//     sum_avg_co2_dose += avg_co2_dose[k];
// }

//     // /* Calculate objective function values */
//     // double worst_freq_TTHM_exceed = max(freq_TTHM_exceed, N_SIMS);
//     // double worst_freq_HAA5_exceed = max(freq_HAA5_exceed, N_SIMS);
//     // double avg_avg_eos_solids = sum_avg_solids / (double)(N_SIMS);
//     // double avg_avg_lime_dose = sum_avg_lime_dose / (double)(N_SIMS);
//     // double avg_avg_co2_dose = sum_avg_co2_dose / (double)(N_SIMS);

//     // objs[0] = worst_freq_TTHM_exceed; // Minimize worst-case frequency of TTHM concentrations greater than MCL
//     // objs[1] = worst_freq_HAA5_exceed; // Minimize worst-case frequency of HAA5 concentrations greater than MCL
//     // objs[2] = avg_avg_eos_solids;     // Minimize expected solids production
//     // objs[3] = avg_avg_lime_dose;      // Minimize expected lime dose
//     // objs[4] = avg_avg_co2_dose;       // Minimize expected carbon dioxide dose

//     // consts[0] = LRAA_TTHM_exceed_count + LRAA_HAA5_exceed_count;

//     /* Free dynamically allocated memory */
//     free(eos_TTHM);
//     free(eos_HAA5);
//     free(eos_solids);
//     free(lime_dose);
//     free(co2_dose);
//     //				free(alum_dose    );
//     //				free(naocl_dose   );

//     return;
// }

// validate that water quality parameter is non-negative
void validate_wq_nonnegative(double value, const char *name)
{
    if (value < 0)
    {
        std::cout << "Error: " << name << " = " << value << ", negative values are non-physical for this parameter." << std::endl;
        exit(EXIT_FAILURE);
    }
}

// validate that water quality parameter is between some specified bounds
void validate_wq_bounds(double value, const char *name, double minimum, double maximum)
{
    if (value < minimum || value > maximum)
    {
        std::cout << "Error: " << name << " = " << value << " which is beyond the typical or allowable bounds for this parameter." << std::endl;
        exit(EXIT_FAILURE);
    }
}

// validate that the decision variable is between some specified bounds
void validate_vars_bounds(double value, const char *name, double minimum, double maximum)
{
    if (value < minimum || value > maximum)
    {
        std::cout << "Error: " << name << " = " << value << " which is beyond the typical or allowable bounds for this decision variable." << std::endl;
        exit(EXIT_FAILURE);
    }
}

// validate that the simulation number, year, and quarter are being read in correctly
void validate_sim_year_quarter(int expected, int actual, const char *name, int column, const char *filename)
{
    if (expected != actual)
    {
        std::cout << "Error: " << name << " being read in " << filename << " is incorrect" << std::endl;
        std::cout << "Expected " << name << ": " << expected << std::endl;
        std::cout << "Acutal " << name << ": " << actual << std::endl;
        exit(EXIT_FAILURE);
    }
}
