/*  wtp.h   Water Treatment Plant Header File -- Sept 21, 1993
*
*     Michael D. Cummins
*       July, 1993
*
*     project modifications done by:
*     Kenneth B. Curry, Malcolm Pirnie, INC.
*             October 12, 1993
*
*     Further modifications done by:
*     Warren J. Swanson, MPI-PHO, 10/98
*     (for WTP Model Modifications Project for U of Col. and U. of Cinn.)
*
*     William J. Raseman, U of Col., 6/17
*/

#ifndef WTP_H
#define WTP_H 1

#define VERSION "Version 2.2"
#define DATE "Aug. 2017"

/* ANSI C header files */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>
#include <ctype.h>
#include <float.h>
#include <time.h>

/* String functions that are Not-ANSI but supplied with
*  Borlands, Lattice, and Manx compilers.
*/
char *strupr(char *string); /* Convert string to upper case */
char *strlwr(char *string); /*   "      "     "  lower  "   */

#define TRUE 1
#define FALSE 0

/* Mathematical constants */
#define PI acos(-1.0)

/****************  WTP Global variables located in Globals.c  *************/

extern int coldflag; /* TRUE=Run model at low temperature  (TRUE/FALSE) */

/*
*  swflag, coagflag and filtflag
*
*  Version 1.4 does not estimate giardia or virus inactivation through
*  the process train.  Instead, a term named 'ct_ratio' is computed.  The
*  ct_ratio is intended to express the level of inactivation relative to
*  the SWTR.  A ct_ratio less than 1.0 indicates that the treatment process
*  is out of compliance, a value greater than 1.0 indicates compliance.
*
*  The 'swflag' is used to indicate if giardia or virus inactivation should
*  be used to determine the required CT.  The flag is needed at the first
*  point of chlorination and does not change through the process train.
*
*  'coagflag' and 'filtflag' are used to indicate that CT credits
*  should be given for coagulation and filtration.  In version 1.4, the
*  flags are cleared and re-set in runmodel() and processed by ct().
*
*  'coagflag' and 'filtflag' could not be implemented in the Effluent data
*  structure because the required CT needs coagflag and filtflag at the
*  first point of chlorination which may be before coagulation and
*  filtration.
*
*  Planned change:
*    1. Track Giardia and virus inactivation/removal through the treatment
*       process and eliminate ct_ratio.
*    2. Re-examine if coagflag and filtflag are needed.
*/
extern int swflag;        // TRUE=Surface Water                                     (TRUE/FALSE) */
extern int gw_virus_flag; // TRUE=Groundwater with virus disinfection required (TRUE/FALSE)*/
extern int coagflag;      //* TRUE=Process train has alum, iron, or lime coagulation (TRUE/FALSE) */
extern int conv_filtflag; //* TRUE=Process train has conv. filtration, coag/flocc/sed/filt, as
                          //    1st filtration stage (TRUE/FALSE) */
extern int filtflag;      //* TRUE=conv_filtflag, dir_filtflag, bagfflag, cartfflag, ssfflag, defflag, OR
                          //	nfflag is TRUE */
extern int filt2flag;     // TRUE= There is a second stage of some kind of filtration
extern int softflag;      //* TRUE=Process train has a lime dose with softening purpose (T/F) */
extern int soft2flag;     //* TRUE=Process train has a 2nd stage of softening (T/F) */
extern int floccflag;     //* TRUE=Process train has flocculation basin      (TRUE/FALSE) */
extern int sedflag;       //* TRUE=Process train has sedimentation           (TRUE/FALSE) */
extern int sedconvflag;   //* TRUE=Process train has sedimentation after coagflag is TRUE before filtflag is TRUE
extern int gacflag;       //* TRUE=Process train has GAC                     (TRUE/FALSE) */
extern int mfufflag;      //* TRUE=Process train has an MF/UF as 1st filter stage (TRUE/FALSE) */
extern int nfflag;        //* TRUE=Process train has an NF as 1st filter stage (TRUE/FALSE) */
extern int ssfflag;       //* TRUE=Process train has a slow sand filter as 1st filter stage (TRUE/FALSE) */
extern int defflag;       //* TRUE=Process train has a diatomaceous earth as 1st filter stage (TRUE/FALSE) */
extern int bagfflag;      //* TRUE=Process train has a bag filter as 1st filter stage (TRUE/FALSE) */
extern int cartfflag;     //* TRUE=Process train has a cartridge filter as 1st filter stage (TRUE/FALSE) */
extern int bankfflag;     //* 0=No bank filtration, 1=bank filtration for small credit, 2=bank filt for large credit */
extern int granf2flag;    //* TRUE=Process train has FILTER as 2nd filter stage (TRUE/FALSE) */
extern int gac2flag;      //* TRUE=Process train has GAC as second filter stage (TRUE/FALSE) */
extern int mfuf2flag;     //* TRUE=Process train has an MF/UF as 2nd filter stage (TRUE/FALSE) */
extern int nf2flag;       //* TRUE=Process train has an NF as 2nd filter stage (TRUE/FALSE) */
extern int ssf2flag;      //* TRUE=Process train has a slow sand filter as 2nd filter stage (TRUE/FALSE) */
extern int def2flag;      //* TRUE=Process train has a diatomaceous earth as 2nd filter stage (TRUE/FALSE) */
extern int bagf2flag;     //* TRUE=Process train has a bag filter as 2nd filter stage (TRUE/FALSE) */
extern int cartf2flag;    //* TRUE=Process train has a cartridge filter as 2nd filter stage (TRUE/FALSE) */
extern int lt2presedflag; //* TRUE=Process train has a lt2 toolbox presed basin with coag. (TRUE/FALSE) */
extern int cfeflag;       //* TRUE=First stage conv./direct filter meets CFE criterion for toolbox (T/F) */
extern int ifeflag;       //* TRUE=First stage conv./direct filter meets IFE criterion for toolbox (T/F) */
extern int uvflag;        //* TRUE=Process train has a UV process            (TRUE/FALSE) */
extern int lt2_wscp_flag; //* TRUE=Credit for 0.5-log Crypto reduction from watershed control granted */
extern int o3flag;        //* TRUE=Process train has ozonation               (TRUE/FALSE) */
extern int pre_o3flag;    //* TRUE=Process train has pre-ozonation           (TRUE/FALSE) */
extern int int_o3flag;    //* TRUE=Process train has int-ozonation           (TRUE/FALSE) */
extern int post_o3flag;   //* TRUE=Process train has post-ozonation          (TRUE/FALSE) */
extern int clo2flag;      //* TRUE=Process train has chlorine dioxide        (TRUE/FALSE) */
extern int dir_filtflag;  //* TRUE=Process train has coag and filtration as 1st filter stage (TRUE/FALSE) */
extern int bio_filtflag;  //* TRUE=Process train has biofiltration           (TRUE/FALSE) */
                          //* These flags determine what DBP formation model is in effect for each
                          //   process in the train as the main loop is executed in runmodel */
extern int rwdbpflag;     //* TRUE=use raw water DBP formation model for current process  */
extern int owdbpflag;     //* TRUE=use ozonated water DBP formation model for current process  */
extern int coagdbpflag;   //* TRUE=use coagulated water DBP formation model for current process  */
extern int modrw1dbpflag; //* TRUE=use raw water DBP formation model
                          //       modified with Pre-RM factor for current process  */
extern int modrw2dbpflag; //* TRUE=use raw water DBP formation model
                          //	   modified with Post-RM factor for current process  */
extern int gacmemdbpflag; //* TRUE=use gac-/nf(membrane)-treated water model for current process */

extern double nonconv_discredit; //sum of the UV inactivation and membrane, bankfilt, bagfilt, cartfilt
                                 // removal credits for LT2 compliance checking
extern double bin34_inactreqd;   //Amount of Crypto CT required based on the 1.0-log requirement for various techs. controlling

/* Global time counter and chlorine variable for modrw2dbp() */
extern double modrw2dbptime; /* =0 if no chlorine before RM; = RM mean det. time (hrs.) if
					chlorine before RM; set in runmodel(), used in
					modrw2dbp()*/
extern double modrw2dbpcl2;

/* Globals added for pathogen removal by filtration/membranes used in ct()*/

//These contain the total amount of disinfection required in logs before any credits for
// treatment, source water protection, etc., are granted
extern double tot_dis_req_g;
extern double tot_dis_req_c;
extern double tot_dis_req_v;

extern double tot_crypto_lr;
extern double tot_giardia_lr;
extern double tot_virus_lr;

/************  Data structures for Water Treatment Plant ***************/

struct Effluent
{                        /******  Data Packet for All Unit Processes *******/
                         /* Operating data:                                */
  double DegK;           /*   Temperature                          (Deg K) */
  double Flow;           /*   Average flow                           (MGD) */
  double Peak;           /*   Max hourly flow                        (MDG) */
                         /* Unit process counters:                         */
  short cl2cnt;          /*   Number of times chlorine added.              */
                         /* Measurable Water Quality Parameters:           */
  double pH;             /*   [H+]=pow(10,-pH)                         (-) */
  double TOC;            /*   Total Organic Carbon                  (mg/L) */
  double UV;             /*   Absorption at 254 nm (internal calcs.)(1/cm) */
  double UV_out;         /*   Absorption at 254 nm (for output)     (1/cm) */
  double SUVA;           /*   TSUVA = UVA/TOC*100                 (L/mg-m) */
  double Br;             /*   Bromide                         (Mole/Liter) */
  double Br_at_last_cl2; /*   Bromide at last point of Cl2 addition (Mole/Liter) */
  double Alk;            /*   Alkalinity                       (equ/Liter) */
  double NH3;            /*   Ammonia                         (Mole/Liter) */
  double Turbidity;      /*   Turbidity                              (NTU) */
  double FreeCl2;        /*   Free chlorine [HOCl]+[OCl-]     (Mole/Liter) */
  double NH2Cl;          /*   Monochloramime                  (Mole/Liter) */
  double NHCl2;          /*   Dichloramine (see breakpt.c)      (??/Liter) */

  /* Disinfection Byproducts:                       */
  double CHCl3;   /*   Chloroform                        (ug/Liter) */
  double CHBrCl2; /*   Bromodichloromethane              (ug/Liter) */
  double CHBr2Cl; /*   Dibromochloromethane              (ug/Liter) */
  double CHBr3;   /*   Bromoform                         (ug/Liter) */
  double TTHM;    /*   Total THM                         (ug/Liter) */

  double MCAA;  /*   Monochloroacetic  acid            (ug/Liter) */
  double DCAA;  /*   Dichloroacetic    acid            (ug/Liter) */
  double TCAA;  /*   Trichloroacetic   acid            (ug/Liter) */
  double MBAA;  /*   Monobromoacetic   acid            (ug/Liter) */
  double DBAA;  /*   Dibromoacetic     acid            (ug/Liter) */
  double BCAA;  /*   Bromochloroacetic acid            (ug/Liter) */
  double BDCAA; /*   Bromodichloroacetic acid          (ug/Liter) */
  double DBCAA; /*   Dibromochloroacetic acid          (ug/Liter) */
  double TBAA;  /*   Tribromoacetic acid               (ug/Liter) */
  double HAA5;  /*   Sum of first 5 HAAS listed above  (ug/Liter) */
  double HAA6;  /*   Sum of above 6 HAAS               (ug/Liter) */
  double HAA9;  /*   Sum of 9 HAAS                     (ug/Liter) */

  /* Internal Parameters:                           */
  double CO2_aq;    /*   Total carbonate                 (Mole/Liter) */
  double Ca_aq;     /*   Dissolved Calcium               (Mole/Liter) */
  double Ca_solid;  /*   CaCO3 in suspension             (Mole/Liter) */
  double Ca_sludge; /*   Sum of Ca_solid removed         (Mole/Liter) */
  double Mg_aq;     /*   Dissolved Magnesium             (Mole/Liter) */
  double Mg_solid;  /*   Mg(OH)2 in suspension           (Mole/Liter) */
  double Mg_sludge; /*   Sum of Mg_solid removed         (Mole/Liter) */
  double solids;    /*   Solids production                 (mg/Liter) */

  /*
    *  AlumDose and FericDose contain the cumulative coagulant dose up to
    *  the point of rapid mixing.  alum_rmv() and fecl_rmv() zero AlumDose
    *  and FericDose after the respective function has performed the
    *  removal calculation.
    */
  double AlumDose;     /*   Cum. Alum dose up to RM           (mg/Liter) */
  double FericDose;    /*   Cum. Feric dose up to RM          (mg/Liter) */
  double LimeDose;     /*   Cumulative lime dose       (mg/L as Ca(OH)2) */
  int floc_carry_cntr; /*=TRUE if floc carry-over occurs; if TRUE,
			     it is not allowed to occur again*/

  /* EquivAlumDose is used only to store the cumulative alum dose at the
  influent of the first RM encountered in the process train.  The value is
  needed in modrw1dbp() and modrw2dbp() and stores cum. alum + cum. ferric
  expressed as alum (mg alum/L). Value is set in runmodel(). */
  double EquivAlumDose;

  /* These sludge components must be tracked now that TOC removal and sludge
     generation occur in different locales in the process train */
  double toc_sludge;
  double alum_to_sludge;
  double iron_to_sludge;

  /*
    *  CBminusCA is the residual electrical charge of all cations and
    *  anions that are NOT tracked by phchange(), see phchange() for
    *  current list of ions that are tracked.  CBminusCA is used in
    *  phchange() to adjust the pH such that a charge balance is obtained.
    *  CBminusCA is initialized by calling phchange(unit,FALSE).
    *  The units are (equ/Liter).
    */
  double CBminusCA; /*                                    (equ/Liter) */

  /*
    *  'hours' is the cumulative contact time from the last point of
    *  chlorination.  'hours' is re-set to zero each time chlorine is
    *  added to the process train.
    */
  double hours;

  /*
    *  ct_ratio contains the ratio of CT achieved to the CT required by
    *  the SWTR for Giardia.  ct_ratio_v is the same, but for viruses
    *  These may be deleted in future versions of WTP model
    *  and replaced with giardia and virus concentrations.
    */
  double ct_ratio;   /*   CT achieved / CT required for Giardia        */
  double ct_ratio_v; /*   CT achieved / CT required for viruses        */
  double ct_ratio_c; /*   CT achieved / CT required for Cryptosporidium*/

  double ct_cl2;   /*   Cumulative CT achieved with free chlorine in the process train */
  double ct_nh2cl; /*   Cumulative CT achieved with chloramines in the process train */
  double ct_clo2;  /*   Cumulative CT achieved with free ClO2 in the process train */
  double ct_o3;    /*   Cumulative CT achieved with ozone in the process train */

  double log_inact_ach_g; /*   Cumulative log inactivation achieved for Giardia in the process train */
  double log_inact_ach_v; /*   Cumulative log inactivation achieved for Viruses in the process train */
  double log_inact_ach_c; /*   Cumulative log inactivation achieved for Crypto. in the process train */

  /*
    *  influent is set by influent() so that other unit process can
    *  access influent water quality parameters.  Example: Giardia in the
    *  source water is needed to determine disinfection requirements.
    */
  struct UnitProcess *influent;

  /*
    *  wtp_eff is set by runmodel() when a UnitProcess.type==WTP_EFFLUENT
    *  is found.  *wtp_eff is used by distribution samples to estimate
    *  chlorine decay and byproduct formation relative to the plant
    *  effluent, specifically 'hours' and other water quality parameters.
    */
  struct UnitProcess *wtp_effluent;

  /* The following variable needs to be documented.  */
  double cl2dose; /* Appears to be total chlorine residual.
		     *  (mg/L) as Cl2.  Check out decay.c     */

  /* The following added by WJS, 10/98 to save the influent WQ
     data to the RAPID MIX for use in determining DBPs for downstream
     processes in the modrw1dbp and modrw2dbp models */
  struct UnitProcess *last_rm_inf;

  double o3_minutes;               /* stores the cumulative time since the last ozone
			   addition point. Value set = 0.0 upon OZONE */
  struct UnitProcess *last_o3_inf; /* stores WQ information at the influent of the
					 last point of ozonation */
  double o3_res;                   /* Ozone Residual        (mg/L) */
  double BrO3;                     /* Bromate residual (ug/L) */
  double K_o3decay;                /* constant for ozone residual decay based on WQ at point of
			   ozonation and calculated in ozonate() */
  int first_o3_chamber;            /* Tells whether or not this u.p. is the first ozone chamber,
                           for the purposes of CT calc. (no CT, if TRUE) */

  /* The following is used to store the type of DBP Model Used in a Given Unit Process
     Assignments occur in influent() and choose_model(). */
  const char *dbpmodel; /* contains the name of the DBP model used in a the unit process */

  /* The following two variables store unit process residence time and cumulative
     train residence time for inclusion in output tables */
  double processtime; /* unit process mean residence time (hours) */
  double traintime;   /* cumulative process train residence time (hours) */

  double clo2_res;     /* chlorine dioxide residual in effluent of unit (mg clo2/L) */
  double clo2dose;     /* Contains the last non-zero ClO2 dose (mg/L) */
  double clo2_avg;     /* Contains average clo2 concentration for a unit process for
			  use in CT calcs. It is determined in basn_dbp and filt_dbp*/
  double clo2_minutes; /* time since last chlorine dioxide addition point (minutes) */
  double chlorite;     /* chlorite concentration in unit process effluent (mg clo2-/L) */
  int clo2_cntr;       /* Counts number of clo2 addition points in the process train */
  int clo2_init_decay; /*TRUE if it is the unit process (with a volume) where initial
			   demand is exerted; FALSE otherwise*/

  double tox;                  /* Total organic halides (TOX) in unit process effluent (ug/L) */
  double log_required_g;       /*Log inactivation req. for Giardia by cl2, chcl3, o3 and/or clo2 (logs)*/
  double log_required_v;       /*Log inactivation req. for Virus by cl2, chcl3, o3 and/or clo2 (logs)*/
  double log_required_c;       /*Log inactivation req. for Crypto by cl2, chcl3, o3 and/or clo2 (logs)*/
  int limesoftening;           /* "TRUE  if a lime addition for softening has occurred in the train
			    FALSE otherwise; if FALSE, phchange() will not consider precip. of
			    CaCO3 or Mg(OH)2.*/
  int ec_meeting_step1;        /* 0 = no, not enough TOC removal; 1 = yes, Step 1 met */
  int ec_exempt;               /* 0 = no, raw water TOC and UV and finished water TOC all >= 2; */
                               /* 1 = yes, one of these is < 2 */
  int pre_re_chlor_flag;       /* 1 = apply pre_re_chlor DBP adjustment factors in coagdbp;
			       0 = do not apply, but possibly apply downstream;
			       2 = do not apply b/c it has already been applied and no longer should be*/
  int pre_chlor_flag;          /* Prechlorination in the process train has occurred */
  double pre_chlor_TTHM;       /* TTHMs formed up to Cl2 point where pre_re_chlor flag becomes TRUE */
  double pre_chlor_HAA6;       /* HAA6 formed up to Cl2 point where pre_re_chlor flag becomes TRUE */
  double pre_chlor_dose_track; /* Cumulative chlorine dose (mg/L) applied prior to a unit process with a volume */
  double pre_chlor_dose_ratio; /* Total pre-chlorine dose (mg/L), applied directly before and/or after the rapid mix,
				  divided by raw water TOC (mg/L) - used to set limits on implementing the
				  pre-/re-chlorination DBP adjustment factors */

}; /************  End of struct Effluent  ************/

/*****************************************/
struct ProcessTrain
{                           /* Control Structure for Process Train   */
  struct UnitProcess *head; /*   First UnitProcess in ProcessTrain   */
  struct UnitProcess *null; /*   Always NULL                         */
  struct UnitProcess *tail; /*   Last UnitProcess in ProcessTrain    */
  char file_name[120];      /*   Full path and extension             */
};                          /*****************************************/

struct UnitProcess
{                           /********** Treatment Process ***************/
  struct UnitProcess *next; /* Double Linked list                       */
  struct UnitProcess *prev; /*   "      "     "                         */
  short type;               /* Defined unit process types               */
  short pad;                /* Maintain 32 bit alinment of pointers     */
  union {                   /* Design and operating parameters:          */
    void *ptr;
    struct Influent *influent;
    struct Alum *alum;
    struct Gac *gac;
    struct Filter *filter;
    struct Basin *basin;
    struct Mfuf *mfuf;
    struct Nf *nf;
    struct Ssf *ssf;
    struct Def *def;
    struct Altf *altf;
    struct Bankf *bankf;
    struct Presed *presed;
    struct Uvdis *uvdis;
    struct Iron *iron;
    struct chemical *chemical;
    struct clo2 *clo2;
    struct lime *lime;
    struct Avg_tap *avg_tap;
    struct End_of_system *end_of_system;
  } data;

  struct Effluent eff;
};

/********** Data structures for simulation parameters **********/
struct SimParam
{                             /* Control structure for simulation/optimization parameters */
  struct Simulation *sim;     /* Simulation mode parameters */
  struct SimOpt *simopt;      /* Simulation-optimization mode parameters */
  struct MonteCarlo *monte;   /* Monte Carlo sampling parameters */
  struct TimeSeries *tseries; /* Time series modeling parameters */
  char file_name_in[120];     /* Parameter input file name (full path and extension) */
};

struct Simulation
{                          /* Simulation mode parameters */
  int problem;             /* Choose among pre-defined problem formulations (1/2/3) */
  int n_func_eval;         /* Number of function evaluations */
  char file_name_in[120];  /* Decision variable input file name (full path and extension) */
  char file_name_out[120]; /* Results output file name (full path and extension) */
};

struct SimOpt
{                          /* Simulation-optimization mode parameters */
  int problem;             /* Choose among pre-defined problem formulations (1/2/3) */
  int n_func_eval;         /* Number of function evaluations */
  char file_name_out[120]; /* Results output file name (full path and extension) */
};

struct MonteCarlo
{                         /* Monte Carlo sampling parameters */
  int mc_flag;            /* Turn Monte Carlo sampling on or off (TRUE/FALSE) */
  int n_samples;          /* Number of Monte Carlo samples for different raw waters */
  char file_name_in[120]; /* Monte Carlo samples input file name (full path and extension) */
};

struct TimeSeries
{                         /* Time series modeling parameters */
  int ts_flag;            /* Turn timeseries mode on or off (TRUE/FALSE) */
  int n_years;            /* Number of years to be modeled */
  int timesteps_per_year; /* Number of timesteps per year */
};

/* define(s) for UnitProcess.type */
#define VACANT 0
#define INFLUENT 1
#define RAPID_MIX 2
#define SLOW_MIX 3
#define SETTLING_BASIN 4
#define FILTER 5
#define BASIN 6
#define CONTACT_TANK 7
#define CLEARWELL 8
#define O3_CONTACTOR 9
#define GAC 10
#define MFUF_UP 11
#define NF_UP 12
#define SLOW_FILTER 13
#define DE_FILTER 14
#define BAG_FILTER 15
#define CART_FILTER 16
#define BANK_FILTER 17
#define PRESED_BASIN 18
#define UV_DIS 19
//#define  MEMBRANE

/* Chemical addition: */
#define ALUM 20
#define IRON 21
#define CHLORINE 22
#define SULFURIC_ACID 23
#define LIME 24
#define SODA_ASH 25
#define AMMONIA 26
#define PERMANGANATE 27
#define CARBON_DIOXIDE 28
#define OZONE 29
#define CHLORINE_DIOXIDE 30
#define SODIUM_HYDROXIDE 31
#define HYPOCHLORITE 32
#define SULFUR_DIOXIDE 33
#define AMMONIUM_SULFATE 34

/* Sample Points: */
#define WTP_EFFLUENT 35
#define AVG_TAP 36
#define END_OF_SYSTEM 37
#define LOCATION_1 38

/* Unit process data packets for UnitProcess.data */

struct Influent
{                          /* Raw Water Data                          */
  double pH;               /* (-)                                     */
  double temp;             /* Average temperature (C)                 */
  double low_temp;         /* Low temperature for disinfection (C)    */
  double toc;              /* (mg/L)                                  */
  double uv254;            /* (1/cm)                                  */
  double bromide;          /* (mg/L)                                  */
  double alkalinity;       /* (mg/L as CaCO3)                         */
  double calcium;          /* Calcium Hardness (mg/L as CaCO3)        */
  double hardness;         /* Total   Hardness (mg/L as CaCO3)        */
  double nh3;              /* Ammonia (mg/L as N)                     */
  double ntu;              /* Turbidity                               */
  double peak_flow;        /* Peak Hourly Flow for disinfection (MGD) */
  double avg_flow;         /* Average Flow (MGD)                      */
  int swflag;              /* TRUE=Surface Water; FALSE=Ground Water  */
  double crypto_conc;      /* Source water Crypto conc used to determine required Crypto. disinfection
			   level based on LT2 Rule */
  int lt2_wscp_flag;       /* There is a watershed control program in place that meets LT2 Rule
			    Toolbox criteria for achieving 0.5-log Crypto LR credit */
  int gw_virus_flag;       /* If a GW system, is virus disinfection required (TRUE/FALSE)*/
  double gw_tot_dis_req_v; //Total virus disinfection required (logs) for a groundwater system
                           //that is required to treat for virus (under proposed GW Rule)
};                         /*******************************************/

struct Alum
{              /* Alum Addition  ********************/
  double dose; /*   (mg/L) as Al2*(SO4)3-14H2O      */
};             /*************************************/

struct Gac
{                     /* Granular Activated Carbon ******************/
  double ebct;        /*   Empty bed contact time         (minutes) */
  double regen;       /*   Regeneration frequency            (days) */
  char config;        /* "S"=Single Contactor, "B"=Blended Parallel
							       System*/
  char toc_calc;      /* "M"=Maximum TOC breakthrough calculated    */
                      /* "A"=Average TOC breakthrough calculated    */
  double crypto_lr_2; /* Crypto log removal for second stage filter  */
  int filt_stage;     //Filtration stage of this filter based on PT position
                      // 0=not 1st or 2nd, 1=1st, 2=2nd
  /***  Additions: ***************************************************/
};

struct Filter
{                         /* Filter design and operating parameters *****/
  double volume;          /*   Filtration volume size for theoretical t */
  double fi_mean;         /*   Ratio mean/theoretical detention time    */
  double fi_t_ten;        /*   Ratio t10/theoretical detention time     */
  int cl2_bkwsh;          /*   TRUE = Chlorinated Backwash Used         */
  char media;             /* "S" = anthracite/sand, "G" = GAC           */
  double giardia_lr_conv; /* Giardia log removal for conventional plant */
  double virus_lr_conv;   /* Virus log removal for conventional plant */
  double crypto_lr_conv;  /* Crypto log removal for conventional plant */
  double giardia_lr_df;   /* Giardia log removal for direct filter plant */
  double virus_lr_df;     /* Virus log removal for direct filter plant */
  double crypto_lr_df;    /* Crypto log removal for direct filter plant */
  int cfe_turb_flag;      /* TRUE = Combined filt. eff. turb. meets LT2 Toolbox
			     criteria, which means eligible for 0.5-log Crypto */
  int ife_turb_flag;      /* TRUE = Individ. filt. eff. turb. meets LT2 Toolbox
			     criteria, which means eligible for 1.0-log Crypto */
  double crypto_lr_2;     /* Crypto log removal for second stage filter  */
  int filt_stage;         //Filtration stage of this filter based on PT position
                          // 0=not 1st or 2nd, 1=1st, 2=2nd
};                        /**********************************************/

struct Basin
{                  /* Settling, contact, storage              */
  double sb_mean;  /*   Ratio mean/theoretical detention time */
  double sb_t_ten; /*   Ratio t10/theoretical detention time  */
  double volume;   /*                                    (MG) */
};                 /*******************************************/
/*
*  Use struct Basin for Rapid_Mix, Slow_Mix, Settling_Basin, Contact_Tank,
*  and Clearwell.
*/

struct Mfuf
{                     /* Microfiltration/Ultrafiltration ***************/
  double recover;     /*   Effluent (permeate) flow/influent flow  (%) */
  double giardia_lr;  /* Log removal of Giardia by mfuf                */
  double virus_lr;    /* Log removal of virus by mfuf                  */
  double crypto_lr_1; /* Log removal of crypto by mfuf as 1st stage filter */
  double crypto_lr_2; /* Log removal of crypto by mfuf as 2nd stage filter */
  int filt_stage;     //Filtration stage of this filter based on PT position
                      // 0=not 1st or 2nd, 1=1st, 2=2nd
};                    /*************************************************/

struct Nf
{                        /* Nanofiltration ********************************/
  double mwc;            /*   Molecular Weight cutoff           (gm/Mole) */
  double recover;        /*   Effluent (permeate) flow/influent flow  (%) */
  double giardia_lr;     /* Log removal of Giardia by membranes           */
  double virus_lr;       /* Log removal of virus by membranes             */
  double crypto_lr;      /* Log removal of crypto by membranes            */
  double toc_rem;        /* TOC removal by nanofiltration (%)             */
  double uva_rem;        /* UVA removal by nanofiltration (%)             */
  double br_rem;         /* Bromide removal by nanofiltration (%)         */
  double treat_fraction; /* Fraction of the flow treated (not by-passed)  */
  int filt_stage;        //Filtration stage of this filter based on PT position
                         // 0=not 1st or 2nd, 1=1st, 2=2nd
};                       /*************************************************/

struct Ssf
{                       /* Slow sand filter design/op. params. *****/
  double volume;        /*   Filtration volume size for theoretical det. time calc. */
  double fi_mean;       /*   Ratio mean/theoretical detention time    */
  double fi_t_ten;      /*   Ratio t10/theoretical detention time     */
  double toc_rem_lowt;  /* TOC removal in low temperature range */
  double toc_rem_midt;  /* TOC removal in mid temperature range */
  double toc_rem_hight; /* TOC removal in high temperature range */
  double giardia_lr;    /* Log removal of giardia (assumes primary filtration stage */
  double virus_lr;      /* Log removal of viruses (assumes primary filtration stage */
  double crypto_lr_1;   /* Log removal of cryptosporidium if primary filter stage  */
  double crypto_lr_2;   /* Log removal of cryptosporidium if secondary filter stage  */
  int filt_stage;       //Filtration stage of this filter based on PT position
                        // 0=not 1st or 2nd, 1=1st, 2=2nd
};

struct Def
{                     /* Diatomaceous earth filter design/op. params. *****/
  double volume;      /*   Filtration volume size for theoretical det. time calc. */
  double fi_mean;     /*   Ratio mean/theoretical detention time    */
  double fi_t_ten;    /*   Ratio t10/theoretical detention time     */
  double giardia_lr;  /* Log removal of giardia (assumes primary filtration stage */
  double virus_lr;    /* Log removal of viruses (assumes primary filtration stage */
  double crypto_lr_1; /* Log removal of cryptosporidium if primary filter stage  */
  double crypto_lr_2; /* Log removal of cryptosporidium if secondary filter stage  */
  int filt_stage;     //Filtration stage of this filter based on PT position
                      // 0=not 1st or 2nd, 1=1st, 2=2nd
};

struct Altf
{                     /* Alternative filter (bags and cartridges) design/op. params. *****/
  double giardia_lr;  /* Log removal of giardia (assumes primary filtration stage */
  double virus_lr;    /* Log removal of viruses (assumes primary filtration stage */
  double crypto_lr_1; /* Log removal of cryptosporidium if primary filter stage  */
  double crypto_lr_2; /* Log removal of cryptosporidium if secondary filter stage  */
  int filt_stage;     //Filtration stage of this filter based on PT position
                      // 0=not 1st or 2nd, 1=1st, 2=2nd
};

struct Bankf
{                         /* Alternative filter (bags and cartridges) design/op. params. *****/
  int eligible_lt2;       /*   TRUE = Bank filter meets criteria to receive Crypto removal credit */
  double distance;        /*   Distance between withdrawal wells and surface water source */
  double crypto_lr_close; /* Log removal of cryptosporidium for distance less than 50 ft  */
  double crypto_lr_far;   /* Log removal of cryptosporidium for distance at least 50 ft */
};

struct Presed
{                   /* Presedimentation basin (primarily included here for LT2 Toolbox */
  double volume;    /* (MG) */
  double sb_mean;   /* Ratio mean/theoretical detention time */
  double sb_t_ten;  /* Ratio t10/theoretical detention time  */
  int eligible_lt2; /* TRUE = Presed is new and meets criteria to receive Crypto removal credit */
  double crypto_lr; /* Log removal of cryptosporidium */
};                  /*******************************************/

struct Uvdis
{                    /* Ultraviolet Disinfection */
  double giardia_li; /* Log inactivation of giardia */
  double virus_li;   /* Log inactivation of viruses */
  double crypto_li;  /* Log inactivation of crypto  */
};

struct Iron
{              /* Feric chloride addition */
  double dose; /*    (mg/L) as FeCl3-6H2O */
};             /***************************/

struct First_tap
{                /* First customer  ****************************/
  double minute; /*   Resident time to first customer (minute) */
};               /**********************************************/

struct Avg_tap
{              /* Also general sample from distribution system */
  double days; /*   Average resident time               (days) */
};             /************************************************/

struct End_of_system
{              /* Last customer *******************/
  double days; /*   Maximum resident time (days)  */
};             /***********************************/

struct chemical
{               /* General chemical addition */
  double chlor; /*   (mg/L) as Cl2     */
  double o3;    /*          as O3      */
  double nh3;   /*          as N       */
  double kmno4; /*          as KMnOH   */
  double naoh;  /*          as NaOH    */
  double h2so4; /*          as H2SO4   */
  double soda;  /*          as Na2CO3  */
  double co2;   /*          as CO2     */
  double naocl; /*          as Cl2     */
  double so2;   /*          as SO2     */
};

struct clo2
{                    /* Addition of Chlorine dioxide */
  double dose;       /* Dose as mg/L of ClO2 */
  double conversion; /* % Coversion of chlorine dioxide to chlorite (%) entered by user */
};

struct lime
{               /* Lime Addition */
  double dose;  /* Dose as mg/L of Ca(OH)2 */
  char purpose; /* Purpose of lime addition: pH adjustment "P" or softening "S" */
};

/* Default Values for Unit Process Data (located in Globals.c) */
extern struct Influent default_influent;
extern struct Alum default_alum;
extern struct Gac default_gac;
extern struct Filter default_filter;
extern struct Basin default_basin;
//extern struct Membrane      default_membrane;
extern struct Mfuf default_mfuf;
extern struct Nf default_nf;
extern struct Ssf default_ssf;
extern struct Def default_def;
extern struct Altf default_bagf;
extern struct Altf default_cartf;
extern struct Bankf default_bankf;
extern struct Presed default_presed;
extern struct Uvdis default_uvdis;
extern struct Iron default_iron;
extern struct chemical default_chemical;
extern struct clo2 default_clo2;
/* extern struct WTP_effluent  default_wtp_effluent;   No longer needed - WJS, 11/98 */
extern struct First_tap default_first_tap;
extern struct Avg_tap default_avg_tap;
extern struct Avg_tap default_location_1;
extern struct End_of_system default_end_of_system;
extern struct Basin default_rapid_mix;
extern struct Basin default_slow_mix;
extern struct Basin default_settling_basin;
extern struct Basin default_contact_tank;
extern struct Basin default_clearwell;
extern struct Basin default_o3_contactor;
extern struct lime default_lime;

/***** DataID contains text information describing a data element.  ********/
struct DataID
{
  const char *id;     /* ID code used in data file.                               */
  const char *name;   /* Full name of data element for use with data entry screen */
  const char *units;  /* units of measurment (mg/L), (MGD), etc                   */
  const char *fmtout; /* output format for data element. Example "%6.2lf" */
  const char *fmtin;  /* input  format for data element. Example "%lf"    */
};

/****  UnitInfo contains information about a unit process.  ***************/
struct UnitInfo
{
  const char *name;          /* Name of unit process used in data file.        */
  struct DataID *data; /* Array of DataID elements in unit process       */
                       /* read() and write() are File IO functions       */
  int (*read)(char *buffer, short i, struct UnitProcess *unit);
  int (*write)(char *buffer, short i, struct UnitProcess *unit);
};

/*
*  UnitProcessTable[] is the controlling structure that contains text
*  describing all unit process known by WTP model and the data elements.
*/
extern struct UnitInfo UnitProcessTable[];

/****  SimParamInfo contains information simulation parameters.  ***************/
struct SimParamInfo
{
  char *name;          /* Name of simulation parameter.                  */
  struct DataID *data; /* Array of DataID elements in unit process       */
                       /* read() and write() are File IO functions       */
  int (*read)(char *buffer, short i, struct SimParam *param);
  int (*write)(char *buffer, short i, struct SimParam *param);
};

/*
*  SimParamTable[] is the controlling structure that contains text
*  describing the simulation parameters. 
*/
extern struct SimParamInfo SimParamTable[];

/******************  WTP Model constants  *****************/

/* Molecular Weights (mg/Mole) */
#define MW_CaCO3 100090L
#define MW_MgOH2 58320L
#define MW_Br 79900L
#define MW_NH3 14010L
#define MW_Cl2 70900L
#define MW_ClO2 67500L

#define MW_CaOH2 74096L
#define MW_CO2 44010L
#define MW_SODA 105990L
#define MW_NaOH 39998L
#define MW_H2SO4 98076L
#define MW_SO2 64058L

#define MW_LIME MW_CaOH2

//DBP molecular weights (mg/mol)
#define MW_CHCL3 119400L
#define MW_BDCM 163800L
#define MW_DBCM 208300L
#define MW_CHBR3 252700L
#define MW_MCAA 93500L
#define MW_MBAA 137900L
#define MW_DBAA 216800L
#define MW_BCAA 172400L
#define MW_DCAA 127900L
#define MW_TCAA 162400L
#define MW_BDCAA 206800L
#define MW_DBCAA 251300L
#define MW_TBAA 295700L

/* Bisection Method constants */
#define ERROR_TOL 1E-2 // error tolerance for search

/***********************   Functions   *******************************/

/* WTP functions:
* The following functions deal with the C language structure of the process
* train and unit process, are ANSI, and support the GUI functions.
*/
int cli_wtp(int argc, char *argv[]);

int new_wtp(struct ProcessTrain *train);
int edit_wtp(struct UnitProcess *unit);

int save_wtp(const char *name, struct ProcessTrain *train, FILE *ferr);
int writewtp(FILE *fout, struct ProcessTrain *train, FILE *ferr);

int open_wtp(const char *file_name, struct ProcessTrain *train, FILE *fout);
int read_wtp(FILE *fin, struct ProcessTrain *train, FILE *fout);

int list_wtp(struct ProcessTrain *train, FILE *fout);
int run_wtp(struct ProcessTrain *train, FILE *fout);
int thm_wtp(struct ProcessTrain *train, FILE *fout);

/* Old function definitions for WTP Model with Borland C++ Builder-based GUI (WJR - 01/2019) */
// int  open_wtp(char *file_name, struct ProcessTrain *train, FILE *fout,
//                               void (*wait)(char *msg), long max_line_count );
// int  read_wtp(FILE *fin, struct ProcessTrain *train, FILE *fout,
//                               void (*wait)(char *msg), long max_line_count );
// int  list_wtp(struct ProcessTrain *train,FILE *fout,
//                               void (*wait)(char *msg), long max_line_count );

// int  run_wtp( struct ProcessTrain *train,FILE *fout,
//                               void (*wait)(char *msg), long max_line_count );
// int  thm_wtp( struct ProcessTrain *train,FILE *fout,
//                               void (*wait)(char *msg), long max_line_count );

/* read_13() is for reading version 1.3 and earlier data files.*/
/* this backward compatibility thing is to be put on hold */
/*int  read_13(FILE *fin, struct ProcessTrain *train, FILE *fout,
                              void (*wait)(char *msg), long max_line_count );  */

/* The following are unit process I/O functions are located in data_wtp.c */
int read_influent(char *buffer, short i, struct UnitProcess *unit);
int read_alum(char *buffer, short i, struct UnitProcess *unit);
int read_gac(char *buffer, short i, struct UnitProcess *unit);
int read_filter(char *buffer, short i, struct UnitProcess *unit);
int read_basin(char *buffer, short i, struct UnitProcess *unit);
int read_mfuf(char *buffer, short i, struct UnitProcess *unit);
int read_nf(char *buffer, short i, struct UnitProcess *unit);
int read_ssf(char *buffer, short i, struct UnitProcess *unit);
int read_def(char *buffer, short i, struct UnitProcess *unit);
int read_altf(char *buffer, short i, struct UnitProcess *unit);
int read_bankf(char *buffer, short i, struct UnitProcess *unit);
int read_presed(char *buffer, short i, struct UnitProcess *unit);
int read_uvdis(char *buffer, short i, struct UnitProcess *unit);
int read_iron(char *buffer, short i, struct UnitProcess *unit);
int read_avg_tap(char *buffer, short i, struct UnitProcess *unit);
int read_end_sys(char *buffer, short i, struct UnitProcess *unit);
int read_chlorine(char *buffer, short i, struct UnitProcess *unit);
int read_h2so4(char *buffer, short i, struct UnitProcess *unit);
int read_naoh(char *buffer, short i, struct UnitProcess *unit);
int read_caoh2(char *buffer, short i, struct UnitProcess *unit);
int read_na2co3(char *buffer, short i, struct UnitProcess *unit);
int read_nh3(char *buffer, short i, struct UnitProcess *unit);
int read_kmno4(char *buffer, short i, struct UnitProcess *unit);
int read_co2(char *buffer, short i, struct UnitProcess *unit);
int read_o3(char *buffer, short i, struct UnitProcess *unit);
int read_clo2(char *buffer, short i, struct UnitProcess *unit);
int read_naocl(char *buffer, short i, struct UnitProcess *unit);
int read_so2(char *buffer, short i, struct UnitProcess *unit);

int write_influent(char *buffer, short i, struct UnitProcess *unit);
int write_alum(char *buffer, short i, struct UnitProcess *unit);
int write_gac(char *buffer, short i, struct UnitProcess *unit);
int write_filter(char *buffer, short i, struct UnitProcess *unit);
int write_basin(char *buffer, short i, struct UnitProcess *unit);
int write_mfuf(char *buffer, short i, struct UnitProcess *unit);
int write_nf(char *buffer, short i, struct UnitProcess *unit);
int write_ssf(char *buffer, short i, struct UnitProcess *unit);
int write_def(char *buffer, short i, struct UnitProcess *unit);
int write_altf(char *buffer, short i, struct UnitProcess *unit);
int write_bankf(char *buffer, short i, struct UnitProcess *unit);
int write_presed(char *buffer, short i, struct UnitProcess *unit);
int write_uvdis(char *buffer, short i, struct UnitProcess *unit);
int write_iron(char *buffer, short i, struct UnitProcess *unit);
int write_avg_tap(char *buffer, short i, struct UnitProcess *unit);
int write_end_sys(char *buffer, short i, struct UnitProcess *unit);
int write_chlorine(char *buffer, short i, struct UnitProcess *unit);
int write_h2so4(char *buffer, short i, struct UnitProcess *unit);
int write_naoh(char *buffer, short i, struct UnitProcess *unit);
int write_caoh2(char *buffer, short i, struct UnitProcess *unit);
int write_na2co3(char *buffer, short i, struct UnitProcess *unit);
int write_nh3(char *buffer, short i, struct UnitProcess *unit);
int write_kmno4(char *buffer, short i, struct UnitProcess *unit);
int write_co2(char *buffer, short i, struct UnitProcess *unit);
int write_o3(char *buffer, short i, struct UnitProcess *unit);
int write_clo2(char *buffer, short i, struct UnitProcess *unit);
int write_naocl(char *buffer, short i, struct UnitProcess *unit);
int write_so2(char *buffer, short i, struct UnitProcess *unit);

/* Supporting I/O functions */
int ftab(FILE *fp, register short tab);
int read_line(FILE *fp, char **id, char **value);
int fdoublef(FILE *fp, const char *fmt, double x);
int strip_buffer(char buffer[], FILE *fout);

/****************  Model functions **************************/

int runmodel(struct ProcessTrain *train);

/* Water Chemistry functions */
double Kw(double DegK);        /* Ionization of water.               */
double K_HCO3(double DegK);    /* 1st ionization of carbonic acid.   */
double K_CO3(double DegK);     /* 2nd ionization of carbonic acid.   */
double K_HOCl(double DegK);    /* Ionization of free chlorine        */
double K_NH3(double DegK);     /* Ionization of ammonia              */
double K_MgOH(double DegK);    /* Ionization of magnesium hydroxide. */
double K_MgOH2(double DegK);   /* Solubility of magnesium.           */
double K_MgOH2aq(double DegK); /* Solubility of magnesium hydorxide. */
double K_CaOH(double DegK);    /* Ionization of Calcium hydroxide.   */
double K_CaOH2aq(double DegK); /* Solubility of Calcium hydroxide.   */
double K_CaCO3(double DegK);   /* Solubility of Calcium carbonate.   */

void phchange(struct UnitProcess *unit, short flag);
void cl2decay(struct UnitProcess *unit, double rxnhours);
void breakpt(struct UnitProcess *unit);

short ncstr(double t_ten_theo);
void ct(struct UnitProcess *unit, double t_ten);

/* Chemical addition: located in Chemical.c */
void chem_reset(struct ProcessTrain *train);
void alumadd(struct UnitProcess *unit);
void fecladd(struct UnitProcess *unit);
void nh3add(struct UnitProcess *unit);
void naohadd(struct UnitProcess *unit);
void limeadd(struct UnitProcess *unit);
void h2so4add(struct UnitProcess *unit);
void kmno4add(struct UnitProcess *unit);
void chloradd(struct UnitProcess *unit);
void clo2add(struct UnitProcess *unit);
void co2_add(struct UnitProcess *unit);
void soda_add(struct UnitProcess *unit);
void so2_add(struct UnitProcess *unit);
void o3add(struct UnitProcess *unit);

/* Unit Process Removal/Production functions: */
void influent(struct UnitProcess *unit);
void alum_rmv(struct UnitProcess *unit);
void fecl_rmv(struct UnitProcess *unit);
void soft_rmv(struct UnitProcess *unit);
void solids_rmv(struct UnitProcess *unit);
void gac_rmv(struct UnitProcess *unit);

void basn_dbp(struct UnitProcess *unit);
void filt_dbp(struct UnitProcess *unit);

/* D/DBP formation equations */
void dist_dbp(struct UnitProcess *unit);

/* New subroutines added by WJS, 10/98 */
void choose_dbpmodel(struct UnitProcess *unit);
void modrw1dbp(struct UnitProcess *unit);
void modrw2dbp(struct UnitProcess *unit);
void rwdbp(struct UnitProcess *unit);
void owdbp(struct UnitProcess *unit);
void coagdbp(struct UnitProcess *unit);
void gacmbdbp(struct UnitProcess *unit);
void biofilt_rmv(struct UnitProcess *unit);
//void biogac_rmv     ( struct UnitProcess *unit );
void ozonate(struct UnitProcess *unit);
void res_time(struct UnitProcess *unit);
void clo2_decay(struct UnitProcess *unit, double cfstr_time);
void nf_rmv(struct UnitProcess *unit);
void mfuf_rmv(struct UnitProcess *unit);
void new_cl2decay(struct UnitProcess *unit, double res_time);
void ec_comply(struct UnitProcess *unit);

/******************* Structure functions ****************************/

struct ProcessTrain *AllocProcessTrain(void);
struct ProcessTrain *FreeProcessTrain(struct ProcessTrain *train);
struct ProcessTrain *InitProcessTrain(struct ProcessTrain *train);

struct UnitProcess *AllocUnitProcess(short type);
struct UnitProcess *FreeUnitProcess(struct UnitProcess *unit);
struct UnitProcess *AddUnitProcess(struct ProcessTrain *train, short type);
struct UnitProcess *RemoveUnitProcess(struct UnitProcess *unit);
struct UnitProcess *MoveUnitProcess(struct UnitProcess *prev,
                                    struct UnitProcess *unit);
struct UnitProcess *FirstUnitProcess(struct ProcessTrain *train);
struct UnitProcess *NextUnitProcess(struct UnitProcess *unit);
struct UnitProcess *LastUnitProcess(struct ProcessTrain *train);
struct UnitProcess *PrevUnitProcess(struct UnitProcess *unit);
struct UnitProcess *GetUnitProcess(struct ProcessTrain *train, int n);

int GetUnitIndex(struct ProcessTrain *train, struct UnitProcess *unit);
int unit_type(char *name);

/* Adjustment functions based on fitting ICR data */
double f_adj_pH(double pH);
double f_adj_toc(double toc, struct UnitProcess *unit);

#endif
