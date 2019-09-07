/* Globals.c -- September 16, 1993
*    Variables for Water Treatment Plant (WTP) model.
*
*  Michael D. Cummins
*     July 1993
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*          Project modifications done by:
*                 Kenneth B. Curry
*                 October 12, 1993

*                Sample Location Points
	       and Volume Calculations for
	       Theoretical Residence times

	       More modifications done by:
	       Warren J. Swanson, MPI-PHO, 10/98
*/
#include "wtp.h"

/* Default values UnitProcess design and operating data: */

struct Influent default_influent = {
    /* Default Influent Data */
    /* pH         */ 8.0,
    /* temp       */ 20.0,
    /* low_temp   */ 5.0,
    /* toc        */ 3.0,
    /* uv254      */ 0.06,
    /* bromide    */ 0.05,
    /* alkalinity */ 100.0,
    /* calcium    */ 100.0,
    /* hardness   */ 120.0,
    /* nh3        */ 0.01,
    /* ntu        */ 5.0,
    //* crypto_req */   3.0,
    //* clo2_crypto_ct_mult */ 7.5,
    /* peak       */ 5.0,
    /* flow       */ 2.0,
    /* swflag     */ TRUE,
    /* crypto_conc*/ 0.000,
    /* lt2_wscp_flag */ FALSE,
    /* gw_virus_flag */ FALSE,
    /* gw_tot_dis_req_v */ 4.0
    //  /* run_name   */  "Run1"
};

struct Alum default_alum = {
    /* dose */ 10};

struct Gac default_gac = {
    /* ebct     */ 15.0,
    /* regen    */ 180.0,
    /* config   */ 'S',
    /* toc_calc */ 'A',
    /* crypto_lr_2*/ 0.5,
    /* filt_stage*/ 0};

struct Filter default_filter = {
    /* detent   */ /*   15.0, */
    /* fi_mean  */ 1.0,
    /* fi_t_ten */ 0.5,
    /* volume   */ 0.028,
    /* cl2_bkwsh*/ TRUE,
    /* media    */ 'S',
    /* giardia_lr_conv */ 2.5,
    /* virus_lr_conv */ 2.0,
    /* crypto_lr_conv */ 3.0,
    /* giardia_lr_df */ 2.0,
    /* virus_lr_df */ 1.0,
    /* crypto_lr_df */ 3.0,
    /* cfe_turb_flag */ FALSE,
    /* ife_turb_flag */ FALSE,
    /* crypto_lr_2 */ 0.5,
    /* filt_stage*/ 0};

struct Basin default_basin = {
    /* detent   */ /*  90.0, */
    /* sb_mean  */ 1.0,
    /* sb_t_ten */ 0.3,
    /* Volume   */ 1.0};

struct Basin default_o3_contactor = {
    /* detent   */ /*   90.0, */
    /* sb_mean  */ 1.0,
    /* sb_t_ten */ 0.5,
    /* Volume   */ 0.007};

//struct Membrane default_membrane = {
//  /* mwc            */  500.0,
//  /* recover        */   90.0,
//  /* giardia_lr     */    2.0,
//  /* virus_lr       */    2.0,
//  /* crypto_lr      */    2.0
//  /* op_psi         *//*160.0*/
//  };

struct Mfuf default_mfuf = {
    /* recover     */ 90.0,
    /* giardia_lr  */ 2.5,
    /* virus_lr    */ 1.0,
    /* crypto_lr_1 */ 3.0,
    /* crypto_lr_2 */ 2.5,
    /* filt_stage*/ 0};

struct Nf default_nf = {
    /* mwc     */ 500.0,
    /* recover */ 90.0,
    /* giardia_lr  */ 2.0,
    /* virus_lr    */ 2.0,
    /* crypto_lr   */ 2.0,
    /* toc_rem     */ 90.0,
    /* uva_rem     */ 90.0,
    /* br_rem      */ 50.0,
    /* treat_fraction*/ 1.0,
    /* filt_stage  */ 0};

struct Ssf default_ssf = {
    /* volume   */ 2.80,
    /* fi_mean  */ 1.0,
    /* fi_t_ten */ 0.5,
    /* toc_rem_lowt*/ 10.0,
    /* toc_rem_midt*/ 15.0,
    /* toc_rem_hight*/ 20.0,
    /* giardia_lr  */ 2.0,
    /* virus_lr    */ 2.0,
    /* crypto_lr_1*/ 3.0,
    /* crypto_lr_2*/ 2.5,
    /* filt_stage*/ 0};

struct Def default_def = {
    /* volume   */ 0.028,
    /* fi_mean  */ 1.0,
    /* fi_t_ten */ 0.5,
    /* giardia_lr */ 2.0,
    /* virus_lr    */ 1.0,
    /* crypto_lr_1*/ 3.0,
    /* crypto_lr_2*/ 0.5,
    /* filt_stage*/ 0};

struct Altf default_bagf = {
    /* giardia_lr  */ 2.0,
    /* virus_lr    */ 0.0,
    /* crypto_lr_1 */ 3.0,
    /* crypto_lr_2 */ 1.0,
    /* filt_stage*/ 0};

struct Altf default_cartf = {
    /* giardia_lr  */ 2.0,
    /* virus_lr    */ 0.0,
    /* crypto_lr_1 */ 3.0,
    /* crypto_lr_2 */ 2.0,
    /* filt_stage*/ 0};

struct Bankf default_bankf = {
    /* eligible_lt2 */ TRUE,
    /* distance */ 30.0,
    /* crypto_lr_close */ 0.5,
    /* crypto_lr_far */ 1.0};

struct Presed default_presed = {
    /* volume */ 0.167,
    /* sb_mean */ 1.0,
    /* sb_t_ten */ 0.3,
    /* eligible_lt2 */ FALSE,
    /* crypto_lr */ 0.5};

struct Uvdis default_uvdis = {
    /* giardia_li */ 3.0,
    /* virus_li */ 0.0,
    /* crypto_li */ 3.0};

struct Iron default_iron = {
    /* dose */ 10.0};

struct chemical default_chemical = {
    /* chlor */ 2.0,
    /* ozone */ 1.0,
    /* clo2  */ /* 1.0, */
    /* nh3   */ 0.4,
    /* kmno4 */ 1.0,
    /* naoh  */ 1.0,
    /* h2so4 */ 1.0,
    /* lime  */ /* 10.0, */
    /* soda  */ 1.0,
    /* co2   */ 1.0,
    /* naocl */ 2.0,
    /* so2   */ 1.0};

struct clo2 default_clo2 = {
    /* dose */ 1.0,
    /* conversion */ 70.0};

struct lime default_lime = {
    /* dose */ 10.0,
    /* purpose */ 'P'};

/* No longer needed - WJS 11/98 
struct WTP_effluent  default_wtp_effluent = { 
   hha_flag   'T' 
  }; */

struct First_tap default_first_tap = {
    /* minute */ 10.0};

struct Avg_tap default_avg_tap = {
    /* days */ 1.0};

struct End_of_system default_end_of_system = {
    /* days */ 3.0};

struct Basin default_rapid_mix = {
    /* detent   */ /*  5.0,  */
    /* sb_mean  */ 1.0,
    /* sb_t_ten */ 0.1,
    /* Volume   */ 0.007};

struct Basin default_slow_mix = {
    /* detent   */ /*   20.0,  */
    /* sb_mean  */ 1.0,
    /* sb_t_ten */ 0.1,
    /* Volume   */ 0.04};

struct Basin default_settling_basin = {
    /* detent   */ /*  90.0,   */
    /* sb_mean  */ 1.0,
    /* sb_t_ten */ 0.3,
    /* Volume   */ 0.167};

struct Basin default_contact_tank = {
    /* detent   */ /*  30.0,  */
    /* sb_mean  */ 1.0,
    /* sb_t_ten */ 0.5,
    /* Volume   */ 1.0};

struct Basin default_clearwell = {
    /* detent   */ /*   60.0,   */
    /* sb_mean  */ 1.0,
    /* sb_t_ten */ 0.3,
    /* Volume   */ 1.0};

struct Avg_tap default_location_1 = {
    /* days */ 1.0};

/* Global flags */
int coldflag = FALSE; /* TRUE=Run model at cold temperature and peak flow */
int coagflag = FALSE;
int conv_filtflag = FALSE;
int filtflag = FALSE;
int filt2flag = FALSE;
int swflag = TRUE;
int gw_virus_flag = FALSE;

/* More global flags added 10/98 by WJS */
int softflag = FALSE;
int soft2flag = FALSE;
int floccflag = FALSE;
int sedflag = FALSE;
int sedconvflag = FALSE;
int gacflag = FALSE;
int mfufflag = FALSE;
int nfflag = FALSE;
int ssfflag = FALSE;
int defflag = FALSE;
int bagfflag = FALSE;
int cartfflag = FALSE;
int bankfflag = 0;
int granf2flag = FALSE;
int gac2flag = FALSE;
int mfuf2flag = FALSE;
int nf2flag = FALSE;
int ssf2flag = FALSE;
int def2flag = FALSE;
int bagf2flag = FALSE;
int cartf2flag = FALSE;
int lt2presedflag = FALSE;
int cfeflag = FALSE;
int ifeflag = FALSE;
int uvflag = FALSE;
int lt2_wscp_flag = FALSE;
int clo2flag = FALSE;
int o3flag = FALSE;
int pre_o3flag = FALSE;
int int_o3flag = FALSE;
int post_o3flag = FALSE;
int dir_filtflag = FALSE;
int bio_filtflag = FALSE;
//int   bio_gacflag   = FALSE;

int rwdbpflag = TRUE;
int owdbpflag = FALSE;
int coagdbpflag = FALSE;
int modrw1dbpflag = FALSE;
int modrw2dbpflag = FALSE;
int gacmemdbpflag = FALSE;

double nonconv_discredit = 0.0;
double bin34_inactreqd = 0.0;

/* Global time counter to be used for modrw1dbp() and modrw2dpb()
   subroutines */
double modrw2dbptime = 0;
double modrw2dbpcl2 = 0;

/* Globals added for pathogen removal by filtration/membranes */

//double fi_crypto_lr  = 0.0;
//double fi_giardia_lr = 0.0;
//double fi_virus_lr   = 0.0;

double tot_dis_req_c = 0.0;
double tot_dis_req_g = 0.0;
double tot_dis_req_v = 0.0;

double tot_crypto_lr = 0.0;
double tot_giardia_lr = 0.0;
double tot_virus_lr = 0.0;