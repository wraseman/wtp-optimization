/* UPTable.c  -- DOS file name for UnitProcessTable[] -- Sept 16, 1993
*
* UnitProcessTable[] is an array that contains information describing
* all the unit processes and data elements known by WTP model.  It is
* intended to be used in input and output of data files, listing
* process train data, and input via data entry screens.  The data
* element of UnitProcessTable[] is a structure called UnitInfo. 
*
* UnitInfo contains information for one unit process.  The information
* includes the name of the unit process, a pointer to two function that
* will input and output each unit process parameter from/to a string,
* and a pointer to a DataID array.
*
* The DataID array contains information about the operating parameters
* for the unit process.  DataID includes the parameter ID used in data
* files, full name of parameter as it appears on the data entry screen,
* units of measurement as appears on the data entry screen, display
* format that specifies number of decimal places, and finally the data
* element type (double, integer, character, or string).
*
* Both UnitProcessTable[] and DataID[] arrays are NULL terminated and can
* be index as follows.
*
* void main(void)
* {
*   struct UnitInfo    *info;
*   struct DataID      *data;
*   int                 i;
*
*   for( info=UnitProcessTable; info->name; info++ )
*     {
*       printf("%s\n", info->name );
*       for( data=info->data,i=0; data->id; data++,i++ )
*         {
*           printf("  %2d %s\n", i,data->id);
*         }
*     }
* }
*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*      Provide for additonal location sample points  Using Location_# type
**/

#include "wtp.h"

struct DataID InfluentID[] = {
    /* id  ****  name  *********************  units  ***********  fmtout  fmtin*/
    {
        "pH",
        "pH",
        "",
        "%5.1lf",
        "%lf",
    },
    {
        "Temp",
        "Influent Temperature",
        "(Celsius)",
        "%5.1lf",
        "%lf",
    },
    {
        "MinTemp",
        "Minimum Temperature",
        "(Celsius)",
        "%5.1lf",
        "%lf",
    },
    {
        "TOC",
        "Total Organic Carbon",
        "(mg/L)",
        "%5.1lf",
        "%lf",
    },
    {
        "UVA",
        "UV Absorbance at 254nm",
        "(1/cm)",
        "%7.3lf",
        "%lf",
    },
    {
        "Br",
        "Bromide",
        "(mg/L)",
        "%7.3lf",
        "%lf",
    },
    {
        "Alk",
        "Alkalinity",
        "(mg/L as CaCO3)",
        "%3.0lf",
        "%lf",
    },
    {
        "Ca",
        "Calcium Hardness",
        "(mg/L as CaCO3)",
        "%3.0lf",
        "%lf",
    },
    {
        "Hard",
        "Total Hardness",
        "(mg/L as CaCO3)",
        "%3.0lf",
        "%lf",
    },
    {
        "NH3",
        "Ammonia",
        "(mg/L as N)",
        "%6.2lf",
        "%lf",
    },
    {
        "NTU",
        "Turbidity",
        "(NTU)",
        "%5.1lf",
        "%lf",
    },
    //{"Crypto"  ,"Cryptosporidium Removal+Inact. Required","(logs)","%5.1lf","%lf",},
    //{"CrypClO2","Multiplier for Crypto. CT by ClO2",""            ,"%5.1lf","%lf",},
    {
        "PeakFlow",
        "Peak Flow",
        "(MGD)",
        "%5.3lf",
        "%lf",
    },
    {
        "Flow",
        "Plant Flow",
        "(MGD)",
        "%5.3lf",
        "%lf",
    },
    {
        "SWFlag",
        "Surface Water by SWTR",
        "(TRUE/FALSE)",
        "%s",
        "%s",
    },
    {
        "CryptoIn",
        "Source Water Crypto. Concentration",
        "(oocysts/Liter)",
        "%6.3lf",
        "%lf",
    },
    {
        "LT2WSCP",
        "LT2 Rule Watershed Control Prog. Credit?",
        "(TRUE/FALSE)",
        "%s",
        "%s",
    },
    {
        "GWvirus?",
        "If GW System, Is Virus Disinfection Req'd?",
        "(TRUE/FALSE)",
        "%s",
        "%s",
    },
    {
        "GWvLogs",
        "Virus Disinfection for GW, if Req'd",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    //{"RunName" ,"Run Description"            ,""                  ,"%c"    ,"c", },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID AlumID[] = {
    {
        "Dose",
        "Alum Dose",
        "(mg/L as Al2(SO4)3*14H2O)",
        "%5.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID GacID[] = {
    {
        "EBCT",
        "Empty Bed Contact Time (at 'Plant Flow')",
        "(minutes)",
        "%3.0lf",
        "%lf",
    },
    {
        "React_Fq",
        "GAC Reactivation Interval",
        "(days)",
        "%3.0lf",
        "%lf",
    },
    {
        "Config",
        "GAC Contacting System (Single/Blended)",
        "(S or B)",
        "%c",
        "c",
    },
    {
        "Eff_TOC",
        "TOC Breakthrough for Single Unit (Max/Avg)",
        "(M or A)",
        "%c",
        "c",
    },
    {
        "LRcrypto2",
        "Crypto. Removal Credit as 2nd Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID FilterID[] = {
    {
        "Vol(liq)",
        "Liquid Volume",
        "(MG)",
        "%7.4lf",
        "%lf",
    },
    {
        "T50",
        "Ratio of T50/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {
        "T10",
        "Ratio of T10/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {
        "Cl2BKWSH",
        "Chlorinated Backwash Water?",
        "(TRUE/FALSE)",
        "%s",
        "%s",
    },
    {
        "Media",
        "Filter Media (Anthracite/Sand or GAC)",
        "(S or G)",
        "%c",
        "%c",
    },
    {
        "LR_G_Conv",
        "Giardia Removal Credit - Conv. Filters",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LR_V_Conv",
        "Virus Removal Credit - Conv. Filters",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LR_C_Conv",
        "Crypto. Removal Credit - Conv. Filters",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LR_G_Df",
        "Giardia Removal Credit - Direct Filters",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LR_V_Df",
        "Virus Removal Credit - Direct Filters",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LR_C_Df",
        "Crypto. Removal Credit - Direct Filters",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "CFE_flag",
        "CFE Turb. Meets LT2 Toolbox Criteria?",
        "(TRUE/FALSE)",
        "%s",
        "%s",
    },
    {
        "IFE_flag",
        "IFE Turb. Meets LT2 Toolbox Criteria?",
        "(TRUE/FALSE)",
        "%s",
        "%s",
    },
    {
        "LRcrypto2",
        "Crypto. Credit as 2nd Stage Filt.",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID BasinID[] = {
    {
        "Vol.",
        "Volume of Basin",
        "(MG)",
        "%7.4lf",
        "%lf",
    },
    {
        "T50",
        "Ratio of T50/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {
        "T10",
        "Ratio of T10/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID NfID[] = {
    {
        "MWC",
        "Molecular Weight Cutoff",
        "(gm/mole)",
        "%5.0lf",
        "%lf",
    },
    {
        "Recovery",
        "Percent Recovery",
        "(%)",
        "%6.0lf",
        "%lf",
    },
    {
        "LRgiardia",
        "Giardia Removal Credit as 1st Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LRvirus",
        "Virus Removal Credit as 1st Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LRcrypto",
        "Crypto. Removal Credit",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "TOCRem",
        "TOC Removal by NF",
        "(%)",
        "%6.0lf",
        "%lf",
    },
    {
        "UVARem",
        "UVA Removal by NF (< or = TOC Rem.)",
        "(%)",
        "%6.0lf",
        "%lf",
    },
    {
        "BrRem",
        "Bromide Removal by NF",
        "(%)",
        "%6.0lf",
        "%lf",
    },
    {
        "Trt_Frac",
        "Fraction of Flow Treated by NF",
        "(--)",
        "%6.2lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL}};

struct DataID MfUfID[] = {
    {
        "Recovery",
        "Percent Recovery",
        "(%)",
        "%6.0lf",
        "%lf",
    },
    {
        "LRgiardia",
        "Giardia Removal Credit as 1st Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LRvirus",
        "Virus Removal Credit as 1st Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LRcrypto1",
        "Crypto. Removal Credit as 1st Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LRcrypto2",
        "Crypto. Removal Credit as 2nd Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL}};

struct DataID IronID[] = {
    {
        "Dose",
        "Ferric Dose",
        "(mg/L as FeCl3*6H2O)",
        "%6.2lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID ChlorineID[] = {
    {
        "Dose",
        "Chlorine Dose",
        "(mg/L as Cl2)",
        "%5.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID OzoneID[] = {
    {
        "Dose",
        "Ozone Dose",
        "(mg/L as O3)",
        "%5.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

/* Added 11/98, WJS */
struct DataID O3CntctrID[] = {
    {
        "Vol.",
        "Volume of Basin",
        "(MG)",
        "%7.4lf",
        "%lf",
    },
    {
        "T50",
        "Ratio of T50/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {
        "T10",
        "Ratio of T10/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID SSFilterID[] = {
    {
        "Vol(liq)",
        "Liquid Volume",
        "(MG)",
        "%7.4lf",
        "%lf",
    },
    {
        "T50",
        "Ratio of T50/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {
        "T10",
        "Ratio of T10/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {
        "TOCrmLowT",
        "TOC Removal: Temp. < 10 C",
        "(%)",
        "%6.0lf",
        "%lf",
    },
    {
        "TOCrmMidT",
        "TOC Removal: 10 C <= Temp. < 20 C",
        "(%)",
        "%6.0lf",
        "%lf",
    },
    {
        "TOCrmHiT",
        "TOC Removal: 20 C <= Temp.",
        "(%)",
        "%6.0lf",
        "%lf",
    },
    {
        "LRgiardia",
        "Giardia Removal Credit as 1st Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LRvirus",
        "Virus Removal Credit as 1st Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LRcrypto1",
        "Crypto. Removal Credit as 1st Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LRcrypto2",
        "Crypto. Removal Credit as 2nd Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID DEFilterID[] = {
    {
        "Vol(liq)",
        "Liquid Volume",
        "(MG)",
        "%7.4lf",
        "%lf",
    },
    {
        "T50",
        "Ratio of T50/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {
        "T10",
        "Ratio of T10/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {
        "LRgiardia",
        "Giardia Removal Credit as 1st Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LRvirus",
        "Virus Removal Credit as 1st Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LRcrypto1",
        "Crypto. Removal Credit as 1st Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LRcrypto2",
        "Crypto. Removal Credit as 2nd Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID BagFilterID[] = {
    {
        "LRgiardia",
        "Giardia Removal Credit as 1st Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LRvirus",
        "Virus Removal Credit as 1st Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LRcrypto1",
        "Crypto. Removal Credit as 1st Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LRcrypto2",
        "Crypto. Removal Credit as 2nd Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID CartFilterID[] = {

    {
        "LRgiardia",
        "Giardia Removal Credit as 1st Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LRvirus",
        "Virus Removal Credit as 1st Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LRcrypto1",
        "Crypto. Removal Credit as 1st Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LRcrypto2",
        "Crypto. Removal Credit as 2nd Stage",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID BankFilterID[] = {
    {
        "Eligible",
        "Eligible for LT2 Toolbox Crypto. Credit?",
        "(TRUE/FALSE)",
        "%s",
        "%s",
    },
    {
        "SW_Dist",
        "Distance from Wells to Surface Water",
        "(feet)",
        "%6.lf",
        "%lf",
    },
    {
        "LRcrypto1",
        "Crypto. Removal Credit: 25 to 49.9 ft",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LRcrypto2",
        "Crypto. Removal Credit: 50.0+ ft",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID PreSedID[] = {
    {
        "Vol.",
        "Volume of Basin",
        "(MG)",
        "%7.4lf",
        "%lf",
    },
    {
        "T50",
        "Ratio of T50/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {
        "T10",
        "Ratio of T10/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {
        "Eligible",
        "Eligible for LT2 Toolbox Crypto. Credit?",
        "(TRUE/FALSE)",
        "%s",
        "%s",
    },
    {
        "LRcrypto",
        "LT2 Toolbox Crypto. Removal Credit",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID UVDisID[] = {
    {
        "LIgiardia",
        "Giardia Inactivation Credit",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LIvirus",
        "Virus Inactivation Credit",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {
        "LIcrypto",
        "Crypto. Inactivation Credit",
        "(logs)",
        "%6.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID ClO2ID[] = {
    {
        "Dose",
        "Chlorine Dioxide Dose",
        "(mg/L as ClO2)",
        "%5.1lf",
        "%lf",
    },
    {
        "ClO2pct",
        "Conversion to Chlorite",
        "(%)",
        "%5.0lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID AmmoniaID[] = {
    {
        "Dose",
        "Ammonia Dose",
        "(mg/L as N)",
        "%5.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID NH42SO4ID[] = {
    {
        "Dose",
        "Ammonium Sulfate Dose",
        "(mg/L as N)",
        "%5.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID KMnO4ID[] = {
    {
        "Dose",
        "Permanganate Dose",
        "(mg/L as KMnO4)",
        "%5.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID NaOHID[] = {
    {
        "Dose",
        "Sodium Hydroxide Dose",
        "(mg/L as NaOH)",
        "%5.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID H2SO4ID[] = {
    {
        "Dose",
        "Sulfuric Acid Dose",
        "(mg/L as H2SO4)",
        "%5.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID CaOH2ID[] = {
    {
        "Dose",
        "Lime Dose",
        "(mg/L as Ca(OH)2)",
        "%5.1lf",
        "%lf",
    },
    {
        "Purpose",
        "For pH adjustment (P) or Softening (S)",
        "(P or S)",
        "%c",
        "%c",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID Na2CO3ID[] = {
    {
        "Dose",
        "Soda Ash Dose",
        "(mg/L as Na2CO3)",
        "%6.2lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID CO2ID[] = {
    {
        "Dose",
        "Carbon Dioxide Dose",
        "(mg/L as CO2)",
        "%5.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID NaOClID[] = {
    {
        "Dose",
        "Sodium Hypochlorite Dose",
        "(mg/L as Cl2)",
        "%5.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID SO2ID[] = {
    {
        "Dose",
        "Sulfur Dioxide Dose",
        "(mg/L as SO2)",
        "%5.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID AvgTapID[] = {
    {
        "DetTime",
        "Average Residence Time (For Average Flow)",
        "(Days)",
        "%5.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID EndSysID[] = {
    {
        "DetTime",
        "Maximum Residence Time (For Average Flow)",
        "(Days)",
        "%5.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID RapidMixID[] = {
    {
        "Vol.",
        "Volume of Basin",
        "(MG)",
        "%7.4lf",
        "%lf",
    },
    {
        "T50",
        "Ratio of T50/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {
        "T10",
        "Ratio of T10/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID SlowMixID[] = {
    {
        "Vol.",
        "Volume of Basin",
        "(MG)",
        "%7.4lf",
        "%lf",
    },
    {
        "T50",
        "Ratio of T50/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {
        "T10",
        "Ratio of T10/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID SettlingID[] = {
    {
        "Vol.",
        "Volume of Basin",
        "(MG)",
        "%7.4lf",
        "%lf",
    },
    {
        "T50",
        "Ratio of T50/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {
        "T10",
        "Ratio of T10/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID ContactID[] = {
    {
        "Vol.",
        "Volume of Basin",
        "(MG)",
        "%7.4lf",
        "%lf",
    },
    {
        "T50",
        "Ratio of T50/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {
        "T10",
        "Ratio of T10/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID ClearwellID[] = {
    {
        "Vol.",
        "Volume of Basin",
        "(MG)",
        "%7.4lf",
        "%lf",
    },
    {
        "T50",
        "Ratio of T50/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {
        "T10",
        "Ratio of T10/Detention Time",
        "(ratio)",
        "%6.2lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID Location_1ID[] = {
    {
        "DetTime",
        "Residence Time (For Plant Flow)",
        "(Days)",
        "%5.1lf",
        "%lf",
    },
    {NULL, NULL, NULL, NULL, NULL}};

struct DataID VacantID[] = {{NULL, NULL, NULL, NULL, NULL}};

/**************  UnitProcessTable  *******************************/
struct UnitInfo UnitProcessTable[] = {
    /* Index     Name              DataIDTable File input     File output   */
    /*   0  */ {
        "Vacant",
        VacantID,
        NULL,
        NULL,
    },
    /*   1  */ {
        "Influent",
        InfluentID,
        read_influent,
        write_influent,
    },
    /*   2  */ {
        "Rapid Mix",
        RapidMixID,
        read_basin,
        write_basin,
    },
    /*   3  */ {
        "Flocculation",
        SlowMixID,
        read_basin,
        write_basin,
    },
    /*   4  */ {
        "Settling Basin",
        SettlingID,
        read_basin,
        write_basin,
    },
    /*   5  */ {
        "Filtration",
        FilterID,
        read_filter,
        write_filter,
    },
    /*   6  */ {
        "Basin",
        BasinID,
        read_basin,
        write_basin,
    },
    /*   7  */ {
        "Contact Tank",
        ContactID,
        read_basin,
        write_basin,
    },
    /*   8  */ {
        "Reservoir",
        ClearwellID,
        read_basin,
        write_basin,
    },
    /*   9  */ {
        "Ozone Chamber",
        O3CntctrID,
        read_basin,
        write_basin,
    },
    /*  10  */ {
        "GAC",
        GacID,
        read_gac,
        write_gac,
    },
    /*  11  */ {
        "MF/UF",
        MfUfID,
        read_mfuf,
        write_mfuf,
    },
    /*  12  */ {
        "Nanofiltration",
        NfID,
        read_nf,
        write_nf,
    },
    /*  13  */ {
        "Slow Sand Filtration",
        SSFilterID,
        read_ssf,
        write_ssf,
    },
    /*  14  */ {
        "D.E. Filtration",
        DEFilterID,
        read_def,
        write_def,
    },
    /*  15  */ {
        "Bag Filtration",
        BagFilterID,
        read_altf,
        write_altf,
    },
    /*  16  */ {
        "Cartridge Filtration",
        CartFilterID,
        read_altf,
        write_altf,
    },
    /*  17  */ {
        "Bank Filtration",
        BankFilterID,
        read_bankf,
        write_bankf,
    },
    /*  18  */ {
        "Presed. Basin",
        PreSedID,
        read_presed,
        write_presed,
    },
    /*  19  */ {
        "UV Disinfection",
        UVDisID,
        read_uvdis,
        write_uvdis,
    },
    /*  20  */ {
        "Alum",
        AlumID,
        read_alum,
        write_alum,
    },
    /*  21  */ {
        "Iron",
        IronID,
        read_iron,
        write_iron,
    },
    /*  22  */ {
        "Chlorine (Gas)",
        ChlorineID,
        read_chlorine,
        write_chlorine,
    },
    /*  23  */ {
        "Sulfuric Acid",
        H2SO4ID,
        read_h2so4,
        write_h2so4,
    },
    /*  24  */ {
        "Lime",
        CaOH2ID,
        read_caoh2,
        write_caoh2,
    },
    /*  25  */ {
        "Soda Ash",
        Na2CO3ID,
        read_na2co3,
        write_na2co3,
    },
    /*  26  */ {
        "Ammonia",
        AmmoniaID,
        read_nh3,
        write_nh3,
    },
    /*  27  */ {
        "Permanganate",
        KMnO4ID,
        read_kmno4,
        write_kmno4,
    },
    /*  28  */ {
        "Carbon Dioxide",
        CO2ID,
        read_co2,
        write_co2,
    },
    /*  29  */ {
        "Ozone",
        OzoneID,
        read_o3,
        write_o3,
    },
    /*  30  */ {
        "Chlorine Dioxide",
        ClO2ID,
        read_clo2,
        write_clo2,
    },
    /*  31  */ {
        "Sod. Hydroxide",
        NaOHID,
        read_naoh,
        write_naoh,
    },
    /*  32  */ {
        "Sod. Hypochlor.",
        NaOClID,
        read_naocl,
        write_naocl,
    },
    /*  33  */ {
        "Sulfur Dioxide",
        SO2ID,
        read_so2,
        write_so2,
    },
    /*  34  */ {
        "Ammon Sulfate",
        NH42SO4ID,
        read_nh3,
        write_nh3,
    },
    /*  35  */ {
        "WTP Effluent",
        VacantID,
        NULL,
        /*read_effluent,*/ NULL /* write_effluent*/,
    },
    /*  36  */ {
        "Average Tap",
        AvgTapID,
        read_avg_tap,
        write_avg_tap,
    },
    /*  37  */ {
        "End of System",
        EndSysID,
        read_end_sys,
        write_end_sys,
    },
    /*  38  */ {
        "Additional Point",
        Location_1ID,
        read_avg_tap,
        write_avg_tap,
    },
    /* End   */ {NULL, NULL, NULL, NULL}};
