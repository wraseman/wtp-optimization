/* pHChange.c -- September 2, 1993 */
#include "wtp.h"

void phchange(struct UnitProcess *unit, short flag)
/*
*  Purpose:
*   1. Adjust pH or CBminusCA such that charge balance is obtained.
*   2. Precipate CaCO3 and Mg(OH)2 for hardness removal.
*   3. Adjust alkalinity based on changes in pH and carbonate concentrations.
*
*  Input:
*    flag=TRUE  adjust pH
*        =FALSE adjust CBminusCA
*
*  Method:
*    The electrical charge of the postive cations must equal the
*    electrical charge of the negative anions.  Version 1.3 of WTP tracks
*    the following ions.
*
*  Cations
*   [H+]
*   [Ca++]
*   [CaOH+]
*   [Mg++]
*   [MgOH+]
*   [NH4+]
*
*  Anions
*   [OH-]
*   [HCO3-]
*   [CO3--]
*   [OCl-]
*
*  CBminusCA is the charge due to all other ions.
*
*  Notes:
*   1. All other ions are tracket by CBmunusCa, cations adding to CBminusCA
*      and anions subtracting from CBminusCA.
*   2. Input and returned values are passed using the UnitProcess data
*      structure.
*   3. CBminusCA is initialized using the same algorithm as in version 1.2.
*      Specifically, CaCO3 and Mg(OH)2 are NOT percipated before adjusting
*      CBminusCA.  I am not sure this logic is valid.
*   4. This routine checks for changes in the inputs before starting the
*      intensive calculations.
*
*  Michael D. Cummins
*    May 13, 1993
*/
{
  //FILE *fptr;
  //iter = 0;
  /* Inputs not changed: *****************************************************/
  double DegK;    /* Temperature                                 (Deg K) */
  double NH3;     /* Ammonia                                    (Mole/L) */
  double FreeCl2; /* [HOCl]+[OCl]                               (Mole/L) */

  /* Inputs changes: *********************************************************/
  double pH;         /*                                                 (-) */
  double Ca_aq;      /* Dissolved calcium   [Ca]+[CaOH]+[Ca(OH)2aq] (Mole/L) */
  double Mg_aq;      /* Dissolved manesium [Mg]+[MgOH]+[Mg(OH)2aq] (Mole/L) */
  double CO3_aq;     /* Dissolved carbonate   [H2CO3]+[HCO3]+[CO3] (Mole/L) */
  double CaCO3_floc; /* CaCO3 floc                                 (Mole/L) */
  double MgOH2_floc; /* Mg(OH)2 floc                               (Mole/L) */
  double CBminusCA;  /* Remainder of charge balance                 (equ/L) */

  /* Internal: *************************************************************/
  double Ca_total;     /* Total calcium including floc              (Mole/L) */
  double Mg_total;     /* Total magnesium including floc            (Mole/L) */
  double CO3_total;    /* Total carbonate including floc            (Mole/L) */
  double lo_pH, hi_pH; /* pH iteration range                             (-) */
                       /* Ions tracked by WTP:                               */
  double H;            /* [H+]                                      (Mole/L) */
  double OH;           /* [OH-]                                     (Mole/L) */
  double CO3;          /* [CO3--]                                   (Mole/L) */
  double HCO3;         /* [HCO3-]                                   (Mole/L) */
  double Ca;           /* [Ca++]                                    (Mole/L) */
  double CaOH;         /* [CaOH+]                                   (Mole/L) */
  double Mg;           /* [Mg++]                                    (Mole/L) */
  double MgOH;         /* [MgOH+]                                   (Mole/L) */
  double NH4;          /* [NH4+]                                    (Mole/L) */
  double OCl;          /* [OCl-]                                    (Mole/L) */
  double cations;      /* Sum of postive ions                        (equ/L) */
  double anions;       /* Sum of negative ions                       (equ/L) */

  double alpha_CO3; /* [CO3--]/[CO3_aq]                        (Fraction) */
  double alpha_Ca;  /* [Ca++]/[Ca_aq]                          (Fraction) */
  double alpha_Mg;  /* [Mg++]/[Mg_aq]                          (Fraction) */

  double error = 0.00001; /* acceptable error between hi and lo guesses to end iterative search */

  double b, c;
  struct Effluent *eff;

  /* Temperature dependent coefficients are static variables and recomputed
  *  only if the temperature changes.
  */
  static double kw;        /* Ionization coeffieient of water               */
  static double k1, k2;    /* Ionization coeffieients of carbonic acid      */
  static double k_hocl;    /* Ionization coeffieient of chlorine            */
  static double k_nh3;     /* Ionization coeffieient of ammonia.            */
  static double k_mgoh;    /* Ionization coeffieient of magnesium hydroxide.*/
  static double k_mgoh2;   /* Solubility of magnesium                       */
  static double k_mgoh2aq; /* Solubility of magnesium hydroxide.            */
  static double k_caoh;    /* Ionization of calcium hydroxide.              */
  static double k_caco3;   /* Solubility of calcium carbonate.              */
  static double k_caoh2aq; /* Solubility of calcium hydroxide.              */

  static struct Effluent old; /* Used to test changes in input parameters */

  /* Get inputs from UnitProcess data structure */
  eff = &unit->eff;
  pH = eff->pH;
  DegK = eff->DegK;
  Ca_aq = eff->Ca_aq;
  CaCO3_floc = eff->Ca_solid;
  Mg_aq = eff->Mg_aq;
  MgOH2_floc = eff->Mg_solid;
  CO3_aq = eff->CO2_aq;
  NH3 = eff->NH3;
  FreeCl2 = eff->FreeCl2;
  CBminusCA = eff->CBminusCA;

  /* Check if any of the inputs have changed.  */
  if (pH == old.pH &&
      DegK == old.DegK &&
      Ca_aq == old.Ca_aq &&
      CaCO3_floc == old.Ca_solid &&
      Mg_aq == old.Mg_aq &&
      MgOH2_floc == old.Mg_solid &&
      CO3_aq == old.CO2_aq &&
      NH3 == old.NH3 &&
      FreeCl2 == old.FreeCl2 &&
      CBminusCA == old.CBminusCA &&
      flag == TRUE)
  {
    return; /* The inputs have not changed. */
  }

  /* Temperature dependent coefficients: */
  if (DegK != old.DegK)
  {
    /* Temperature has changed... Recompute coefficients. */
    kw = Kw(DegK);
    k1 = K_HCO3(DegK);
    k2 = K_CO3(DegK);
    k_hocl = K_HOCl(DegK);
    k_nh3 = K_NH3(DegK);
    k_mgoh = K_MgOH(DegK);
    k_mgoh2 = K_MgOH2(DegK);
    k_mgoh2aq = K_MgOH2aq(DegK);
    k_caoh = K_CaOH(DegK);
    k_caco3 = K_CaCO3(DegK);
    k_caoh2aq = K_CaOH2aq(DegK);
  }

  /* Initial values */
  lo_pH = 3.0;
  hi_pH = 13.0;

  /*
  *  Ca_total, Mg_total, and CO3_total are mass balances and do not change
  *  by adjusting pH.
  */
  Ca_total = Ca_aq + CaCO3_floc;
  Mg_total = Mg_aq + MgOH2_floc;
  CO3_total = CO3_aq + CaCO3_floc;

#define CHECK FALSE

#if CHECK
  fptr = fopen("pHdata.dat", "a+");
  fprintf(fptr, "Module: %d  bottom \n", unit->type);
  fprintf(fptr, "phchange():\n");
  fprintf(fptr, "flag      :%d\n", flag);
  fprintf(fptr, "DegK      :%f\n", DegK);
  fprintf(fptr, "NH3       :%e\n", NH3);
  fprintf(fptr, "FreeCl2   :%e\n", FreeCl2);
  fprintf(fptr, "CBminusCA :%e\n", CBminusCA);
  fprintf(fptr, "Alk       :%e\n", eff->Alk);
  fprintf(fptr, "kw        :%e\n", kw);
  fprintf(fptr, "k_hco3    :%e\n", k1);
  fprintf(fptr, "k_co3     :%e\n", k2);
  fprintf(fptr, "k_hocl    :%e\n", k_hocl);
  fprintf(fptr, "k_nh3     :%e\n", k_nh3);
  fprintf(fptr, "k_mgoh    :%e\n", k_mgoh);
  fprintf(fptr, "k_mgoh2   :%e\n", k_mgoh2);
  fprintf(fptr, "k_mgoh2aq :%e\n", k_mgoh2aq);
  fprintf(fptr, "k_caoh    :%e\n", k_caoh);
  fprintf(fptr, "k_caco3   :%e\n", k_caco3);
  fprintf(fptr, "k_caoh2aq :%e\n", k_caoh2aq);
  fprintf(fptr, "Ca_total  :%e\n", Ca_total);
  fprintf(fptr, "Mg_total  :%e\n", Mg_total);
  fprintf(fptr, "CO3_total :%e\n", CO3_total);
#endif

  /* pH iteration loop */
  while ((hi_pH - lo_pH) > error) // changed by WJR to get more accurate estimate (6/17)
  {
    //   iter++;
    if (flag == TRUE)
      pH = (lo_pH + hi_pH) / 2;
    H = pow(10.0, -(pH));
    OH = kw / H;

    /*Assume any floc is dissolved then check solubility products.
      */
    CO3_aq = CO3_total;
    Ca_aq = Ca_total;
    Mg_aq = Mg_total;
    CaCO3_floc = 0.0;
    MgOH2_floc = 0.0;

    /* These are a function of pH and temperature.. Not concentrations. */
    alpha_CO3 = (k1 * k2) / ((H * H) + (k1 * H) + (k1 * k2));    /* A-30 */
    alpha_Ca = 1.0 / (1 + (k_caoh / H) + (k_caoh2aq / (H * H))); /* A-48 */
    alpha_Mg = 1.0 / (1 + (k_mgoh / H) + (k_mgoh2aq / (H * H))); /* Eq A-52 */
#if CHECK
    fprintf(fptr, "\n");
    fprintf(fptr, "pH        :%f\n", pH);
    fprintf(fptr, "H         :%e\n", H);
    fprintf(fptr, "OH        :%e\n", OH);
    fprintf(fptr, "alpha_CO3 :%e\n", alpha_CO3);
    fprintf(fptr, "alpha_Ca  :%e\n", alpha_Ca);
    fprintf(fptr, "alpha_Mg  :%e\n", alpha_Mg);
#endif

    /*
      * Determine if CaCO3 controls solubility of [Ca++] and [CO3--].
      * If so, calculate new calcium and carbonate distributions.
      *
      * The following [Ca++] and [CO3--] include any floc.
      */
    CO3 = CO3_aq * alpha_CO3;
    Ca = Ca_aq * alpha_Ca;
#if CHECK
    fprintf(fptr, "Before precipitation\n");
    fprintf(fptr, "CO3       :%e\n", CO3);
    fprintf(fptr, "Ca        :%e\n", Ca);
#endif

    /*if (unit->type == INFLUENT && flag == FALSE)
{
fptr=fopen("pHdata.dat","a+");
fprintf(fptr,"COLDFLAG: %d  \n",coldflag);
fprintf(fptr,"Module: %d  \n",unit->type);
fprintf(fptr,"phchange():\n");
fprintf(fptr,"Ca*CO3      :%e\n", Ca*CO3      );
fprintf(fptr,"k_caco3      :%e\n", k_caco3      );
fprintf(fptr,"ratio      :%f\n", Ca*CO3/k_caco3      );
fclose(fptr);
}   */

    if ((Ca * CO3) > (k_caco3) && flag == TRUE && unit->eff.limesoftening == TRUE)
    { /*
          * Solubility product of CaCO3 is exceded. Estimate CO3_aq
          * using eq A-70 through A-73.
          */
      b = Ca_total - CO3_total;
      c = -k_caco3 / (alpha_CO3 * alpha_Ca);
      CO3_aq = (-b + sqrt((b * b) - (4 * c))) / 2.0;
      CaCO3_floc = CO3_total - CO3_aq;
      Ca_aq = Ca_total - CaCO3_floc;

      /* Update [CO3] and [Ca] */
      CO3 = CO3_aq * alpha_CO3;
      Ca = Ca_aq * alpha_Ca;

      /*
if (unit->type == INFLUENT)
{
fptr=fopen("pHdata.dat","a+");
fprintf(fptr,"Module: %d  \n",unit->type);
fprintf(fptr,"in precip statement:\n");
fprintf(fptr,"ratio      :%e\n", Ca*CO3/k_caco3      );
fprintf(fptr,"pH      :%f\n", pH      );
fprintf(fptr,"lo_pH      :%f\n", lo_pH      );
fprintf(fptr,"hi_pH      :%f\n", hi_pH      );
fclose(fptr);
}
 */
    }

    /*
      * Determine if Mg(OH)2 controls solubility of magnesium.
      * If so, calculate new magnesium dissolved and floc concentrations.
      */
    Mg = Mg_aq * alpha_Mg;
#if CHECK
    fprintf(fptr, "Mg        :%e\n", Mg);
#endif

    if ((Mg / (H * H)) > k_mgoh2 && flag == TRUE && unit->eff.limesoftening == TRUE) /* Eq A-75 */
    {
      Mg_aq = k_mgoh2 * ((H * H) + (H * k_mgoh) + k_mgoh2aq); /* Eq A-78 */
      MgOH2_floc = Mg_total - Mg_aq;

      /* Update [Mg++] */
      Mg = Mg_aq * alpha_Mg;
    }

    /* Other Ions tracked by WTP */
    HCO3 = CO3_aq * (k1 * H) / ((H * H) + (k1 * H) + (k1 * k2)); /* A-30 */
    CaOH = Ca * k_caoh / H;                                      /* A-40 */
    MgOH = Mg * k_mgoh / H;                                      /* A-41 */
    OCl = FreeCl2 / (1 + (H / k_hocl));                          /* A-56 */
    NH4 = NH3 / (1 + (k_nh3 / H));                               /* A-60 */

    cations = CBminusCA + H + 2 * Ca + CaOH + 2 * Mg + MgOH + NH4; /* equ/L */
    anions = OH + HCO3 + 2 * CO3 + OCl;                            /* equ/L */

#if CHECK
    fprintf(fptr, "After precipitation\n");
    fprintf(fptr, "CO3       :%e\n", CO3);
    fprintf(fptr, "Ca        :%e\n", Ca);
    fprintf(fptr, "Mg        :%e\n", Mg);
    fprintf(fptr, "CO3_aq    :%e\n", CO3_aq);
    fprintf(fptr, "Ca_aq     :%e\n", Ca_aq);
    fprintf(fptr, "CaCO3_floc:%e\n", CaCO3_floc);
    fprintf(fptr, "Mg_aq     :%e\n", Mg_aq);
    fprintf(fptr, "MgOH2_floc:%e\n", MgOH2_floc);

    fprintf(fptr, "HCO3      :%e\n", HCO3);
    fprintf(fptr, "CaOH      :%e\n", CaOH);
    fprintf(fptr, "MgOH      :%e\n", MgOH);
    fprintf(fptr, "OCl       :%e\n", OCl);
    fprintf(fptr, "NH4       :%e\n", NH4);

    fprintf(fptr, "cations   :%e\n", cations);
    fprintf(fptr, "anions    :%e\n", anions);
#endif

    if (flag == TRUE)
    { /* Adjust pH */
      if ((cations - anions) > 0)
        lo_pH = pH;
      else
        hi_pH = pH;
    }
    else
    { /* Adjust CBminusCA such that cations==anions */
      CBminusCA -= (cations - anions);
      lo_pH = pH;
      hi_pH = pH;
    }
  }
  /*
if (unit->type == LIME && coldflag == FALSE)
{
fptr=fopen("pHdata.dat","a+");
fprintf(fptr,"Module: %d  \n",unit->type);
fprintf(fptr,"iter:            :%d\n",iter);
fprintf(fptr,"CO3_aq (CT)      :%e\n", CO3_aq);
fprintf(fptr,"CO3              :%e\n", CO3);
fprintf(fptr,"HCO3             :%e\n", HCO3);
fprintf(fptr,"OH	       :%e\n", OH);
fprintf(fptr,"kw	       :%e\n", kw);
fprintf(fptr,"Alk              :%e\n", (HCO3 + 2*CO3 + OH - H));
fprintf(fptr,"CT based on CO3  :%e\n", CO3/alpha_CO3);
fprintf(fptr,"CT based on Alk  :%e\n", CO3);
fprintf(fptr,"Ca               :%e\n", Ca);
fprintf(fptr,"Ca_aq            :%e\n", Ca_aq);
fprintf(fptr,"pH         :%f\n", pH      );
fprintf(fptr,"lo_pH      :%f\n", lo_pH      );
fprintf(fptr,"hi_pH      :%f\n\n", hi_pH      );
fclose(fptr);
}
*/

  /* Copy outputs to UnitProcess data structure */
  eff->Alk = HCO3 + 2 * CO3 + OH - H;
  eff->pH = f_adj_pH(pH);
  eff->Ca_aq = Ca_aq;
  eff->Ca_solid = CaCO3_floc;
  eff->Mg_aq = Mg_aq;
  eff->Mg_solid = MgOH2_floc;
  eff->CO2_aq = CO3_aq;
  eff->CBminusCA = CBminusCA;

#if CHECK
  fprintf(fptr, "CBminusCA :%e\n", CBminusCA);
  fprintf(fptr, "end phchange()\n");
  fclose(fptr);
#endif

  /* Update 'old' */
  old = *eff;
}
