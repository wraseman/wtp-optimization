/* gacmbdbp.c */
#include "wtp.h"

void gacmbdbp(struct UnitProcess *unit)
/*
* Purpose: Estimate trihalomethane and haloacetic acid formation in the
*          current UnitProcess based upon chlorination of GAC-treated or
*          membrane-treated water.  These equations are based upon the ICR
*          Treatment Study Data for GAC-treated water with TOC < 2.0 mg/L
*
* Notes:
*  1. gacmbdbp() is called by choose_dbpmodel when "gacmemdbpflag == TRUE"
*  2. The calling routine must update unit.effluent.hours to include the
*     contact time of this unit process.
*  3. inf_hours and eff_hours are cumulative from the last point of
*     chlorine addition.
*
* Documentation and code by WJS, 12/98
*/
{
  /* Inputs to model: */
  double pH;         /*                         (-) */
  double DegC;       /* Temperature         (Deg C) */
  double doc;        /* Really, its TOC      (mg/L) */
  double uv;         /* UV254 absorbance     (1/cm) */
  double doc_uv;     /* product of the above two    */
  double br;         /* Bromide              (ug/L) */
  double cl2dose;    /* Free chlorine (mg/L) as Cl2 */
  double eff_hours;  /* Contact time        (hours) */
  double inf_FreeCl; /* Process influent Free Chlorine      (Mole/L) */
  double inf_nh2cl;  /* Process influent Monochloramine     (Mole/L) */
  double inf_hours;  /* Cumulative to influent of UnitProcess */

  /* Outputs: */
  double delta_chcl3 = 0.0;   /* Chloroform           (ug/L) */
  double delta_chbrcl2 = 0.0; /* Bromodichloromethane (ug/L) */
  double delta_chbr2cl = 0.0; /* Dibromochloromethane (ug/L) */
  double delta_chbr3 = 0.0;   /* Bromoform            (ug/L) */
  double delta_tthm = 0.0;    /* Total THM            (ug/L) */

  /* Outputs: */
  double delta_mcaa = 0.0; /* Monochloroacetic acid (ug/L) */
  double delta_dcaa = 0.0; /* Dichloroacetic   acid (ug/L) */
  double delta_tcaa = 0.0; /* Trichloroacetic  acid (ug/L) */
  double delta_mbaa = 0.0; /* Monobromoacetic  acid (ug/L) */
  double delta_dbaa = 0.0; /* Dibromoacetic    acid (ug/L) */
  double delta_haa5 = 0.0; /* sum of above 5 HAAs   (ug/L) */
  double delta_bcaa = 0.0; /* Bromochloroacetic acid (ug/L)*/
  double delta_haa6 = 0.0; /* sum of above 6 HAAs   (ug/L) */
  double delta_tbaa = 0.0;
  double delta_dbcaa = 0.0;
  double delta_bdcaa = 0.0;
  double delta_haa9 = 0.0;
  double delta_haa6_for_haa9 = 0.0;
  double delta_haa9_haa6 = 0.0;
  double delta_tox = 0.0;

  /* Internal: */
  double factor, term1, term2, term3, term4, term5, term6;
  double eff_Br, moles_br_incorporated;
  register struct Effluent *eff;

  /* Get inputs from UnitProcess data structure */
  eff = &unit->eff;
  pH = eff->pH;
  DegC = eff->DegK - 273.15;
  //FreeCl    = eff->FreeCl2;
  //nh2cl     = eff->NH2Cl;
  doc = eff->TOC;
  uv = eff->UV;
  br = eff->Br_at_last_cl2 * MW_Br * 1000.0; /* to convert to ug/L */
  cl2dose = eff->cl2dose;
  eff_hours = eff->hours;

  /* Self-protection - Set reasonable minimum values */
  if (DegC < 1.0)
    DegC = 1.0;
  if (br < 1.0)
    br = 1.0;
  if (uv <= 0.001)
    uv = 0.001;

  if (doc > 0.0 && doc < 0.1)
    doc = 0.1; /* Min. non-zero doc */

  doc_uv = doc * uv;

  /*  Get Cumulative time to influent of UnitProcess.    */
  if (eff->wtp_effluent != NULL)
  {
    /* 'unit' is a distribution sample -- get influent time from 
    *  WTP_EFFLUENT.
    */
    inf_hours = eff->wtp_effluent->eff.hours;
    inf_FreeCl = eff->wtp_effluent->eff.FreeCl2;
    inf_nh2cl = eff->wtp_effluent->eff.NH2Cl;
  }
  else
  {
    if (PrevUnitProcess(unit) != NULL)
    {
      inf_hours = PrevUnitProcess(unit)->eff.hours;
      inf_FreeCl = PrevUnitProcess(unit)->eff.FreeCl2;
      inf_nh2cl = PrevUnitProcess(unit)->eff.NH2Cl;
    }
    else
    { /* 'unit' is first unit process in process train. */
      inf_hours = 0.0;
      inf_FreeCl = 0.0;
      inf_nh2cl = 0.0;
    }
  }

  if ((inf_FreeCl > 0.0 || inf_nh2cl > 0.0) && doc > 0.0 && cl2dose > 0.0 && eff_hours > 0.0)
  { /* THMs and/or HAAs should be forming */

    if (doc <= 2.0)
    { /* Use the GAC-treated water DBP Formation Equations */

      /* Estimate TTHM (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.474);
      term2 = pow(cl2dose, 0.173);
      term3 = pow(br, 0.245);
      term4 = pow(1.036, (DegC - 20));
      term5 = pow(1.316, (pH - 8.0));
      term6 = pow(eff_hours, 0.365) - pow(inf_hours, 0.365);
      delta_tthm = 17.73 * term1 * term2 * term3 * term4 * term5 * term6;

      /* Estimate CHCl3 (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.614);
      term2 = pow(cl2dose, 0.698);
      term3 = pow(br, -0.468);
      term4 = pow(1.035, (DegC - 20));
      term5 = pow(1.099, (pH - 8.0));
      term6 = pow(eff_hours, 0.336) - pow(inf_hours, 0.336);
      delta_chcl3 = 101.0 * term1 * term2 * term3 * term4 * term5 * term6;

      /* Estimate CHBrCl2 (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.430);
      term2 = pow(cl2dose, 0.555);
      term3 = pow(br, 0.075);
      term4 = pow(1.030, (DegC - 20));
      term5 = pow(1.349, (pH - 8.0));
      term6 = pow(eff_hours, 0.278) - pow(inf_hours, 0.278);
      delta_chbrcl2 = 7.37 * term1 * term2 * term3 * term4 * term5 * term6;

      /* Estimate CHBr2Cl (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.525);
      term2 = pow(cl2dose, 0.125);
      term3 = pow(br, 0.362);
      term4 = pow(1.027, (DegC - 20));
      term5 = pow(1.432, (pH - 8.0));
      term6 = pow(eff_hours, 0.319) - pow(inf_hours, 0.319);
      delta_chbr2cl = 3.987 * term1 * term2 * term3 * term4 * term5 * term6;

      /* Estimate CHBr3  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.405);
      term2 = pow(cl2dose, -0.118);
      term3 = pow(br, 0.950);
      term4 = pow(1.048, (DegC - 20));
      term5 = pow(1.430, (pH - 8.0));
      term6 = pow(eff_hours, 0.323) - pow(inf_hours, 0.323);
      delta_chbr3 = 0.158 * term1 * term2 * term3 * term4 * term5 * term6;

      /* Estimate HAAs */

      /* Estimate HAA6  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.502);
      term2 = pow(cl2dose, 0.371);
      term3 = pow(br, -0.078);
      term4 = pow(1.021, (DegC - 20));
      term5 = pow(0.913, (pH - 8.0));
      term6 = pow(eff_hours, 0.279) - pow(inf_hours, 0.279);
      delta_haa6 = 36.9 * term1 * term2 * term3 * term4 * term5 * term6;

      /* Estimate HAA5  (ug/liter) ------------------------ Equation XX  */
      /*  term1 = pow( doc_uv   , 0.488  );
      term2 = pow( cl2dose  , 0.385  );
      term3 = pow( br	    ,-0.155  );
      term4 = pow( 1.021    , (DegC-20));
      term5 = pow( 0.867    , (pH-8.0));
      term6 = pow( eff_hours, 0.262 ) - pow( inf_hours, 0.262 );
      delta_haa5 = 40.0*term1*term2*term3*term4*term5*term6;    */

      /* Estimate MCAA  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.043);
      term2 = pow(cl2dose, 0.126);
      term3 = pow(br, -0.135);
      term4 = pow(1.009, (DegC - 20));
      term5 = pow(0.873, (pH - 8.0));
      term6 = pow(eff_hours, 0.087) - pow(inf_hours, 0.087);
      delta_mcaa = 1.75 * term1 * term2 * term3 * term4 * term5 * term6;
      if (delta_mcaa > 8.0)
        delta_mcaa = 8.0; //based on max. value in development database

      /* Estimate MBAA  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.067);
      term2 = pow(cl2dose, 0.253);
      term3 = pow(br, -0.044);
      term4 = pow(1.017, (DegC - 20));
      term5 = pow(0.845, (pH - 8.0));
      term6 = pow(eff_hours, 0.097) - pow(inf_hours, 0.097);
      delta_mbaa = 0.666 * term1 * term2 * term3 * term4 * term5 * term6;
      if (delta_mbaa > 6.0)
        delta_mbaa = 6.0; //based on max. value in development database

      /* Estimate DCAA  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.473);
      term2 = pow(cl2dose, 0.415);
      term3 = pow(br, -0.380);
      term4 = pow(1.019, (DegC - 20));
      term5 = pow(0.870, (pH - 8.0));
      term6 = pow(eff_hours, 0.288) - pow(inf_hours, 0.288);
      delta_dcaa = 34.6 * term1 * term2 * term3 * term4 * term5 * term6;

      /* Estimate TCAA  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.557);
      term2 = pow(cl2dose, 0.692);
      term3 = pow(br, -0.394);
      term4 = pow(1.010, (DegC - 20));
      term5 = pow(0.619, (pH - 8.0));
      term6 = pow(eff_hours, 0.167) - pow(inf_hours, 0.167);
      delta_tcaa = 37.4 * term1 * term2 * term3 * term4 * term5 * term6;

      /* Estimate DBAA  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.467);
      term2 = pow(cl2dose, -0.237);
      term3 = pow(br, 0.640);
      term4 = pow(1.019, (DegC - 20));
      term5 = pow(1.286, (pH - 8.0));
      term6 = pow(eff_hours, 0.295) - pow(inf_hours, 0.295);
      delta_dbaa = 0.467 * term1 * term2 * term3 * term4 * term5 * term6;

      /* Estimate BCAA  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.484);
      term2 = pow(cl2dose, 0.245);
      term3 = pow(br, 0.110);
      term4 = pow(1.018, (DegC - 20));
      term5 = pow(1.064, (pH - 8.0));
      term6 = pow(eff_hours, 0.313) - pow(inf_hours, 0.313);
      delta_bcaa = 3.43 * term1 * term2 * term3 * term4 * term5 * term6;

      /*  Estimate TBAA  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.385);
      term2 = pow(cl2dose, -0.844);
      term3 = pow(br, 0.532);
      term4 = pow(1.003, (DegC - 20.0));
      term5 = pow(0.829, (pH - 8.0));
      term6 = pow(eff_hours, 0.464) - pow(inf_hours, 0.464);
      delta_tbaa = 1.83 * pow(10, -1.0) * term1 * term2 * term3 * term4 * term5 * term6;
      if (delta_tbaa > 10.0)
        delta_tbaa = 10.0; //based on max. value in development database

      /*  Estimate DBCAA  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.316);
      term2 = pow(cl2dose, -0.081);
      term3 = pow(br, 0.368);
      term4 = pow(1.008, (DegC - 20.0));
      term5 = pow(0.716, (pH - 8.0));
      term6 = pow(eff_hours, 0.363) - pow(inf_hours, 0.363);
      delta_dbcaa = 4.35 * pow(10, -1.0) * term1 * term2 * term3 * term4 * term5 * term6;

      /*  Estimate BDCAA  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.518);
      term2 = pow(cl2dose, 0.248);
      term3 = pow(br, 0.190);
      term4 = pow(0.987, (DegC - 20.0));
      term5 = pow(0.613, (pH - 8.0));
      term6 = pow(eff_hours, 0.317) - pow(inf_hours, 0.317);
      delta_bdcaa = 2.0 * term1 * term2 * term3 * term4 * term5 * term6;

      /*  Estimate HAA9  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.509);
      term2 = pow(cl2dose, 0.253);
      term3 = pow(br, 0.053);
      term4 = pow(1.019, (DegC - 20.0));
      term5 = pow(0.828, (pH - 8.0));
      term6 = pow(eff_hours, 0.425) - pow(inf_hours, 0.425);
      delta_haa9 = 20.6 * term1 * term2 * term3 * term4 * term5 * term6;

      /*  Estimate HAA6 for HAA9  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.320);
      term2 = pow(cl2dose, 0.510);
      term3 = pow(br, -0.106);
      term4 = pow(1.017, (DegC - 20.0));
      term5 = pow(0.938, (pH - 8.0));
      term6 = pow(eff_hours, 0.353) - pow(inf_hours, 0.353);
      delta_haa6_for_haa9 = 18.6 * term1 * term2 * term3 * term4 * term5 * term6;

      /*  Estimate TOX  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.529);
      term2 = pow(cl2dose, 0.349);
      term4 = pow(1.009, (DegC - 20.0));
      term6 = pow(eff_hours, 0.239) - pow(inf_hours, 0.239);
      delta_tox = 168.0 * term1 * term2 * term4 * term6;

      /* No correction based on ICR data for the GAC-treated water equations*/

    } /* End "if(doc < 2.0)"*/
    else
    { /* doc > 2.0 --> Use coagulated water DBP-formation equations with
                         ICR correction factors*/

      /* Estimate TTHM (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.403);
      term2 = pow(cl2dose, 0.225);
      term3 = pow(br, 0.141);
      term4 = pow(1.026, (DegC - 20));
      term5 = pow(1.156, (pH - 7.5));
      term6 = pow(eff_hours, 0.264) - pow(inf_hours, 0.264);
      delta_tthm = 23.9 * term1 * term2 * term3 * term4 * term5 * term6;

      /* Estimate CHCl3 (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.483);
      term2 = pow(cl2dose, 0.424);
      term3 = pow(br, -0.679);
      term4 = pow(1.018, (DegC - 20));
      term5 = pow(1.132, (pH - 7.5));
      term6 = pow(eff_hours, 0.333) - pow(inf_hours, 0.333);
      delta_chcl3 = 266.0 * term1 * term2 * term3 * term4 * term5 * term6;

      /* Estimate CHBrCl2 (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.260);
      term2 = pow(cl2dose, 0.114);
      term3 = pow(br, 0.462);
      term4 = pow(1.026, (DegC - 20));
      term5 = pow(1.098, (pH - 7.5));
      term6 = pow(eff_hours, 0.196) - pow(inf_hours, 0.196);
      delta_chbrcl2 = 1.68 * term1 * term2 * term3 * term4 * term5 * term6;

      /* Estimate CHBr2Cl (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, -0.056);
      term2 = pow(cl2dose, -0.157);
      term3 = pow(br, 1.425);
      term4 = pow(1.021, (DegC - 20));
      term5 = pow(1.127, (pH - 7.5));
      term6 = pow(eff_hours, 0.148) - pow(inf_hours, 0.148);
      delta_chbr2cl = 0.0080 * term1 * term2 * term3 * term4 * term5 * term6;

      /* Estimate CHBr3  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, -0.300);
      term2 = pow(cl2dose, -0.221);
      term3 = pow(br, 2.134);
      term4 = pow(1.037, (DegC - 20));
      term5 = pow(1.391, (pH - 7.5));
      term6 = pow(eff_hours, 0.143) - pow(inf_hours, 0.143);
      delta_chbr3 = 4.4 * pow(10, -5.0) * term1 * term2 * term3 * term4 * term5 * term6;

      /* Estimate HAAs */

      /* Estimate HAA5  (ug/liter) ------------------------ Equation XX */
      /*    term1 = pow( doc_uv   , 0.302 );
      term2 = pow( cl2dose  , 0.541 );
      term3 = pow( br	    ,-0.012 );
      term4 = pow( 1.022    , (DegC-20));
      term5 = pow( 0.922    , (pH-7.5));
      term6 = pow( eff_hours, 0.161 ) - pow( inf_hours, 0.161 );
      delta_haa6 = 30.7*term1*term2*term3*term4*term5*term6;  */

      /* WJR COMMENT OUT BEGIN 
     * Estimate HAA6  (ug/liter) ------------------------ Equation XX 
      term1 = pow( doc_uv   , 0.328  );
      term2 = pow( cl2dose  , 0.585  );
      term3 = pow( br	    ,-0.121  );
      term4 = pow( 1.021    , (DegC-20));
      term5 = pow( 0.932    , (pH-7.5));
      term6 = pow( eff_hours, 0.150 ) - pow( inf_hours, 0.150 );
      delta_haa5 = 41.6*term1*term2*term3*term4*term5*term6;  
       * WJR COMMENT OUT END, 4/10 */

      // BEGIN WJR ADD 4/10
      /* Estimate HAA6  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.328);
      term2 = pow(cl2dose, 0.585);
      term3 = pow(br, -0.121);
      term4 = pow(1.022, (DegC - 20));
      term5 = pow(0.922, (pH - 7.5));
      term6 = pow(eff_hours, 0.150) - pow(inf_hours, 0.150);
      delta_haa6 = 30.7 * term1 * term2 * term3 * term4 * term5 * term6;

      // END WJR ADD 4/10

      /* Estimate MCAA  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, -0.090);
      term2 = pow(cl2dose, 0.662);
      term3 = pow(br, -0.224);
      term4 = pow(1.024, (DegC - 20));
      term5 = pow(1.042, (pH - 7.5));
      if (inf_hours > 0.0)
        term6 = pow(eff_hours, -0.043) - pow(inf_hours, -0.043);
      else
        term6 = pow(eff_hours, -0.043);
      delta_mcaa = 4.58 * term1 * term2 * term3 * term4 * term5 * term6;
      if (delta_mcaa > 33.0)
        delta_mcaa = 33.0; //based on model development database inputs and eqn.

      /* Estimate MBAA  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.358);
      term2 = pow(cl2dose, -0.101);
      term3 = pow(br, 0.812);
      term4 = pow(1.162, (DegC - 20));
      term5 = pow(0.653, (pH - 7.5));
      term6 = pow(eff_hours, 0.088) - pow(inf_hours, 0.088);
      delta_mbaa = 0.0206 * term1 * term2 * term3 * term4 * term5 * term6;

      /* Estimate DCAA  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.397);
      term2 = pow(cl2dose, 0.665);
      term3 = pow(br, -0.558);
      term4 = pow(1.017, (DegC - 20));
      term5 = pow(1.034, (pH - 7.5));
      term6 = pow(eff_hours, 0.222) - pow(inf_hours, 0.222);
      delta_dcaa = 60.4 * term1 * term2 * term3 * term4 * term5 * term6;

      /* Estimate TCAA  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.403);
      term2 = pow(cl2dose, 0.749);
      term3 = pow(br, -0.416);
      term4 = pow(1.014, (DegC - 20));
      term5 = pow(0.874, (pH - 7.5));
      term6 = pow(eff_hours, 0.163) - pow(inf_hours, 0.163);
      delta_tcaa = 52.6 * term1 * term2 * term3 * term4 * term5 * term6;

      /* Estimate DBAA  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.059);
      term2 = pow(cl2dose, 0.182);
      term3 = pow(br, 2.109);
      term4 = pow(1.007, (DegC - 20));
      term5 = pow(1.21, (pH - 7.5));
      term6 = pow(eff_hours, 0.070) - pow(inf_hours, 0.070);
      delta_dbaa = 9.42 * pow(10, -5.0) * term1 * term2 * term3 * term4 * term5 * term6;

      /* Estimate BCAA  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.153);
      term2 = pow(cl2dose, 0.257);
      term3 = pow(br, 0.586);
      term4 = pow(1.042, (DegC - 20));
      term5 = pow(1.181, (pH - 7.5));
      term6 = pow(eff_hours, 0.201) - pow(inf_hours, 0.201);
      delta_bcaa = 0.323 * term1 * term2 * term3 * term4 * term5 * term6;

      //  The following are coag.-water HAA9 species and TOX models developed by MPI in 4/2001

      /*  Estimate TBAA  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.0657);
      term2 = pow(cl2dose, -2.51);
      term3 = pow(br, 2.32);
      term4 = pow(1.059, (DegC - 20.0));
      term5 = pow(0.555, (pH - 8.0));
      term6 = pow(eff_hours, 1.26) - pow(inf_hours, 1.26);
      delta_tbaa = 5.59 * pow(10, -6.0) * term1 * term2 * term3 * term4 * term5 * term6;

      /*  Estimate DBCAA  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, -0.0162);
      term2 = pow(cl2dose, -0.170);
      term3 = pow(br, 0.972);
      term4 = pow(1.054, (DegC - 20.0));
      term5 = pow(0.839, (pH - 8.0));
      term6 = pow(eff_hours, 0.685) - pow(inf_hours, 0.685);
      delta_dbcaa = 3.7 * pow(10, -3.0) * term1 * term2 * term3 * term4 * term5 * term6;

      /*  Estimate BDCAA  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.230);
      term2 = pow(cl2dose, 0.140);
      term3 = pow(br, 0.301);
      term4 = pow(1.022, (DegC - 20.0));
      term5 = pow(0.700, (pH - 8.0));
      term6 = pow(eff_hours, 0.422) - pow(inf_hours, 0.422);
      delta_bdcaa = 0.589 * term1 * term2 * term3 * term4 * term5 * term6;

      /*  Estimate HAA9  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.250);
      term2 = pow(cl2dose, 0.500);
      term3 = pow(br, 0.054);
      term4 = pow(1.015, (DegC - 20.0));
      term5 = pow(0.894, (pH - 8.0));
      term6 = pow(eff_hours, 0.348) - pow(inf_hours, 0.348);
      delta_haa9 = 10.783 * term1 * term2 * term3 * term4 * term5 * term6;

      /*  Estimate HAA6 for HAA9  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.320);
      term2 = pow(cl2dose, 0.510);
      term3 = pow(br, -0.106);
      term4 = pow(1.017, (DegC - 20.0));
      term5 = pow(0.938, (pH - 8.0));
      term6 = pow(eff_hours, 0.353) - pow(inf_hours, 0.353);
      delta_haa6_for_haa9 = 18.6 * term1 * term2 * term3 * term4 * term5 * term6;

      /*  Estimate TOX  (ug/liter) ------------------------ Equation XX */
      term1 = pow(doc_uv, 0.362);
      term2 = pow(cl2dose, 0.129);
      term4 = pow(DegC, 0.211);
      term6 = pow(eff_hours, 0.182) - pow(inf_hours, 0.182);
      delta_tox = 109.0 * term1 * term2 * term4 * term6;

      /********CORRECTION FACTORS BASED ON ICR DATA CALIBRATION WORK (WJS, 5/2001)********/

      /* Correction factors to apply for free chlorine based on comparing
       predicted vs. observed in months 7-18 in AUX8 at the FINISH water
       location for plantmonths with Validation mode = VALID,
       WTP_DIS = CL2, predicted alk > 0 at FINISH, and predicted nh2cl = 0
       at FINISH */
      delta_tthm /= 0.77;
      delta_chcl3 /= 1.00;
      delta_chbr2cl /= 0.50;
      delta_chbrcl2 /= 0.86;
      delta_chbr3 /= 1.00;
      //    delta_haa5    /= 1.17;
      delta_haa6 /= 1.0;
      delta_mcaa /= 1.00;
      delta_dcaa /= 0.71;
      delta_tcaa /= 1.3;
      delta_mbaa /= 1.00;
      delta_dbaa /= 1.00;
      delta_bcaa /= 0.82;

      /* In Distribution system, need several different correction factors */
      if (unit->type == AVG_TAP || unit->type == END_OF_SYSTEM ||
          unit->type == LOCATION_1)
      {
        delta_dcaa *= (0.71 / 1.3);
        delta_bcaa *= (0.82 / 2.0);
        delta_chcl3 /= 1.1;
        delta_chbr2cl /= 0.80;
      }

      /******************************END OF FREE CL2 CORRECTION FACTORS*******************************/
      /*(FOR COAG.-TREATED WATER DBPs ONLY)*/

    } /* End "if-else(doc < 2.0)"*/

    /* If chloramines exist, use appropriate factors for the ratio of formation
       with chloramines to that with free chlorine (20% is the default, the
       others were set based on calibration with AUX8 ICR plants with CLM as
       the DS disinfectant) */

    if (inf_FreeCl <= 0.0 && inf_nh2cl > 0.0)
    {
      delta_tthm *= 0.3;
      delta_chcl3 *= 0.3;
      delta_chbr2cl *= 0.3;
      delta_chbrcl2 *= 0.3;
      delta_chbr3 *= 0.3;
    }

    if (inf_nh2cl > 0.0)
    {
      //     delta_haa5 *= 0.05;
      delta_haa6 *= 0.35;
      delta_haa9 *= 0.2;
      delta_haa6_for_haa9 *= 0.2;
      delta_mcaa *= 0.2;
      delta_dcaa *= 0.35;
      delta_tcaa *= 0.05;
      delta_mbaa *= 0.2;
      delta_dbaa *= 0.2;
      delta_bcaa *= 0.3;
      delta_tbaa *= 0.2;
      delta_dbcaa *= 0.2;
      delta_bdcaa *= 0.2;
      delta_tox *= 0.325;
    }

    // MCAA should only be able to go to zero
    if ((delta_mcaa + eff->MCAA) < 0.0)
      delta_mcaa = -eff->MCAA;

    /* Proportion the species based on the bulk TTHM */
    if (delta_chcl3 > 0.0 || delta_chbrcl2 > 0.0 || delta_chbr2cl > 0.0 || delta_chbr3 > 0.0)
    {
      factor = delta_tthm / (delta_chcl3 + delta_chbrcl2 +
                             delta_chbr2cl + delta_chbr3);
      delta_chcl3 *= factor;
      delta_chbrcl2 *= factor;
      delta_chbr2cl *= factor;
      delta_chbr3 *= factor;
    }
    else
      delta_tthm = 0.0;

    /* Proportion-out 6 HAA species based on HAA6 equation */
    if (delta_mcaa > 0.0 || delta_mbaa > 0.0 || delta_bcaa > 0.0 ||
        delta_dcaa > 0.0 || delta_dbaa > 0.0 || delta_tcaa > 0.0)
    {
      factor = delta_haa6 / (delta_mcaa + delta_mbaa + delta_bcaa +
                             delta_dcaa + delta_dbaa + delta_tcaa);
      delta_mcaa *= factor;
      delta_mbaa *= factor;
      delta_dcaa *= factor;
      delta_dbaa *= factor;
      delta_tcaa *= factor;
      delta_bcaa *= factor;
    }
    else
      delta_haa6 = 0.0;

    /* Calc. HAA5 based on HAA6 equation, minus proportioned BCAA */
    delta_haa5 = delta_haa6 - delta_bcaa;

    /* Proportion the HAA9 species based on HAA9-HAA6 */
    if (delta_haa9 > delta_haa6_for_haa9)
    {
      delta_haa9_haa6 = delta_haa9 - delta_haa6_for_haa9;
    }
    else
    {
      delta_haa9_haa6 = 0.0;
      delta_dbcaa = 0.0;
      delta_bdcaa = 0.0;
      delta_tbaa = 0.0;
    }

    if (delta_dbcaa > 0.0 || delta_bdcaa > 0.0 || delta_tbaa > 0.0)
    {
      factor = delta_haa9_haa6 / (delta_dbcaa + delta_bdcaa + delta_tbaa);
      delta_dbcaa *= factor;
      delta_bdcaa *= factor;
      delta_tbaa *= factor;
    }

    /* Calc. HAA9 based on HAA5, plus species */
    delta_haa9 = delta_haa6 + delta_tbaa + delta_dbcaa + delta_bdcaa;

    /* Add incremental outputs in UnitProcess data structure. */
    /* THMs */
    eff->CHCl3 += delta_chcl3;
    eff->CHBrCl2 += delta_chbrcl2;
    eff->CHBr2Cl += delta_chbr2cl;
    eff->CHBr3 += delta_chbr3;
    eff->TTHM += delta_tthm;

    /* HAAs */
    eff->MCAA += delta_mcaa;
    eff->DCAA += delta_dcaa;
    eff->TCAA += delta_tcaa;
    eff->MBAA += delta_mbaa;
    eff->DBAA += delta_dbaa;
    eff->BCAA += delta_bcaa;
    eff->TBAA += delta_tbaa;
    eff->DBCAA += delta_dbcaa;
    eff->BDCAA += delta_bdcaa;
    eff->HAA5 += delta_haa5;
    eff->HAA6 += delta_haa6;
    eff->HAA9 += delta_haa9;

    /*TOX */
    eff->tox += delta_tox;

    /*Stoichiometric Bromide incorporation (non-THM or HAA TOX not considered)*/
    moles_br_incorporated =
        (delta_chbrcl2 / 1000.0) / MW_BDCM * 1.0 +
        (delta_chbr2cl / 1000.0) / MW_DBCM * 2.0 +
        (delta_chbr3 / 1000.0) / MW_CHBR3 * 3.0 +
        (delta_mbaa / 1000.0) / MW_MBAA * 1.0 +
        (delta_bcaa / 1000.0) / MW_BCAA * 1.0 +
        (delta_dbaa / 1000.0) / MW_DBAA * 2.0 +
        (delta_bdcaa / 1000.0) / MW_BDCAA * 1.0 +
        (delta_dbcaa / 1000.0) / MW_DBCAA * 2.0 +
        (delta_tbaa / 1000.0) / MW_TBAA * 3.0;

    eff_Br = eff->Br - moles_br_incorporated;

    if (eff_Br > eff->Br)
      eff_Br = eff->Br; /* No increasing Br allowed */
    if (eff_Br < 0.0)
      eff_Br = 0.0; /* No negative Br allowed   */

    eff->Br = eff_Br; /*Update UnitProcess data structure */

    /***********************End of Bromide Incorporation**********************/

    /* Because of possibility of a negative delta_mcaa, the HAA5 species
       could all possibly be negative due to proportioning*/
    if (eff->MCAA < 0.0)
      eff->MCAA = 0.0;
    if (eff->DCAA < 0.0)
      eff->DCAA = 0.0;
    if (eff->TCAA < 0.0)
      eff->TCAA = 0.0;
    if (eff->MBAA < 0.0)
      eff->MBAA = 0.0;
    if (eff->DBAA < 0.0)
      eff->DBAA = 0.0;
    if (eff->BCAA < 0.0)
      eff->BCAA = 0.0;
  }
}
