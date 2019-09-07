/* Nf_Rmv.c */
#include "wtp.h"

void nf_rmv(struct UnitProcess *unit)
/*
*   This subroutine calculates changes in alkalinity, pH, hardness,
*   TOC, UV-254 and bromide through NF_UP.  pH is assumed
*   to remain the same and the new alkalinity is used to calculate a new
*   total carbonate concentration.  CBminusCA is also then re-calculated.
*   Solids are removed.
*
* Note: 1. Developed for SWAT
*/
{
  //FILE *fptr;

  /* Inputs: */
  double mwc;            /* Molecular weight cutoff (gm/mole = daltons) */
  double recovery;       /* % of flow through the NF membrane recovered */
  double toc_rem;        /* % TOC removal across the NF membrane */
  double uva_rem;        /* % UV-254 removal across the NF membrane */
  double br_rem;         /* % Bromide removal across the NF membrane */
  double treat_fraction; /* Fraction of flow treated */

  /* Inputs changed: */
  double ca;         /* Total Calcium               (Mole/L) */
  double mg;         /* Total Magnesium             (Mole/L) */
  double alk_equ;    /* Alkalinity                   (equ/L) */
  double bromide;    /* Bromide			  (ug/L)  */
  double eff_doc;    /* DOC                           (mg/L) */
  double eff_uv;     /* UV for internal calcs.        (1/cm) */
  double eff_uv_out; /* UV for outputting             (1/cm) */

  /* Outputs: */
  double co3; /* Total carbonate             (Mole/L) */
  double pH;  /*                                  (-) */
  double CBminusCA;

  /* Internal: */
  double kw;        /* Ionization of water                */
  double k1, k2;    /* Ionization of carbonic acid.       */
  double k_hocl;    /* Ionization of chlorine.            */
  double k_nh3;     /* Ionization of ammonia.             */
  double k_mgoh;    /* Ionization of magnesium hydorxide. */
  double k_mgoh2aq; /* Solubility of magnesium hydroxide. */
  double k_caoh;    /* Ionization of calcium hydroxide.   */
  double k_caoh2aq; /* Solubility of calcium hydroxide.   */
  double H;         /* [H+] (Mole/Liter)                  */
  double term;
  double removal;
  double hardness; /* Total hardness     (mg/L as CaCO3) */
  double alk;      /* Alkalinity         (mg/L as CaCO3) */
  double co3_term, ca_term, mg_term, nh3_term, ocl_term;

  double flow_fraction; /* The fraction of the flow in the effluent of the
                               previous unit process that flows to the next
                               unit process */
  double effluent_term;
  double bypass_term;

  struct Nf *nf;
  struct Effluent *eff;

  if (unit == NULL || unit->data.ptr == NULL || unit->type != NF_UP)
    return;

  /* Get inputs */
  nf = unit->data.nf;
  eff = &unit->eff;

  mwc = nf->mwc;
  recovery = nf->recover;
  toc_rem = nf->toc_rem;
  uva_rem = nf->uva_rem;
  br_rem = nf->br_rem;
  treat_fraction = nf->treat_fraction;

  eff_doc = eff->TOC;
  eff_uv = eff->UV;
  eff_uv_out = eff->UV_out;
  ca = eff->Ca_aq;
  mg = eff->Mg_aq;
  alk_equ = eff->Alk;
  bromide = eff->Br * MW_Br * 1000.0; /* converts from mol/L to ug/L*/
  pH = eff->pH;

  /* Temperature and PChem Coeffieients: */
  kw = Kw(eff->DegK);
  k1 = K_HCO3(eff->DegK);
  k2 = K_CO3(eff->DegK);
  k_hocl = K_HOCl(eff->DegK);
  k_nh3 = K_NH3(eff->DegK);
  k_mgoh = K_MgOH(eff->DegK);
  k_mgoh2aq = K_MgOH2aq(eff->DegK);
  k_caoh = K_CaOH(eff->DegK);
  k_caoh2aq = K_CaOH2aq(eff->DegK);

  /* Some Protection/Range-Checking */
  if (recovery < 1.0)
    recovery = 1.0;
  if (recovery > 100)
    recovery = 100;
  if (treat_fraction < 0.0)
    treat_fraction = 0.0;
  if (treat_fraction > 1.0)
    treat_fraction = 1.0;
  if (toc_rem < 0.0)
    toc_rem = 1.0;
  if (toc_rem > 100)
    toc_rem = 100;
  if (mwc < 1.0)
    mwc = 1.0;

  /* First, update eff->Flow */
  flow_fraction = 1.0 + treat_fraction * (recovery / 100.0 - 1.0);
  eff->Flow *= flow_fraction;

  /* Estimate TOC removal */
  if (eff_doc > 0.0)
  {
    effluent_term = treat_fraction * recovery / 100.0 *
                    eff_doc * (100.0 - toc_rem) / 100.0;
    bypass_term = (1.0 - treat_fraction) * eff_doc;
    eff_doc = (effluent_term + bypass_term) / flow_fraction;
    if (eff_doc > eff->TOC)
      eff_doc = eff->TOC;
  }

  /* Estimate UVA removal */
  if (eff_uv > 0.0)
  {
    effluent_term = treat_fraction * recovery / 100.0 *
                    eff_uv * (100.0 - uva_rem) / 100.0;
    bypass_term = (1.0 - treat_fraction) * eff_uv;
    eff_uv = (effluent_term + bypass_term) / flow_fraction;
    if (eff_uv > eff->UV)
      eff_uv = eff->UV;
  }
  if (eff_uv_out > 0.0)
  {
    effluent_term = treat_fraction * recovery / 100.0 *
                    eff_uv_out * (100.0 - uva_rem) / 100.0;
    bypass_term = (1.0 - treat_fraction) * eff_uv_out;
    eff_uv_out = (effluent_term + bypass_term) / flow_fraction;
    if (eff_uv_out > eff->UV_out)
      eff_uv_out = eff->UV_out;
  }

  /* Estimate Bromide removal */
  if (bromide > 0.0)
  {
    effluent_term = treat_fraction * recovery / 100.0 *
                    bromide * (100.0 - br_rem) / 100.0;
    bypass_term = (1.0 - treat_fraction) * bromide;
    bromide = (effluent_term + bypass_term) / flow_fraction;
    if (bromide > (eff->Br * MW_Br * 1000.0))
      bromide = eff->Br * MW_Br * 1000.0;
  }

  /* Estimate Hardness removal */
  hardness = ca * MW_CaCO3 + mg * MW_CaCO3;
  if (hardness > 0.0)
  {
    term = exp(36.801 - 3.327 * log(hardness) - 6.787 * log(mwc) - 0.027 * recovery * log(hardness) + 0.229 * log(recovery) * log(hardness) * log(mwc));
    removal = 1 / (term + 1.0);

    /* Calculate Ca */
    effluent_term = treat_fraction * recovery / 100.0 * ca * (1.0 - removal);
    bypass_term = (1.0 - treat_fraction) * ca;
    ca = (effluent_term + bypass_term) / flow_fraction;
    if (ca > eff->Ca_aq)
      ca = eff->Ca_aq;

    /* Calculate Mg */
    effluent_term = treat_fraction * recovery / 100.0 * mg * (1.0 - removal);
    bypass_term = (1.0 - treat_fraction) * mg;
    mg = (effluent_term + bypass_term) / flow_fraction;
    if (mg > eff->Mg_aq)
      mg = eff->Mg_aq;
  }

  /* Estimate alkalinity removal */
  alk = alk_equ * MW_CaCO3 / 2; /* (alk is mg/L as CaCO3 ) */
  if (alk > 0.0)
  {
    term = exp(14.602 - 1.667 * log(mwc) - 0.054 * log(alk) * log(mwc) - 0.203 * log(recovery) * log(mwc));
    removal = 1 / (term + 1.0);

    effluent_term = treat_fraction * recovery / 100.0 * alk * (1.0 - removal);
    bypass_term = (1.0 - treat_fraction) * alk;
    alk = (effluent_term + bypass_term) / flow_fraction;
    if (alk > (alk_equ * MW_CaCO3 / 2))
      alk = (alk_equ * MW_CaCO3 / 2);

    alk_equ = alk / (MW_CaCO3 / 2); /* (equ/L) */
  }

  /*  Assume pH does not change through membrane. Use new alk. and inf. pH to
      calculate new total carbonate concentration */
  H = pow(10.0, -(pH));
  co3_term = ((k1 * H) + (2 * k1 * k2)) / ((H * H) + (k1 * H) + (k1 * k2));
  co3 = (alk_equ - kw / H + H) / co3_term;

  /*  Calculate new CBminusCA  */
  ca_term = ca * (2 + (k_caoh / H)) / (1 + (k_caoh / H) + (k_caoh2aq / (H * H)));
  mg_term = mg * (2 + (k_mgoh / H)) / (1 + (k_mgoh / H) + (k_mgoh2aq / (H * H)));
  nh3_term = eff->NH3 / (1 + (k_nh3 / H));
  ocl_term = eff->FreeCl2 / (1 + (H / k_hocl));
  CBminusCA = alk_equ - (ca_term + mg_term + nh3_term - ocl_term);

  /* Copy outputs to UnitProcess data structure */
  eff->TOC = eff_doc;
  eff->UV = eff_uv;
  eff->UV_out = eff_uv_out;
  eff->Br = bromide / 1000.0 / MW_Br; /* Convert back to mol/L from ug/L */
  eff->Ca_aq = ca;
  eff->Mg_aq = mg;
  eff->CO2_aq = co3;
  eff->Alk = alk_equ;
  eff->CBminusCA = CBminusCA;

  eff->Ca_solid = 0.0;
  eff->Mg_solid = 0.0;
  /* Note: runmodel() will call phchange() following nf_rmv() */
}