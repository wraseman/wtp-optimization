/* Coef.c --  Temperature dependent ionization and solubility products.
*
* Notation:
*   1. The input to all functions is temperature in (Deg K).
*   2. Numerical equations are from Appendix A of 1.21 WTP manual.
*
* The following functions are in this file.
*/
#include "wtp.h"

double Kw(double DegK) /* Ionization of water. */
/*
*  Reaction: H2O <=> [H+] + [OH-]
*  Returned: Kw = [H+][OH-]  range 1.0E-13 to 1.0E-16
*/
{
  /* Self protection */
  if (DegK < 273.15)
  {
    printf("Temperature must be in Deg K in Kw()\n");
    return (0.0);
  }

  /* From Stumm and Morgan (1981) */
  /* Equation A-36 in 1.21 version of WTP manual Appendix A */
  return (pow(10.0, (6.0875 - (4470.99 / DegK) - (0.01706 * DegK)))); /* A-36 */
}

double K_HCO3(double DegK)
/*
*  First ionization coefficient for carbonic acid.
*
*  From Harned, H. S. and R. Davies, Jr., Am. Chem. Soc. J. 65, 2030, 1943.
*  Temperature range 0 to 50 deg C.
*
*  Equation 15 on page 2036 of Harned, 1943:
*   log10(K1) = -(3404.71/T) + 14.8435 -(0.032786*T)
*
*  Michael D. Cummins
*    May 12, 1993
*/
{
  return (pow(10.0, (-(3404.71 / DegK) + 14.8435 - (0.032786 * DegK))));
}

double K_CO3(double DegK)
/*
*  Second ionization coefficient of carbonic acid.
*  Temperature range is from 0 to 50 deg C.
*
*  Reaction:
*    [HCO3-] <=> [H+] + [CO3--]
*    K2 = [H+] [CO3--] / [HCO3-]
*
*  Equation 9 on page 1709 from Harned, H. S. and Scholes, S. R.,
*  Am. Chem. Soc. J. 63, 1706, 1941. 
*    log10(K2) = -2902.39/T + 6.4980 -0.02379*T
*
*  Michael D. Cummins
*     May, 12 1993
*/
{
  return (pow(10.0, (-(2902.39 / DegK) + 6.4980 - (0.02379 * DegK))));
}

double K_HOCl(double DegK) /* Ionization of free chlorine. */
/*
*  Reaction:  [HOCl] <=> [OCl=] + [H+]
*  Returned:  K_HOCl = [H+][OCl-]/[HOCl]
*/
{
  double f_temp;

  f_temp = (1 / 293.15) - (1 / DegK);
  return (exp((13800 * f_temp / 8.31441) - 17.500)); /* A-55 */
}

double K_NH3(double DegK) /* Ionization of ammonia. */
/*
*  Reaction:  [NH4+] <=> [H+] + [NH3]
*  Returned:  K_NH3 = [H+][NH3]/[NH4+]
*/
{
  double f_temp;

  f_temp = (1 / 293.15) - (1 / DegK);
  return (exp((52210.0 * f_temp / 8.31441) - 21.414)); /* A-59 */
}

double K_CaOH(double DegK) /* Ionization of Calcium hydroxide. */
/*
*  Reaction: [Ca++] + [H2O] <=> [CaOH+] + [H+]
*  Returned: k_CaOH = [CaOH+][H+]/[Ca++]
*/
{
  return (exp(-72320.0 / (8.31441 * DegK))); /* A-42 */
}

double K_CaOH2aq(double DegK) /* Solubility of Calcium hydroxied. */
/*
*  Reaction: [Ca++] + 2[H2O] <=> 2[H+] + [Ca(OH)2] (Solid)
*  Returned: k_CaOH2aq = [Ca(OH)2][H+][H+]/[Ca++]
*/
{
  return (exp(-159800.0 / (8.31441 * DegK))); /* A-47 */
}

double K_CaCO3(double DegK) /* Solubility of Calcium carbonate. */
/*
*  Reaction: [Ca++] + [CO3--] => [CaCO3] (Solid)
*  Returned: k_CaCO3 = [Ca++][CO3--]
*/
{
  double f_temp;

  f_temp = (1 / 293.15) - (1 / DegK);
  return (exp((-12530 * f_temp / 8.31441) - 19.111)); /* A-65 */
}

double K_MgOH2(double DegK) /* Solubility of magnesium hydroxide. */
/*
*  Reaction: [Mg++] + 2[H2O] => 2[H+] + Mg(OH)2 (Solid)
*  Returned: k_MgOH2 = [Mg++]/[H+][H+]  
*/
{
  double f_temp;

  f_temp = (1 / 293.15) - (1 / DegK);
  return (exp(38.78 + (-113960.0 * f_temp / 8.31441))); /* A-76 */
}

double K_MgOH(double DegK) /* Ionization of magnesium hydroxide. */
/*
*  Reaction:  [Mg++] + [H2O] <=> [MgOH-] + [H+]
*  Returned:  k_MgOH = [H+][MgOH+]/[Mg++]
*
*/
{
  return (exp(-65180.0 / (8.31441 * DegK))); /* A-43 */
}

double K_MgOH2aq(double DegK) /* Solubility of magnesium hydorxide */
/*
*
*  Reaction: [Mg++ ] + 2[H2O] <=> 2[H+] + [Mg(OH)2aq]
*  Returned: k_MgOH2aq = [Mg(OH)2aq][H+][H+]/[Mg++]
*/
{
  return (exp(-159760.0 / (8.31441 * DegK))); /* A-51 */
}
