#include "wtp.h"

/* ---------------------------------------------------------
* Functions to adjust pH and TOC predictions
* of the model based on fitting of ICR AUX8 database.
* -----------------------------------------------------------*/

//******************************************************************

/* Define pH adjustment function */
double f_adj_pH(double pH)
{
  double Adj_pH;
  double m, b;

  if (softflag == TRUE)
  { //Use correction factors for Softening Plants (at all locations
    //in the plant - thus the use of 'softflag').  Based on analysis of ph
    // at unit processes with a volume, excluding plant-months with negative
    // finished water alkalinity for ICR plants

    b = 1.86;
    m = 0.712;

    Adj_pH = (pH - b) / m;

    // Correction may lead to unreasonable values in some cases - do not let
    // out of 1.0 to 14.0 range
    if (Adj_pH > 14.0)
      Adj_pH = 14.0;
    if (Adj_pH < 1.0)
      Adj_pH = 1.0;
  }
  else
  { //Correction factors for Conventional Plants

    Adj_pH = pH; //No adjustment for conventional plants
  };

  return (Adj_pH);
};

/* Define toc adjustment function */
double f_adj_toc(double toc, struct UnitProcess *unit)
{
  double Adj_toc;
  double m;

  if (unit->eff.limesoftening == TRUE)
  { //Use correction factors for Softening Plants - based on analysis of toc
    // at unit processes with a volume, excluding plant-months with negative
    // finished water alkalinity, using data points representing < 90th
    // percentile toc error (all for softening plants only)

    m = 0.874; /* Slope-only adjustment */

    Adj_toc = toc / m;
  }
  else
  { //Correction factors for Conventional Plants

    Adj_toc = toc; //No adjustment for conventional plants
  };

  return (Adj_toc);
};
