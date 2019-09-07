/* NCSTR.c */
#include "wtp.h"

short ncstr(double t_ten_theo)

/* -------------------------------------------------------------------- */
/*                                                                      */
/*  Determines the number of continuous flow, completely mixed          */
/*  reactors in series that will simulate the process unit calling      */
/*  this subroutine.                                                    */
/*                                                                      */
/*  ------------------------------------------------------------------  */

{
  short ncstrs;

  if (t_ten_theo <= 0.186)
    ncstrs = 1;
  if (t_ten_theo > 0.186 && t_ten_theo <= 0.317)
    ncstrs = 2;
  if (t_ten_theo > 0.317 && t_ten_theo <= 0.402)
    ncstrs = 3;
  if (t_ten_theo > 0.402 && t_ten_theo <= 0.461)
    ncstrs = 4;
  if (t_ten_theo > 0.461 && t_ten_theo <= 0.506)
    ncstrs = 5;
  if (t_ten_theo > 0.506 && t_ten_theo <= 0.540)
    ncstrs = 6;
  if (t_ten_theo > 0.540 && t_ten_theo <= 0.569)
    ncstrs = 7;
  if (t_ten_theo > 0.569 && t_ten_theo <= 0.593)
    ncstrs = 8;
  if (t_ten_theo > 0.593 && t_ten_theo <= 0.613)
    ncstrs = 9;
  if (t_ten_theo > 0.613 && t_ten_theo <= 0.630)
    ncstrs = 10;
  if (t_ten_theo > 0.630 && t_ten_theo <= 0.645)
    ncstrs = 11;
  if (t_ten_theo > 0.645 && t_ten_theo <= 0.659)
    ncstrs = 12;
  if (t_ten_theo > 0.659 && t_ten_theo <= 0.671)
    ncstrs = 13;
  if (t_ten_theo > 0.671 && t_ten_theo <= 0.682)
    ncstrs = 14;
  if (t_ten_theo > 0.682 && t_ten_theo <= 0.691)
    ncstrs = 15;
  if (t_ten_theo > 0.691 && t_ten_theo <= 0.700)
    ncstrs = 16;
  if (t_ten_theo > 0.700 && t_ten_theo <= 0.708)
    ncstrs = 17;
  if (t_ten_theo > 0.708 && t_ten_theo <= 0.716)
    ncstrs = 18;
  if (t_ten_theo > 0.716 && t_ten_theo <= 0.723)
    ncstrs = 19;
  if (t_ten_theo > 0.723 && t_ten_theo <= 0.729)
    ncstrs = 20;
  if (t_ten_theo > 0.729 && t_ten_theo <= 0.735)
    ncstrs = 21;
  if (t_ten_theo > 0.735 && t_ten_theo <= 0.741)
    ncstrs = 22;
  if (t_ten_theo > 0.741 && t_ten_theo <= 0.746)
    ncstrs = 23;
  if (t_ten_theo > 0.746 && t_ten_theo <= 0.751)
    ncstrs = 24;
  if (t_ten_theo > 0.751 && t_ten_theo <= 1.0)
    ncstrs = 25;

  /*
* Tests were run for ncstr=25, 100, and 1000 with no significant
* differences  observed in chlorine concentrations. 
*/

  return (ncstrs);
}
