/* dbpslect.c */
#include "wtp.h"

void choose_dbpmodel(struct UnitProcess *unit)
/*
* Purpose:
*   To call the appropriate DBP formation subroutine based upon the settings
*   of the following global variables:
*                   	rwdbpflag
* 			owdbpflag
*			coagdbpflag
*		      	modrw1dbpflag
*			modrw2dbpflag
*			gacmemdbpflag
*
* Inputs and outputs:
*   The UnitProcess data structure is passed to supporting routines which
*   estimate DBP Formation.
*
* Notes:
*  1. choose_dbp() calls rwdbp(), owdbp (), coagdbp(), modrw1dbp(),
*		   modrw2dbp(), or gacmemdbp().
*  2. choose_dbp() is called from basn_dbp, filt_dbp, or dist_dbp
*  
*
* Documention and code by WJS, 10/98 */

{
	//FILE *fptr;

	/* Make sure that if this is a groundwater (swflag==FALSE), then we are using
   coagdbp() function */

	/* Call the appropriate support function to estimate DBP formation */
	if (swflag == FALSE)
	{ /* Make sure that if this is a groundwater (swflag==FALSE), then we are using
              coagulated water model for DBPs */
		coagdbp(unit);
	}

	else if (rwdbpflag == TRUE)
	{
		rwdbp(unit);
	}

	else if (owdbpflag == TRUE)
	{
		owdbp(unit);
	}

	else if (coagdbpflag == TRUE)
	{
		coagdbp(unit);
	}
	else if (gacmemdbpflag == TRUE)
	{
		gacmbdbp(unit);
	}

	else if (modrw1dbpflag == TRUE)
	{
		modrw1dbp(unit);
	}

	else if (modrw2dbpflag == TRUE && unit->type != RAPID_MIX)
	{
		modrw2dbp(unit);
	}

	else if (modrw2dbpflag == TRUE && unit->type == RAPID_MIX)
	{
		modrw1dbp(unit);
	}
	else /* The Do Nothing Alternative */
		;
}
