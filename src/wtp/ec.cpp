/* ec.c */
#include "wtp.h"

void ec_comply(struct UnitProcess *unit)
/*
* Purpose:
*   Return an integer value indicating the status of compliance 
*   with enhanced coagulation (Step 1) requirements or SUVA
*   and or TOC exemptions
*
* Documentation and code by WJS 5/2001
*/

{
	double eff_toc;
	double raw_uva;
	double raw_toc;
	double raw_alk;
	double toc_rem_reqd;
	double toc_rem;
	double raw_tsuva;
	struct Effluent *eff; /* Effluent data  */

	eff = &unit->eff;

	eff_toc = eff->TOC;
	raw_toc = eff->influent->eff.TOC;
	raw_uva = eff->influent->eff.UV;
	raw_alk = eff->influent->eff.Alk;

	//First determine compliance status with Step 1
	if (raw_toc > 0.0)
		toc_rem = (raw_toc - eff_toc) / raw_toc * 100.0;
	else
		toc_rem = 0.0;

	if (raw_toc <= 2.0)
		eff->ec_meeting_step1 = TRUE;

	else if (raw_toc > 2.0 && raw_toc <= 4.0)
	{
		if (raw_alk <= 60.0)
			toc_rem_reqd = 35.0;
		else if (raw_alk > 60.0 && raw_alk <= 120.0)
			toc_rem_reqd = 25.0;
		else
			toc_rem_reqd = 15.0;

		if (toc_rem >= toc_rem_reqd)
			eff->ec_meeting_step1 = TRUE;
	}

	else if (raw_toc > 4.0 && raw_toc <= 8.0)
	{
		if (raw_alk <= 60.0)
			toc_rem_reqd = 45.0;
		else if (raw_alk > 60.0 && raw_alk <= 120.0)
			toc_rem_reqd = 35.0;
		else
			toc_rem_reqd = 25.0;

		if (toc_rem >= toc_rem_reqd)
			eff->ec_meeting_step1 = TRUE;
	}

	else //if (raw_toc > 8.0 && raw_toc <= 4.0)
	{
		if (raw_alk <= 60.0)
			toc_rem_reqd = 50.0;
		else if (raw_alk > 60.0 && raw_alk <= 120.0)
			toc_rem_reqd = 40.0;
		else
			toc_rem_reqd = 30.0;

		if (toc_rem >= toc_rem_reqd)
			eff->ec_meeting_step1 = TRUE;
	}

	//Next determine exemption status with respect to raw water TOC, raw water UVA, and
	//finished water TOC

	if (raw_toc > 0.0)
		raw_tsuva = raw_uva / raw_toc * 100.0;
	else
		raw_tsuva = 999999.9;

	if (eff_toc <= 2.0 || raw_toc <= 2.0 || raw_tsuva <= 2.0)
		eff->ec_exempt = TRUE;

} /* End of ec_comply() */