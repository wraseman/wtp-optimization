/* res_time.c */
#include "wtp.h"

void res_time(struct UnitProcess *unit)
/*
* Purpose:
*   Update process residence time and cumulative process train residence time
*   variables in Unit Process data structure
*/
{
	register struct Effluent *eff;

	eff = &unit->eff;

	double theo_res_time;	 /* This is the EBCT of the unit process basin (calc.- hrs.)   */
	double mean_theo;		  /* This is the ratio of the mean residence time to the
			       theoretical res. time                                      */
	double peak_factor = 1.0; /* Used when coldflag == TRUE to implement peak flow scenario */

	if (coldflag == TRUE)
		peak_factor = eff->Peak / eff->influent->data.influent->avg_flow;

	if (unit->type == GAC)
	{
		if (unit->data.gac->config == 'S' && coldflag == TRUE)
			eff->processtime = unit->data.gac->ebct / 60.0 / peak_factor;
		else
			eff->processtime = unit->data.gac->ebct / 60.0;
	}

	else if (unit->type == BASIN ||
			 unit->type == O3_CONTACTOR ||
			 unit->type == RAPID_MIX ||
			 unit->type == SLOW_MIX ||
			 unit->type == SETTLING_BASIN ||
			 unit->type == CONTACT_TANK ||
			 unit->type == CLEARWELL)
	{
		theo_res_time = unit->data.basin->volume / eff->Flow * 24.0;
		mean_theo = unit->data.basin->sb_mean;
		eff->processtime = theo_res_time * mean_theo;
	}

	else if (unit->type == FILTER)
	{
		theo_res_time = unit->data.filter->volume / eff->Flow * 24.0;
		mean_theo = unit->data.filter->fi_mean;
		eff->processtime = theo_res_time * mean_theo;
	}

	else if (unit->type == AVG_TAP ||
			 unit->type == LOCATION_1)
	{
		eff->processtime = unit->data.avg_tap->days * 24.0 / peak_factor;
	}

	else if (unit->type == END_OF_SYSTEM)
	{
		eff->processtime = unit->data.end_of_system->days * 24.0 / peak_factor;
	}

	else
	{
		eff->processtime = 0.0;
	}

	eff->traintime += eff->processtime;
}
