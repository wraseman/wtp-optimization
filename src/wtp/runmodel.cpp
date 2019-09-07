/* RunModel.c --  September 16, 1993*/
#include "wtp.h"

int runmodel(struct ProcessTrain *train)
/*
*  Purpose: Run WTP Model. This function will update all UnitProcess data
*           elements in the process train.
*
*  Inputs:
*    *train = The process train controlling structure.
*
*  Return:
*    TRUE/FALSE for success/fail. 
*
*  Notes:
*    1. runmodel() deals with the internal process train structure of WTP
*       and calls the Model functions to estimate water quality parameters.
*       Note that run_wtp() and thm_wtp() deal with the format of the
*       printed table whereas runmodel() deals with the internal working of
*       the model.
*    2. runmodel() is called from run_wtp() and thm_wtp().
*    3. runmodel() is ANSI.
*
*  Michael D. Cummins
*    July 1993
*/
{
	//  FILE *fptr2;
	//  FILE *fptr;

	struct UnitProcess *unit;
	struct UnitProcess *prev = NULL;
	int success = TRUE;
	int sed_cntr = 0;			 /* counter for sed basins in process train      */
	int soft_cntr = 0;		 /* approx. counter for coagulant or lime soft. doses
                                        that are in unique softening stages in train */
	int coag_cntr = 0;		 /* counter for coag. add. pts. in process train */
	int gac_cntr = 0;			 /* counter for gac units in process train       */
	int rm_cntr = 0;			 /* counter for rapid mix units in process train */
	int o3_cntr = 0;			 /* counter for ozone app. pts. in process train */
	int nf_cntr = 0;			 /* counter for nanofilter units in process train  */
	int o3_chamber_cntr = 0; /* counter for O3 contactors in process train  */
	int cl2uvox = FALSE;		 /* Determines if uv has been reduced by cl2*/
	int modrw2cl2_cntr = 0;  /* counter for number of chlorine points
					      in process train after the RM when calc.
					      DBPs by the "modrw2dbp" method */
	int presedflag = FALSE;  //Presed with no coag in front - helps to deter-
									 //mine whether or not pre-ozonation exists

	/* Initialize Globals: I don't think these should be in globals.c - WJS, 11/98*/
	gw_virus_flag = FALSE;
	coagflag = FALSE;
	conv_filtflag = FALSE;
	filtflag = FALSE;
	filt2flag = FALSE;
	softflag = FALSE;
	soft2flag = FALSE;
	rwdbpflag = TRUE;
	owdbpflag = FALSE;
	modrw1dbpflag = FALSE;
	modrw2dbpflag = FALSE;
	coagdbpflag = FALSE;
	gacmemdbpflag = FALSE;
	floccflag = FALSE;
	sedflag = FALSE;
	sedconvflag = FALSE;
	gacflag = FALSE;
	mfufflag = FALSE;
	nfflag = FALSE;
	ssfflag = FALSE;
	defflag = FALSE;
	bagfflag = FALSE;
	cartfflag = FALSE;
	bankfflag = FALSE;
	granf2flag = FALSE;
	gac2flag = FALSE;
	mfuf2flag = FALSE;
	nf2flag = FALSE;
	ssf2flag = FALSE;
	def2flag = FALSE;
	bagf2flag = FALSE;
	cartf2flag = FALSE;
	cfeflag = FALSE;
	ifeflag = FALSE;
	uvflag = FALSE;
	lt2_wscp_flag = FALSE;
	o3flag = FALSE;
	clo2flag = FALSE;
	pre_o3flag = FALSE;
	int_o3flag = FALSE;
	post_o3flag = FALSE;
	dir_filtflag = FALSE;
	bio_filtflag = FALSE;
	lt2presedflag = FALSE;
	nonconv_discredit = 0.0;
	bin34_inactreqd = 0.0;
	modrw2dbptime = 0.0;
	modrw2dbpcl2 = 0.0;
	tot_dis_req_g = 0.0;
	tot_dis_req_v = 0.0;
	tot_dis_req_c = 0.0;
	tot_crypto_lr = 0.0;
	tot_giardia_lr = 0.0;
	tot_virus_lr = 0.0;

	/* SET PRIMARY GLOBAL FLAGS
   (and global variables related to pathogen log removal */

	for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
	{
		switch (unit->type)
		{
		case ALUM:
			if (unit->data.alum->dose > 0.0)
				coagflag = TRUE;
			if (soft_cntr > 0 && sedconvflag == TRUE)
				soft_cntr += 1;
			break;

		case FILTER:
			if (filtflag == TRUE && filt2flag == FALSE) //2ND STAGE FILTRATION - GRANULAR MEDIA
			{
				filt2flag = TRUE;
				granf2flag = TRUE;
				tot_crypto_lr += unit->data.filter->crypto_lr_2;
				unit->data.filter->filt_stage = 2;
			}
			if (filtflag == FALSE && coagflag == TRUE && floccflag == TRUE && sedconvflag == TRUE)
			{ //PRIMARY FILTRATION - CONVENTIONAL
				filtflag = TRUE;
				conv_filtflag = TRUE;
				tot_giardia_lr += unit->data.filter->giardia_lr_conv;
				tot_virus_lr += unit->data.filter->virus_lr_conv;
				tot_crypto_lr += unit->data.filter->crypto_lr_conv;
				if (unit->data.filter->cfe_turb_flag == TRUE)
					tot_crypto_lr += 0.5;
				if (unit->data.filter->ife_turb_flag == TRUE &&
					 unit->data.filter->cfe_turb_flag == TRUE)
					tot_crypto_lr += 0.5;
				if (unit->data.filter->ife_turb_flag == TRUE &&
					 unit->data.filter->cfe_turb_flag == FALSE)
					tot_crypto_lr += 1.0;
				unit->data.filter->filt_stage = 1;
			}
			if (filtflag == FALSE && coagflag == TRUE && (sedconvflag == FALSE || (sedconvflag == TRUE && floccflag == FALSE)))
			{ //PRIMARY FILTRATION - DIRECT
				filtflag = TRUE;
				dir_filtflag = TRUE;
				tot_giardia_lr += unit->data.filter->giardia_lr_df;
				tot_virus_lr += unit->data.filter->virus_lr_df;
				tot_crypto_lr += unit->data.filter->crypto_lr_df;
				if (unit->data.filter->cfe_turb_flag == TRUE)
					tot_crypto_lr += 0.5;
				if (unit->data.filter->ife_turb_flag == TRUE &&
					 unit->data.filter->cfe_turb_flag == TRUE)
					tot_crypto_lr += 0.5;
				if (unit->data.filter->ife_turb_flag == TRUE &&
					 unit->data.filter->cfe_turb_flag == FALSE)
					tot_crypto_lr += 1.0;
				unit->data.filter->filt_stage = 1;
			}
			break;

		case IRON:
			if (unit->data.iron->dose > 0.0)
				coagflag = TRUE;
			if (soft_cntr > 0 && sedconvflag == TRUE)
				soft_cntr += 1;
			break;

		case LIME:
			if (unit->data.lime->purpose == 'S')
			{
				coagflag = TRUE;
				softflag = TRUE;
				if (soft_cntr == 0)
					soft_cntr = 1;
				else if (sedconvflag == TRUE)
					soft_cntr += 1;
				else
					;
			}
			break;

		case SLOW_MIX:
			floccflag = TRUE;
			break;

		case SETTLING_BASIN:
			sedflag = TRUE;
			if (sedconvflag == TRUE && soft_cntr > 1) //2-STAGE PRECIP. SOFTENING
			{
				soft2flag = TRUE;
				tot_crypto_lr += 0.5; //This assumes that the user will put filters downstream
			}
			if (coagflag == TRUE && filtflag == FALSE)
				sedconvflag = TRUE;
			break;

		case GAC:
			gacflag = TRUE;
			if (filtflag == TRUE && filt2flag == FALSE) //2ND STAGE FILTRATION - GAC MEDIA
			{
				filt2flag = TRUE;
				gac2flag = TRUE;
				tot_crypto_lr += unit->data.gac->crypto_lr_2;
				unit->data.gac->filt_stage = 2;
			}
			break;

		case MFUF_UP:
			if (filtflag == TRUE && filt2flag == FALSE) //2ND STAGE FILTRATION - MFUF
			{
				filt2flag = TRUE;
				mfuf2flag = TRUE;
				tot_crypto_lr += unit->data.mfuf->crypto_lr_2;
				unit->data.mfuf->filt_stage = 2;
				nonconv_discredit += unit->data.mfuf->crypto_lr_2;
			}
			if (filtflag == FALSE)
			{ //PRIMARY FILTRATION - MFUF
				filtflag = TRUE;
				mfufflag = TRUE;
				tot_giardia_lr += unit->data.mfuf->giardia_lr;
				tot_virus_lr += unit->data.mfuf->virus_lr;
				tot_crypto_lr += unit->data.mfuf->crypto_lr_1;
				unit->data.mfuf->filt_stage = 1;
				nonconv_discredit += unit->data.mfuf->crypto_lr_1;
			}
			break;

		case NF_UP:
			if (filtflag == TRUE && filt2flag == FALSE && unit->data.nf->treat_fraction >= 1.0)
			{ //2ND STAGE FILTRATION - NF
				filt2flag = TRUE;
				nf2flag = TRUE;
				tot_crypto_lr += unit->data.nf->crypto_lr;
				unit->data.nf->filt_stage = 2;
				nonconv_discredit += unit->data.nf->crypto_lr;
			}
			if (filtflag == FALSE && unit->data.nf->treat_fraction >= 1.0)
			{ //PRIMARY FILTRATION - NF
				filtflag = TRUE;
				nfflag = TRUE;
				tot_giardia_lr += unit->data.nf->giardia_lr;
				tot_virus_lr += unit->data.nf->virus_lr;
				tot_crypto_lr += unit->data.nf->crypto_lr;
				unit->data.nf->filt_stage = 1;
				nonconv_discredit += unit->data.nf->crypto_lr;
			}
			break;

		case SLOW_FILTER:
			if (filtflag == TRUE && filt2flag == FALSE)
			{ //2ND STAGE FILTRATION - SLOW SAND FILTER
				filt2flag = TRUE;
				ssf2flag = TRUE;
				tot_crypto_lr += unit->data.ssf->crypto_lr_2;
				unit->data.ssf->filt_stage = 2;
			}
			if (filtflag == FALSE)
			{ //PRIMARY FILTRATION - SLOW SAND FILTER
				filtflag = TRUE;
				ssfflag = TRUE;
				tot_giardia_lr += unit->data.ssf->giardia_lr;
				tot_virus_lr += unit->data.ssf->virus_lr;
				tot_crypto_lr += unit->data.ssf->crypto_lr_1;
				unit->data.ssf->filt_stage = 1;
			}
			break;

		case DE_FILTER:
			if (filtflag == TRUE && filt2flag == FALSE)
			{ //2ND STAGE FILTRATION - DE FILTER
				filt2flag = TRUE;
				def2flag = TRUE;
				tot_crypto_lr += unit->data.def->crypto_lr_2;
				unit->data.def->filt_stage = 1;
			}
			if (filtflag == FALSE)
			{ //PRIMARY FILTRATION - DE FILTER
				filtflag = TRUE;
				defflag = TRUE;
				tot_giardia_lr += unit->data.def->giardia_lr;
				tot_virus_lr += unit->data.def->virus_lr;
				tot_crypto_lr += unit->data.def->crypto_lr_1;
				unit->data.def->filt_stage = 1;
			}
			break;

		case BAG_FILTER:
			if (filtflag == TRUE && filt2flag == FALSE)
			{ //2ND STAGE FILTRATION - BAG FILTER
				filt2flag = TRUE;
				bagf2flag = TRUE;
				tot_crypto_lr += unit->data.altf->crypto_lr_2;
				unit->data.altf->filt_stage = 2;
				nonconv_discredit += unit->data.altf->crypto_lr_2;
			}
			if (filtflag == FALSE)
			{ //PRIMARY FILTRATION - BAG FILTER
				filtflag = TRUE;
				bagfflag = TRUE;
				tot_giardia_lr += unit->data.altf->giardia_lr;
				tot_virus_lr += unit->data.altf->virus_lr;
				tot_crypto_lr += unit->data.altf->crypto_lr_1;
				unit->data.altf->filt_stage = 1;
				nonconv_discredit += unit->data.altf->crypto_lr_1;
			}
			break;

		case CART_FILTER:
			if (filtflag == TRUE && filt2flag == FALSE)
			{ //2ND STAGE FILTRATION - CARTRIDGE FILTER
				filt2flag = TRUE;
				cartf2flag = TRUE;
				tot_crypto_lr += unit->data.altf->crypto_lr_2;
				unit->data.altf->filt_stage = 2;
				nonconv_discredit += unit->data.altf->crypto_lr_2;
			}
			if (filtflag == FALSE)
			{ //PRIMARY FILTRATION - CARTRIDGE FILTER
				filtflag = TRUE;
				cartfflag = TRUE;
				tot_giardia_lr += unit->data.altf->giardia_lr;
				tot_virus_lr += unit->data.altf->virus_lr;
				tot_crypto_lr += unit->data.altf->crypto_lr_1;
				unit->data.altf->filt_stage = 1;
				nonconv_discredit += unit->data.altf->crypto_lr_1;
			}
			break;

		case BANK_FILTER:
			if (bankfflag == 0 && unit->data.bankf->eligible_lt2 == TRUE)
			{ //BANK FILTRATION TOOLBOX OPTION
				if (unit->data.bankf->distance >= 25.0 && unit->data.bankf->distance < 50.0)
				{
					bankfflag = 1;
					tot_crypto_lr += unit->data.bankf->crypto_lr_close;
					nonconv_discredit += unit->data.bankf->crypto_lr_close;
				}
				if (unit->data.bankf->distance >= 50.0)
				{
					bankfflag = 2;
					tot_crypto_lr += unit->data.bankf->crypto_lr_far;
					nonconv_discredit += unit->data.bankf->crypto_lr_far;
				}
			}
			break;

		case PRESED_BASIN:
			if (lt2presedflag == FALSE && coagflag == TRUE &&
				 filtflag == FALSE && unit->data.presed->eligible_lt2 == TRUE)
			{ //PRESED TOOLBOX OPTION
				lt2presedflag = TRUE;
				tot_crypto_lr += unit->data.presed->crypto_lr;
			}
			break;

		case UV_DIS:
			uvflag = TRUE;
			tot_giardia_lr += unit->data.uvdis->giardia_li;
			tot_virus_lr += unit->data.uvdis->virus_li;
			tot_crypto_lr += unit->data.uvdis->crypto_li;
			nonconv_discredit += unit->data.uvdis->crypto_li;
			break;

		case OZONE:
			o3flag = TRUE;
			pre_o3flag = TRUE; // Pre-ozonation is assumed true until location
			break;				 // of ozonation set in next loop */

		case CHLORINE_DIOXIDE:
			clo2flag = TRUE;
			break;

			/* No longer needed - WJS, 10/98
          case WTP_EFFLUENT:
	    break; */

		default:
			break;
		}

		/*DEBUGGING CODE*/
		/*      fptr2=fopen("debug.txt","a+");
      fprintf(fptr2,"Module:                %d\n",unit->type);
      fprintf(fptr2,"filtflag:              %d\n",filtflag);
      fprintf(fptr2,"dir_filtflag:          %d\n",dir_filtflag);
      fprintf(fptr2,"conv_filtflag:         %d\n",conv_filtflag);
      fprintf(fptr2,"sedconvflag:           %d\n",sedconvflag);
      fprintf(fptr2,"coagflag:              %d\n",coagflag);
      fprintf(fptr2,"floccflag:             %f\n",floccflag);
      fclose(fptr2);  */
		/*DEBUGGING CODE*/
	}

	// Some protection against giving 2nd stage softening credit with no filters in the train
	//  (not that this protects against there being 2nd stage softening after 1st stage filters
	//   and no 2nd stage filters in the train)
	if (soft2flag == TRUE && filtflag == FALSE)
		tot_crypto_lr -= 0.5;

	/* SET SECONDARY GLOBAL FLAGS (based on process order and Primary Global Flag values */

	/* This code determines if pre-sedimentation exists by checking if sedimentation
   has occurred before any coagulation */
	for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
	{
		if (unit->type == ALUM && unit->data.alum->dose > 0.0)
			coag_cntr += 1;
		if (unit->type == IRON && unit->data.iron->dose > 0.0)
			coag_cntr += 1;
		if (unit->type == LIME && unit->data.lime->purpose == 'S')
			coag_cntr += 1;
		if (unit->type == SETTLING_BASIN && coag_cntr == 0)
			presedflag = TRUE;
	} /* End pre_sedflag determination loop */

	/* The following code section sets the ozonation flags by seeing what processes
   precede ozonation - WJS 10/98 */
	if (o3flag == TRUE)
	{ /* ozone exists, so update flag values */
		if (presedflag == FALSE)
		{ /* presedimentation does not exist, so step thru train
	               			     to find a sed basin or filter before ozone*/
			for (unit = FirstUnitProcess(train); unit->type != OZONE; unit = NextUnitProcess(unit))
			{
				if (unit->type == SETTLING_BASIN)
				{ /* a sed basin precedes ozone, so int-ozone exists */
					pre_o3flag = FALSE;
					int_o3flag = TRUE;
				}
				if (unit->type == FILTER || unit->type == GAC)
				{ /* filtration (or GAC) precedes ozone,
		  			   so post-ozone exists */
					pre_o3flag = FALSE;
					int_o3flag = FALSE;
					post_o3flag = TRUE;
				}
			} /* END of FOR loop through process train */
		}
		else
		{ /* Pre-sedimentation exists, so look for more than one
	   				    sed basin, or a filter before ozone*/
			for (unit = FirstUnitProcess(train); unit->type != OZONE; unit = NextUnitProcess(unit))
			{
				if (unit->type == SETTLING_BASIN)
				{
					sed_cntr += 1;
					if (sed_cntr > 1)
					{ /* If only one basin precedes ozonation, then it must
					    be the pre-sed basin, so pre-ozonation exists.  If
					    more than one basin precedes ozonation, then
					    intermediate-ozonation might exist */
						pre_o3flag = FALSE;
						int_o3flag = TRUE;
					}
				}
				if (unit->type == FILTER || unit->type == GAC)
				{ /* If a filter (or GAC) precedes ozonation,
		                            then post-ozonation exists */
					pre_o3flag = FALSE;
					int_o3flag = FALSE;
					post_o3flag = TRUE;
				}
			} /* End of FOR loop through process train */

		} /* End of IF-ELSE for pre-sed condition */
	}
	/* End of IF for o3flag == TRUE */

	/*********************************************************************************************/
	/* THIS IS THE NEW MAIN LOOP TO RUN THE MODEL. - WJS, 10/98 */
	for (unit = FirstUnitProcess(train); unit; unit = NextUnitProcess(unit))
	{
		/* Copy Effluent data from previous unit process */
		if (prev != NULL)
		{
			unit->eff = prev->eff;
		}
		prev = unit;

		switch (unit->type) /* Compute effluent */
		{

		case INFLUENT:
			influent(unit);

			break;

		case RAPID_MIX:
			rm_cntr += 1; /* Update this process counter */

			/* Must save RM influent WQ data for use in modrw1dbp()
			   and modrw2dbp() */
			unit->eff.last_rm_inf = PrevUnitProcess(unit);

			/* Must know the equivalent alum dose at the point of the
			   first RM for use in modrw1dbp() and modrw2dpb() */
			if (rm_cntr == 1)
			{
				unit->eff.EquivAlumDose = unit->eff.AlumDose + unit->eff.FericDose / 270.0 * 297.0;
				/*270 = MW of FerricCl and 594 = MW of Alum (but, there are
				 two meq of metal per mol for Alum) */
			}

			/* For modrw2dbp(), need to know if chlorination has occurred before
			   the RM so that we know whether or not the TOC being chlorinated
			   at the post-RM chlorine add. pt. has already been chlorinated or not.
			   If it has, then we don't want to start out at the beginning of the
			   DBP vs. time curve in modrw2dbp() */
			if (unit->eff.cl2cnt > 0)
			{
				modrw2dbptime /*(hrs)*/ = (unit->data.basin->sb_mean) * (unit->data.basin->volume) / (unit->eff.Flow / 24.0);
			}
			else
			{
				modrw2dbptime = 0.0;
			}

			/* Once an RM is hit, these two DBP models definitely do not apply*/
			rwdbpflag = FALSE;
			owdbpflag = FALSE;

			/* modrw1dbp() model applies if chlorine has been added and the
			   next process is not chlorine, in which case modrw2dbp() applies.
			   The modrw2cl2cntr is needed to determine when the 2nd
			   chlorination pt after the RM has occurred, so that modrw2dbpflag
			   will become FALSE*/
			if (unit->eff.cl2cnt > 0 && NextUnitProcess(unit)->type != CHLORINE && NextUnitProcess(unit)->type != HYPOCHLORITE)
			{
				modrw1dbpflag = TRUE;
			}
			if (NextUnitProcess(unit)->type == CHLORINE ||
				 NextUnitProcess(unit)->type == HYPOCHLORITE)
			{
				modrw2dbpflag = TRUE;
				modrw2cl2_cntr = 0; /* initialize chlorine add. pt. counter */
			}

			/* If this is the first RM and no chlorine has been added, nor is
			   CHLORINE the following unit process, we need to be in CoagDBP mode hereafter */
			if (unit->eff.cl2cnt == 0 && rm_cntr == 1 && NextUnitProcess(unit)->type != CHLORINE && NextUnitProcess(unit)->type != HYPOCHLORITE)
			{
				coagdbpflag = TRUE;
			}

			/* If any of the following TOC removal processes have occurred, the coagulated water
			   DBP formation model will be in effect for remaining DBP calcs until GAC or
			   Membranes are hit. */
			if (bio_filtflag == TRUE || rm_cntr > 1 || o3_cntr > 0)
			{
				coagdbpflag = TRUE;
				modrw1dbpflag = FALSE;
				modrw2dbpflag = FALSE;
			}
			/*Of course, if GAC or NFUF occurred, that DBP model will be used*/
			if (gac_cntr > 0 || nf_cntr > 0)
			{
				coagdbpflag = FALSE;
				modrw1dbpflag = FALSE;
				modrw2dbpflag = FALSE;
			}

			/* TOC removal will now occur at the RM */
			if (unit->eff.limesoftening == TRUE)
				soft_rmv(unit);
			else if (unit->eff.AlumDose > 0.0)
				alum_rmv(unit);
			else if (unit->eff.FericDose > 0.0)
				fecl_rmv(unit);
			else /* Do nothing - no coagulants have been added */
				;

			//basn_rmv(unit);

			/* Calculate DBPs (now that proper flag is set) */
			basn_dbp(unit);

			if (NextUnitProcess(unit)->type == CHLORINE)
				modrw2dbpcl2 = NextUnitProcess(unit)->data.chemical->chlor + (unit->eff.FreeCl2 + unit->eff.NH2Cl) * MW_Cl2;

			if (NextUnitProcess(unit)->type == HYPOCHLORITE)
				modrw2dbpcl2 = NextUnitProcess(unit)->data.chemical->naocl + (unit->eff.FreeCl2 + unit->eff.NH2Cl) * MW_Cl2;

			/* If there is true prechlorination (simultaneous addition of chlorine and coagulant
			   as in Miguel Arias' study) do the following */
			if (unit->eff.pre_chlor_flag == FALSE && rwdbpflag == FALSE &&
				 (unit->eff.pre_chlor_dose_track > 0 || NextUnitProcess(unit)->type == CHLORINE || NextUnitProcess(unit)->type == HYPOCHLORITE))
			{
				// Set this flag - used in implementing pre-/re-chlor adjustment factor
				unit->eff.pre_chlor_flag = TRUE;

				// Calc pre-Cl2 dose to raw TOC ratio here
				unit->eff.pre_chlor_dose_ratio = unit->eff.pre_chlor_dose_track;
				if (NextUnitProcess(unit)->type == CHLORINE)
					unit->eff.pre_chlor_dose_ratio += NextUnitProcess(unit)->data.chemical->chlor;
				if (NextUnitProcess(unit)->type == HYPOCHLORITE)
					unit->eff.pre_chlor_dose_ratio += NextUnitProcess(unit)->data.chemical->naocl;
				if (unit->data.influent->toc > 0.0)
					unit->eff.pre_chlor_dose_ratio /= unit->data.influent->toc;
				else													 // raw TOC = 0.0
					unit->eff.pre_chlor_dose_ratio = 999.0; //high number - doesn't matter b/c no DBPs will form
			}															 // end "if(unit->eff.pre_chlor_flag == FALSE &&..."

			break;

		case GAC:
			gac_cntr += 1; /* Update this process counter */

			/* GAC requires that the gac/mem DBP formation
			   model be in effect throughout the rest of the process train */
			gacmemdbpflag = TRUE;
			coagdbpflag = FALSE;
			rwdbpflag = FALSE;
			owdbpflag = FALSE;
			modrw1dbpflag = FALSE;
			modrw2dbpflag = FALSE;

			/*Set ozone residual to zero*/
			unit->eff.o3_res = 0.0;

			/* No DBP formation will be considered in the GAC bed itself
			   since disinfect residuals go to zero */

			gac_rmv(unit);

			break;

		case NF_UP:
			nf_cntr += 1; /* Update this process counter */

			/* NF_UP requires that the gac/mem DBP formation
			   model be in effect throughout the rest of the process train */
			gacmemdbpflag = TRUE;
			coagdbpflag = FALSE;
			rwdbpflag = FALSE;
			owdbpflag = FALSE;
			modrw1dbpflag = FALSE;
			modrw2dbpflag = FALSE;

			nf_rmv(unit);

			break;

		case OZONE:
			o3_cntr += 1; /* Update this process counter */

			o3add(unit);

			/* When an ozone unit occurs in the process train, DBPs will only be
			   calculated with owdbp() if it is ozonation of water in which TOC
			   has not yet been removed (consistent with the experimental conditions
			   used to develop the algorithms*/
			rwdbpflag = FALSE;
			if (rm_cntr == 0 && nf_cntr == 0 && gac_cntr == 0 && bio_filtflag == FALSE)
			{ /* If no TOC removal yet */
				owdbpflag = TRUE;
			}
			if (modrw1dbpflag == TRUE || modrw2dbpflag == TRUE)
			{
				modrw1dbpflag = FALSE;
				modrw2dbpflag = FALSE;
				coagdbpflag = TRUE;
			}

			/* ELSE, coagdbpflag or gacmemdbpflag will 
			   remain TRUE and owdbpflag will remain FALSE */

			break;

		case FILTER:
			if (gac_cntr == 0 && nf_cntr == 0 && rm_cntr > 0)
				coagdbpflag = TRUE;
			/* ELSE, gacmemdbpflag remains TRUE and coagdbpflag FALSE*/

			/* If water has been ozonated and cl2 residual (Free + Combined) is
			   below the semi-arbitrary 0.1 mg/L threshhold, then filter is
			   bio-active, and TOC removal will be calculated with biofilt_rmv()
			   and all further processes will use coagdbpflag */

			if (o3_cntr > 0 && (unit->eff.FreeCl2 * MW_Cl2 + unit->eff.NH2Cl * MW_Cl2) < 0.1
				 /*&&  unit->eff.clo2_res          < 0.1*/
				 /*	&&  unit->eff.o3_res            < 0.1*/
				 && unit->data.filter->cl2_bkwsh == FALSE)
			{
				bio_filtflag = TRUE;
				biofilt_rmv(unit);
			}
			/* If it is a biofilter, then Cl2 residual must be low enough
			       that DBP formation is low enough so filt_dbp() does not need
			       to be called; also, disinfectant residuals are set to zero
			       in 'biofilt_rmv()' */
			else
			{
				filt_dbp(unit);
			}

			solids_rmv(unit);
			break;

		case SLOW_FILTER:
			if (gac_cntr == 0 && nf_cntr == 0)
				coagdbpflag = TRUE;
			/* ELSE, gacmemdbpflag remains TRUE and coagdbpflag FALSE*/

			/* If disinfectant residuals are below the semi-arbitrary 0.1 mg/L
			   threshhold, then the slow sand filter is bio-active, and TOC
			   removal will be calculated with biofilt_rmv()
			   and all further processes will use coagdbpflag */

			if ((unit->eff.FreeCl2 * MW_Cl2 + unit->eff.NH2Cl * MW_Cl2) < 0.1
				 /*&&  unit->eff.clo2_res     < 0.1*/
				 && unit->eff.o3_res < 0.1)
			{
				bio_filtflag = TRUE; //just to signal TOC removal for purpose
											//of selecting proper DBP routine
				biofilt_rmv(unit);
			}

			/* If it is a biofilter, then Cl2 residual must be low enough
			       that DBP formation is low enough so filt_dbp() does not need
			       to be called; also, disinfectant residuals are set to zero
			       in 'biofilt_rmv()' */
			else
			{
				filt_dbp(unit);
			}

			solids_rmv(unit);
			break;

		case MFUF_UP:
		case DE_FILTER:
		case BAG_FILTER:
		case CART_FILTER:
			if (gac_cntr == 0 && nf_cntr == 0 && rm_cntr > 0)
				coagdbpflag = TRUE;
			/* ELSE, gacmemdbpflag remains TRUE and coagdbpflag FALSE*/

			if (unit->type == DE_FILTER)
				filt_dbp(unit);

			if (unit->type == MFUF_UP)
				mfuf_rmv(unit);

			solids_rmv(unit);
			break;

			/* These do not change the DBP model in effect, nor do they
                           remove TOC, they just form DBPs */
		case BASIN:
		case SLOW_MIX:
		case SETTLING_BASIN:
		case CONTACT_TANK:
		case CLEARWELL:
		case PRESED_BASIN:

			basn_dbp(unit);

			/*DEBUGGING CODE*/
			//      fptr2=fopen("debug.txt","a+");
			//      fprintf(fptr2,"Module: %d  after  'if'  \n",unit->type);
			//      fprintf(fptr2,"Free_Cl2:              %f\n",unit->eff.FreeCl2);
			//      fclose(fptr2);
			/*DEBUGGING CODE*/

			if (unit->type == SETTLING_BASIN || unit->type == PRESED_BASIN)
			{
				solids_rmv(unit);
			}

			break;

		case O3_CONTACTOR:
			/* Take care of setting the first chamber variable to make
                           sure that CT is not calculated for it when "basn_dbp()" is called */
			o3_chamber_cntr++;
			if (o3_chamber_cntr == 1)
				unit->eff.first_o3_chamber = TRUE;

			/* Calc. o3_residual, bromate, bromide as necessary*/
			if (unit->eff.last_o3_inf != NULL)
				ozonate(unit);

			basn_dbp(unit);

			/* Once CT has been calculated, this should be put back to FALSE*/
			unit->eff.first_o3_chamber = FALSE;

			break;

		case CHLORINE:
		case HYPOCHLORITE:
			/*The purpose of the pre_re_chlor_flag is to track the condition
			  where we are 1st time rechlorinating water that was prechlorinated*/
			/* We don't want to apply the adjustment for follow-on re_chlor doses*/
			if (unit->eff.pre_re_chlor_flag == 1)
				unit->eff.pre_re_chlor_flag = 2;

			/* This is to cover the case where there was pre-chlorination, but we
			   are rechlorinating late enough in the train that the mdrwdbpflags
			   are no longer TRUE */
			if (coagdbpflag == TRUE && unit->eff.pre_chlor_flag == TRUE &&
				 unit->eff.pre_chlor_dose_ratio >= 0.2 && unit->eff.pre_re_chlor_flag == 0)
			{
				unit->eff.pre_re_chlor_flag = 1;
				unit->eff.pre_chlor_TTHM = unit->eff.TTHM;
				unit->eff.pre_chlor_HAA6 = unit->eff.HAA6;
			}

			/* When chlorination occurs beyond what the modrw1dbp and modrw2dbp
			   models inherently consider, the coagulated water DBP formation
			   model takes effect */
			if (modrw1dbpflag == TRUE)
			{
				modrw1dbpflag = FALSE;
				coagdbpflag = TRUE;
				if (unit->eff.pre_chlor_dose_ratio >= 0.2)
				{
					unit->eff.pre_re_chlor_flag = 1;
					unit->eff.pre_chlor_TTHM = unit->eff.TTHM;
					unit->eff.pre_chlor_HAA6 = unit->eff.HAA6;
				}
			}
			if (modrw2dbpflag == TRUE)
			{
				modrw2cl2_cntr += 1;
				if (modrw2cl2_cntr > 1)
				{
					modrw2dbpflag = FALSE;
					coagdbpflag = TRUE;
					if (unit->eff.pre_chlor_dose_ratio >= 0.2)
					{
						unit->eff.pre_re_chlor_flag = 1;
						unit->eff.pre_chlor_TTHM = unit->eff.TTHM;
						unit->eff.pre_chlor_HAA6 = unit->eff.HAA6;
					}
				}
			}
			chloradd(unit);

			/*For internal calculation purposes, we only want to decrease
			  UV by chlorination if it is at least the second point of
			  chlorination and an RM has already occurred because, DOC
			  removal by coag. is not typically negatively affected by
			  pre-chlorination, and all DBP models inherently consider
			  UV reduction by chlorine (they are based on unchlorinated
			  UV); The only time we want UV to go down is for multiple
			  chlorination points after an RM */
			if (rm_cntr > 0 && unit->eff.cl2cnt > 1 && cl2uvox == FALSE)
			{
				unit->eff.UV *= 0.7;
				cl2uvox = TRUE;
			}
			/*For output purposes, we want UV term to decrease at first
			  point of chlorination only*/
			if (unit->eff.cl2cnt == 1)
				unit->eff.UV_out *= 0.7;

			break;

			/* Chemical addition: */
		case ALUM:
			alumadd(unit);
			break;
		case IRON:
			fecladd(unit);
			break;
		case SULFURIC_ACID:
			h2so4add(unit);
			break;
		case LIME:
			limeadd(unit);
			break;
		case SODA_ASH:
			soda_add(unit);
			break;
		case AMMONIA:
			nh3add(unit);
			break;
		case AMMONIUM_SULFATE:
			nh3add(unit);
			break;
		case PERMANGANATE:
			kmno4add(unit);
			break;
		case CARBON_DIOXIDE:
			co2_add(unit);
			break;
		case SODIUM_HYDROXIDE:
			naohadd(unit);
			break;
		case SULFUR_DIOXIDE:
			so2_add(unit);
			break;
		case CHLORINE_DIOXIDE:
			clo2add(unit);
			break;

			/* Sample points: */
		case WTP_EFFLUENT:

			/*Set ozone residual to zero*/
			unit->eff.o3_res = 0.0;

			/* Have to save the WTP_Effluent WQ data for distribution system
	       samples to use */
			unit->eff.wtp_effluent = unit;

			/* Take care of EC status */
			ec_comply(unit);

			break;

		case AVG_TAP:
		case LOCATION_1:
		case END_OF_SYSTEM:

			dist_dbp(unit);
			break;

		case UV_DIS:
		case BANK_FILTER: //nothing to do with these here
			break;

		default:
#if DEBUG == TRUE
			printf("Unknown unit process %d in runmodel()\n", unit->type);
#endif
			success = FALSE;
			break;

		} /* end switch(unit->type) ** Compute effluent */

		/* Update dbpmodel variable in Unit Process Effluent data structure
         (for Debugging only)*/
		/*      if     ( rwdbpflag     == TRUE )
  	      unit->eff.dbpmodel = "Raw Water Model";
      else if( owdbpflag     == TRUE )
              unit->eff.dbpmodel = "Ozonated Water Model";
      else if( coagdbpflag   == TRUE )
	      unit->eff.dbpmodel = "Coagulated Water Model";
      else if( gacmemdbpflag   == TRUE )
  	      unit->eff.dbpmodel = "GAC/Membrane-Treated Water Model";
      else if( modrw1dbpflag == TRUE )
	      unit->eff.dbpmodel = "Raw Water Model (PRE-RM chlorination modified)";
      else if( modrw2dbpflag == TRUE && unit->type != RAPID_MIX )
  	      unit->eff.dbpmodel = "Raw Water Model (POST-RM chlorination modified)";
      else if( modrw2dbpflag == TRUE && unit->type == RAPID_MIX )
              unit->eff.dbpmodel = "Raw Water Model (PRE-RM chlorination modified)";
      else;  */

		breakpt(unit);

		phchange(unit, TRUE);

		//Update SUVA here by definition (actually TSUVA) ; SUVA variable in
		//Effluent data structure not used in any calcs.
		if (unit->eff.TOC > 0.0)
			unit->eff.SUVA = unit->eff.UV_out / unit->eff.TOC * 100.0;
		else
			unit->eff.SUVA = 999999.9;

		/* Update residence time variables in Unit Process Effluent data structure */
		res_time(unit);

		// Update this variable.  It gets zeroed out if any res time has elapsed...
		if (unit->eff.processtime > 0.0)
			unit->eff.pre_chlor_dose_track = 0.0;

	} /* end for(unit=...) loop */

	return (success);
}
