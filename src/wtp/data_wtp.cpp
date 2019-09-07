/* Data_WTP.c -- September 17, 1993
*
*  Purpose:
*    This file contains read() and write() functions for each unit process.
*    The functions are accessed via;
*      UnitProcessTable[unit->type].read(buffer,i,unit)
*    and
*      UnitProcessTable[unit->type].write(buffer,i,unit)
*
*    These read() and write() functions have two intended purposed.  One
*    purpose is disk IO specificaly, saving the process train to a data
*    file and reading the data file into a unit process data structure.
*    This is implemented in read_wtp() and write_wtp().
*
*    The second purpose is to support user interface data entry screens.
*    The read() functions transfer ASCII data from a buffer into the unit
*    process data structure while the write() functions transfer data from
*    the unit process data structure into a buffer.  The user interface must
*    display and allow the user to edit the ASCII buffer.
*
*    For both disk IO and user interface purposes, the read() and write()
*    functions obtain the format from;
*      UnitProcessTable[unit->type].data[i].fmtin for the read()
*    and ...fmtout for write().  The fmtin can be used to indicates the data
*    type (%lf==double, %d==integer, %c==character, %s==string) thus allowing
*    the user interface code to generate an approate edit control.  The fmtout
*    will round the value to the approate number of significant places.
*
*    Error detection is weak in the current implementation and needs to be
*    re-thought.  At the present error detection is limited to read()
*    insuring that 'buffer' can be converted into the proper data type.
*    
*  Inputs/Outputs:
*    buffer: A conventional null terminated string.  For the read() function,
*            buffer should contain a value to be transfered into the unit
*            process data structure.  For the write() function, data is
*            transfered from the unit process data structure to buffer.
*            The string in buffer will be terminate with a '\0'  and not
*            contain a '\n'.  The the calling routine must insure that buffer
*            is adequately dimensioned, char buffer[20] is satifactory.
*    i     : Index to data element in the unit process data structure.  The
*            UnitProcessTable[].data[i]... implements an array accessed to
*            each data element.
*    unit  : Pointer to unit process data structure.
*
*  Returned: Error code.
*    In the present implementation, an fprintf() error code is returned.
*    Negative values indicate an error, zero and postives values indicate
*    no error.
*
*  Planded changes:
*   1.  The read() and write() functions will reorganized into files for each
*       unit process.
*   2.  An error and warning convention needs to be developed for checking
*       acceptable limits.
*   3.  Rethink fmtin.  The Information Collection Requirement (ICR) will
*       need a varity of data types.
*
*  Michale D. Cummins
*  September 20, 1993
*/
#include "wtp.h"

int write_influent(register char *buffer, short i, struct UnitProcess *unit)
/*
*  Purpose: Format data element 'i' and put in 'buffer'.  On return 
*           buffer is a C conventional null terminated string without
*           a '\n'.  The format is specified in UnitProcessTable.
*
*  Retrun:
*    error code conforming to fprintf() convention. Negative value
*    indicates an error, postive value indicates no error.
*
*  Notes:
*   1. The sprintf() function returns the number of characters written.
*      A return value of zero from sprintf() indicates an error.  The
*      sprintf() return value is checked and if zero then this function
*      returns a negative value to comply with the fprintf() return value
*      convention. See writewtp() for more information of fprintf()
*      return value.
*/
{
  register int e;
  register const char *fmt;
  register struct Influent *inf;

  inf = unit->data.influent;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, inf->pH);
    break;
  case 1:
    e = sprintf(buffer, fmt, inf->temp);
    break;
  case 2:
    e = sprintf(buffer, fmt, inf->low_temp);
    break;
  case 3:
    e = sprintf(buffer, fmt, inf->toc);
    break;
  case 4:
    e = sprintf(buffer, fmt, inf->uv254);
    break;
  case 5:
    e = sprintf(buffer, fmt, inf->bromide);
    break;
  case 6:
    e = sprintf(buffer, fmt, inf->alkalinity);
    break;
  case 7:
    e = sprintf(buffer, fmt, inf->calcium);
    break;
  case 8:
    e = sprintf(buffer, fmt, inf->hardness);
    break;
  case 9:
    e = sprintf(buffer, fmt, inf->nh3);
    break;
  case 10:
    e = sprintf(buffer, fmt, inf->ntu);
    break;
    //   case 11: e = sprintf( buffer, fmt, inf->crypto_req); break;
    //   case 12: e = sprintf( buffer, fmt, inf->clo2_crypto_ct_mult ); break;
  case 11:
    e = sprintf(buffer, fmt, inf->peak_flow);
    break;
  case 12:
    e = sprintf(buffer, fmt, inf->avg_flow);
    break;
  case 13:
    if (inf->swflag == TRUE)
      e = sprintf(buffer, fmt, "TRUE");
    else
      e = sprintf(buffer, fmt, "FALSE");
    break;
  case 14:
    e = sprintf(buffer, fmt, inf->crypto_conc);
    break;
  case 15:
    if (inf->lt2_wscp_flag == TRUE)
      e = sprintf(buffer, fmt, "TRUE");
    else
      e = sprintf(buffer, fmt, "FALSE");
    break;
  case 16:
    if (inf->gw_virus_flag == TRUE)
      e = sprintf(buffer, fmt, "TRUE");
    else
      e = sprintf(buffer, fmt, "FALSE");
    break;
  case 17:
    e = sprintf(buffer, fmt, inf->gw_tot_dis_req_v);
    break;
    //    case 16: e = sprintf(buffer, inf->run_name); break;
  default:
    e = -1;
    break;
  }

  if (e <= 0)
    e = -1; /* An error occoured... change 'e' fprintf() error */
  return (e);
}

int read_influent(register char *buffer, short i, struct UnitProcess *unit)
{
  //FILE *fptr2;
  register int e;
  register const char *fmt;
  register struct Influent *inf;
  char test[80];
  //  char			    run_name[128];

  inf = unit->data.influent; /* Read in default values for influent */
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &inf->pH);
    break;
  case 1:
    e = sscanf(buffer, fmt, &inf->temp);
    break;
  case 2:
    e = sscanf(buffer, fmt, &inf->low_temp);
    break;
  case 3:
    e = sscanf(buffer, fmt, &inf->toc);
    break;
  case 4:
    e = sscanf(buffer, fmt, &inf->uv254);
    break;
  case 5:
    e = sscanf(buffer, fmt, &inf->bromide);
    break;
  case 6:
    e = sscanf(buffer, fmt, &inf->alkalinity);
    break;
  case 7:
    e = sscanf(buffer, fmt, &inf->calcium);
    break;
  case 8:
    e = sscanf(buffer, fmt, &inf->hardness);
    break;
  case 9:
    e = sscanf(buffer, fmt, &inf->nh3);
    break;
  case 10:
    e = sscanf(buffer, fmt, &inf->ntu);
    break;
    //  case 11: e = sscanf( buffer, fmt, &inf->crypto_req); break;
    //  case 12: e = sscanf( buffer, fmt, &inf->clo2_crypto_ct_mult); break;
  case 11:
    e = sscanf(buffer, fmt, &inf->peak_flow);
    break;
  case 12:
    e = sscanf(buffer, fmt, &inf->avg_flow);
    break;
  case 13:
    strncpy(test, buffer, sizeof(test));
    strupr(test);
    if (strstr(test, "T") != NULL || strstr(test, "Y") != NULL)
      inf->swflag = TRUE;
    else
      inf->swflag = FALSE;
    break;
  case 14:
    e = sscanf(buffer, fmt, &inf->crypto_conc);
    break;
  case 15:
    strncpy(test, buffer, sizeof(test));
    strupr(test);
    if (strstr(test, "T") != NULL || strstr(test, "Y") != NULL)
      inf->lt2_wscp_flag = TRUE;
    else
      inf->lt2_wscp_flag = FALSE;
    break;
  case 16:
    strncpy(test, buffer, sizeof(test));
    strupr(test);
    if (strstr(test, "T") != NULL || strstr(test, "Y") != NULL)
      inf->gw_virus_flag = TRUE;
    else
      inf->gw_virus_flag = FALSE;
    break;
  case 17:
    e = sscanf(buffer, fmt, &inf->gw_tot_dis_req_v);
    break;

    /*    case 16: 	strncpy(run_name,buffer,sizeof(run_name));
		inf->run_name = run_name;
                break;  */
  default:
    e = -1;
    break;
  }

  /*DEBUGGING CODE*/
  //      fptr2=fopen("debug.txt","a+");
  //      fprintf(fptr2,"%s",inf->run_name);
  //      fprintf(fptr2,"Module: %d  after  'if'  \n",unit->type);
  //      fprintf(fptr2,"Free_Cl2:              %f\n",unit->eff.FreeCl2);
  //      fclose(fptr2);
  /*DEBUGGING CODE*/

  //Input range limit checking
  if (inf->pH < 0.0)
    inf->pH = 0.0;
  if (inf->pH > 14.0)
    inf->pH = 14.0;
  if (inf->temp < 0.0)
    inf->temp = 0.0;
  if (inf->low_temp < 0.0)
    inf->low_temp = 0.0;
  if (inf->toc < 0.1)
    inf->toc = 0.1;
  if (inf->uv254 < 0.0)
    inf->uv254 = 0.0;
  if (inf->bromide < 0.0)
    inf->bromide = 0.0;
  if (inf->alkalinity < 0.0)
    inf->alkalinity = 0.0;
  if (inf->calcium < 0.0)
    inf->calcium = 0.0;
  if (inf->hardness < 0.0)
    inf->hardness = 0.0;
  if (inf->nh3 < 0.0)
    inf->nh3 = 0.0;
  if (inf->ntu < 0.0)
    inf->ntu = 0.0;
  //  if(inf->crypto_req < 0.0) inf->crypto_req = 0.0;
  //  if(inf->clo2_crypto_ct_mult < 0.0) inf->clo2_crypto_ct_mult = 0.0;
  if (inf->peak_flow < 0.001)
    inf->peak_flow = 0.001;
  if (inf->avg_flow < 0.001)
    inf->avg_flow = 0.001;
  if (inf->crypto_conc < 0.0)
    inf->crypto_conc = 0.0;
  if (inf->gw_tot_dis_req_v < 0.0)
    inf->gw_tot_dis_req_v = 0.0;

  if (e == 0)
    e = -1; /* convert to fprintf() error convention */
  return (e);
}

int write_alum(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Alum *alum;

  alum = unit->data.alum;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, alum->dose);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_alum(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Alum *alum;

  alum = unit->data.alum;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &alum->dose);
    break;
  default:
    e = -1;
    break;
  }

  //Input Check
  if (alum->dose < 0.0)
    alum->dose = 0.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_gac(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Gac *gac;

  gac = unit->data.gac;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, gac->ebct);
    break;
  case 1:
    e = sprintf(buffer, fmt, gac->regen);
    break;
  case 2:
    switch (gac->config)
    {
    case 'S':
      e = sprintf(buffer, "Single");
      break;
    case 'B':
      e = sprintf(buffer, "Blended");
      break;
    default:
      e = sprintf(buffer, "Error");
      break;
    }
    break;
  case 3:
    switch (gac->toc_calc)
    {
    case 'M':
      e = sprintf(buffer, "Max_TOC");
      break;
    case 'A':
      e = sprintf(buffer, "Avg_TOC");
      break;
    default:
      e = sprintf(buffer, "Error");
      break;
    }
    break;
  case 4:
    e = sprintf(buffer, fmt, gac->crypto_lr_2);
    break;

  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_gac(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Gac *gac;
  char test[80];

  gac = unit->data.gac;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &gac->ebct);
    break;
  case 1:
    e = sscanf(buffer, fmt, &gac->regen);
    break;
  case 2:
    strncpy(test, buffer, sizeof(test));
    strupr(test);
    if (strchr(test, 'S'))
    {
      gac->config = 'S';
      e = 1;
    }
    else if (strchr(test, 'B'))
    {
      gac->config = 'B';
      e = 1;
    }
    else
    {
      e = -1;
    }
    break;
  case 3:
    strncpy(test, buffer, sizeof(test));
    strupr(test);
    if (strchr(test, 'M'))
    {
      gac->toc_calc = 'M';
      e = 1;
    }
    else if (strchr(test, 'A'))
    {
      gac->toc_calc = 'A';
      e = 1;
    }
    else
    {
      e = -1;
    }
    break;
  case 4:
    e = sscanf(buffer, fmt, &gac->crypto_lr_2);
    break;

  default:
    e = -1;
    break;
  }

  //Input Check
  if (gac->ebct < 0.0)
    gac->ebct = 0.0;
  if (gac->regen < 0.0)
    gac->regen = 0.0;
  if (gac->crypto_lr_2 < 0.0)
    gac->crypto_lr_2 = 0.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_filter(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Filter *filter;

  filter = unit->data.filter;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, filter->volume);
    break;
  case 1:
    e = sprintf(buffer, fmt, filter->fi_mean);
    break;
  case 2:
    e = sprintf(buffer, fmt, filter->fi_t_ten);
    break;
  case 3:
    if (filter->cl2_bkwsh == TRUE)
      e = sprintf(buffer, fmt, "TRUE");
    else
      e = sprintf(buffer, fmt, "FALSE");
    break;
  case 4:
    switch (filter->media)
    {
    case 'S':
      e = sprintf(buffer, "A/S");
      break;
    case 'G':
      e = sprintf(buffer, "GAC");
      break;
    default:
      e = sprintf(buffer, "Error");
      break;
    }
    break;
  case 5:
    e = sprintf(buffer, fmt, filter->giardia_lr_conv);
    break;
  case 6:
    e = sprintf(buffer, fmt, filter->virus_lr_conv);
    break;
  case 7:
    e = sprintf(buffer, fmt, filter->crypto_lr_conv);
    break;
  case 8:
    e = sprintf(buffer, fmt, filter->giardia_lr_df);
    break;
  case 9:
    e = sprintf(buffer, fmt, filter->virus_lr_df);
    break;
  case 10:
    e = sprintf(buffer, fmt, filter->crypto_lr_df);
    break;
  case 11:
    if (filter->cfe_turb_flag == TRUE)
      e = sprintf(buffer, fmt, "TRUE");
    else
      e = sprintf(buffer, fmt, "FALSE");
    break;
  case 12:
    if (filter->ife_turb_flag == TRUE)
      e = sprintf(buffer, fmt, "TRUE");
    else
      e = sprintf(buffer, fmt, "FALSE");
    break;
  case 13:
    e = sprintf(buffer, fmt, filter->crypto_lr_2);
    break;

  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_filter(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Filter *filter;
  char test[80];

  filter = unit->data.filter;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &filter->volume);
    break;
  case 1:
    e = sscanf(buffer, fmt, &filter->fi_mean);
    break;
  case 2:
    e = sscanf(buffer, fmt, &filter->fi_t_ten);
    break;
  case 3:
    strncpy(test, buffer, sizeof(test));
    strupr(test);
    if (strstr(test, "T") != NULL || strstr(test, "Y") != NULL)
    {
      filter->cl2_bkwsh = TRUE;
      e = 1;
    }
    else if (strstr(test, "F") != NULL || strstr(test, "N") != NULL)
    {
      filter->cl2_bkwsh = FALSE;
      e = 1;
    }
    else
    {
      e = -1;
    }
    break;
  case 4:
    strncpy(test, buffer, sizeof(test));
    strupr(test);
    if (strchr(test, 'S'))
    {
      filter->media = 'S';
      e = 1;
    }
    else if (strchr(test, 'G'))
    {
      filter->media = 'G';
      e = 1;
    }
    else
    {
      e = -1;
    }
    break;
  case 5:
    e = sscanf(buffer, fmt, &filter->giardia_lr_conv);
    break;
  case 6:
    e = sscanf(buffer, fmt, &filter->virus_lr_conv);
    break;
  case 7:
    e = sscanf(buffer, fmt, &filter->crypto_lr_conv);
    break;
  case 8:
    e = sscanf(buffer, fmt, &filter->giardia_lr_df);
    break;
  case 9:
    e = sscanf(buffer, fmt, &filter->virus_lr_df);
    break;
  case 10:
    e = sscanf(buffer, fmt, &filter->crypto_lr_df);
    break;
  case 11:
    strncpy(test, buffer, sizeof(test));
    strupr(test);
    if (strstr(test, "T") != NULL || strstr(test, "Y") != NULL)
    {
      filter->cfe_turb_flag = TRUE;
      e = 1;
    }
    else if (strstr(test, "F") != NULL || strstr(test, "N") != NULL)
    {
      filter->cfe_turb_flag = FALSE;
      e = 1;
    }
    else
    {
      e = -1;
    }
    break;
  case 12:
    strncpy(test, buffer, sizeof(test));
    strupr(test);
    if (strstr(test, "T") != NULL || strstr(test, "Y") != NULL)
    {
      filter->ife_turb_flag = TRUE;
      e = 1;
    }
    else if (strstr(test, "F") != NULL || strstr(test, "N") != NULL)
    {
      filter->ife_turb_flag = FALSE;
      e = 1;
    }
    else
    {
      e = -1;
    }
    break;
  case 13:
    e = sscanf(buffer, fmt, &filter->crypto_lr_2);
    break;

  default:
    e = -1;
    break;
  }

  //Input Check
  if (filter->volume < 0.0)
    filter->volume = 0.0;
  if (filter->fi_mean < 0.0)
    filter->fi_mean = 0.0;
  if (filter->fi_mean > 1.0)
    filter->fi_mean = 1.0;
  if (filter->fi_t_ten < 0.0)
    filter->fi_t_ten = 0.0;
  if (filter->fi_t_ten > 1.0)
    filter->fi_t_ten = 1.0;
  if (filter->giardia_lr_conv < 0.0)
    filter->giardia_lr_conv = 0.0;
  if (filter->virus_lr_conv < 0.0)
    filter->virus_lr_conv = 0.0;
  if (filter->crypto_lr_conv < 0.0)
    filter->crypto_lr_conv = 0.0;
  if (filter->giardia_lr_df < 0.0)
    filter->giardia_lr_df = 0.0;
  if (filter->virus_lr_df < 0.0)
    filter->virus_lr_df = 0.0;
  if (filter->crypto_lr_df < 0.0)
    filter->crypto_lr_df = 0.0;
  if (filter->crypto_lr_2 < 0.0)
    filter->crypto_lr_2 = 0.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_ssf(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Ssf *ssf;

  ssf = unit->data.ssf;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, ssf->volume);
    break;
  case 1:
    e = sprintf(buffer, fmt, ssf->fi_mean);
    break;
  case 2:
    e = sprintf(buffer, fmt, ssf->fi_t_ten);
    break;
  case 3:
    e = sprintf(buffer, fmt, ssf->toc_rem_lowt);
    break;
  case 4:
    e = sprintf(buffer, fmt, ssf->toc_rem_midt);
    break;
  case 5:
    e = sprintf(buffer, fmt, ssf->toc_rem_hight);
    break;
  case 6:
    e = sprintf(buffer, fmt, ssf->giardia_lr);
    break;
  case 7:
    e = sprintf(buffer, fmt, ssf->virus_lr);
    break;
  case 8:
    e = sprintf(buffer, fmt, ssf->crypto_lr_1);
    break;
  case 9:
    e = sprintf(buffer, fmt, ssf->crypto_lr_2);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_ssf(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Ssf *ssf;

  ssf = unit->data.ssf;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &ssf->volume);
    break;
  case 1:
    e = sscanf(buffer, fmt, &ssf->fi_mean);
    break;
  case 2:
    e = sscanf(buffer, fmt, &ssf->fi_t_ten);
    break;
  case 3:
    e = sscanf(buffer, fmt, &ssf->toc_rem_lowt);
    break;
  case 4:
    e = sscanf(buffer, fmt, &ssf->toc_rem_midt);
    break;
  case 5:
    e = sscanf(buffer, fmt, &ssf->toc_rem_hight);
    break;
  case 6:
    e = sscanf(buffer, fmt, &ssf->giardia_lr);
    break;
  case 7:
    e = sscanf(buffer, fmt, &ssf->virus_lr);
    break;
  case 8:
    e = sscanf(buffer, fmt, &ssf->crypto_lr_1);
    break;
  case 9:
    e = sscanf(buffer, fmt, &ssf->crypto_lr_2);
    break;
  default:
    e = -1;
    break;
  }

  //Input Check
  if (ssf->volume < 0.0)
    ssf->volume = 0.0;
  if (ssf->fi_mean < 0.0)
    ssf->fi_mean = 0.0;
  if (ssf->fi_mean > 1.0)
    ssf->fi_mean = 1.0;
  if (ssf->fi_t_ten < 0.0)
    ssf->fi_t_ten = 0.0;
  if (ssf->fi_t_ten > 1.0)
    ssf->fi_t_ten = 1.0;
  if (ssf->toc_rem_lowt < 0.0)
    ssf->toc_rem_lowt = 0.0;
  if (ssf->toc_rem_lowt > 100.0)
    ssf->toc_rem_lowt = 100.0;
  if (ssf->toc_rem_midt < 0.0)
    ssf->toc_rem_midt = 0.0;
  if (ssf->toc_rem_midt > 100.0)
    ssf->toc_rem_midt = 100.0;
  if (ssf->toc_rem_hight < 0.0)
    ssf->toc_rem_hight = 0.0;
  if (ssf->toc_rem_hight > 100.0)
    ssf->toc_rem_hight = 100.0;
  if (ssf->giardia_lr < 0.0)
    ssf->giardia_lr = 0.0;
  if (ssf->virus_lr < 0.0)
    ssf->virus_lr = 0.0;
  if (ssf->crypto_lr_1 < 0.0)
    ssf->crypto_lr_1 = 0.0;
  if (ssf->crypto_lr_2 < 0.0)
    ssf->crypto_lr_2 = 0.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_def(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Def *def;

  def = unit->data.def;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, def->volume);
    break;
  case 1:
    e = sprintf(buffer, fmt, def->fi_mean);
    break;
  case 2:
    e = sprintf(buffer, fmt, def->fi_t_ten);
    break;
  case 3:
    e = sprintf(buffer, fmt, def->giardia_lr);
    break;
  case 4:
    e = sprintf(buffer, fmt, def->virus_lr);
    break;
  case 5:
    e = sprintf(buffer, fmt, def->crypto_lr_1);
    break;
  case 6:
    e = sprintf(buffer, fmt, def->crypto_lr_2);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_def(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Def *def;

  def = unit->data.def;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &def->volume);
    break;
  case 1:
    e = sscanf(buffer, fmt, &def->fi_mean);
    break;
  case 2:
    e = sscanf(buffer, fmt, &def->fi_t_ten);
    break;
  case 3:
    e = sscanf(buffer, fmt, &def->giardia_lr);
    break;
  case 4:
    e = sscanf(buffer, fmt, &def->virus_lr);
    break;
  case 5:
    e = sscanf(buffer, fmt, &def->crypto_lr_1);
    break;
  case 6:
    e = sscanf(buffer, fmt, &def->crypto_lr_2);
    break;
  default:
    e = -1;
    break;
  }

  //Input Check
  if (def->volume < 0.0)
    def->volume = 0.0;
  if (def->fi_mean < 0.0)
    def->fi_mean = 0.0;
  if (def->fi_mean > 1.0)
    def->fi_mean = 1.0;
  if (def->fi_t_ten < 0.0)
    def->fi_t_ten = 0.0;
  if (def->fi_t_ten > 1.0)
    def->fi_t_ten = 1.0;
  if (def->giardia_lr < 0.0)
    def->giardia_lr = 0.0;
  if (def->virus_lr < 0.0)
    def->virus_lr = 0.0;
  if (def->crypto_lr_1 < 0.0)
    def->crypto_lr_1 = 0.0;
  if (def->crypto_lr_2 < 0.0)
    def->crypto_lr_2 = 0.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_altf(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Altf *altf;

  altf = unit->data.altf;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, altf->giardia_lr);
    break;
  case 1:
    e = sprintf(buffer, fmt, altf->virus_lr);
    break;
  case 2:
    e = sprintf(buffer, fmt, altf->crypto_lr_1);
    break;
  case 3:
    e = sprintf(buffer, fmt, altf->crypto_lr_2);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_altf(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Altf *altf;

  altf = unit->data.altf;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &altf->giardia_lr);
    break;
  case 1:
    e = sscanf(buffer, fmt, &altf->virus_lr);
    break;
  case 2:
    e = sscanf(buffer, fmt, &altf->crypto_lr_1);
    break;
  case 3:
    e = sscanf(buffer, fmt, &altf->crypto_lr_2);
    break;
  default:
    e = -1;
    break;
  }

  //Input Check
  if (altf->giardia_lr < 0.0)
    altf->giardia_lr = 0.0;
  if (altf->virus_lr < 0.0)
    altf->virus_lr = 0.0;
  if (altf->crypto_lr_1 < 0.0)
    altf->crypto_lr_1 = 0.0;
  if (altf->crypto_lr_2 < 0.0)
    altf->crypto_lr_2 = 0.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_bankf(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Bankf *bankf;

  bankf = unit->data.bankf;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    if (bankf->eligible_lt2 == TRUE)
      e = sprintf(buffer, fmt, "TRUE");
    else
      e = sprintf(buffer, fmt, "FALSE");
    break;
  case 1:
    e = sprintf(buffer, fmt, bankf->distance);
    break;
  case 2:
    e = sprintf(buffer, fmt, bankf->crypto_lr_close);
    break;
  case 3:
    e = sprintf(buffer, fmt, bankf->crypto_lr_far);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_bankf(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Bankf *bankf;
  char test[80];

  bankf = unit->data.bankf;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    strncpy(test, buffer, sizeof(test));
    strupr(test);
    if (strstr(test, "T") != NULL || strstr(test, "Y") != NULL)
    {
      bankf->eligible_lt2 = TRUE;
      e = 1;
    }
    else if (strstr(test, "F") != NULL || strstr(test, "N") != NULL)
    {
      bankf->eligible_lt2 = FALSE;
      e = 1;
    }
    else
    {
      e = -1;
    }
    break;
  case 1:
    e = sscanf(buffer, fmt, &bankf->distance);
    break;
  case 2:
    e = sscanf(buffer, fmt, &bankf->crypto_lr_close);
    break;
  case 3:
    e = sscanf(buffer, fmt, &bankf->crypto_lr_far);
    break;
  default:
    e = -1;
    break;
  }

  //Input Check
  if (bankf->distance < 0.0)
    bankf->distance = 0.0;
  if (bankf->crypto_lr_close < 0.0)
    bankf->crypto_lr_close = 0.0;
  if (bankf->crypto_lr_far < 0.0)
    bankf->crypto_lr_far = 0.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_presed(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Presed *presed;

  presed = unit->data.presed;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, presed->volume);
    break;
  case 1:
    e = sprintf(buffer, fmt, presed->sb_mean);
    break;
  case 2:
    e = sprintf(buffer, fmt, presed->sb_t_ten);
    break;
  case 3:
    if (presed->eligible_lt2 == TRUE)
      e = sprintf(buffer, fmt, "TRUE");
    else
      e = sprintf(buffer, fmt, "FALSE");
    break;
  case 4:
    e = sprintf(buffer, fmt, presed->crypto_lr);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_presed(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Presed *presed;
  char test[80];

  presed = unit->data.presed;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &presed->volume);
    break;
  case 1:
    e = sscanf(buffer, fmt, &presed->sb_mean);
    break;
  case 2:
    e = sscanf(buffer, fmt, &presed->sb_t_ten);
    break;
  case 3:
    strncpy(test, buffer, sizeof(test));
    strupr(test);
    if (strstr(test, "T") != NULL || strstr(test, "Y") != NULL)
    {
      presed->eligible_lt2 = TRUE;
      e = 1;
    }
    else if (strstr(test, "F") != NULL || strstr(test, "N") != NULL)
    {
      presed->eligible_lt2 = FALSE;
      e = 1;
    }
    else
    {
      e = -1;
    }
    break;
  case 4:
    e = sscanf(buffer, fmt, &presed->crypto_lr);
    break;
  default:
    e = -1;
    break;
  }

  //Input Check
  if (presed->volume < 0.0)
    presed->volume = 0.0;
  if (presed->sb_mean < 0.0)
    presed->sb_mean = 0.0;
  if (presed->sb_mean > 1.0)
    presed->sb_mean = 1.0;
  if (presed->sb_t_ten < 0.0)
    presed->sb_t_ten = 0.0;
  if (presed->sb_t_ten > 1.0)
    presed->sb_t_ten = 1.0;
  if (presed->crypto_lr < 0.0)
    presed->crypto_lr = 0.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_uvdis(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Uvdis *uvdis;

  uvdis = unit->data.uvdis;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, uvdis->giardia_li);
    break;
  case 1:
    e = sprintf(buffer, fmt, uvdis->virus_li);
    break;
  case 2:
    e = sprintf(buffer, fmt, uvdis->crypto_li);
    break;

  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_uvdis(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Uvdis *uvdis;

  uvdis = unit->data.uvdis;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &uvdis->giardia_li);
    break;
  case 1:
    e = sscanf(buffer, fmt, &uvdis->virus_li);
    break;
  case 2:
    e = sscanf(buffer, fmt, &uvdis->crypto_li);
    break;

  default:
    e = -1;
    break;
  }

  //Input Check
  if (uvdis->giardia_li < 0.0)
    uvdis->giardia_li = 0.0;
  if (uvdis->virus_li < 0.0)
    uvdis->virus_li = 0.0;
  if (uvdis->crypto_li < 0.0)
    uvdis->crypto_li = 0.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_basin(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Basin *basin;

  basin = unit->data.basin;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, basin->volume);
    break;
  case 1:
    e = sprintf(buffer, fmt, basin->sb_mean);
    break;
  case 2:
    e = sprintf(buffer, fmt, basin->sb_t_ten);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_basin(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Basin *basin;

  basin = unit->data.basin;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &basin->volume);
    break;
  case 1:
    e = sscanf(buffer, fmt, &basin->sb_mean);
    break;
  case 2:
    e = sscanf(buffer, fmt, &basin->sb_t_ten);
    break;
  default:
    e = -1;
    break;
  }

  //Input Check
  if (basin->volume < 0.0)
    basin->volume = 0.0;
  if (basin->sb_mean < 0.0)
    basin->sb_mean = 0.0;
  if (basin->sb_mean > 1.0)
    basin->sb_mean = 1.0;
  if (basin->sb_t_ten < 0.0)
    basin->sb_t_ten = 0.0;
  if (basin->sb_t_ten > 1.0)
    basin->sb_t_ten = 1.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_mfuf(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Mfuf *mfuf;

  mfuf = unit->data.mfuf;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, mfuf->recover);
    break;
  case 1:
    e = sprintf(buffer, fmt, mfuf->giardia_lr);
    break;
  case 2:
    e = sprintf(buffer, fmt, mfuf->virus_lr);
    break;
  case 3:
    e = sprintf(buffer, fmt, mfuf->crypto_lr_1);
    break;
  case 4:
    e = sprintf(buffer, fmt, mfuf->crypto_lr_2);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_mfuf(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Mfuf *mfuf;

  mfuf = unit->data.mfuf;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &mfuf->recover);
    break;
  case 1:
    e = sscanf(buffer, fmt, &mfuf->giardia_lr);
    break;
  case 2:
    e = sscanf(buffer, fmt, &mfuf->virus_lr);
    break;
  case 3:
    e = sscanf(buffer, fmt, &mfuf->crypto_lr_1);
    break;
  case 4:
    e = sscanf(buffer, fmt, &mfuf->crypto_lr_2);
    break;
  default:
    e = -1;
    break;
  }

  //Input Check
  if (mfuf->recover < 0.0)
    mfuf->recover = 0.0;
  if (mfuf->recover > 100.0)
    mfuf->recover = 100.0;
  if (mfuf->giardia_lr < 0.0)
    mfuf->giardia_lr = 0.0;
  if (mfuf->virus_lr < 0.0)
    mfuf->virus_lr = 0.0;
  if (mfuf->crypto_lr_1 < 0.0)
    mfuf->crypto_lr_1 = 0.0;
  if (mfuf->crypto_lr_2 < 0.0)
    mfuf->crypto_lr_2 = 0.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_nf(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Nf *nf;

  nf = unit->data.nf;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, nf->mwc);
    break;
  case 1:
    e = sprintf(buffer, fmt, nf->recover);
    break;
  case 2:
    e = sprintf(buffer, fmt, nf->giardia_lr);
    break;
  case 3:
    e = sprintf(buffer, fmt, nf->virus_lr);
    break;
  case 4:
    e = sprintf(buffer, fmt, nf->crypto_lr);
    break;
  case 5:
    e = sprintf(buffer, fmt, nf->toc_rem);
    break;
  case 6:
    e = sprintf(buffer, fmt, nf->uva_rem);
    break;
  case 7:
    e = sprintf(buffer, fmt, nf->br_rem);
    break;
  case 8:
    e = sprintf(buffer, fmt, nf->treat_fraction);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_nf(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Nf *nf;

  nf = unit->data.nf;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &nf->mwc);
    break;
  case 1:
    e = sscanf(buffer, fmt, &nf->recover);
    break;
  case 2:
    e = sscanf(buffer, fmt, &nf->giardia_lr);
    break;
  case 3:
    e = sscanf(buffer, fmt, &nf->virus_lr);
    break;
  case 4:
    e = sscanf(buffer, fmt, &nf->crypto_lr);
    break;
  case 5:
    e = sscanf(buffer, fmt, &nf->toc_rem);
    break;
  case 6:
    e = sscanf(buffer, fmt, &nf->uva_rem);
    break;
  case 7:
    e = sscanf(buffer, fmt, &nf->br_rem);
    break;
  case 8:
    e = sscanf(buffer, fmt, &nf->treat_fraction);
    break;
  default:
    e = -1;
    break;
  }

  //Input checking
  if (nf->mwc < 0.0)
    nf->mwc = 0.0;
  if (nf->recover < 0.0)
    nf->recover = 0.0;
  if (nf->recover > 100.0)
    nf->recover = 100.0;
  if (nf->crypto_lr < 0.0)
    nf->crypto_lr = 0.0;
  if (nf->giardia_lr < 0.0)
    nf->giardia_lr = 0.0;
  if (nf->virus_lr < 0.0)
    nf->virus_lr = 0.0;
  if (nf->toc_rem < 0.0)
    nf->toc_rem = 0.0;
  if (nf->toc_rem > 100.0)
    nf->toc_rem = 100.0;
  if (nf->uva_rem < 0.0)
    nf->uva_rem = 0.0;
  if (nf->uva_rem > nf->toc_rem)
    nf->uva_rem = nf->toc_rem;
  if (nf->br_rem < 0.0)
    nf->br_rem = 0.0;
  if (nf->br_rem > 100.0)
    nf->br_rem = 100.0;
  if (nf->treat_fraction < 0.0)
    nf->treat_fraction = 0.0;
  if (nf->treat_fraction > 1.0)
    nf->treat_fraction = 1.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_iron(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Iron *iron;

  iron = unit->data.iron;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, iron->dose);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_iron(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Iron *iron;

  iron = unit->data.iron;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &iron->dose);
    break;
  default:
    e = -1;
    break;
  }
  //Input checking
  if (iron->dose < 0.0)
    iron->dose = 0.0;

  if (e == 0)
    e = -1;
  return (e);
}

/*
int write_effluent(register char *buffer, short i, struct UnitProcess *unit)
{
  register int                  e; 
  register char                *fmt;
  register struct WTP_effluent *eff;  

  eff = unit->data.wtp_effluent;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch( i ) {
    case  0:
     e = 1 ;

      switch( eff->haa_flag ) {
        case 'H': e = sprintf(buffer,"Haas" ); break;
        case 'T': e = sprintf(buffer,"TAW"  ); break;
        default : e = sprintf(buffer,"Error"); break;
        } 
      break;
    default: e = -1;break;
    }
  if( e<=0 ) e = -1;
  return(e);
}  */
/*
int read_effluent(register char *buffer,short i,struct UnitProcess *unit)
{
  register int                  e;
  register char                *fmt;
  register struct WTP_effluent *eff; 
  char test[80];

  eff = unit->data.wtp_effluent; 
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch( i ) {
    case 0:
      strncpy(test,buffer,sizeof(test));
      strupr(test);

      if(      strchr(test,'H') ) {eff->haa_flag = 'H'; e=1;}
      else if( strchr(test,'T') ) {eff->haa_flag = 'T'; e=1;}  
      else {e=-1;}
      break;
    default:
      e = -1;
      break;
    }

  return(e);
}  */

int write_avg_tap(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Avg_tap *tap;

  tap = unit->data.avg_tap;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:

    e = sprintf(buffer, fmt, tap->days);
    break;

  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_avg_tap(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct Avg_tap *tap;

  tap = unit->data.avg_tap;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &tap->days);

    break;
  default:
    e = -1;
    break;
  }

  //Input checking
  if (tap->days < 0.0)
    tap->days = 0.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_end_sys(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct End_of_system *end;

  end = unit->data.end_of_system;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, end->days);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_end_sys(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct End_of_system *end;

  end = unit->data.end_of_system;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &end->days);
    break;
  default:
    e = -1;
    break;
  }

  //Input checking
  if (end->days < 0.0)
    end->days = 0.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_chlorine(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, chem->chlor);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_chlorine(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &chem->chlor);
    break;
  default:
    e = -1;
    break;
  }

  //Input checking
  if (chem->chlor < 0.0)
    chem->chlor = 0.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_h2so4(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, chem->h2so4);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_h2so4(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &chem->h2so4);
    break;
  default:
    e = -1;
    break;
  }

  //Input checking
  if (chem->h2so4 < 0.0)
    chem->h2so4 = 0.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_naoh(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, chem->naoh);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_naoh(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &chem->naoh);
    break;
  default:
    e = -1;
    break;
  }

  //Input checking
  if (chem->naoh < 0.0)
    chem->naoh = 0.0;
  if (e == 0)
    e = -1;
  return (e);
}

int write_caoh2(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct lime *lime;

  lime = unit->data.lime;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, lime->dose);
    break;
  case 1:
    switch (lime->purpose)
    {
    case 'P':
      e = sprintf(buffer, "PH_ADJ.");
      break;
    case 'S':
      e = sprintf(buffer, "SOFTEN");
      break;
    default:
      e = sprintf(buffer, "Error");
      break;
    }
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_caoh2(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct lime *lime;
  char test[80];

  lime = unit->data.lime;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &lime->dose);
    break;
  case 1:
    strncpy(test, buffer, sizeof(test));
    strupr(test);
    if (strchr(test, 'P'))
    {
      lime->purpose = 'P';
      e = 1;
    }
    else if (strchr(test, 'S'))
    {
      lime->purpose = 'S';
      e = 1;
    }
    else
    {
      e = -1;
    }
    break;
  default:
    e = -1;
    break;
  }

  //Input checking
  if (lime->dose < 0.0)
    lime->dose = 0.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_na2co3(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, chem->soda);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_na2co3(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &chem->soda);
    break;
  default:
    e = -1;
    break;
  }

  //Input checking
  if (chem->soda < 0.0)
    chem->soda = 0.0;
  if (e == 0)
    e = -1;
  return (e);
}

int write_nh3(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    /*	if (unit->type == AMMONIA)*/ e = sprintf(buffer, fmt, chem->nh3);
    /*	else			   e = sprintf( buffer, fmt, chem->nh42so4 ); */
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_nh3(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &chem->nh3);
    break;
  default:
    e = -1;
    break;
  }

  //Input checking
  if (chem->nh3 < 0.0)
    chem->nh3 = 0.0;
  if (e == 0)
    e = -1;
  return (e);
}

int write_kmno4(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, chem->kmno4);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_kmno4(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &chem->kmno4);
    break;
  default:
    e = -1;
    break;
  }

  //Input checking
  if (chem->kmno4 < 0.0)
    chem->kmno4 = 0.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_co2(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, chem->co2);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_co2(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &chem->co2);
    break;
  default:
    e = -1;
    break;
  }

  //Input checking
  if (chem->co2 < 0.0)
    chem->co2 = 0.0;
  if (e == 0)
    e = -1;
  return (e);
}

int write_o3(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, chem->o3);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_o3(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &chem->o3);
    break;
  default:
    e = -1;
    break;
  }

  //Input checking
  if (chem->o3 < 0.0)
    chem->o3 = 0.0;
  if (e == 0)
    e = -1;
  return (e);
}

int write_clo2(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct clo2 *clo2;

  clo2 = unit->data.clo2;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, clo2->dose);
    break;
  case 1:
    e = sprintf(buffer, fmt, clo2->conversion);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_clo2(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct clo2 *clo2;

  clo2 = unit->data.clo2;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &clo2->dose);
    break;
  case 1:
    e = sscanf(buffer, fmt, &clo2->conversion);
    break;
  default:
    e = -1;
    break;
  }

  //Input checking
  if (clo2->dose < 0.0)
    clo2->dose = 0.0;
  if (clo2->conversion < 0.0)
    clo2->conversion = 0.0;
  if (clo2->conversion > 100.0)
    clo2->conversion = 100.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_naocl(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, chem->naocl);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_naocl(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &chem->naocl);
    break;
  default:
    e = -1;
    break;
  }

  //Input Check
  if (chem->naocl < 0.0)
    chem->naocl = 0.0;

  if (e == 0)
    e = -1;
  return (e);
}

int write_so2(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtout;
  switch (i)
  {
  case 0:
    e = sprintf(buffer, fmt, chem->so2);
    break;
  default:
    e = -1;
    break;
  }
  if (e <= 0)
    e = -1;
  return (e);
}

int read_so2(register char *buffer, short i, struct UnitProcess *unit)
{
  register int e;
  register const char *fmt;
  register struct chemical *chem;

  chem = unit->data.chemical;
  fmt = UnitProcessTable[unit->type].data[i].fmtin;
  switch (i)
  {
  case 0:
    e = sscanf(buffer, fmt, &chem->so2);
    break;
  default:
    e = -1;
    break;
  }

  //Input Check
  if (chem->so2 < 0.0)
    chem->so2 = 0.0;

  if (e == 0)
    e = -1;
  return (e);
}
