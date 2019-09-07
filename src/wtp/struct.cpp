/* Struct.c  -- September 16, 1993
*
*  Functions in this file support management of a linked list of
*  unit processes.
*
*  Memory management functions for struct ProcessTrain:
*    AllocProcessTrain()
*    FreeProcessTrain()
*
*  Memory management functions for struct UnitProcess:
*    AllocUnitProcess()
*    FreeUnitProcess()
*
*  Functions which support building a process train:
*    InitProcessTrain()
*    AddUnitProcess()
*    RemoveUnitProcess()
*    MoveUnitProcess()
*
*  Functions which support indexing through a process train:
*    FirstUnitProcess()
*    LastUnitProcess()
*    NextUnitProcess()
*    PrevUnitProcess()
*    GetUnitProcess()
*    GetUnitIndex()
*  
*  Michael D. Cummins
*    December 1992
*/
#include "wtp.h"

/******************   AllocProcessTrain  ********************************/
struct ProcessTrain *AllocProcessTrain(void)
/*
*  Purpose:  Allocate memory and return an initialized process train.
*/
{
  register struct ProcessTrain *train;

  if ((train = (struct ProcessTrain *)calloc(1, sizeof(struct ProcessTrain))) != NULL)
  {
    /* Initialize list as empty */
    train->head = (struct UnitProcess *)&train->null;
    train->null = NULL;
    train->tail = (struct UnitProcess *)&train->head;

    /*clear memory*/
    memset(train->file_name, 0, sizeof(train->file_name));

    /* Let InitProcessTrain() do its thing. */
    InitProcessTrain(train);
  }

  return (train);
}

/******************   FreeProcessTrain  ***********************************/
struct ProcessTrain *FreeProcessTrain(register struct ProcessTrain *train)
/*
*  Purpose:Deallocate all memory associated with a process train.
*/
{
  register struct UnitProcess *unit;

  if (train)
  {
    /* Deallocate all unit processes */
    while ((unit = FirstUnitProcess(train)) != NULL)
    {
      MoveUnitProcess(NULL, unit);
      FreeUnitProcess(unit);
    }

    free(train);
  }

  return (NULL);
}

/******************   InitProcessTrain  ************************************/
struct ProcessTrain *InitProcessTrain(register struct ProcessTrain *train)
/*
*  Purpose:
*    1. Remove and deallocate all unit processes from process train.
*    2. Add 'influent' unit process.
*    3. train->file_name is not changes.
*/
{
  register struct UnitProcess *unit;

  if (train)
  {

    /* Deallocate all unit processes */
    while ((unit = FirstUnitProcess(train)) != NULL)
    {
      MoveUnitProcess(NULL, unit); /* Delink unit */
      FreeUnitProcess(unit);
    }

    /* First unit process must be raw water */
    AddUnitProcess(train, INFLUENT);
  }

  return (train);
}

/*****************   AllocUnitProcess  ****************************************/
struct UnitProcess *AllocUnitProcess(register short type)
/*
*  Purpose: Allocate and initialize a UnitProcess and UnitProcess.data
*           data structure.  The unit process data is initialized to
*           default values.
*
*  Input:
*    type = defined UnitProcess type
*
*  Return:
*    Pointer to UnitProcess data structure or NULL.
*/
{
  register struct UnitProcess *unit;

  if ((unit = (struct UnitProcess *)calloc(1, sizeof(struct UnitProcess))) != NULL)
  {
    unit->next = NULL; /* Flag as not linked into process train */
    unit->prev = NULL;
    unit->type = type;
    unit->pad = 0;
    switch (type)
    {

    case INFLUENT:
      if ((unit->data.influent = (struct Influent *)
               calloc(1, sizeof(struct Influent))) != NULL)
      {
        *(unit->data.influent) = default_influent;
      }
      break;

    case ALUM:
      if ((unit->data.alum = (struct Alum *)
               calloc(1, sizeof(struct Alum))) != NULL)
      {
        *(unit->data.alum) = default_alum;
      }
      break;

    case GAC:
      if ((unit->data.gac = (struct Gac *)
               calloc(1, sizeof(struct Gac))) != NULL)
      {
        *(unit->data.gac) = default_gac;
      }
      break;

    case FILTER:
      if ((unit->data.filter = (struct Filter *)
               calloc(1, sizeof(struct Filter))) != NULL)
      {
        *(unit->data.filter) = default_filter;
      }
      break;

    case BASIN:
      if ((unit->data.basin = (struct Basin *)
               calloc(1, sizeof(struct Basin))) != NULL)
      {
        *(unit->data.basin) = default_basin;
      }
      break;

    case O3_CONTACTOR:
      if ((unit->data.basin = (struct Basin *)
               calloc(1, sizeof(struct Basin))) != NULL)
      {
        *(unit->data.basin) = default_o3_contactor;
      }
      break;

    case MFUF_UP:
      if ((unit->data.mfuf = (struct Mfuf *)
               calloc(1, sizeof(struct Mfuf))) != NULL)
      {
        *(unit->data.mfuf) = default_mfuf;
      }
      break;

    case NF_UP:
      if ((unit->data.nf = (struct Nf *)
               calloc(1, sizeof(struct Nf))) != NULL)
      {
        *(unit->data.nf) = default_nf;
      }
      break;

    case BANK_FILTER:
      if ((unit->data.bankf = (struct Bankf *)
               calloc(1, sizeof(struct Bankf))) != NULL)
      {
        *(unit->data.bankf) = default_bankf;
      }
      break;

    case PRESED_BASIN:
      if ((unit->data.presed = (struct Presed *)
               calloc(1, sizeof(struct Presed))) != NULL)
      {
        *(unit->data.presed) = default_presed;
      }
      break;

    case UV_DIS:
      if ((unit->data.uvdis = (struct Uvdis *)
               calloc(1, sizeof(struct Uvdis))) != NULL)
      {
        *(unit->data.uvdis) = default_uvdis;
      }
      break;

    case SLOW_FILTER:
      if ((unit->data.ssf = (struct Ssf *)
               calloc(1, sizeof(struct Ssf))) != NULL)
      {
        *(unit->data.ssf) = default_ssf;
      }
      break;

    case DE_FILTER:
      if ((unit->data.def = (struct Def *)
               calloc(1, sizeof(struct Def))) != NULL)
      {
        *(unit->data.def) = default_def;
      }
      break;

    case BAG_FILTER:
      if ((unit->data.altf = (struct Altf *)
               calloc(1, sizeof(struct Altf))) != NULL)
      {
        *(unit->data.altf) = default_bagf;
      }
      break;

    case CART_FILTER:
      if ((unit->data.altf = (struct Altf *)
               calloc(1, sizeof(struct Altf))) != NULL)
      {
        *(unit->data.altf) = default_cartf;
      }
      break;

    case IRON:
      if ((unit->data.iron = (struct Iron *)
               calloc(1, sizeof(struct Iron))) != NULL)
      {
        *(unit->data.iron) = default_iron;
      }
      break;

    case CHLORINE_DIOXIDE:
      if ((unit->data.clo2 = (struct clo2 *)
               calloc(1, sizeof(struct clo2))) != NULL)
      {
        *(unit->data.clo2) = default_clo2;
      }
      break;

    case LIME:
      if ((unit->data.lime = (struct lime *)
               calloc(1, sizeof(struct lime))) != NULL)
      {
        *(unit->data.lime) = default_lime;
      }
      break;

    /* Chemical Additions */
    case CHLORINE:
    case SULFURIC_ACID:
      /*  case LIME:  */
    case SODA_ASH:
    case AMMONIA:
    case AMMONIUM_SULFATE:
    case PERMANGANATE:
    case CARBON_DIOXIDE:
    case OZONE:
    case SODIUM_HYDROXIDE:
    case HYPOCHLORITE:
    case SULFUR_DIOXIDE:
      if ((unit->data.chemical = (struct chemical *)
               calloc(1, sizeof(struct chemical))) != NULL)
      {
        *(unit->data.chemical) = default_chemical;
      }
      break;

    case WTP_EFFLUENT:
      if ((unit->data.avg_tap = (struct Avg_tap *)
               calloc(1, sizeof(struct Avg_tap))) != NULL)
      {
        *(unit->data.avg_tap) = default_avg_tap;
      }
      /*  if( (unit->data.wtp_effluent = (struct WTP_effluent *)
                            calloc(1,sizeof(struct WTP_effluent))) != NULL )
            {
	      *(unit->data.wtp_effluent) = default_wtp_effluent;
              DO NOTHING - WJS, 11/98 
            }      */
      break;

    case AVG_TAP:
      if ((unit->data.avg_tap = (struct Avg_tap *)
               calloc(1, sizeof(struct Avg_tap))) != NULL)
      {
        *(unit->data.avg_tap) = default_avg_tap;
      }
      break;

    case END_OF_SYSTEM:
      if ((unit->data.end_of_system = (struct End_of_system *)
               calloc(1, sizeof(struct End_of_system))) != NULL)
      {
        *(unit->data.end_of_system) = default_end_of_system;
      }
      break;

    case RAPID_MIX:
      if ((unit->data.basin = (struct Basin *)
               calloc(1, sizeof(struct Basin))) != NULL)
      {
        *(unit->data.basin) = default_rapid_mix;
      }
      break;

    case SLOW_MIX:
      if ((unit->data.basin = (struct Basin *)
               calloc(1, sizeof(struct Basin))) != NULL)
      {
        *(unit->data.basin) = default_slow_mix;
      }
      break;

    case SETTLING_BASIN:
      if ((unit->data.basin = (struct Basin *)
               calloc(1, sizeof(struct Basin))) != NULL)
      {
        *(unit->data.basin) = default_settling_basin;
      }
      break;

    case CONTACT_TANK:
      if ((unit->data.basin = (struct Basin *)
               calloc(1, sizeof(struct Basin))) != NULL)
      {
        *(unit->data.basin) = default_contact_tank;
      }
      break;

    case CLEARWELL:
      if ((unit->data.basin = (struct Basin *)
               calloc(1, sizeof(struct Basin))) != NULL)
      {
        *(unit->data.basin) = default_clearwell;
      }
      break;

    case LOCATION_1:
      if ((unit->data.avg_tap = (struct Avg_tap *)
               calloc(1, sizeof(struct Avg_tap))) != NULL)
      {
        *(unit->data.avg_tap) = default_avg_tap;
      }
      break;

      /*	case LOCATION_2:
	  if( (unit->data.avg_tap = (struct Avg_tap *)
                                 calloc(1,sizeof(struct Avg_tap))) != NULL)
            {
              *(unit->data.avg_tap) = default_avg_tap;
            }
	  break;

	case LOCATION_3:
	  if( (unit->data.avg_tap = (struct Avg_tap *)
                                 calloc(1,sizeof(struct Avg_tap))) != NULL)
            {
              *(unit->data.avg_tap) = default_avg_tap;
            }
	  break;

	case LOCATION_4:
	  if( (unit->data.avg_tap = (struct Avg_tap *)
                                 calloc(1,sizeof(struct Avg_tap))) != NULL)
            {
              *(unit->data.avg_tap) = default_avg_tap;
            }
	  break;

	case LOCATION_5:
	  if( (unit->data.avg_tap = (struct Avg_tap *)
                                 calloc(1,sizeof(struct Avg_tap))) != NULL)
            {
              *(unit->data.avg_tap) = default_avg_tap;
            }
	  break; */

    default:
      unit->data.ptr = NULL;
      printf("Unknown unit process type in AllocUnitProcess:%d\n", type);
      break;
    }

    if (unit->data.ptr == NULL)
    {
      free(unit);
      unit = NULL;
    }
  }

  return (unit);
}

/*****************   FreeUnitProcess  ********************************/
struct UnitProcess *FreeUnitProcess(register struct UnitProcess *unit)
/*
*  Purpose: Deallocate memory allocated by AllocUnitProcess.
*           NULL is always returned.
*/
{
  if (unit != NULL)
  {
    if (unit->next && unit->prev)
      MoveUnitProcess(NULL, unit);
    if (unit->data.ptr != NULL)
      free(unit->data.ptr);
    free(unit);
  }

  return (NULL);
}

/*****************   AddUnitProcess  ***********************************/
struct UnitProcess *AddUnitProcess(
    register struct ProcessTrain *train,
    register short type)
/*
*  Purpose: Allocate, initialize, and add a unit process to the end of
*           the process train.  A unit process design and operating
*           parameter data packet is also allocated and initialized with
*           default values.
*
*  Notes:
*   1.  AllocUnitProcess() allocates the memory and initializes the data
*       fields.  AddUnitProcess() links the unit process to the tail 
*       of the process train.
*/
{
  register struct UnitProcess *unit;

  if ((unit = AllocUnitProcess(type)) != NULL)
  {
    /* Link 'unit' to tail of process train */
    MoveUnitProcess(train->tail, unit);
  }

  return (unit);
}

/*****************   RemoveUnitProcess  ***********************************/
struct UnitProcess *RemoveUnitProcess(register struct UnitProcess *unit)
/*
*  Purpose: Delink and dealloc 'unit' from the process train.  The 'unit'
*           does not need to be at the end of the process train.
*
*  Return:
*    NULL is always returned.
*/
{
  if (unit != NULL)
  {
    MoveUnitProcess(NULL, unit); /* Delink unit */
    FreeUnitProcess(unit);
  }

  return (NULL);
}

/*****************   MoveUnitProcess  *************************************/
struct UnitProcess *MoveUnitProcess(
    register struct UnitProcess *prev, /* May be NULL, train->head or tail */
    register struct UnitProcess *unit  /* May be linked or un-linked       */
)
/*
*  Purpose: Move 'unit' to follow 'prev'.
*
*           MoveUnitProcess() manipulates the double linked list for many
*           structure functions.  It does not allocate or deallocate
*           memory or change values in the unit process data packet.
*
*  Inputs:
*    prev = If 'prev' is a pointer to a linked UnitProcess then 'unit'
*           is moved to follow 'prev'.
*         . If 'prev' is NULL or un-linked then 'unit' is removed from its 
*           list and not re-linked.
*         . If 'prev' is train->head then 'unit' is inserted as the first
*           unit process in the process train.
*         . If 'prev' is train->tail then 'unit' is inserted as the last
*           unit process in the process train.
* 
*    unit = Unit process to be operated on.
*         . If 'unit' is linked into a list then it is first removed from
*           the list.
*
*  Return:
*   unit is always returned.
*
*  Note:
*   1. MoveUnitProcess() is the workhorse for manipulating the double
*      linked list.  Many other functions call on this function.
*/
{

  if (unit != NULL && unit != prev)
  {

    /* De-link: */
    if (unit->prev && unit->next)
    {
      /* 'unit' is in a linked list -- remove unit from list. */
      unit->prev->next = unit->next;
      unit->next->prev = unit->prev;

      unit->next = NULL; /* Flag 'unit' as removed */
      unit->prev = NULL;
    }
    else
    {
      /* 'unit' is not in a linked list -- do not delink 'unit'. */
    }

    /* Re-Link: */
    if (prev && (prev->next || prev->prev))
    {
      /* Note on (prev->next || prev->prev) logic:                  */
      /*   prev->next!=NULL && prev->prev!=NULL; prev is in a list  */
      /*   prev->next==NULL && prev->prev!=NULL; prev is train.head (this seems reversed!) */
      /*   prev->next!=NULL && prev->prev==NULL; prev is train.tail (seems reversed!)*/
      /* In all 3 of the above cases, insert unit after prev        */
      /* The above 3 cases reduce to (prev->next || prev->prev)     */

      /* Insert 'unit' after 'prev' */
      unit->next = prev->next;
      unit->prev = prev;

      prev->next->prev = unit;
      prev->next = unit;
    }
    else
    {
      /* prev->next==NULL && prev->prev==NULL; prev not in a list or */
      /* prev==NULL --- do not link 'unit'. */
    }
  }

  return (unit);
}

/*****************   FirstUnitProcess  **********************************/
struct UnitProcess *FirstUnitProcess(struct ProcessTrain *train)
/*
*  Purpose: Return first unit process in process train or NULL if the
*           process train list is empty.
*/
{
  struct UnitProcess *unit;
  if (train)
  {
    unit = train->head;
    if (unit->next)
      return (unit); /* Valid unit process */
    else
      return (NULL); /* Empty list */
  }
  else
  {
    return (NULL); /* Bad input */
  }
}

/*****************   LastUnitProcess  **********************************/
struct UnitProcess *LastUnitProcess(struct ProcessTrain *train)
/*
*  Purpose: Return last unit process in process train or NULL if the
*           process train list is empty.
*/
{
  struct UnitProcess *unit;
  if (train)
  {
    unit = train->tail;
    if (unit->prev)
      return (unit); /* Valid unit process */
    else
      return (NULL); /* Empty list */
  }
  else
  {
    return (NULL); /* Bad input */
  }
}

/*****************   NextUnitProcess  ************************************/
struct UnitProcess *NextUnitProcess(register struct UnitProcess *unit)
/*
*  Purpose: Return next unit process in process train or NULL if 'unit'
*           is the last unit process.
*/
{
  if (unit->next->next)
    return (unit->next);
  else
    return (NULL);
}

/*****************   PrevUnitProcess  ***********************************/
struct UnitProcess *PrevUnitProcess(register struct UnitProcess *unit)
/*
*  Purpose: Return previous unit process in process train of NULL if 'unit'
*           is the first unit process.
*/
{
  if (unit->prev->prev)
    return (unit->prev);
  else
    return (NULL);
}

/**************     GetUnitProcess   ************************************/
struct UnitProcess *GetUnitProcess(struct ProcessTrain *train, int n)
/*
*  Purpose: Return the 'n'th unit process in the process train.
*
*  Notes:
*    1. n=0 will return a pointer to Influent.
*    2. If n is greater than the number of unit processes in the process
*       train then NULL is returned.
*/
{
  int i;
  struct UnitProcess *unit = NULL;

  if (n >= 0)
  {
    for (unit = FirstUnitProcess(train), i = 0;
         unit != NULL && i < n;
         unit = NextUnitProcess(unit), i++)
      ;
  }

  return (unit);
}

/** GetUnitIndex  ******************************************************/
int GetUnitIndex(struct ProcessTrain *train, struct UnitProcess *unit)
/*
*  Purpoes: Return the location (index) of 'unit' in the process train.
*
*  Inputs:
*   train = Process train control structure.
*   unit  = Unit process in 'train'
*
*  Return: The index of 'unit' in the process train or -1 if 'unit' is
*          not in 'train'
*   
*  Notes:
*   1) GetUnitIndex() is the opposit of GetUnitProcess().
*/
{
  int index = -1;
  struct UnitProcess *test = NULL;

  if (unit != NULL)
  {
    for (test = FirstUnitProcess(train), index = 0;
         test != NULL && test != unit;
         test = NextUnitProcess(test), index++)
      ;
    if (test != unit)
      index = -1;
  }

  return (index);
}
