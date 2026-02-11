// Copyright (C) Jonathan Lawder 2001-2011

#ifdef __MSDOS__
	#include "db.h"
	#include "..\hilbert\hilbert.h"
	#include "..\utils\utils.h"
#else
	#include "db.h"
	#include "../hilbert/hilbert.h"
	#include "../utils/utils.h"
#endif

#define		ACTIVE			1
#define		PARTIAL_MATCH		2
#define		RANGE_QUERY		4

using namespace std;

/*============================================================================*/
/*                            RET_SET::RET_SET                                */
/*============================================================================*/
RET_SET::RET_SET( int dims )
{
	flags = 0;
	Qsaf = 0;
	buffslot = numspec = pos = 0;

	LB = new PU_int[dims];
	UB = new PU_int[dims];
	memset( LB, 0, sizeof(PU_int) * dims );
	memset( UB, 0, sizeof(PU_int) * dims );
}

/*============================================================================*/
/*                            RET_SET::~RET_SET                                */
/*============================================================================*/
RET_SET::~RET_SET()
{
	delete [] LB;
	delete [] UB;
}

/*============================================================================*/
/***                   DBASE::db_open_set				    ***/
/*============================================================================*/
//	FOR PARTIAL MATCH QUERIES
bool DBASE::db_open_set( PU_int *point, int *set_id )
{
	int	i;
	int	lpage;
	HU_int	*minmatch = new U_int[dimensions];
	HU_int	*key = new U_int[dimensions];

	// find an inactive set...
	if (FreeRet_setList.size() > 0)
	{
		*set_id = FreeRet_setList.top();
		FreeRet_setList.pop();
	}
	else
	// ...or create a new one
	{
		*set_id = Ret_set.size();
		RET_SET *r = new RET_SET( dimensions );
		Ret_set.push_back( r );
	}


	// initialise the RET_SET
	for (i = 0; i < dimensions; i++)
	{
// add code in here to check values are in range
		if (point[i] != _UNSPECIFIED_)
		{
			Ret_set[*set_id]->Qsaf |= ((U_int)1 << (dimensions-1-i));
			Ret_set[*set_id]->numspec++;
			Ret_set[*set_id]->LB[i] = point[i];
		}
		else
		{
			Ret_set[*set_id]->LB[i] = MINTOKEN;
		}
	}

	if (Ret_set[*set_id]->numspec == 0)
	{
		cout << "ERROR 1 in dbi_open_set(): not programmed to answer "
		     << "fully unspecified queries\n";
		delete [] minmatch;
		delete [] key;
		return false;
	}

       	// find the minimum match to the query
	memset( minmatch, 0, sizeof(PU_int) * dimensions );
	memset( key, 0, sizeof(PU_int) * dimensions );

       	if (false == H_nextmatch_PM( Ret_set[*set_id]->LB, minmatch, key,
       			Ret_set[*set_id]->Qsaf, dimensions ))
       	{
       		errorexit("ERROR 2 in db_open_set(): lowest match not found\n");
       	}
	
	// find the page that may contain the minimum match       		
	lpage = BT.idx_search( minmatch );
	
	Ret_set[*set_id]->flags = ACTIVE | PARTIAL_MATCH;
	// bring in the first page to search
	Ret_set[*set_id]->buffslot = Buffer.b_page_retrieve( lpage );

	// find the position on the page from which to start the search
	if (Ret_set[*set_id]->Qsaf >= ((U_int)1 << (dimensions-1)))
	{
		// binary search possible
		i = Buffer.BSlot[Ret_set[*set_id]->buffslot]->BPage.p_find_pageslot( Ret_set[*set_id]->LB );
		if (i < 0)
			i = -i;
		Ret_set[*set_id]->pos = i;
	}
	else
		Ret_set[*set_id]->pos = 1;

	// no need to tag the buffslot as 'FIXED' - done in b_page_retrieve()
	Buffer.BSlot[Ret_set[*set_id]->buffslot]->query = true; // QUERY;

	delete [] minmatch;
	delete [] key;
	return true;
}

/*============================================================================*/
/***                   DBASE::db_range_open_set				    ***/
/*============================================================================*/
//	FOR RANGE QUERIES
bool DBASE::db_range_open_set( PU_int* LB, PU_int *UB, int *set_id )
{
	int	i;
	int	lpage;
	HU_int	*minmatch = new U_int[dimensions];
	HU_int	*key = new U_int[dimensions];

	// find an inactive set...
	if (FreeRet_setList.size() > 0)
	{
		*set_id = FreeRet_setList.top();
		FreeRet_setList.pop();
	}
	else
	// ...or create a new one
	{
		*set_id = Ret_set.size();
		RET_SET *r = new RET_SET( dimensions );
		Ret_set.push_back( r );
	}
	
	for (i = 0, Ret_set[*set_id]->numspec = dimensions; i < dimensions; i++)
	{
// add code in here to check values are in range
		if (LB[i] != _UNSPECIFIED_)
			Ret_set[*set_id]->LB[i] = LB[i];
		else
			Ret_set[*set_id]->LB[i] = MINTOKEN;

		if (UB[i] != _UNSPECIFIED_)
			Ret_set[*set_id]->UB[i] = UB[i];
		else
			Ret_set[*set_id]->UB[i] = MAXTOKEN;

		if (Ret_set[*set_id]->LB[i] == MINTOKEN &&
			Ret_set[*set_id]->UB[i] == MAXTOKEN)
			Ret_set[*set_id]->numspec--;
	}
	if (Ret_set[*set_id]->numspec == 0)
	{
		cout << "ERROR 1 in dbi_range_open_set(): not programmed to answer "
		     << "unspecified queries\n";
		delete [] minmatch;
		delete [] key;
		return false;
	}

	// find lowest match to query
	memset( minmatch, 0, sizeof(PU_int) * dimensions );
	memset( key, 0, sizeof(PU_int) * dimensions );

     	if (false == H_nextmatch_RQ( Ret_set[*set_id]->LB,
       			Ret_set[*set_id]->UB,
       			minmatch, key, dimensions ))
        {
		errorexit("ERROR 2 in dbi_range_open_set(): lowest match not found\n");
	}
	
	// find the page that may contain the minimum match       		
	lpage = BT.idx_search( minmatch );

	Ret_set[*set_id]->flags = ACTIVE | RANGE_QUERY;
	// bring in the first page to search
	Ret_set[*set_id]->buffslot = Buffer.b_page_retrieve( lpage );

	i = Buffer.BSlot[Ret_set[*set_id]->buffslot]->BPage.p_find_pageslot( Ret_set[*set_id]->LB );
	if (i < 0)
		i = -i;
	Ret_set[*set_id]->pos = i;

	Buffer.BSlot[Ret_set[*set_id]->buffslot]->query = true; // QUERY;

	delete [] minmatch;
	delete [] key;
	return true;
}

/*============================================================================*/
/***                   DBASE::db_close_set				    ***/
/*============================================================================*/
// this function will do for both partial match and range queries
// we NEVER want to remove a RET_SET from the vector as we rely on the set_id
// of a RET_SET remaining constant while it is in use
bool DBASE::db_close_set( int set_id )
{
	int i;

	if (! (Ret_set[set_id]->flags & ACTIVE))
	{
		cerr << "ERROR - Ret_set " << set_id << " is inactive - unexpected\n";
		return false;
	}

	if (Ret_set.size() == FreeRet_setList.size())
		errorexit( "ERROR in db_close_set - no active RET_SETs\n" );

	if (Ret_set.size() - FreeRet_setList.size() == 1)
	{
		// this is the only ACTIVE Ret_set
		i = -1;
	}
	else
	{
	       	// check all of the Ret_sets to see if any is accessing the same
       		// page that Ret_set[set_id] last accessed
	       	// (could use an iterator here but this is simple)
       		for (i = Ret_set.size() - 1; i >= 0; i--)
		{
			if (i == set_id || !(Ret_set[i]->flags & ACTIVE))
				continue;
			if (Ret_set[i]->buffslot == Ret_set[set_id]->buffslot)
			// another RET_SET is accessing the same page as Ret_set[set_id]
				break;
		}
	}

       	if (i < 0)
       	{
       		// no other Ret_set is using the page
       		Buffer.BSlot[Ret_set[set_id]->buffslot]->fix =
       		Buffer.BSlot[Ret_set[set_id]->buffslot]->query = false;
       	}

	// free-up the RET_SET
	Ret_set[set_id]->flags = 0;
	Ret_set[set_id]->Qsaf = 0;
	Ret_set[set_id]->buffslot = Ret_set[set_id]->numspec = Ret_set[set_id]->pos = 0;
	memset( Ret_set[set_id]->LB, 0, sizeof(PU_int) * dimensions );
	memset( Ret_set[set_id]->UB, 0, sizeof(PU_int) * dimensions );
	FreeRet_setList.push( set_id );

	return true;
}

/*============================================================================*/
/***                   DBASE::db_fetch_another				    ***/
/*============================================================================*/
// memory must already have been allocated to retval
// retval holds a match to the query (if this function returns true)
// this function returns as soon as it finds a match to the query
bool DBASE::db_fetch_another( int set_id, PU_int *retval )
{
	int	i, buffslot;
	HU_int	*next_pagekey = new U_int[dimensions], // the key of the next page after the current one
		*next_match = new U_int[dimensions];
	PU_int	*query = Ret_set[set_id]->LB,
		*data;
	int	lpage, pos, end;
	U_int	Qsaf = Ret_set[set_id]->Qsaf;

	// while there are still pages left to search
	for (;;)
	{
		buffslot = Ret_set[set_id]->buffslot;
		pos = Ret_set[set_id]->pos;
		end = Buffer.BSlot[buffslot]->BPage.page_hdr->size;
		data = Buffer.BSlot[buffslot]->BPage.data[pos];

		// while not at end of page (currently being searched)
		// NB incrementing data moves it along one PU_int, not one record
		for ( ;pos <= end; pos++, data+= (dimensions + NO_EXTRA_TOKENS) )
		{
			i = 0;

			if (Qsaf >= ((U_int)1 << (dimensions-1)))
			// 'high' attribute(s) specified
			{
				for (; i < dimensions && Qsaf & ((U_int)1 << (dimensions-1-i)); i++)
					if (query[i] != data[i])
						break;
				if (i < dimensions && Qsaf & ((U_int)1 << (dimensions-1-i)))
					break; /* no more matches on this page */
			}

			for (; i < dimensions; i++)
			{
				if (Qsaf & ((U_int)1 << (dimensions-1-i)))
				{
					if (query[i] != data[i])
					{
						break;
					}
				}
			}
			if (i == dimensions) /* next_match found */
			{
				keycopy( retval, data, dimensions );
				Ret_set[set_id]->pos = pos + 1;
				delete [] next_pagekey;
				delete [] next_match;
				return true;
			}
		}

		if (Qsaf == ((U_int)((1 << dimensions)-1)))
		{
			// query is fully specified - to get this far, either we've already
			// found the (single) next_match or searched the page that might contain it
			delete [] next_pagekey;
			delete [] next_match;
			return false;
		}

		// have we just searched the last page?
		if (Buffer.BSlot[buffslot]->BPage.page_hdr->lpage == LastPage)
		{
			// there can be no higher match
			delete [] next_pagekey;
			delete [] next_match;
			return false;
		}

		// find key of next page - there will be one
		keycopy( next_pagekey,
			 BT.idx_get_next_key(
				Buffer.BSlot[buffslot]->BPage.index,
				Buffer.BSlot[buffslot]->BPage.page_hdr->lpage ),
			 dimensions );

		// find next match above this key
	 	memset( next_match, 0, sizeof(HU_int) * dimensions );

		if (false == H_nextmatch_PM( query, next_match, next_pagekey,
					Qsaf, dimensions ))
		{
			delete [] next_pagekey;
			delete [] next_match;
			return false; // no higher matching hilbert codes
		}

		// clear the current page's flags if not in use by another Ret_set
        	if (Ret_set.size() - FreeRet_setList.size() == 1)
        	{
        		// this is the only ACTIVE Ret_set
        		i = -1;
        	}
        	else
        	{
        	       	// check all of the Ret_sets to see if any is accessing the same
               		// page that Ret_set[set_id] last accessed
        	       	// (could use an iterator here but this is simple)
               		for (i = Ret_set.size() - 1; i >= 0; i--)
        		{
        			if (i == set_id || !(Ret_set[i]->flags & ACTIVE))
        				continue;
        			if (Ret_set[i]->buffslot == Ret_set[set_id]->buffslot)
        			// another RET_SET is accessing the same page as Ret_set[set_id]
        				break;
        		}
        	}
            	if (i < 0)
          	{
          		// no other Ret_set is using the page
          		Buffer.BSlot[buffslot]->fix =
          		Buffer.BSlot[buffslot]->query = false;
          	}

		// find the page that may contain the match
		lpage = BT.idx_search( next_match );

		Ret_set[set_id]->buffslot = Buffer.b_page_retrieve( lpage );
		Buffer.BSlot[Ret_set[set_id]->buffslot]->query = true; // QUERY;

		if (Qsaf >= ((U_int)1 << (dimensions-1)))
		{
			// binary search possible
			i = Buffer.BSlot[Ret_set[set_id]->buffslot]->BPage.p_find_pageslot( query );
			if (i < 0)
				i = -i;
			Ret_set[set_id]->pos = i;
		}
		else
			Ret_set[set_id]->pos = 1;
	}
}

/*============================================================================*/
/***                   DBASE::db_range_fetch_another			    ***/
/*============================================================================*/
// memory must already have been allocated to retval
bool DBASE::db_range_fetch_another( int set_id, PU_int *retval )
{
	int	i, buffslot;
	bool	maybe, no_more;
	HU_int	*next_pagekey = new U_int[dimensions], // the key of the next page after the current one
		*next_match = new U_int[dimensions];
	PU_int	*Lobound = Ret_set[set_id]->LB,
		*Upbound = Ret_set[set_id]->UB,
		*data;
	int	lpage, pos, end;

	for (;;)
	{
		buffslot = Ret_set[set_id]->buffslot;
		pos = Ret_set[set_id]->pos;
		end = Buffer.BSlot[buffslot]->BPage.page_hdr->size;
		data = Buffer.BSlot[buffslot]->BPage.data[pos];

		// (see note in db_fetch_another())
		for ( ;pos <= end; pos++, data+= (dimensions + NO_EXTRA_TOKENS) )
		{
			for (maybe = no_more = false, i = 0; i < dimensions; i++)
			{
				if (data[i] > Upbound[i])
				{
					if (!maybe)
						no_more = true;
					break;
				}
				if (data[i] < Lobound[i])
					break;
				maybe = true;
			}

			if (i == dimensions) // next_match found
			{
				keycopy( retval, data, dimensions);
				Ret_set[set_id]->pos = pos + 1;
				delete [] next_pagekey;
				delete [] next_match;
				return true;
			}
			if (no_more)
				break;
		}

		// have we just searched the last page?
		if (Buffer.BSlot[buffslot]->BPage.page_hdr->lpage == LastPage)
		{
			// there can be no higher match
			delete [] next_pagekey;
			delete [] next_match;
			return false;
		}

		// find key of next page - there will be one
		keycopy( next_pagekey,
			 BT.idx_get_next_key(
				Buffer.BSlot[buffslot]->BPage.index,
				Buffer.BSlot[buffslot]->BPage.page_hdr->lpage ),
			 dimensions );

		// find next match above this next_pagekey, place result in next_match
	 	memset( next_match, 0, sizeof(PU_int) * dimensions );

		if (false == H_nextmatch_RQ( Ret_set[set_id]->LB,
			Ret_set[set_id]->UB,
			next_match, next_pagekey,
			dimensions ))
		{
			delete [] next_pagekey;
			delete [] next_match;
			return false; // no higher matching hilbert codes
		}

		// clear the current page's flags if not in use by another Ret_set
        	if (Ret_set.size() - FreeRet_setList.size() == 1)
        	{
        		// this is the only ACTIVE Ret_set
        		i = -1;
        	}
        	else
        	{
        	       	// check all of the Ret_sets to see if any is accessing the same
               		// page that Ret_set[set_id] last accessed
        	       	// (could use an iterator here but this is simple)
               		for (i = Ret_set.size() - 1; i >= 0; i--)
        		{
        			if (i == set_id || !(Ret_set[i]->flags & ACTIVE))
        				continue;
        			if (Ret_set[i]->buffslot == Ret_set[set_id]->buffslot)
        			// another RET_SET is accessing the same page as Ret_set[set_id]
        				break;
        		}
        	}
            	if (i < 0)
          	{
          		// no other Ret_set is using the page
          		Buffer.BSlot[buffslot]->fix =
          		Buffer.BSlot[buffslot]->query = false;
          	}

		// find the page that may contain the match
		lpage = BT.idx_search( next_match );

		Ret_set[set_id]->buffslot = Buffer.b_page_retrieve( lpage );
		Buffer.BSlot[Ret_set[set_id]->buffslot]->query = true; // QUERY;

		i = Buffer.BSlot[Ret_set[set_id]->buffslot]->BPage.p_find_pageslot( Ret_set[set_id]->LB );
		if (i < 0)
			i = -i;
		Ret_set[set_id]->pos = i;
	}
}
