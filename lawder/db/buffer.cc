// Copyright (C) Jonathan Lawder 2001-2011

#include "buffer.h"
#include "db.h"
#ifdef __MSDOS__
	#include "..\utils\utils.h"
#else
	#include "../utils/utils.h"
#endif
#include <stdio.h>

using namespace std;


#define debug 1

// db flags
#define		ALREADY_PRESENT		-1
#define		NOT_PRESENT		-1

// work needs to be done to allow this to be set to 1
// BUFF_PAGE's would need access to DBASE* (or possibly just Ret_sets)
// Ret_sets would ideally want to know lpages - but they do know buffslots
#define		ALLOW_UPDATES		0

/*============================================================================*/
/*                     BUFF_PAGE::BUFF_PAGE				    */
/*============================================================================*/
BUFF_PAGE::BUFF_PAGE( int dims, int p_entries, MED *m )
	: BPage( dims, p_entries, m )
{
	dimensions = dims;
	bp_page_entry_bytes = sizeof(U_int) * dimensions;
	mod = fix = query = false;	// not needed?
	lru = 0;				// not needed?
}

/*============================================================================*/
/*                     BUFF_PAGE::~BUFF_PAGE				    */
/*============================================================================*/
// BUFF_PAGE::~BUFF_PAGE(){}

/*============================================================================*/
/***                   BUFF_PAGE::bp_insert_on_page			    ***/
/*============================================================================*/
//int BUFF_PAGE::bp_insert_on_page(DBASE & DB, Point& DATA)
int BUFF_PAGE::bp_insert_on_page( const PU_int * const DATA )
{
	int pageslot, nobj;

	pageslot = BPage.p_find_pageslot(DATA);  //@@@
	if (pageslot > 0)
		return ALREADY_PRESENT;
	if (pageslot == 0)
		errorexit("ERROR in bp_insert_on_page\n");
	pageslot *= -1;

	nobj = bp_page_entry_bytes * (BPage.page_hdr->size - pageslot + 1);
	memmove(BPage.data[pageslot + 1], BPage.data[pageslot], nobj);
	keycopy ( BPage.data[pageslot], DATA, dimensions );

	BPage.page_hdr->size++;

	mod = true; // CHANGED;

#if ALLOW_UPDATES
// ..................need to make 'query' boolean and Ret_set visiblee
	if (true == query)
	{
		for (int i = DB->Ret_set.size() - 1; i >= 0; i--)
		{
			if (DB->Ret_set[i]->flags & ACTIVE &&
				BPage.page_hdr->lpage == DB->Ret_set[i]->lpage &&
				pageslot < DB->Ret_set[i]->pos)
			{
				DB->Ret_set[i]->pos++;
			}
		}
	}
#endif

	return static_cast<int>(BPage.page_hdr->size);
}

/*============================================================================*/
/***                   BUFF_PAGE::bp_delete_from_page			    ***/
/*============================================================================*/
/* Returns NOT_PRESENT or the size of the page.
   Serially searches page for data to delete. */
//int BUFF_PAGE::bp_delete_from_page(DBASE & DB, Point& DATA)
int BUFF_PAGE::bp_delete_from_page( PU_int *DATA )
{
	int pageslot, nobj;

	pageslot = BPage.p_find_pageslot(DATA);  //@@@
	if (pageslot < 1)
		return NOT_PRESENT;

	nobj = bp_page_entry_bytes * (BPage.page_hdr->size - pageslot);
	memmove(BPage.data[pageslot], BPage.data[pageslot + 1], nobj);
	BPage.page_hdr->size--;
	mod = true; // CHANGED;
	fix = query; // finished with it

#if ALLOW_UPDATES
// ..................need to make 'query' boolean and Ret_set visiblee
	if (true == query)
	{
		for (int i = DB->Ret_set.size() - 1; i >= 0; i--)
		{
			if (DB->Ret_set[i]->flags & ACTIVE &&
				BPage.page_hdr->lpage == DB->Ret_set[i]->lpage &&
				pageslot < DB->Ret_set[i]->pos)
			{
				DB->Ret_set[i]->pos--;
			}
		}
	}
#endif

	return static_cast<int>(BPage.page_hdr->size);
}

/*============================================================================*/
/*                    BUFFER::BUFFER				    */
/*============================================================================*/
BUFFER::BUFFER( int dims, int b_slots, int p_entries, DBASE *db )
{
	DB =	db;

	if ( THRESHOLD + EXTRA_RECORDS > p_entries )
	{
		p_entries = THRESHOLD + EXTRA_RECORDS;
	}

	for ( int i = 0; i < b_slots; i++ )
	{
		BSlot.push_back( new BUFF_PAGE( dims, p_entries, &DB->dbMED ) );
	}

	dimensions = dims;
	num_Bslots = b_slots;
	// (b_slots - 1) because num_Bslots are numbered in the range [ 0 .. num_Bslots-1 ]
	free_Bslots = b_slots-1;
	LRU = 0;
	b_page_bytes = sizeof(pageheader_t) + p_entries * sizeof(U_int) * dimensions;
}

/*============================================================================*/
/*                    BUFFER::~BUFFER				    */
/*============================================================================*/
BUFFER::~BUFFER()
{
	for (int i = BSlot.size() - 1; i >= 0; i--)
		delete BSlot[i];

/* not needed  ????

	BSlot.erase( BSlot.begin(), BSlot.end() );

	FreeBufferList...........
	LRU_idx.erase( LRU_idx.begin(), LRU_idx.end() );
	Buff_idx.erase( Buff_idx.begin(), Buff_idx.end() );
*/
}

/*============================================================================*/
/***                   BUFFER::inc_LRU				    ***/
/*============================================================================*/
// usage: BSlot[x]->lru = inc_LRU(x)
U_int BUFFER::inc_LRU(int buffslot)
{
		// cout << "calling inc_LRU - LRU: " << LRU << " BSlot[buffslot]->lru: "
			// << BSlot[buffslot]->lru << endl; // Debugging

	if (! LRU_idx.erase(BSlot[buffslot]->lru) && LRU > 0)
		errorexit("ERROR in BUFFER::inc_LRU - missing LRU entry in LRU_idx map\n");
	LRU++;
	LRU_idx.insert(valTypeB1(LRU, buffslot));
	return LRU;

#if debug
	if (Buff_idx.size() != LRU_idx.size())
	{
		cout << "Buff_idx.size(): " << Buff_idx.size()
			<< " LRU_idx.size(): " << LRU_idx.size() << endl;
		errorexit("ERROR in inc_LRU()\n");
	}
#endif
}

/*============================================================================*/
/***                   BUFFER::in_Buffer				    ***/
/*============================================================================*/
// returns the buffslot number of a page or -1 if it's not in the buffer
inline int BUFFER::in_Buffer(int lpage)
{
	map<int, int>::iterator iter = Buff_idx.find(lpage);
	if (iter != Buff_idx.end())
		return (*iter).second;
	return -1;
}

/*============================================================================*/
/***                   BUFFER::Buff_idx_insert				    ***/
/*============================================================================*/
void BUFFER::Buff_idx_insert(int lpage, int buffslot)
{

		// cout << "!!!!!!! in Buff_idx_insert !!! Buff_idx.size(): " << Buff_idx.size()
			// << " LRU_idx.size(): " << LRU_idx.size() << endl; // Debugging
		// cout<< "------ " << lpage << " ----- " << buffslot <<endl; // Debugging

//#if debug
	if (Buff_idx.count(lpage))
		errorexit("ERROR 1 in Buff_idx_insert - trying to insert a key that's already in index\n");
//#endif

	Buff_idx.insert(valTypeB1(lpage, buffslot));
	LRU++;
	LRU_idx.insert(valTypeB1(LRU, buffslot));
	BSlot[buffslot]->lru = LRU;

//#if debug
//	int nbsize = Buff_idx.size();
//	int nlsize = LRU_idx.size();

	if (Buff_idx.size() != LRU_idx.size())
	{
		cout << "Buff_idx.size(): " << Buff_idx.size()
			<< " LRU_idx.size(): " << LRU_idx.size() << endl;
		errorexit("ERROR 2 in Buff_idx_insert() - sizes of Buff_idx and LRU_idx are different\n");
	}
//#endif
}

/*============================================================================*/
/***                   BUFFER::Buff_idx_erase				    ***/
/*============================================================================*/
void BUFFER::Buff_idx_erase(int lpage, int buffslot)
{
	Buff_idx.erase(lpage);
	LRU_idx.erase(BSlot[buffslot]->lru);

#if debug
	if (Buff_idx.size() != LRU_idx.size())
	{
		cout << "Buff_idx.size(): " << Buff_idx.size()
			<< " LRU_idx.size(): " << LRU_idx.size() << endl;
		errorexit("ERROR in Buff_idx_erase()\n");
	}
#endif
}

/*============================================================================*/
/***                   BUFFER::b_get_buffer_slot	  		    ***/
/*============================================================================*/
/* returns a free buffer slot, swapping a page out if necessary */
inline int BUFFER::b_get_buffer_slot()
{
	int buffslot;
#ifdef JKLDEBUGxxx
	if (! FreeBufferList.empty())
		cout << "buffer list not empty - FreeBufferList.top : " << FreeBufferList.top() << endl;
	else
		cout << "buffer list empty ";
	cout << " free_Bslot : " << free_Bslots << " Buff_idx.size : " << Buff_idx.size() << endl;
#endif

	if (! FreeBufferList.empty())
	{
#ifdef JKLDEBUGxxx
		cout << "\tusing buffslot from FreeBufferList\n\n";
#endif
		buffslot = FreeBufferList.top();
		FreeBufferList.pop();
		return buffslot;
	}

// don't do this comaprison because Buff_idx may be out of date (awaiting updating)
// when this function is called when it is called twice in succession to get 2
// buffslots at the same time
//	if (Buff_idx.size() < num_Bslots)
	if (free_Bslots >= 0)
	{
#ifdef JKLDEBUGxxx
		cout << "\tusing buffslot not previously used\n\n";
#endif
//		if (free_Bslots < 0)
//			errorexit( "ERROR in b_get_buffer_slot()\n" );
		buffslot = free_Bslots;
		free_Bslots--;
		return buffslot;
	}
#ifdef JKLDEBUGxxx
	cout << "Calling b_swapout" << endl;
#endif

	return b_swapout();
}

/*============================================================================*/
/***                   BUFFER::b_swapout				    ***/
/*============================================================================*/
int BUFFER::b_swapout()
{
	map<U_int, int>::iterator iter = LRU_idx.begin();

	for ( ; iter != LRU_idx.end(); iter++)
		if (false == BSlot[iter->second]->fix)
			break;
	if (iter == LRU_idx.end())
		errorexit("ERROR 1 in b_swapout() - can't find a page to swapout\n");

	int buffslot = iter->second;
	if (true == BSlot[buffslot]->mod)
	{
		int offset = BSlot[buffslot]->BPage.page_hdr->lpage * b_page_bytes;

		DB->fDB.seekp( offset, ios::beg );
		DB->fDB.write( reinterpret_cast<char*>(BSlot[buffslot]->BPage.raw_data), b_page_bytes );
		if (! DB->fDB)
			errorexit("ERROR 2 in b_swapout(): writing to database\n");

//		BSlot[buffslot]->mod = BSlot[buffslot]->query = 0; - do in page_retrieve()
	}
	Buff_idx_erase(BSlot[buffslot]->BPage.page_hdr->lpage, buffslot);

// no need to bother putting the free bufferslot on the stack - it's going to be used immediately

	return buffslot;
}

/*============================================================================*/
/***                   BUFFER::b_page_retrieve				    ***/
/*============================================================================*/
int BUFFER::b_page_retrieve( int lpage )
{
	int offset, buffslot = in_Buffer(lpage);

	if (buffslot > -1)
		BSlot[buffslot]->lru = inc_LRU(buffslot);
	else	// not in buffer
	{
		buffslot = b_get_buffer_slot();

		/* read the page in and insert in buffer index */
		offset = lpage * b_page_bytes;

		DB->fDB.seekg( offset, ios::beg );
		DB->fDB.read( reinterpret_cast<char*>(BSlot[buffslot]->BPage.raw_data),
			b_page_bytes );
		if (! DB->fDB)
			errorexit("ERROR 1 in page_retrieve(): "
				"reading database\n");

		Buff_idx_insert( lpage, buffslot );

			// cout << "in b_page_retrieve - Buff_idx.size: "
			// 	<< Buff_idx.size() << endl; // Debugging

		BSlot[buffslot]->BPage.page_hdr->lpage = lpage;
		BSlot[buffslot]->mod = BSlot[buffslot]->query = false;

//		Disk_reads++;
	}
	BSlot[buffslot]->fix = true; // FIXED;

//	Pages_retrieved++;

	return buffslot;
}

/*============================================================================*/
/***                   BUFFER::b_data_insert				    ***/
/*============================================================================*/
int BUFFER::b_data_insert( PU_int *data, int lpage )
{
	int i, buffslot = b_page_retrieve( lpage );

	if (lpage != BSlot[buffslot]->BPage.page_hdr->lpage)
		errorexit("ERROR 1 in b_data_insert(): index inconsistent\n");

			// debug: cout << "lpage " << lpage << " buffslot " << buffslot << endl;

#if ALLOW_UPDATES
	// don't insert data on a page being queried if this will lead to overflow
//	if (BSlot[buffslot]->query && BSlot[buffslot]->BPage.page_hdr->size >= MAX_DATA - 1)
	if (true == BSlot[buffslot]->query &&
		BSlot[buffslot]->BPage.page_hdr->size >=
		BSlot[buffslot]->BPage.p_page_entries - 2)
	{
		cout << "WARNING in b_data_insert(): attempting to\n" <<
			"insert data on a retrieval set's current page which is " <<
			"full\n - insertion abandoned\n";
		return -1;
	}
#endif
#if !ALLOW_UPDATES
	if (true == BSlot[buffslot]->query)
	{
		cout << "Cannot insert data as page is in use by a query\n";
		return -1;
	}
#endif

	i = BSlot[buffslot]->bp_insert_on_page( data );

//	if (i == MAX_DATA)
	if (i == BSlot[buffslot]->BPage.p_page_entries - 1)
	{
//		printf( "i = %i, BSlot[buffslot]->BPage.p_page_entries - 1 = %i\n", i, BSlot[buffslot]->BPage.p_page_entries - 1);
		i = b_process_overflow( buffslot );  // deals with flags
	}
	else
		BSlot[buffslot]->fix = BSlot[buffslot]->query; // finished with it

#if debug
	if (Buff_idx.size() != num_Bslots - free_Bslots - 1 - FreeBufferList.size())
	{
		cout << "Buff_idx.size(): " << Buff_idx.size()
			<< " num_Bslots: " << num_Bslots
			<< " free_Bslots: " << free_Bslots
			<< " FreeBufferList.size(): " << FreeBufferList.size() << endl;
		errorexit("ERROR 2 in b_data_insert(): buffer index size error\n");
	}
#endif

	return i;
}

/*============================================================================*/
/***                   BUFFER::b_data_delete				    ***/
/*============================================================================*/
int BUFFER::b_data_delete( PU_int *data, int lpage )
{
	int buffslot = b_page_retrieve( lpage );
	int MIN_DAT = (int)((double)(BSlot[0]->BPage.p_page_entries - 1) * 4 / 10);

#if ALLOW_UPDATES
	if (true == BSlot[buffslot]->query && BSlot[buffslot]->BPage.page_hdr->size <= MIN_DAT)
	{
		cout << "WARNING in b_data_delete(): attempting to\n" <<
			"delete data from a retrieval set's current page which is " <<
			"at minimum occupancy\n - deletion abandoned\n";
		return -1;
	}
#endif
#if !ALLOW_UPDATES
	if (true == BSlot[buffslot]->query)
	{
		cout << "Cannot delete data as page is in use by a query\n";
		return -1;
	}
#endif

	int i = BSlot[buffslot]->bp_delete_from_page( data );

	if (i > 0 && i <= MIN_DAT)
	/* Logically, the test should be (i == MIN_DATA) but it's safer to say
	   (i <= MIN_DATA) as pages could conceivably be smaller than MIN_DATA
	   when data is unsorted or ordered by attribute if a suitable median
	   was not found during an earlier attempt to process underflow.
	   We also need to bear in mind that bp_delete_from_page will return
	   FALSE (0) if there was an attempt to delete from a minimal sized
	   page where data is unsorted and where the page is 'locked' by a
	   retrieval set in which case the deletion will have failed and we
	   certainly don't want to process underflow.
	   Also, more importantly, bp_delete_from_page will return
	   -1 if the data was not present, regardless of data order/retrieval
	   sets etc
	*/
		b_process_underflow( buffslot );        /* deals with flags */

//#if debug
	if ((int)Buff_idx.size() != num_Bslots - free_Bslots - 1)
	{
		cout << "Buff_idx.size(): " << Buff_idx.size()
			<< " num_Bslots: " << num_Bslots
			<< " free_Bslots: " << free_Bslots << endl;
		errorexit("ERROR 2 in b_data_delete(): buffer index size error\n");
	}
//#endif

	return i;
}

/*============================================================================*/
/***                   BUFFER::b_process_overflow	  		    ***/
/*============================================================================*/
/* the median goes on the new page and becomes its key
   returns new page's buffslot
   the full page is split between 2 buffer slots (therefore 3 slots are used)
   this would not get called for a page which was being queried */

int BUFFER::b_process_overflow( int oflowslot )
{
	int	newleft, newright, newlpage;

	if (true == BSlot[oflowslot]->query)
		errorexit ("ERROR in b_process_overflow() - about to split a query page\n");

	newleft = b_get_buffer_slot();
	newright = b_get_buffer_slot();

	/* no need to fix newright */
	if (newleft == newright || newleft == oflowslot || newright == oflowslot)
		errorexit("ERROR 1 in dbi_process_overflow():"
				  "\n\tover-populated page has been swapped out\n");

	/* get a lpage number for new page - it comes from the free page list or is added to db */
	newlpage = DB->dbi_get_new_page();

// 	split oflowslot between newleft & newright - this also assigns with page nos. & page keys
	BSlot[oflowslot]->BPage.p_split_page( BSlot[newleft]->BPage, BSlot[newright]->BPage, newlpage );

//	release oflowslot - MUST be done before insertions into the indexes
	Buff_idx_erase( BSlot[oflowslot]->BPage.page_hdr->lpage, oflowslot ); // deals with LRU_idx too
	BSlot[oflowslot]->mod = BSlot[oflowslot]->fix = BSlot[oflowslot]->query = false;
	FreeBufferList.push( oflowslot );

//	also deals with LRU values
	Buff_idx_insert( BSlot[newleft]->BPage.page_hdr->lpage, newleft );
	Buff_idx_insert( newlpage, newright );

//	put new page in index
	DB->BT.idx_insert_key( BSlot[newright]->BPage.index, newlpage );

//	deal with flags
	BSlot[newleft]->fix = BSlot[newright]->fix = false;
	BSlot[newleft]->query = BSlot[newright]->query = false;
	BSlot[newleft]->mod = BSlot[newright]->mod = true; // CHANGED;

//	if the page just split was the last page, update DB_LastPage
	if (DB->LastPage == BSlot[newleft]->BPage.page_hdr->lpage)
		DB->LastPage = BSlot[newright]->BPage.page_hdr->lpage;

	return newright;
}

/*============================================================================*/
/***                   BUFFER::b_process_underflow			    ***/
/*============================================================================*/
// this method feels a bit long and clumsy
// - I don't like the use of flags 'left' and 'right'
// - calls to midx_search() are effectively repeated within any subsequent call to b_page_retrieve()
int BUFFER::b_process_underflow( int uflowslot )
{
	int		left = 0, right = 0, PageRight, PageLeft, buffleft, buffright;
	int	uflowpage = BSlot[uflowslot]->BPage.page_hdr->lpage;

// top priority: MERGE with a page that's ALREADY in the buffer
	/* get next page's lpage */
	PageRight = DB->BT.idx_get_next( BSlot[uflowslot]->BPage.index, uflowpage );
	if (PageRight > 0)
	{
		right++;  // a right hand page exist

		map<int, int>::iterator iter = Buff_idx.find( PageRight );
		if (iter != Buff_idx.end())
		{
			buffright = (*iter).second;

			right++;  // page is in buffer
			if (false == BSlot[buffright]->query)
			{
				right++;  //  page not in use by query
				if (BSlot[uflowslot]->BPage.page_hdr->size + BSlot[buffright]->BPage.page_hdr->size <=
					(int)((double)(BSlot[0]->BPage.p_page_entries-1)*0.9))
				{
					b_merge_pages( uflowslot, buffright );
					return 0;
				}
			}
		}
	}

	/* get prev page's lpage */
	PageLeft = DB->BT.idx_get_prev( BSlot[uflowslot]->BPage.index, uflowpage );
	if (PageLeft > 0)
	{
		left++;  // a left hand page exists

		map<int, int>::iterator iter = Buff_idx.find( PageLeft );
		if (iter != Buff_idx.end())
		{
			buffleft = (*iter).second;

			left++;  // page is in buffer
			if (false == BSlot[buffleft]->query)
			{
				left++;  //  page not in use by query
				if (BSlot[uflowslot]->BPage.page_hdr->size + BSlot[buffleft]->BPage.page_hdr->size <=
					(int)((double)(BSlot[0]->BPage.p_page_entries-1)*0.9))
				{
					b_merge_pages( buffleft, uflowslot );
					return 0;
				}
			}
		}
	}

	if (left == 0 && right == 0) // this is the only page in the database
		if (DB->nextPID - DB->NumFreePages > 1)
			errorexit("ERROR 1 in b_process_underflow(): cannot find page to merge with\n");


// second priority: SHIFT from a page that's ALREADY in the buffer
	if (right == 3)
	{
		b_shift_from_right( uflowslot, buffright );
		return 1;
	}
	if (left == 3)
	{
		b_shift_from_left( uflowslot, buffleft );
		return 1;
	}
// third priority: bring a page into memory - it won't be in the middle of being searched!!!
	if (right == 1)
	{
		buffright = b_page_retrieve( PageRight );
		if (BSlot[uflowslot]->BPage.page_hdr->size + BSlot[buffright]->BPage.page_hdr->size <=
			(int)((double)(BSlot[0]->BPage.p_page_entries-1)*0.9))
		{
			b_merge_pages( uflowslot, buffright );
			return 0;
		}
		else
		{
			b_shift_from_right( uflowslot, buffright );
			return 1;
		}

	}
	if (left == 1)
	{
		buffleft = b_page_retrieve( PageLeft);
		if (BSlot[uflowslot]->BPage.page_hdr->size + BSlot[buffleft]->BPage.page_hdr->size <=
			(int)((double)(BSlot[0]->BPage.p_page_entries-1)*0.9))
		{
			b_merge_pages( buffleft, uflowslot );
			return 0;
		}
		else
		{
			b_shift_from_left( uflowslot, buffleft );
			return 1;
		}

	}
	// therefore left == 2, right == 2; ie left AND right pages are in buffer AND
	// both are being queried

	cout << "Unable to process underflow in page number " <<  uflowpage << endl;
	return -1;   // NB return values not currently used
}

/*============================================================================*/
/***                   BUFFER::b_merge_pages	  			    ***/
/*============================================================================*/
/*
   The right page's contents are moved to the left.
   It doesn't matter which one is underpopulated.
   The right page is deleted.
   left and right are buffer slots. */
int BUFFER::b_merge_pages( int left, int right )
{
	int	newleft;

	newleft = b_get_buffer_slot();

	/* no need to fix newleft */
	if (newleft == left || newleft == right)
		errorexit("ERROR 1 in b_merge_pages(): page taking "
				  "part in\nmerge swapped out.of buffer\n");

	BSlot[newleft]->BPage.p_merge_pages( BSlot[left]->BPage, BSlot[right]->BPage );

// NB currently the following situation can't happen as it's prevented within b_process_underflow()
#if ALLOW_UPDATES  // make this into a separate function at some point
	if (true == BSlot[left]->query || true == BSlot[right]->query)
	{
		int x;
		for (i = DB->Ret_set.size() - 1; i >= 0; i--)
		{
			if (!(DB->Ret_set[i]->flags & ACTIVE))
				continue;
			if (DB->Ret_set[i]->buffslot == left)
			{
				DB->Ret_set[i]->buffslot = newleft;
				x = BSlot[left]->BPage.p_find_pageslot(BSlot[left]->BPage.data[DB->Ret_set[i]->pos]);  //@@@

				if (x < 1)
					errorexit("ERROR 3 in b_merge_pages() : "
								"retrieval set current data not found\n");
				DB->Ret_set[i]->pos = x;
			}
			else
				if (DB->Ret_set[i]->buffslot == right)
				{
					DB->Ret_set[i]->buffslot = newleft;
					DB->Ret_set[i]->lpage = BSlot[left]->BPage.page_hdr->lpage;
					x = BSlot[left]->BPage.p_find_pageslot(BSlot[left]->BPage.data[DB->Ret_set[i]->pos]);  //@@@

					if (x < 1)
						errorexit("ERROR 3 in b_merge_pages() : "
									"retrieval set current data not found\n");
					DB->Ret_set[i]->pos = x;
				}
		}
	}
#endif

//	adjust flags
	BSlot[newleft]->mod = true; // CHANGED;
	BSlot[newleft]->fix = BSlot[left]->query || BSlot[right]->query;
	BSlot[newleft]->query = BSlot[left]->query || BSlot[right]->query;

//	release left & right pages' buffer slots
	Buff_idx_erase( BSlot[left]->BPage.page_hdr->lpage, left );
	Buff_idx_erase( BSlot[right]->BPage.page_hdr->lpage, right );
	BSlot[left]->mod = BSlot[left]->fix = BSlot[left]->query = false;
	BSlot[right]->mod = BSlot[right]->fix = BSlot[right]->query = false;
	FreeBufferList.push( left );
	FreeBufferList.push( right );

//	MUST be done AFTER releasing left & right buffer slots, also deals with LRU values
	Buff_idx_insert(BSlot[newleft]->BPage.page_hdr->lpage, newleft);

//	add right lpage to free page list
	DB->FreePageList.push( BSlot[right]->BPage.page_hdr->lpage );
	DB->NumFreePages++;

	DB->BT.idx_delete_key( BSlot[right]->BPage.index, BSlot[right]->BPage.page_hdr->lpage );

//	update LastPage if necessary
	if (DB->LastPage == BSlot[right]->BPage.page_hdr->lpage)
		DB->LastPage = BSlot[newleft]->BPage.page_hdr->lpage;

	return 1;
}

/*============================================================================*/
/***                   BUFFER::b_shift_from_left			    ***/
/*============================================================================*/
/* left and right are bufferslots, right is the under-populated one.
   carry the median to the right page (it becomes its new key).  */
int BUFFER::b_shift_from_left( int left, int right )
{
	int	newleft, newright;

#if ALLOW_UPDATES
	/* bp_delete_from_page should ensure that the right page is not part
	   of this shifting process */
	if (true == BSlot[left]->query)
	{
		printf("WARNING in dbi_shift_from_left(): attempting to "
			"move data from a\nretrieval set's current page to an "
			"under-populated page\n - not allowed! - page re-balance "
			"aborted\n");
		return -1;
	}
#endif

	newleft = b_get_buffer_slot();

//	no need to fix newleft
	if (newleft == left || newleft == right)
		errorexit("ERROR 1 in b_shift_from_left(): page taking "
				  "part in\nmerge swapped out.of buffer\n");

	newright = b_get_buffer_slot();
//	no need to fix newright
	if (newright == left || newright == right || newright == newleft)
		errorexit("ERROR 2 in b_shift_from_left(): page taking "
				  "part in\nmerge swapped out.of buffer\n");

	if (BSlot[left]->BPage.p_shift_from_left( BSlot[right]->BPage,
			BSlot[newleft]->BPage, BSlot[newright]->BPage ))
		return 1; // the median was found on the right page instead of the left

	return b_admin( left, right, newleft, newright ); // update buffer, index, etc
}

/*============================================================================*/
/***                   BUFFER::b_shift_from_right	  		    ***/
/*============================================================================*/
/* left and right are bufferslots, left is the under-populated one.
   leave the median on the right page (it becomes its new key). */
int BUFFER::b_shift_from_right( int left, int right )
{
	int	newleft, newright;

#if ALLOW_UPDATES
	/* bp_delete_from_page should ensure that the left page is not part
	   of this shifting process */
	if (true == BSlot[right]->query)
	{
		printf("WARNING in b_shift_from_right(): attempting to "
			"move data from a\nretrieval set's current page to an "
			"under-populated page\n - not allowed! - page re-balance "
			"aborted\n");
		return -1;
	}
#endif

	newright = b_get_buffer_slot();

//	no need to fix newright
	if (newright == left || newright == right)
		errorexit("ERROR 1 in b_shift_from_right(): page taking "
				  "part in\nmerge swapped out.of buffer\n");
	newleft = b_get_buffer_slot();

//	no need to fix newright
	if (newleft == left || newleft == right || newleft == newright)
		errorexit("ERROR 2 in b_shift_from_right(): page taking "
				  "part in\nmerge swapped out.of buffer\n");

	if (BSlot[left]->BPage.p_shift_from_right(BSlot[right]->BPage,
			BSlot[newleft]->BPage, BSlot[newright]->BPage))
		return 1; // the median was found on the left page instead of the right

	return b_admin( left, right, newleft, newright ); // update buffer, index, etc
}

/*============================================================================*/
/***                   BUFFER::b_admin					    ***/
/*============================================================================*/
// tidy up after shift_from_left or shift_from_right
int BUFFER::b_admin( int left, int right, int newleft, int newright )
{
//	release vacated buffer slots
	Buff_idx_erase( BSlot[left]->BPage.page_hdr->lpage, left );
	Buff_idx_erase( BSlot[right]->BPage.page_hdr->lpage, right );
	BSlot[left]->mod = BSlot[left]->fix = BSlot[left]->query = false;
	BSlot[right]->mod = BSlot[right]->fix = BSlot[right]->query = false;
	FreeBufferList.push( left );
	FreeBufferList.push( right );

//	update buffer index - AFTER erasing old slots -  also deals with LRU
	Buff_idx_insert( BSlot[newleft]->BPage.page_hdr->lpage, newleft );
	Buff_idx_insert( BSlot[newright]->BPage.page_hdr->lpage, newright );

//	deal with flags
	BSlot[newleft]->fix = BSlot[newright]->fix = false;
	BSlot[newleft]->query = BSlot[newright]->query = false;
	BSlot[newleft]->mod = BSlot[newright]->mod = true; // CHANGED;

	DB->BT.idx_delete_key( BSlot[right]->BPage.index, BSlot[right]->BPage.page_hdr->lpage );
	DB->BT.idx_insert_key( BSlot[newright]->BPage.index, BSlot[newright]->BPage.page_hdr->lpage );

	return 0;
}