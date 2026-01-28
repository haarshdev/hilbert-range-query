// Copyright (C) Jonathan Lawder 2001-2011

#include <iomanip>
#include "db.h"
#ifdef __MSDOS__
	#include "..\hilbert\hilbert.h"
	#include "..\utils\utils.h"
#else
	#include "../hilbert/hilbert.h"
	#include "../utils/utils.h"
#endif
#include <stdio.h> // for db_getquery()
#include <stdlib.h> // for db_getquery()

#define		MAX_PAGES		UINT_MAX

using namespace std;

#ifdef JKLDEBUG
fstream fjunk1;
fstream fjunk2;
fstream fjunk3;
fstream fjunk4;
#endif

#ifdef xJKLDEBUGxxxx
fstream fjunk3;
#endif


int		Pages_retrieved, Disk_reads;
unsigned short SEED[] = {3000,1000,2000};


/*============================================================================*/
/***                   MED::constructor					    ***/
/*============================================================================*/
MED::MED( int dims, int p_entries )
{
	dimensions = dims;

	if ( THRESHOLD + EXTRA_RECORDS > p_entries )
	{
		p_entries = THRESHOLD + EXTRA_RECORDS;
	}

	#define MAX_DAT		(p_entries - 1)
	#define MIN_DAT		MAX_DAT / 2 - MAX_DAT / 10

	int i = p_median_calc_workspace( MAX_DAT + MIN_DAT + 10 ); // +10 for safety

	for ( int j = 0; j < i; j++ )
	{
		Hcode h( dimensions );	// can this go outside the for loop?
		MEDdata.push_back( h );
	}

#ifdef xJKLDEBUGxxxx
	cout << "(in MED::MED) Contents of MEDdata :\n";
	for ( int j = 0; j < i; j++ )
	{
		for (int i = 0; i < dims; i++)
			cout << MEDdata[j].hcode[i]<< " ";
		cout << j << endl;
	}
#endif

	MEDmap.resize( MAX_DAT );
	p_shuffle_medmap( MAX_DAT );


#ifdef xJKLDEBUGxxxx
	cout << "RAND_MAX = " << RAND_MAX << endl;
	cout << "(in MED::MED) Contents of MEDmap :\n";
	for (int i = 0; i < MEDmap.size(); i++)
		cout << MEDmap[i] << " ";
	cout << endl;
#endif

	#undef MAX_DAT
	#undef MIN_DAT
}

/*============================================================================*/
/***                   MED::p_median_calc_workspace			    ***/
/*============================================================================*/
/* How to calculate how much space is required for the original data and
	 for the intermediate and final results when calculating the median
	 of medians */
int MED::p_median_calc_workspace(int numentries)
{
	int num;

	if (numentries == 0)
		return 0;
	if (numentries <= MEDIAN)
		return numentries + 1;
	num = numentries / MEDIAN;
	if (numentries % MEDIAN)
		num++;
	return numentries + p_median_calc_workspace(num);
}

/*============================================================================*/
/***                   MED:p_shuffle_medmap				    ***/
/*============================================================================*/
/* this sets up an array with numbers 0 - MAX_DATA-1 which is then shuffled
   it's supposed to help with finding medians of medians which are more
   consistently close to the true median */
void MED::p_shuffle_medmap( int size )
{
	int i, a, b, temp;
	double da, db;

	seed48(SEED);

	for (i = 0, a = size; i < a; i++)
		MEDmap[i] = i;

	for (i = 0; i < 10000; i++)
	{
// not all random number generators produce the same numbers on different platforms

//		da = static_cast<double>(rand()); //		da = da / RAND_MAX * size;
//		db = static_cast<double>(rand());
//		db = db / RAND_MAX * size;

//		da = drand48() * size;
//		cout << da << " ";
//		db = drand48() * size;

		da = static_cast<double>(lrand48()) / (U_int)(1 << 31 ) * size;
		db = static_cast<double>(lrand48()) / (U_int)(1 << 31 ) * size;

		a  = static_cast<int>(da);
		b  = static_cast<int>(db);
		temp = MEDmap[a];
		MEDmap[a] = MEDmap[b];
		MEDmap[b] = temp;
	}
//	cout << endl;
}

/*============================================================================*/
/*                            MED::~MED	                          	      */
/*============================================================================*/
MED::~MED()
{
// not needed?????
	MEDmap.erase( MEDmap.begin(), MEDmap.end() );
	MEDdata.erase( MEDdata.begin(), MEDdata.end() );
}

/*============================================================================*/
/*                            DBASE::DBASE	                          	      */
/*============================================================================*/
DBASE::DBASE( string db_name, int dims, int bt_n_entries, int b_slots, int p_entries )
	:
	BT( db_name, dims, bt_n_entries ),
	dbMED( dims, p_entries ),
	Buffer( dims, b_slots, p_entries, this )
{
/*	BT = BTree( db_name, dims, bt_n_entries );*/
/*	BT = BTree(  );*/
	if ( THRESHOLD + EXTRA_RECORDS > p_entries )
	{
		cout << "WARNING: Specified page size of " << p_entries << " records is too small\n";
		cout << "using minimum default of " << THRESHOLD + EXTRA_RECORDS << endl;
		p_entries = THRESHOLD + EXTRA_RECORDS;
	}
	dbname			= db_name;
	dimensions		= dims;
	page_entries	= p_entries;
	bt_node_entries = bt_n_entries;
	num_Bslots		= b_slots;
}

/*============================================================================*/
/***                   DBASE::dbi_freepagelist_setup                        ***/
/*============================================================================*/
/* on opening a db:
   set up a linked list of free physical pages in the db, reading the
   list from a file */
void DBASE::dbi_freepagelist_setup()
{
	int	i, n = 0;
	int	page;
	string	fname = dbname + ".fpl";
	fstream f;

	f.open( fname.c_str(), ios::in | ios::binary );

	if (! f)
		errorexit("ERROR 1 reading .fpl in dbi_freepagelist_setup()\n");

	for (i = 0; i < NumFreePages; i++)
	{
		f.read( reinterpret_cast<char*>(&page), sizeof page );
		if (! f)
			errorexit("ERROR 2 reading .fpl in dbi_freepagelist_setup()\n");
		n += f.gcount();
		FreePageList.push( page );
	}
	if (n != NumFreePages * (int)(sizeof page))
		errorexit("ERROR 3 reading .fpl in dbi_freepagelist_setup()\n");

	/* make sure there's nothing more to read */
	(void) f.get();
	if (! f.eof())
		errorexit("ERROR 4 reading .fpl in dbi_freepagelist_setup()\n");
	f.close();
}

/*============================================================================*/
/***                   DBASE::dbi_freepagelist_save     	  	    ***/
/*============================================================================*/
/* on closing a db:
   write the free physical db page linked list to file */
void DBASE::dbi_freepagelist_save()
{
	int	i;
	string	fname = dbname + ".fpl";
	fstream f;

	if (FreePageList.empty())
	{   /* there are no free pages */
		f.close(); /* empties the file */
		return;
	}

	f.open(fname.c_str(), ios::out | ios::binary);

	if (! f)
		errorexit("ERROR 1 in dbi_freepagelist_save()\n");

	if ((int)FreePageList.size() != NumFreePages)
	{
		cerr << "FreePageList.size: " << FreePageList.size()
			<< endl << "NumFreePages: " << NumFreePages
			<< endl;
		errorexit("ERROR in dbi_freepagelist_save(): "
			"mis-match between FreePageList.size "
			"and NumFreePages\n");
	}

	int temp;

	for (i = NumFreePages - 1; i >= 0; i--)
	{
		temp = FreePageList.top();
		f.write(reinterpret_cast<char*>(&temp), sizeof temp);
		if (! fDB)
			errorexit("ERROR 3 in dbi_freepagelist_save()\n");
		FreePageList.pop();
	}

	f.close();
}

/*============================================================================*/
/***                   DBASE::db_freepagelist_dump			    ***/
/*============================================================================*/
/*  dumps the contents of the free page list and the value of NumFreePages:
 the count of the no. of entries in the list should equal NumFreePages
unless there's something wrong!
*/
void DBASE::db_freepagelist_dump()
{
	if (FreePageList.empty())
	{
		cout << "\nFree Page List is empty" << endl;
		cout << "\nNumFreePages: " << NumFreePages << endl;
		return;
	}

	int	i = FreePageList.size();
	vector<int> temp(i);

	cout << "Free page list contents:\n";
	for ( ; !FreePageList.empty(); FreePageList.pop(), i--)
	{
		cout << setw(6) << FreePageList.top();
		temp[i] = FreePageList.top();
	}

// now stack has been emptied into the vector; rebuild the stack - clumsy! should have
// written a copy constructor;
	vector<int>::iterator iter = temp.begin();
	vector<int>::iterator iter_end = temp.end();
	for ( ; iter != iter_end; iter++)
		FreePageList.push( *iter );

	cout << "\nFree page list size: " << FreePageList.size() << endl;
	cout << "\nNumber of free pages: " << NumFreePages << endl;
}

/*============================================================================*/
/***                   DBASE::dbi_get_last_value			    ***/
/*============================================================================*/
/*
this was used by the double indexing experiment ?
finds the maximum Hcode on a page
HU_int* DBASE::dbi_get_last_value( int lpage )
{
	int i, j, buffslot;
	HU_int *temp = new HU_int[dimensions], *max = new HU_int[dimensions];

	buffslot = Buffer.b_page_retrieve( lpage );
	max = ENCODE( max, Buffer.BSlot[buffslot]->BPage.data[1], dimensions );
	for (i = 2; i <= Buffer.BSlot[buffslot]->BPage.page_hdr->size; i++)
	{
		temp = ENCODE( temp, Buffer.BSlot[buffslot]->BPage.data[i], dimensions );
		for (j = dimensions - 1; j >=0; j--)
		{
			if (temp[j] > max[j])
			{
				keycopy( max, temp, dimensions );
				break;
			}
			if (temp[j] < max[j])
				break;
		}
	}
	Buffer.BSlot[buffslot]->fix = Buffer.BSlot[buffslot]->query;

	delete [] temp;

	return max;
}
*/

/*============================================================================*/
/*---                          dbi_create_info 				   ---*/
/*============================================================================*/
/* stores information about the database configuration */
bool DBASE::dbi_create_info()
{
	int info[INF_SIZE] = {0};

//	int	i;
	string	fname;

	info[0]  =  1;	// nextPID
	info[1]  =  0;	// NumFreePages
	info[2]  =  0;	// LastPage
	info[3]  =  dimensions;
	info[4]  =  page_entries;
	info[5]  =  bt_node_entries;

	if ( INF_SIZE != 6 )
		errorexit("ERROR 1 in dbi_create_info()\n");
//	info[15] = CURVE;
//	info[16] = ORDER;
//	info[17] = NON_KEY_INFO;
//	info[18] = NO_EXTRA_TOKENS;
//	info[19] = LTC;

	fname = dbname + ".inf";

	fstream f;
	f.open(fname.c_str(), ios::out | ios::binary);

	if (! f)
		errorexit("ERROR 2 in dbi_create_info()\n");

	f.write(reinterpret_cast<char*>(info), sizeof (info[0]) * INF_SIZE);

	if (! f)
		errorexit("ERROR 3 in dbi_create_info(): writing to .inf file\n");
	f.close();

	return true;
}

/*============================================================================*/
/***                   DBASE::dbi_open_info 				    ***/
/*============================================================================*/
/* checks whether an executable is consistent with an existing database */
bool DBASE::dbi_open_info()
{
	int info[INF_SIZE] = {0};

	string	fname;
	fstream	f;
	int	errors = 0;

	fname = dbname + ".inf";
	f.open (fname.c_str(), ios::in | ios::binary);
	if (! f)
		errorexit("ERROR 1 in dbi_open_info(): opening .inf file\n");

	f.read(reinterpret_cast<char*>(info), sizeof (info[0]) * INF_SIZE);
	if (! f)
		errorexit("ERROR 2 in dbi_open_info(): .inf file inconsistent\n");

	/* make sure there's nothing more to read */
	(void) f.get();
	if (! f.eof())
		errorexit("ERROR 3 in dbi_open_info(): .inf file inconsistent\n");
	f.close();

	nextPID = info[0];
	NumFreePages = info[1];

	if (nextPID < NumFreePages)	/* check for consistency */
		errorexit("ERROR 4 in dbi_open_info(): .inf file inconsistent\n");

	LastPage = info[2];

/* elements 3 - 11 available for future use */

	if (info[3] != dimensions)
	{
		cout << "No. of Dimensions: Database: " << info[3] << ", Executable: "
			<< dimensions << "\n";
		cout << "ERROR in dbi_open_info(): incompatible no. of dimensions\n";
		errors = 1;
	}
	if (info[4] != page_entries)
	{
		cout << "Database Page Size (bytes): Database: " << info[4]
			<< ", Executable: " << page_entries << "\n";
		cout << "ERROR in dbi_open_info(): incompatible no. of page_entries\n";
		errors = 1;
	}
	if (info[5] != bt_node_entries)
	{
		cout << "Index Node Size (bytes): Database: " << info[5]
			<< ", Executable: " << bt_node_entries << "\n";
		cout << "ERROR in dbi_open_info(): incompatible no. of bt_node_entries\n";
		errors = 1;
	}
/*
	if (info[17] != NON_KEY_INFO)
	{
		cout << "Storage of non-key data: Database: " << whether[info[17]]
			<< ", Executable: " << whether[NON_KEY_INFO] << "\n";
		cout << "ERROR in dbi_open_info(): "
				"incompatible #define: NON_KEY_INFO\n";
		errors = 1;
	}
	if (info[18] != NO_EXTRA_TOKENS)
	{
		cout << "Size of non-key data (words): Database: " << info[18]
			<< ", Executable: " << NO_EXTRA_TOKENS << "\n";
		cout << "ERROR in dbi_open_info(): incompatible #define: NO_EXTRA_TOKENS\n";
		errors = 1;
	}
 */
	if (errors)
		errorexit("Aborting program...\n");
	return true;
}

/*============================================================================*/
/*---                           db_create 				   ---*/
/*============================================================================*/
/* When a database is created, the first page is inserted. It has a 'dummy'
	 value of zero in the first slot. The page is put in the index. The
	 reason for this is that if an empty database was created and values
	 were subsequently inserted in descending key order, each value would
	 otherwise be placed on a separate page!
   Retruns 1, if successful, 0 otherwise.
   RETURN VALUE MUST BE TESTED. */
bool DBASE::db_create()
{
//	int		i;
	string		fname;

	// create database file
	// check a database with dbname doesn't aready exist:
	// create new database if it doesn't
	fname = dbname + ".db";

	fDB.open(fname.c_str(), ios::in);
	if (fDB)
	{
		fDB.close();
		cerr << "ERROR 1 in db_create(), "
			<< dbname << ".db already exists\n";
		return false;
	}

	fDB.open(fname.c_str(), ios::out | ios::binary);
	if (! fDB)
		errorexit("ERROR 2 in db_create(), can't create db\n");

	// create first page - as it's local it'll be destroyed when this func. finishes
	PAGE	page1( dimensions, page_entries );
	int		page_size = sizeof(pageheader_t) + page_entries * sizeof(U_int) * dimensions;

	fDB.write( reinterpret_cast<char*>(page1.raw_data), page_size );
	if (! fDB)
		errorexit("ERROR 3 in db_create(): writing to .db file\n");

	fDB.close();

	// insert first page into index: key = 0, lpage = 0
	HU_int *key = new HU_int[dimensions];
	// initialise key
	memset( key, 0, sizeof(U_int) * dimensions );
	BT.idx_insert_key( key, 0 );
	delete [] key;
	BT.idx_write();
	// delete the BTree since we'll read it from file when we open the database
	BT.free_root();

	// write info to .inf file
	dbi_create_info();

	// create free page list file: empty
	fname = dbname + ".fpl";

	fstream f;
	f.open( fname.c_str(), ios::out | ios::binary );
	if (! f)
		errorexit("ERROR 4 in db_create(), creating .fpl file\n");
	f.close();

	return true;
}

/*============================================================================*/
/***                   DBASE::db_open   				    ***/
/*============================================================================*/
bool DBASE::db_open()
{
	fstream	f;
	string 	fname;
	int		i;

#ifdef xJKLDEBUGxxxx
fjunk3.open("fjunk3.txt", ios::out);
#endif


	// open db
	fname = dbname + ".db";
	fDB.open (fname.c_str(), ios::in | ios::out | ios::binary);
	if (! fDB)
	{
		cerr << "ERROR 2 in db_open() - can't open " << fname << endl;
		return false;
	}

	// read data from .inf file
	dbi_open_info();

	// read btree index into memory
	if (nextPID > NumFreePages)
		i = BT.idx_read();
	// i is the number of pages indexed
	if (i != nextPID - NumFreePages)
		errorexit("ERROR 3 in db_open(): index not consistent with db\n");

	// set up the free page list
	if (NumFreePages > 0)
		dbi_freepagelist_setup();

/*	This is now dealt with by dbi_open_info()
	LastPage = (u2BYTES)idx_get_last_page(BT);*/

	Pages_retrieved = 0;
	Disk_reads = 0;

	return true;
}

/*============================================================================*/
/***                   DBASE::db_close  				    ***/
/*============================================================================*/
bool DBASE::db_close()
{
	int		i, offset;
	fstream	f;
	string	fname;
	int		page_size = sizeof(pageheader_t) + page_entries * sizeof(U_int) * dimensions;


#ifdef xJKLDEBUGxxxx
fjunk3.close();
#endif

//printf("sizeof(PAGE) = %i    page_size = %i\n",sizeof(PAGE), page_size);

	// flush changed pages in buffer to db
	for (i = 0; i < Buffer.num_Bslots; i++)
		if (Buffer.BSlot[i]->mod)
		{
//			offset = Buffer.BSlot[i]->BPage.page_hdr->lpage * sizeof(PAGE);
			offset = Buffer.BSlot[i]->BPage.page_hdr->lpage * page_size;
//printf("offset = %i\n",offset);

			fDB.seekp(offset, ios::beg);
			fDB.write(reinterpret_cast<char*>(Buffer.BSlot[i]->BPage.raw_data),
				page_size);
			if (! fDB)
				errorexit("ERROR in db_close(): writing to database\n");
		}

	// write out info : overwrite existing values of nextPID & NumFreePages
	fname = dbname + ".inf";
	f.open (fname.c_str(), ios::in | ios::out | ios::binary);
	if (! f)
		errorexit("ERROR in db_close(): in .inf file in db_close\n");

	f.write(reinterpret_cast<char*>(&nextPID), sizeof nextPID);
	f.write(reinterpret_cast<char*>(&NumFreePages), sizeof NumFreePages);
	f.write(reinterpret_cast<char*>(&LastPage), sizeof LastPage);

	if (! f)
		errorexit("ERROR in db_close(): re-writing to .inf in db_close\n");
	f.close();

	// write out index
	BT.idx_write();

	// write out free page list (this also frees storage)
	dbi_freepagelist_save();

	fDB.close();

	// don't free index, buffer and MED - this is done by DBASE destructor

	return true;
}

/*============================================================================*/
/***                   DBASE::db_data_present				    ***/
/*============================================================================*/
/* For finding whether a FULLY specified point exists.
   Returns true or false, depending on whether 'data' exists on a page. */
bool DBASE::db_data_present( PU_int *data )
{
//	long	offset;
	int	buffslot, pageslot;
	int	lpage;
	HU_int	*key = new HU_int[dimensions];

	key = ENCODE( key, data, dimensions );
	lpage = BT.idx_search( key );
	buffslot = Buffer.b_page_retrieve( lpage );
	if (lpage != Buffer.BSlot[buffslot]->BPage.page_hdr->lpage)
	{
/*
		BT.idx_dump( );
		cout << "\ndb_data_present : Key :\n";
		for (int j = 0; j < dimensions; j++)
			cout << key[j];
		cout << endl << "Page nos." << endl;
		cout << lpage << "   "  << Buffer.BSlot[buffslot]->BPage.page_hdr->lpage << endl;
		cout << "buffslot : " << buffslot << endl;
*/
		errorexit("ERROR 1 in db_data_present(): index inconsistent\n");
	}
	pageslot = Buffer.BSlot[buffslot]->BPage.p_find_pageslot( data ); //@@@
	delete [] key;
	if (pageslot > 0) // data is present
	{
		return true;
	}
	return false;
}

/*============================================================================*/
/***                   DBASE::db_data_insert				    ***/
/*============================================================================*/
// Returns page size (occupancy) or ALREADY_PRESENT (or possibly MAX_DATA?)
int DBASE::db_data_insert( PU_int* data )
{
//	int		buffslot;
	int	lpage;
	HU_int*	key = new HU_int[dimensions];

	for (int i = 0; i < dimensions; i++)
		if (data[i] == _UNSPECIFIED_)
		{
			cout << "Coordinate " << i << " is unspecified: not allowed!" << endl;
			return false;
		}

	key = ENCODE( key, data, dimensions );
	lpage = BT.idx_search( key );
	delete [] key;

	return Buffer.b_data_insert( data, lpage );
}

/*============================================================================*/
/***                   DBASE::db_data_delete				    ***/
/*============================================================================*/
// Returns MIN_DATA or NOT_PRESENT or the size of the page
int DBASE::db_data_delete( PU_int* data )
{
//	int		buffslot;
	int	lpage;
	HU_int*	key = new HU_int[dimensions];

	if (NumFreePages == nextPID)
		errorexit("ERROR 1 in db_data_delete(): database is empty\n");

	key = ENCODE( key, data, dimensions );
	lpage = BT.idx_search( key );
	delete [] key;

	return Buffer.b_data_delete( data, lpage );
}

/*============================================================================*/
/***                   DBASE::dbi_get_new_page	 	 		    ***/
/*============================================================================*/
/* returns an empty page's number, adding one to the end of the db file if
 page not taken from the free page list */
// inline
int DBASE::dbi_get_new_page()
{
	int	newpage;

	if (FreePageList.size() > 0)
	{   /* re-use an old page */
		newpage = FreePageList.top();
		FreePageList.pop();
		NumFreePages--;
	}
	else /* create a new page */
	{
		if (nextPID == (int)MAX_PAGES)
			errorexit("ERROR 1 in dbi_get_new_page(): database full\n");
		newpage = nextPID;
		nextPID++;
		/* write a page-sized block of memory to the end of the file */
		PAGE	emptypage( dimensions, page_entries );
		int		page_size = sizeof(pageheader_t) + page_entries * sizeof(U_int) * dimensions;

		fDB.seekp(0, ios::end);
		if (! fDB)
			errorexit("ERROR 2 in dbi_get_new_page(): writing to database\n");
		fDB.write( reinterpret_cast<char*>(emptypage.raw_data), page_size );
		if (! fDB)
			errorexit("ERROR 3 in dbi_get_new_page(): writing to database\n");
	}
	return newpage;
}

/*============================================================================*/
/*                            		                          	      */
/*                            Diagnostic functions                    	      */
/*                            		                          	      */
/*============================================================================*/

/*============================================================================*/
/*                            db_info					      */
/*============================================================================*/
void DBASE::db_info()
{
  cout << "\nnumber of dimensions : " << dimensions << "\n";
  cout << "number of B-Tree node entries : " << bt_node_entries << "\n";
  cout << "number of pages in buffer : " << num_Bslots << "\n";
  cout << "number of records on page : " << page_entries << "\n";
}

/*============================================================================*/
/*                            db_key_dump				      */
/*============================================================================*/
// similar to idx_dump but it outputs page no. - key pairs read from the
// database pages themselves rather than from the index.
// should produce the same output as idx_dump unless there's something wrong!
void DBASE::db_key_dump( string fname )
{
	int		i, j, buffslot;
	BTnode		*p = BT.root;
	fstream		f;
	U_int*		idx;	// no need to allocate storage

	f.open( fname.c_str(), ios::out | ios::binary );

	f << "From the database\n\n";

	if (!p)
	{
		f << "Index empty\n";
		f.close();
		return;
	}
	while (!(p->in_HDR->flags & isLEAF))
		p = p->in_HDR->firstptr;

	do
	{
		for (i = 1; i <= p->lf_HDR->size; i++)
		{
			buffslot = Buffer.b_page_retrieve( p->lf_ENTRY[i]->lpage );
			f << "PAGE NO.: " << Buffer.BSlot[buffslot]->BPage.page_hdr->lpage << endl;
			f << "KEY : ";
			idx =  Buffer.BSlot[buffslot]->BPage.index;
			for (j = dimensions - 1; j >= 0; j--)
				f << setw(15) << idx[j];
			f << endl;
			Buffer.BSlot[buffslot]->fix = 0;
		}
		p = p->lf_HDR->nextptr;
	} while (p);
	f.close();
}

/*============================================================================*/
/*                            db_data_dump				      */
/*============================================================================*/
// dumps the entire contents of a database to a file, DATA ONLY
// database must already be open
void DBASE::db_data_dump( string fname, char type )
{
	int 		i, j, width, buffslot;
	HU_int		*key = new U_int[dimensions];
	BTnode		*p;
	fstream		f;
	U_int		*temp;

	// initialise hcode
	memset( key, 0, sizeof(U_int) * dimensions );

	f.open( fname.c_str(), ios::out | ios::binary );

	if (type == 'd')
		width = 11;
	else
		width = 34;

	p = BT.root->idxi_find_leaf( key );
	while (p)
	{
		for (int k = 1; k <= p->lf_HDR->size; k++)
		{
			buffslot = Buffer.b_page_retrieve( p->lf_ENTRY[k]->lpage );
			f << "PAGE NO. : " << setw(5) << Buffer.BSlot[buffslot]->BPage.page_hdr->lpage
				<< "  INDEX : ";
			temp = Buffer.BSlot[buffslot]->BPage.index;
			for (i = dimensions-1; i >= 0; i--)
				f << setw(width) << temp[i];
			f << "  NUM RECORDS : " << Buffer.BSlot[buffslot]->BPage.page_hdr->size;
			f << endl;

			// dump the contents of the page
			for (i = 1; i <= Buffer.BSlot[buffslot]->BPage.page_hdr->size; i++)
			{
				temp = Buffer.BSlot[buffslot]->BPage.data[i];
				for (j = 0; j < dimensions; j++)
					if (type == 'd')
						f << setw(width) << temp[j];
					else
						f << setw(width) << int2bins( temp[j], WORDBITS );
			#if NON_KEY_INFO  // not translated to C++ yet
				for (j = 0; j < NO_EXTRA_TOKENS; j++)
					if (type == 'd')
						fprintf(f, "%-*u", width,
							DB.DB_Buffer.BSlot[buffslot].BPage.data[i].extra[j]);
					else
						fprintf(f, "%-*s", width,
							int2bins(DB.DB_Buffer.BSlot[buffslot].BPage.data[i].extra[j],
							WORDBITS));
			#endif
				f << endl;
			}
			Buffer.BSlot[buffslot]->fix = 0;
		}
		p = p->lf_HDR->nextptr;
	}
	f.close();
}

/*============================================================================*/
/*                            db_batch_update				      */
/*============================================================================*/
void DBASE::db_batch_update( char type, int num )
{
	int	i, j = 0;
	PU_int	*k = new PU_int[dimensions];

	if (type == 'i')
	{
		seed48(SEED);

		cout << num << " record insertion\n";
		for (i = 0; i < num; i++)
		{
			for (j = 0; j < dimensions; j++)
			{
				k[j] = lrand48();
			}

			db_data_insert( k );
			if (!(i%50000))
			{
				cout << setw(6) << i / 50000 << " " << nextPID << endl;
			}
		}
		cout << setw(6) << i / 50000 << " " << nextPID << endl;
	}
	delete [] k;
}

/*============================================================================*/
/*                            db_getquery				      */
/*============================================================================*/
// storage MUST be allocated to 'p' before this is called
bool DBASE::db_getquery( PU_int *p, string message)
{
	char input[256], junk[20];
	int count, i, n;
	unsigned int temp;

	for(;;)
	{
		cout << message;
		fgets(input,sizeof input,stdin);
		if (input[0] == 'q')
		{
			break;
		}

		for ( i = count = 0; count < dimensions; i += n, count++)
		{
			if (sscanf(&input[i], "%s%n", junk, &n) != 1)
				break;
			if (sscanf(junk,"%u", &temp) == 1)
				p[count] = (U_int)temp;
			else
				p[count] = _UNSPECIFIED_;
		}

		if (count != dimensions)
			cout << "Incorrect number of elements entered\n";
		else
			return true;
	}
	return false;
}

/*============================================================================*/
/*                            db_getrange				      */
/*============================================================================*/
bool DBASE::db_getrange( PU_int *LB, PU_int *UB )
{
	int i;

	for (;;)
	{
		if (false == db_getquery( LB, "Enter LOWER bounds of the range query; 'q' to abort\n" ))
		{
			break;
		}
		if (false == db_getquery( UB, "Enter UPPER bounds of the range query: 'q' to abort\n" ))
		{
			break;
		}

		for (i = 0; i < dimensions; i++)
		{
			if (LB[i] > UB[i])
			{
				break;
			}
		}

		if (i < dimensions)
			printf("Invalid range: all lower bounds must be less than upper bounds\n");
		else
			return true;
	}
	return false;
}

/*============================================================================*/
/*                            db_querytest				      */
/*============================================================================*/
void DBASE::db_querytest(void)
{
	char	c;
	string	junk;
	int	i, j;
	PU_int	*LB = new U_int[dimensions];
	PU_int	*UB = new U_int[dimensions];
	PU_int	*result = new U_int[dimensions];
	HU_int	*hresult = new U_int[dimensions];

	for (;;)
	{
		cout << "\nWhat do you want to do? (choose a letter - `q' to quit)\n"
		     <<	"\tp : for partial match query\n"
		     <<	"\tr : for range query\n\n> ";

		cin >> c;
		getline( cin, junk );

		switch (c)
		{
		case 'p':

			if (false == db_getquery( LB,
				"Enter query (any letter except 'q' for 'unspecified'), 'q' to abort\n" ))
			{
				goto end;
			}
			if (true == db_open_set( LB, &i ))
			{
 				cout << setw(12 * dimensions) << "Matching point coordinates";
 				cout << "     "
				     << setw(12 * dimensions) << "Matching point sequence number\n";

				while (true == db_fetch_another( i, result ))
				{
					for (j = 0; j < dimensions; j++)
					{
						cout << setw(12) << result[j];
					}
					cout << "    "; // endl;
					// output the point's hcode
					hresult = ENCODE( hresult, result, dimensions );
					for (j = dimensions - 1; j >= 0; j--)
					{
						cout << setw(12) << hresult[j];
					}
					cout << endl;
				}
				db_close_set( i );
			}
			else
			{
				cout << "ERROR in partial match query\n";
			}

			break;

		case 'r':

			if (false == db_getrange( LB, UB ))
			{
				goto end;
			}
			if (true == db_range_open_set( LB, UB, &i ))
			{
 				cout << setw(12 * dimensions) << "Matching point coordinates";
 				cout << "     "
				     << setw(12 * dimensions) << "Matching point sequence number\n";

				while(true == db_range_fetch_another( i, result ))
				{
					// output the point
					for (j = 0; j < dimensions; j++)
					{
						cout << setw(12) << result[j];
					}
					cout << "    "; // endl;
					// output the point's hcode
					hresult = ENCODE( hresult, result, dimensions );
					for (j = dimensions - 1; j >= 0; j--)
					{
						cout << setw(12) << hresult[j];
					}
					cout << endl;
				}
				db_close_set( i );
			}
			else
			{
				cout << "ERROR in range query\n";
			}

			break;

		default:

			goto end;
		}
	}
end:
	delete [] LB;
	delete [] UB;
	delete [] result;
	delete [] hresult;
}




#if 0

/*============================================================================*/
/***                   DBASE::db_PM_delete				    ***/
/*============================================================================*/
/* similar to db_data_delete except it accepts Point P with _UNSPECIFIED_
   values, in which a SET of Points may be deleted
   returns the number of 'records' deleted */
int DBASE::db_PM_delete (Point P)
{
	int count = 0, set_id;
	Point retp;

#if ! ALLOW_UPDATES
	printf("WARNING! attempting to use db_PM_delete() when "
		"'ALLOW_UPDATES' not defined as 'true'\n"
		"Deletion not allowed\n");
	return 0;
#endif
	if (! db_open_set(P, &set_id))
		return 0;

	while(db_fetch_another(set_id, retp))
		if (db_data_delete(retp) != NOT_PRESENT)
			count++;

	db_close_set(set_id);

	return count;
}

/*============================================================================*/
/***                   DBASE::db_PM_data_present			    ***/
/*============================================================================*/
/* For finding whether any match to a PARTIALLY specified point exists.
   Returns true or false, depending on whether 'data' exists on a page. */
int DBASE::db_PM_data_present( PU_int* ) )
{
	int set_id, any_there;

	if (! db_open_set(P, &set_id))
		return false;

	any_there = db_fetch_another(set_id, P);

	db_close_set(set_id);

	return any_there;
}

#endif // #if 0
