// Copyright (C) Jonathan Lawder 2001-2011

#ifdef DEV
#ifdef __MSDOS__
	#include "..\db\db.h"
#else
	#include "../db/db.h"
#endif
#else
	#include "db.h"
#endif

#include <iomanip>


using namespace std;

/*============================================================================*/
/*                            randomi					      */
/*============================================================================*/
/* this is for inserting a lot of randomly generated records at the same time */

int randomi( DBASE *DB )
{
	string	junk;
	int	num;

	cout << "How many records to insert? : ";
	cin >> num;
	getline(cin, junk);

	DB->db_batch_update( 'i', num );
	return 0;
}

/*============================================================================*/
/*                            main					      */
/*============================================================================*/
int main(void)
{
	DBASE	*DB;
	char	c;
	string	dbname, junk;
	int	DIMS, bt_node_entries, n_bslots, p_entries;

	cout << "Enter database name : ";
	cin >> dbname;
	getline(cin, junk);

// have we already created the db?
	string	fname;
	fstream	f;

	fname = dbname + ".inf";
	f.open (fname.c_str(), ios::in | ios::binary);
	if (! f)
	{
		cout << "Creating a new database ...\n";
		cout << "\nEnter no. of dimensions (between 3 and 31): ";
		cin >> DIMS;
		if( DIMS < 3 || DIMS > 31 )
		{
			cout << "Invalid number of dimensions : " << DIMS << "\n";
			return 0;
		}
		getline(cin, junk);
		cout << "\nEnter no. of BTree node entries (suggested value is 10) : ";
		cin >> bt_node_entries;
		getline(cin, junk);
		cout << "\nEnter size of buffer (pages) (suggested value is 10) : ";
		cin >> n_bslots;
		getline(cin, junk);
		cout << "\nEnter size of page (records) (suggested value is 35 - 1500) : ";
		cin >> p_entries;
		getline(cin, junk);

		DB = new DBASE( dbname, DIMS, bt_node_entries, n_bslots, p_entries );

		DB->db_info();

		DB->db_create();
	}
	else
	{
		cout << "Opening an existing database ...\n";
		int info[INF_SIZE] = {0};
		f.read(reinterpret_cast<char*>(info), sizeof (info[0]) * INF_SIZE);
		f.close();
		cout << "\nEnter size of buffer (pages) (suggested value is 10) : ";
		cin >> n_bslots;
		getline(cin, junk);
		DB = new DBASE( dbname, info[3], info[5], n_bslots, info[4] );

		DB->db_info();

	}
 	string idxdump = dbname + "_idx.dump";
 	string keydump = dbname + "_key.dump";
 	string datadump = dbname + "_data.dump";

	DB->db_open();

	for(;;)
	{
		cout << "\n---- TOP LEVEL MENU ---- " <<
			"What do you want to do? " <<
			"(enter a number or `q' to quit)\n\n" <<
			"  1 : insert randomly generated data\n" <<
			"  2 : ad hoc updating\n" <<
			"  3 : ad hoc querying\n" <<
			"  4 : dump index to file : " << idxdump << endl <<
			"  5 : dump database keys to file : " << keydump << endl <<
			"  6 : dump database to file : " << datadump << endl <<
			"  > ";


		cin >> c;
		getline(cin, junk);

		switch (c)
		{
			case '1':  randomi( DB );   break;	// insert random data
			case '2':  cout << "ad hoc updating - not yet implemented!!\n";    break;
			case '3':  DB->db_querytest();    break;
			case '4':  DB->BT.idx_dump( idxdump );     break;	//
			case '5':  DB->db_key_dump( keydump );     break;	//
			case '6':  DB->db_data_dump( datadump, 'd' );     break;	//

			default :  cout << "closing database " << dbname << "\n"; DB->db_close();  return 0;
		}
	}

	return 0;
}




