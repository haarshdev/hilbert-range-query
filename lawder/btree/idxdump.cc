// Copyright (C) Jonathan Lawder 2001-2011

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __MSDOS__
	#include "..\gendefs.h"
	#include "..\BTREE\btree.h"
#else
	#include "../gendefs.h"
	#include "../BTREE/btree.h"
#endif
#include "db.h"
#include "dbtest.h"


#ifdef __MSDOS__
	#define		DATADIR		"..\\data\\"
#else
	#define		DATADIR		"../DATA/"
#endif

extern	DBASE		*DB;			/* database page index */

int main (void)
{
	char dbname[15];
	FILE *f1, *f2;

	printf("This program dumps the contents of a database index to file\n"
		"and dumps the page no. - index values read from the "
		"database\npages themselves\n");
	printf("Enter the name of the database (excluding file extension):  ");
	fgets(dbname, sizeof(dbname),stdin);
	dbname[strlen(dbname)-1] = '\0';
	printf("\nOutput will be sent to index.dump and key.dump\n");

	db_open(dbname);
	idx_dump(DB->BT, f1);
	db_key_dump(DB->BT, f2);
	db_freepagelist_dump();

	return 0;
}