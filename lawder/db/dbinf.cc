// Copyright (C) Jonathan Lawder 2001-2011

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include<fstream>
#include<iostream>

#ifdef __MSDOS__
	#include "..\gendefs.h"
	#include "..\UTILS\utils.h"
#else
	#include "../gendefs.h"
	#include "../UTILS/utils.h"
#endif

//.................

DBASE		DB;

//.................

/*============================================================================*/
/*                            db_info	 				      */
/*============================================================================*/
/* outputs the contents of a database's '.inf' file */
static int db_info(string dbname)
{
	u2BYTES info[INF_SIZE] = {0};

	fstream	f;
	int	i, errors = 0;
	char	*curve[] = {"HILBERT", "MOORE", "ZORDER", "TRUE_Z", "SNAKE", "GRAYCODEF", "GRAYCODEA", "GRAYCODEB"};
	char	*order[] = {"BYATTRIBUTE", "BYHCODE", "UNSORTED"};
	char	*whether[] = {"FALSE", "TRUE"};
	string	fname = dbname + ".inf";


	f.open (fname.c_str(), ios::in | ios::binary);

	if (! f)
	{
		printf("ERROR 1 in db_info(): .inf file inconsistent\n");
		return 0;
	}

	f.seekg(0);
	f.read(reinterpret_cast<char*>(info), sizeof (info[0]) * INF_SIZE);

	if (! f)
	{
		printf("ERROR 2 in db_info(): .inf file inconsistent\n");
		return 0;
	}
	/* make sure there's nothing more to read */
	(void) f.get();
	if (! f.eof())
	{
		printf("ERROR 3 in db_info(): .inf file inconsistent\n");
		return 0;
	}
	f.close();

	printf("NEXT PAGE ID:..... %u\n", info[0]);
	printf("NO. OF FREE PAGES: %u\n", info[1]);
	printf("LAST PAGE:........ %u\n", info[2]);
	printf("NO. OF DIMS:...... %u\n", info[12]);
	printf("DB PAGE SIZE:..... %u\n", info[13]);
	printf("INDEX PAGE SIZE:.. %u\n", info[14]);
	printf("CURVE TYPE:....... %s\n", curve[info[15]]);
	printf("RECORD ORDER:..... %s\n", order[info[16]]);
	printf("NON KEY INFO:..... %s\n", whether[info[17]]);
	printf("NO_EXTRA_TOKENS:.. %u\n", info[18]);
	printf("LTC PROVIDED:..... %s\n", whether[info[19]]);

	return 1;
}

int main (void)
{
	char dbname[15];

	printf("This program reads and outputs the contents\n"
		"of a database's '.inf' file\n\n");
	printf("Enter the name of the database (excluding file extension):  ");
	fgets(dbname, sizeof(dbname), stdin);
	dbname[strlen(dbname)-1] = '\0';
	printf("\n\n");
	db_info(dbname);

	return 0;
}
