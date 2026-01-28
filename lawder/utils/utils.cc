// Copyright (C) Jonathan Lawder 2001-2011

#include <string>
#include <cassert>

#ifdef __MSDOS__
	#include "..\gendefs.h"
#else
	#include "../gendefs.h"
#endif


/*============================================================================*/
/*                            Hcode::Hcode                         	      */
/*============================================================================*/
Hcode::Hcode( int dims )
{
	for ( int i = 0; i < dims; i++ )
		hcode.push_back( (U_int)0 );
/*	for ( int i = 0; i < dims; i++ )
		cout << hcode[i];
	cout << endl;*/
}

/*============================================================================*/
/*                            keycopy                          	      */
/*============================================================================*/
void keycopy ( HU_int *destination, const HU_int * const source, const int dim )
{
	for ( int i = 0; i < dim; i++ )
	{
		destination[i] = source[i];
	}
}

//*============================================================================*/
/*                            keycopy                          	      */
/*============================================================================*/
void keycopy ( HU_int *destination, const Hcode& source )
{
	for ( int i = source.hcode.size() - 1; i >= 0; i-- )
	{
		destination[i] = source.hcode[i];
	}
}

/*============================================================================*/
/*                            keycopy                          	      */
/*============================================================================*/
void keycopy ( Hcode& destination, const HU_int * const source )
{
	for ( int i = destination.hcode.size() - 1; i >= 0; i-- )
	{
		destination.hcode[i] = source[i];
	}
}
 #if 0
        /*============================================================================*/
        /*                            getstorage				      */
        /*============================================================================*/
        void *getstorage(int size)
        {
        	void *p;
        	assert ((p=malloc(size)) != NULL);
        /*	memset (p, NULL, size);*/
        	return p;
        }
 #endif
	
/*============================================================================*/
/*                            errorexit					      */
/*============================================================================*/
void errorexit(string s)
{
	cerr << s;
#ifdef JKLDEBUG
	assert(0);	
#else	
	exit(EXIT_FAILURE);
#endif	
}

/*============================================================================*/
/*                            int2bins			 		      */
/*============================================================================*/
char* int2bins(unsigned int num, int bits)
{
	int i;
	static char bin[WORDBITS+1];

	for (i = 0; i < (int)(sizeof bin); i++)
		bin[i] = '\0';

	if (bits > 32)
	{
		cout << "error in int2bins()\n";
		exit(1);
	}

	for (i = 0; i < bits ; i++)
	{
		if (num % 2)
			bin[bits-1-i] = '1';
		else
			bin[bits-1-i] = '0';
		num = num /2;
	}
	return bin;
}

