// Copyright (C) Jonathan Lawder 2001-2011

#include <algorithm>	// for sort()
#include <string>
	
#include "page.h"
#include "db.h"

#ifdef __MSDOS__
	#include "..\hilbert\hilbert.h"
	#include "..\utils\utils.h"
#else
	#include "../hilbert/hilbert.h"
	#include "../utils/utils.h"
#endif

#ifdef xJKLDEBUGxxxx
#include <iomanip>
#endif

#define 	MAX_DATA		(p_page_entries - 1)

using namespace std;


#ifdef xJKLDEBUGxxxx
extern fstream fjunk3;
#endif

/*============================================================================*/
/*                            PAGE::PAGE	                          	      */
/*============================================================================*/
// p_page_entries - includes the index entry
PAGE::PAGE( int dims, int p_entries, MED *m ) {
	dimensions = dims;
	p_page_entries = p_entries;
	p_page_entry_size = sizeof(U_int) * dimensions;
	p_page_bytes = sizeof(pageheader_t) + p_page_entries * p_page_entry_size;
	M = m;

 	// the data block that makes up a page (inc. page header and index entry)
	// +1 for the index
/*	raw_data = (unsigned char*)malloc( sizeof(pageheader_t) + p_page_entries * p_page_entry_size );*/
	raw_data = new unsigned char[p_page_bytes];
	memset( raw_data, '\0', p_page_bytes );

	// create the array of pointers to hcodes in a page
/*	data = (U_int**)malloc( sizeof(U_int*) * p_page_entries );*/
	data = new HU_int*[p_page_entries];

	// setup the pointers into raw_data
	// NB data[0] is not used; it overlap with 'index'
	// NB the BTnodehdr 'occupies' the same amount of space as a <Hkey,XX>
	page_hdr = (pageheader_t*)raw_data;

	// CHECK CHECK CHECK !!!
	unsigned char *base_ptr = raw_data + sizeof(pageheader_t);
	index = (U_int*)base_ptr;
	for ( int i = 0; i < p_page_entries; i++ )
	{
		data[i] = (U_int*)base_ptr + i * dimensions;
	}
}

/*============================================================================*/
/*                            PAGE::PAGE	                          	      */
/*============================================================================*/
// slightly more efficient than the normal constructor - used in db creation
// p_page_entries - includes the index entry
PAGE::PAGE( int dims, int p_entries ) {
	dimensions = dims;
	p_page_entries = p_entries;
	p_page_entry_size = sizeof(U_int) * dimensions;
	p_page_bytes = sizeof(pageheader_t) + p_page_entries * p_page_entry_size;

 	// the data block that makes up a page (inc. page header and index entry)
	// +1 for the index
/*	raw_data = (unsigned char*)malloc( sizeof(pageheader_t) + p_page_entries * p_page_entry_size );*/
	raw_data = new unsigned char[p_page_bytes];
	memset( raw_data, '\0', p_page_bytes );

	// create the array of pointers to hcodes in a page
/*	data = (U_int**)malloc( sizeof(U_int*) * p_page_entries );*/
	data = new HU_int*[1];

	// setup the pointers into raw_data
	// NB data[0] is not used; it overlap with 'index'
	// NB the BTnodehdr 'occupies' the same amount of space as a <Hkey,XX>
	page_hdr = (pageheader_t*)raw_data;

	// CHECK CHECK CHECK !!!
	unsigned char *base_ptr = raw_data + sizeof(pageheader_t);
	index = (U_int*)base_ptr;
//	data[0] = (U_int*)base_ptr + dimensions;
	data[0] = (U_int*)base_ptr;
}

/*============================================================================*/
/*                            PAGE::~PAGE	                          	      */
/*============================================================================*/
PAGE::~PAGE() {
	delete [] data;
	delete raw_data;
}

/*============================================================================*/
/***                   PAGE::p_find_pageslot				    ***/
/*============================================================================*/
int PAGE::p_find_pageslot( const PU_int * const search_point )
{
	int	low = 1, mid, high = page_hdr->size, i, j;
	PU_int *page_point;	// no storage required

	if (high == 0) /* page is empty */
		return -1;

	while (high >= low)
	{
		mid = (high + low) / 2;
		page_point = data[mid];
		for (i = j = 0; i < dimensions; i++)
		{
			if (search_point[i] < page_point[i])
			{
				high = mid - 1;
				break;
			}
			if (search_point[i] > page_point[i])
			{
				j = 1;
				low  = mid + 1;
				break;
			}
		}
		if (i == dimensions)
			return mid;
	}

	/* not found - return a neg. no. which is the slot data would go in */
	if (j == 1)			/* key > temp */
		mid++;
	return -mid;
}

/*============================================================================*/
/***                   PAGE::p_split_page	  				    ***/
/*============================================================================*/
// split the page which calls this between left and right pages
// (left and right are empty pages in the buffer
void PAGE::p_split_page(PAGE& left, PAGE& right, int newlpage)
{
	// allocate storage to median, pass it to median funcs, then delete it

	int		i, j, k, lcount, rcount;

	Hcode	median = p_find_median();

	/* re-distribute the data */
	for (i = lcount = rcount = 0; i < page_hdr->size; i++)
	{
		for (j = 0, k = dimensions - 1; k >= 0; k--)
		/* we're dealing with CODES, not ATTRIBUTES, here:
			therefore init k to DIM - 1 */
		{
			if (M->MEDdata[i].hcode[k] == median.hcode[k])
				continue;
			if (M->MEDdata[i].hcode[k] < median.hcode[k])
			{	j = -1; break;}  /* a < b */
			else
			{	j = 1; break;} /* a > b */
		}

		if (j == -1)
		{
			lcount++;
			keycopy( left.data[lcount], data[i + 1], dimensions );
		}
		else
		{
			rcount++;
			keycopy( right.data[rcount], data[i + 1], dimensions );
		}
	}

	left.page_hdr->lpage = page_hdr->lpage;
	left.page_hdr->size = lcount;

	keycopy( left.index, index, dimensions );

	right.page_hdr->lpage = newlpage;
	right.page_hdr->size = rcount;

	keycopy( right.index, median );
}

/*============================================================================*/
/***                   PAGE::p_merge_pages				    ***/
/*============================================================================*/
// merge the contents of left and right pages onto the page which calls this
void PAGE::p_merge_pages(PAGE& left, PAGE& right)
{
	int	L, R, NL, i, j;
	PU_int *ldata, *rdata;	// no storage required

	for (L = NL = R = 1; L <= left.page_hdr->size && R <= right.page_hdr->size; NL++)
	{
		for (j = 0, i = 0; i < dimensions; i++)
		/* we're dealing with ATTRIBUTES, not CODES, here:
			therefore init i to 0 */
		{
			ldata = left.data[L];
			rdata = right.data[R];
			if (ldata[i] < rdata[i])
			{j = -1;  break;  }
			if (ldata[i] > rdata[i])
			{	j = 1;  break;  }
		}
		if (j < 0)	/* L < R */
		{
			keycopy( data[NL], left.data[L], dimensions );
			L++;
		}
		else
		{
			keycopy( data[NL], right.data[R], dimensions );
			R++;
		}
	}

	for (; L <= left.page_hdr->size; L++, NL++)
		keycopy( data[NL], left.data[L], dimensions );

	for (; R <= right.page_hdr->size; R++, NL++)
		keycopy( data[NL], right.data[R], dimensions );

	if (NL - 1 != left.page_hdr->size + right.page_hdr->size)
		errorexit("ERROR 1 in merge_pages()\n");

	/* adjust details of new left page */
	page_hdr->size = left.page_hdr->size + right.page_hdr->size;
	page_hdr->lpage = left.page_hdr->lpage;
	keycopy( index, left.index, dimensions );
}

/*============================================================================*/
/***                   PAGE::p_shift_from_left				    ***/
/*============================================================================*/
// merge the contents of left and right pages onto the page which calls this
// called by the left-hand page
int PAGE::p_shift_from_left(PAGE& right, PAGE& newleft, PAGE& newright )
{
	int		L, NL, R, NR, i, j;
	PU_int	*ldata, *rdata;	// no storage required

	Hcode	median = p_find_median_left( right );	// returns an hcode

	for (i = dimensions - 1; i >= 0; i--)
	{
		if (median.hcode[i] == UINT_MAX)
			continue;
		break;
	}
	if (i == -1)
		return 1; /* median is in the right page - do nothing */

	/* move data from left and right to newleft or newright
		while (MEDmap[L-1] < median)
			DB_BuffPage[newleft].data[NL] = DB_BuffPage[left].data[L];
		while (DB_BuffPage[right].data[R] < DB_BuffPage[left].data[L])
			DB_BuffPage[newright].data[NR] = DB_BuffPage[right].data[R];
		DB_BuffPage[newright].data[NR] = DB_BuffPage[left].data[L]; */

	for (L = NL = R = NR = 1; L <= page_hdr->size && R <= right.page_hdr->size; L++, NR++)
	{
		for (; L <= page_hdr->size; L++, NL++)
		{
			for (j = 0, i = dimensions - 1; i >= 0; i--)
			/* we're dealing with CODES, not ATTRIBUTES, here:
				therefore init i to dimensions - 1 */
			{
				if (M->MEDdata[L-1].hcode[i] < median.hcode[i])
				{	j = -1;  break;  }
				if (M->MEDdata[L-1].hcode[i] > median.hcode[i])
				{	j = 1;  break;  }
			}
			if (j < 0)	/* L < median */
				keycopy( newleft.data[NL], data[L], dimensions );
			else
				break;
		}

		if (L > page_hdr->size)
			break;

		for (; R <= right.page_hdr->size; R++, NR++)
		{
			for (j = 0, i = 0; i < dimensions; i++)
			/* we're dealing with ATTRIBUTES, not CODES, here:
				therefore init i to 0 */
			{
				rdata = right.data[R];
				ldata = data[L];
				if (rdata[i] < ldata[i])
				{	j = -1;  break;  }
				if (rdata[i] > ldata[i])
				{	j = 1;  break;  }
			}
			if (j < 0)	/* R < L */
				keycopy( newright.data[NR], right.data[R], dimensions );
			else
				break;
		}
		keycopy( newright.data[NR], data[L], dimensions );
	}

	for (; L <= page_hdr->size; L++)
	{
		for (j = 0, i = dimensions - 1; i >= 0; i--)
		/* we're dealing with CODES, not ATTRIBUTES, here:
			therefore init i to dimensions - 1 */
		{
			if (M->MEDdata[L - 1].hcode[i] < median.hcode[i])
			{	j = -1;  break;  }
			if (M->MEDdata[L - 1].hcode[i] > median.hcode[i])
			{	j = 1;  break;  }
		}
		if (j < 0)
		{
			keycopy( newleft.data[NL], data[L], dimensions );
			NL++;
		}
		else
		{
			keycopy( newright.data[NR], data[L], dimensions );
			NR++;
		}
	}

	for (; R <= right.page_hdr->size; R++, NR++)
		keycopy( newright.data[NR], right.data[R], dimensions );

	if (NL + NR - 2 != page_hdr->size + right.page_hdr->size)
		errorexit("ERROR 1 in p_shift_from_left()\n");

	newleft.page_hdr->lpage = page_hdr->lpage;
	newleft.page_hdr->size = NL - 1;
	keycopy( newleft.index, index, dimensions );

	newright.page_hdr->lpage = right.page_hdr->lpage;
	newright.page_hdr->size = NR - 1;
	keycopy( newright.index, median );

	return 0;
}

/*============================================================================*/
/***                   PAGE::p_shift_from_right				    ***/
/*============================================================================*/
// merge the contents of left and right pages onto the page which calls this
// called by the left-hand page
int PAGE::p_shift_from_right(PAGE& right, PAGE& newleft, PAGE& newright )
{
	int		L, NL, R, NR, i, j;
	PU_int	*ldata, *rdata;	// no storage required

	Hcode	median = p_find_median_right( right );	// returns an hcode

	for (i = dimensions - 1; i >= 0; i--)
	{
		if (median.hcode[i] == 0)
			continue;
		break;
	}
	if (i == -1)
		return 1; /* median is in the left page - do nothing */

	for (L = NL = R = NR = 1; L <= page_hdr->size && R <= right.page_hdr->size; R++, NL++)
	{
		for (; R <= right.page_hdr->size; R++, NR++)
		{
			for (j = 0, i = dimensions - 1; i >= 0; i--)
			/* we're dealing with CODES, not ATTRIBUTES, here:
				therefore init i to dimensions - 1 */
			{
				if (M->MEDdata[page_hdr->size + R - 1].hcode[i] < median.hcode[i])
				{	j = -1;  break;  }
				if (M->MEDdata[page_hdr->size + R - 1].hcode[i] > median.hcode[i])
				{	j = 1;  break;  }
			}
			if (j >= 0)	/* R >= median */
				keycopy( newright.data[NR], right.data[R], dimensions );
			else
				break;
		}

		if (R > right.page_hdr->size)
			break;

		for (; L <= page_hdr->size; L++, NL++)
		{
			for (j = 0, i = 0; i < dimensions; i++)
			/* we're dealing with ATTRIBUTES, not CODES, here:
				therefore init i to 0 */
			{
				ldata = data[L];
				rdata = right.data[R];
				if (ldata[i] < rdata[i])
				{	j = -1;  break;  }
				if (ldata[i] > rdata[i])
				{	j = 1;  break;  }
			}
			if (j < 0)	/* L < R */
				keycopy( newleft.data[NL], data[L], dimensions );
			else
				break;
		}
		keycopy( newleft.data[NL], right.data[R], dimensions );
	}

	for (; R <= right.page_hdr->size; R++)
	{
		for (j = 0, i = dimensions - 1; i >= 0; i--)
		/* we're dealing with CODES, not ATTRIBUTES, here:
			therefore init i to dimensions - 1 */
		{
			if (M->MEDdata[page_hdr->size + R - 1].hcode[i] < median.hcode[i])
			{	j = -1;  break;  }
			if (M->MEDdata[page_hdr->size + R - 1].hcode[i] > median.hcode[i])
			{	j = 1;  break;  }
		}
		if (j < 0)
		{
			keycopy( newleft.data[NL], data[R], dimensions );
			NL++;
		}
		else
		{
			keycopy( newright.data[NR], data[R], dimensions );
			NR++;
		}
	}

	for (; L <= page_hdr->size; L++, NL++)
		keycopy( newleft.data[NL], data[L], dimensions );

	if (NL + NR - 2 != page_hdr->size + right.page_hdr->size)
		errorexit("ERROR 1 in shift_from_right()\n");

	newleft.page_hdr->lpage = page_hdr->lpage;
	newleft.page_hdr->size = NL - 1;
	keycopy( newleft.index, index, dimensions );

	newright.page_hdr->lpage = right.page_hdr->lpage;
	newright.page_hdr->size = NR - 1;
	keycopy( newright.index, median );

	return 0;
}

/*============================================================================*/
/***                   median_of_3	  				      */
/*============================================================================*/
/* two or more parameters may point to addresses containing the same values */
static Hcode median_of_3( Hcode& first, Hcode& second, Hcode& third,
							int dimensions )
{
	int		j, k;
	Hcode	min = second, max = first, median(dimensions);

	/* find max and min of first 2 */
	for (k = dimensions - 1; k >= 0; k--)
	{
		if (first.hcode[k] > second.hcode[k])
			break;
		if (first.hcode[k] < second.hcode[k])
		{
			min = first;
			max = second;
			break;
		}
	}

	/* find the median of the first 3 */
	for (j = 0, k = dimensions - 1; k >= 0; k--)
	{
		if (third.hcode[k] > max.hcode[k])
		{
			median = max;
			max = third;
			j = 1;
			break;
		}
		if (third.hcode[k] < max.hcode[k])
			break;
	}


	if (j == 0)	/* third < max */
	{
		median = third;
		for (k = dimensions - 1; k >= 0; k--)
		{
			if (third.hcode[k] > min.hcode[k])
				break;
			if (third.hcode[k] < min.hcode[k])
			{
				median = min;
				min = third;
				break;
			}
		}
	}
	return median;
}

/*============================================================================*/
/***                   median_of_5	  				      */
/*============================================================================*/
/* two or more parameters may point to addresses containing the same values */
Hcode median_of_5( Hcode& first, Hcode& second, Hcode& third,
							Hcode& fourth, Hcode& fifth, int dimensions)
{
	int		j, k, mincount = 0, maxcount = 0;
	Hcode	min = second, max = first, median(dimensions);

	/* find max and min of first 2 */
	for (k = dimensions - 1; k >= 0; k--)
	{
		if (first.hcode[k] > second.hcode[k])
			break;
		if (first.hcode[k] < second.hcode[k])
		{
			min = first;
			max = second;
			break;
		}
	}

	/* find the median of the first 3 */
	for (j = 0, k = dimensions - 1; k >= 0; k--)
	{
		if (third.hcode[k] > max.hcode[k])
		{
			median = max;
			max = third;
			j = 1;
			break;
		}
		if (third.hcode[k] < max.hcode[k])
			break;
	}


	if (j == 0)	/* third < max */
	{
		median = third;
		for (k = dimensions - 1; k >= 0; k--)
		{
			if (third.hcode[k] > min.hcode[k])
				break;
			if (third.hcode[k] < min.hcode[k])
			{
				median = min;
				min = third;
				break;
			}
		}
	}

	/* see if fourth < or > median */
	for (k = dimensions - 1; k >= 0; k--)
	{
		if (fourth.hcode[k] > median.hcode[k])
		{
			maxcount++;
			break;
		}
		if (fourth.hcode[k] < median.hcode[k])
		{
			mincount++;
			break;
		}
	}

	/* see if fifth < or > median */
	for (k = dimensions - 1; k >= 0; k--)
	{
		if (fifth.hcode[k] > median.hcode[k])
		{
			maxcount++;
			break;
		}
		if (fifth.hcode[k] < median.hcode[k])
		{
			mincount++;
			break;
		}
	}

	/* both fourth & fifth are < median */
	if (mincount > maxcount)
	{
		median = min;
		for (k = dimensions - 1; k >= 0; k--)
		{
			if (min.hcode[k] > fourth.hcode[k])
				break;
			if (min.hcode[k] < fourth.hcode[k])
			{
				median = fourth;
				break;
			}
		}

		for (k = dimensions - 1; k >= 0; k--)
		{
			if (median.hcode[k] > fifth.hcode[k])
				break;
			if (median.hcode[k] < fifth.hcode[k])
			{
				median = fifth;
				break;
			}
		}
	}

	/* both fourth & fifth are > median */
	if (mincount < maxcount)
	{
		median = fourth;
		for (k = dimensions - 1; k >= 0; k--)
		{
			if (max.hcode[k] > fourth.hcode[k])
				break;
			if (max.hcode[k] < fourth.hcode[k])
			{
				median = max;
				break;
			}
		}

		for (k = dimensions - 1; k >= 0; k--)
		{
			if (median.hcode[k] > fifth.hcode[k])
			{
				median = fifth;
				break;
			}
			if (median.hcode[k] < fifth.hcode[k])
				break;
		}
	}
	return median;
}

/*============================================================================*/
/***                   median_sortfun	  			      */
/*============================================================================*/
bool median_sortfun( Hcode a, Hcode b )
{
	for ( int i = a.hcode.size() - 1; i >= 0; i-- )
		if ( a.hcode[i] < b.hcode[i] )
			return true;
		else
			if ( a.hcode[i] > b.hcode[i] )
				return false;
	// should never get here - a and b should always be different
	errorexit("ERROR 1 in median_sortfn()\n");
	return false;
}


#ifdef xJKLDEBUGxxxx
void print_hcode(string str, Hcode * hval, int dims)
{
	fjunk3 << "\t\t" << str << "\t";
//	for (int z = dims - 1; z >= 0; z--)
	for (int z = 0; z < dims; z++)
	{
		fjunk3 << setw(15) << hval->hcode[z];
	}
	fjunk3 << "\n";
}
#endif


/*============================================================================*/
/***                   PAGE::p_find_median  (of 5)			    ***/
/*============================================================================*/
// the result is placed in 'median'
Hcode PAGE::p_find_median()
{
	int		start = 0, end = page_hdr->size - 1, i, j, k, n, m, nobj, count;
	HU_int	*key = new HU_int[dimensions];
 	Hcode	*first, *second, *third, *fourth, *fifth, *min, *min2, *min3,
			*max, *max2, *max3, *median;
#ifdef xJKLDEBUGxxxx
	fjunk3 << "Outputting page no. "<< page_hdr->lpage <<"'s record's hcodes followed by data values\n";
#endif
	/* place page in oflowslot's data's hilbert codes in MEDdata array */
	for (i = 1; i <= page_hdr->size; i++)
	{
		// copy to an Hcode, an encoding of a PU_int*
		// ENCODE returns 'key' (a HU_int*)
		keycopy( M->MEDdata[i - 1],  ENCODE( key, data[i], dimensions ) );
#ifdef xJKLDEBUGxxxx
		for (int j = 0; j < dimensions; j++)
			fjunk3 << setw(15) << key[j];
		fjunk3 << " " << i << " ";
		PU_int *D = data[i];
		for (int j = 0; j < dimensions; j++)
			fjunk3 << setw(15) << D[j];
		fjunk3 << endl;

#endif
	}
	delete [] key;

#ifdef xJKLDEBUGxxxx
	fjunk3 << "\nOutputting MEDdata's hcodes\n";
	for ( int j = 0; j < page_hdr->size; j++ )
	{
		for (int z = 0; z < dimensions; z++)
			fjunk3 << setw(15) << M->MEDdata[j].hcode[z];
		fjunk3 << " " << j << endl;
	}
#endif

	/* the first iteration of the outer for loop should NOT ALTER THE POSITIONS
	of the hcodes in MEDdata (ie sort them) because the positions of the
	hcodes originally placed in MEDdata correspond to positions of records
	(stored as attribute values) on the origial data page and, later on, we want
	to compare the 1st (etc) MEDdata hcode with the median and do something
	with the 1st (etc) data page's record according to the results of the
	comparison (eg move the record to the left or right page of a pair of pages):
	that is why we check the value of end-start before entering the loop in the
	first place */

	if (end - start <= THRESHOLD)
	{
		cout << "start = " << start << " end = " << end << " THRESHOLD = " << THRESHOLD << endl;
		errorexit("ERROR in p_find_median(): data pages are too small!\n");
	}

	for (n = 0; end != start; n++)
	{
		// I found the sort function couldn't cope with values of about 30 for THRESHOLD
		if (end - start <= 10)//THRESHOLD)
		{
			nobj = end - start + 1;

			sort( &(M->MEDdata[start]), &(M->MEDdata[start + nobj]), median_sortfun );
//			sort( (M->MEDdata[start]), (M->MEDdata[start + nobj]), median_sortfun );

			return M->MEDdata[start + nobj / 2];
		}

		for (i = start, count = 0; i <= end; i += 5)
		// find median of 3, 4 or 5 - otherwise ignore the values
		{
			if ((m = end + 1 - i) > 2)
			{
				// first, find the median of the first three
				// first = &MEDdata[i]; second = first + 1; third = first + 2;
				if (n == 0)
				{
//  NEVER EXECUTED? - yes! first time (only) the outer for loop iterates - hcodes are
// taken 'randomly' from MEDdata according to sequnce of numbers held in MEDmap

 					first = &(M->MEDdata[M->MEDmap[i]]);
					second = &(M->MEDdata[M->MEDmap[i + 1]]);
					third = &(M->MEDdata[M->MEDmap[i + 2]]);
				}
				else
				{
// CHECK VALUES OF SECOND, THIRD
					first = &(M->MEDdata[i]);
					second = first + 1;
					third = first + 2;
				}
#ifdef xJKLDEBUGxxxx
fjunk3 << "step 1 - " << n <<endl;
#endif

				for (k = dimensions - 1; k >= 0; k--)
				{
					if (first->hcode[k] > second->hcode[k])
					{
						max = first;
						min = second;
						break;
					}
					if (first->hcode[k] < second->hcode[k])
					{
						min = first;
						max = second;
						break;
					}
				}
#ifdef xJKLDEBUGxxxx
fjunk3 << "step 2" <<endl;
print_hcode( "max   ",max  , dimensions );
print_hcode(  "min   ",min , dimensions );
#endif

				
				for (j = 0, k = dimensions - 1; k >= 0; k--)
				{
					if (third->hcode[k] > max->hcode[k])
					{
						median = max;
						max = third;
						j = 1;
						break;
					}
					if (third->hcode[k] < max->hcode[k])
						break;
				}
#ifdef xJKLDEBUGxxxx
fjunk3 << "step 3" <<endl;
if ( j==1)
{	
	print_hcode( "median",median  , dimensions );
	print_hcode(  "max   ",max , dimensions );
}
#endif

				if (j == 0)
				{
					for (k = dimensions - 1; k >= 0; k--)
					{
						if (third->hcode[k] > min->hcode[k])
						{
							median = third;
							break;
						}
						if (third->hcode[k] < min->hcode[k])
						{
							median = min;
							min = third;
							break;
						}
					}
#ifdef xJKLDEBUGxxxx
fjunk3 << "step 4" <<endl;
print_hcode( "median",median  , dimensions );
print_hcode(  "min   ",min , dimensions );
#endif
				}

				if (m > 4) /* we are finding the median of 5 */
				{
/*					Hcode *fourth = first + 3, *fifth = first + 4,*/
					int mincount = 0, maxcount = 0;

					if (n == 0)
					{
						fourth = &(M->MEDdata[M->MEDmap[i + 3]]);
						fifth = &(M->MEDdata[M->MEDmap[i + 4]]);
					}
					else
					{
						fourth = first + 3;
						fifth = first + 4;
					}
#ifdef xJKLDEBUGxxxx
fjunk3 << "step 5 - " << n <<endl;
#endif
					for (k = dimensions - 1; k >= 0; k--)
					{
						if (fourth->hcode[k] > median->hcode[k])
						{
							max2 = fourth;
							maxcount++;
#ifdef xJKLDEBUGxxxx
fjunk3 << "step 6a" <<endl;
print_hcode( "max2  ",max2  , dimensions );
#endif
							break;
						}
						if (fourth->hcode[k] < median->hcode[k])
						{
							min2 = fourth;
							mincount++;
#ifdef xJKLDEBUGxxxx
fjunk3 << "step 6b" <<endl;
print_hcode(  "min2  ",min2 , dimensions );
#endif
							break;
						}
					}

					for (k = dimensions - 1; k >= 0; k--)
					{
						if (fifth->hcode[k] > median->hcode[k])
						{
							max3 = fifth;
							maxcount++;
#ifdef xJKLDEBUGxxxx
fjunk3 << "step 7a" <<endl;
print_hcode(  "max3  ",max3 , dimensions );
#endif
							break;
						}
						if (fifth->hcode[k] < median->hcode[k])
						{
							min3 = fifth;
							mincount++;
#ifdef xJKLDEBUGxxxx
fjunk3 << "step 7b" <<endl;
print_hcode( "min3  ",min3  , dimensions );
#endif
							break;
						}
					}

					if (mincount > maxcount)
					{
						for (k = dimensions - 1; k >= 0; k--)
						{
							if (min->hcode[k] > min2->hcode[k])
							{
								median = min;
								break;
							}
							if (min->hcode[k] < min2->hcode[k])
							{
								median = min2;
								break;
							}
						}

						for (k = dimensions - 1; k >= 0; k--)
						{
							if (median->hcode[k] > min3->hcode[k])
								break;
							if (median->hcode[k] < min3->hcode[k])
							{
								median = min3;
								break;
							}
						}
#ifdef xJKLDEBUGxxxx
fjunk3 << "step 8" <<endl;
print_hcode( "median",median  , dimensions );
#endif
					}

					if (mincount < maxcount)
					{
						for (k = dimensions - 1; k >= 0; k--)
						{
							if (max->hcode[k] > max2->hcode[k])
							{
								median = max2;
								break;
							}
							if (max->hcode[k] < max2->hcode[k])
							{
								median = max;
								break;
							}
						}

						for (k = dimensions - 1; k >= 0; k--)
						{
							if (median->hcode[k] > max3->hcode[k])
							{
								median = max3;
								break;
							}
							if (median->hcode[k] < max3->hcode[k])
								break;
						}
#ifdef xJKLDEBUGxxxx
fjunk3 << "step 9" <<endl;
print_hcode( "median",median  , dimensions );
#endif
					}
				}	/* end if (m > 4) */
				count++;
				M->MEDdata[end + count] = *median;

				#ifdef xJKLDEBUGxxxx
					fjunk3 << "\n....\n";
 					for (int z = 0; z < dimensions; z++)
//					for (int z = dimensions - 1; z >= 0; z--)
						fjunk3 << setw(15) << M->MEDdata[end + count].hcode[z];
					fjunk3 << " " << end + count << endl;

				#endif

			}	/* end if ((m = end + 1 - i) > 2) */
		}	/* end inner for */
		start = end + 1;
		end += count;
	}	/* end outer for */
}

/*============================================================================*/
/***                   PAGE::p_find_median_left	  		    ***/
/*============================================================================*/
// called by the left-hand page
Hcode PAGE::p_find_median_left( PAGE& right )
{
	int i, n, m, temp, start, end, count, rstart, nobj, modifytemp, size = page_hdr->size;
	HU_int	*key = new HU_int[dimensions];
	Hcode dummy( dimensions );

	for (i = 0; i < dimensions; i++)
		dummy.hcode[i] = UINT_MAX;

	/* place left page's data's hilbert codes in MEDdata array */
	for (i = 1; i <= size; i++)
		// copy to an Hcode, an encoding of a PU_int*
		// ENCODE returns 'key' (a HU_int*)
		keycopy( M->MEDdata[i - 1],  ENCODE( key, data[i], dimensions ) );
	delete [] key;

	/* fill the rest of left page with dummy data (max hcode) */
	for (; i <= MAX_DATA; i++)    // CHECK ????
		M->MEDdata[i - 1] = dummy;

	/* add extra dummy values to make a multiple of MEDIAN values */
	temp = MAX_DATA % MEDIAN;    // CHECK ????
	for (; temp != 0 && temp < MEDIAN; temp++, i++)
		M->MEDdata[i - 1] = dummy;

	rstart = i - 1;
	start  = 0;
	end    = size + right.page_hdr->size - 1;

/* see lengthy comment at beginning of the MEDIAN of 5 version of
	p_find_median ()*/

	if (end - start <= THRESHOLD)
		errorexit("ERROR in p_find_median_left(): data pages are too small!\n");

	for (n = 0; end != start; n++)
	{
		if (end - start <= THRESHOLD)
		{
			nobj = end - start + 1;

			sort( &(M->MEDdata[start]), &(M->MEDdata[start + nobj]), median_sortfun );

			return M->MEDdata[start + nobj / 2];
		}

		for (i = start, count = 0, modifytemp = true; i <= end; i += MEDIAN)
		/* find median of 3, 4 or 5 - otherwise ignore the values */
		{
			if ((m = end + 1 - i) > 2) /* at least 3 hcodes left */
			{
				count++;
				if (i >= rstart)
				{
					if (modifytemp)
					{
						temp = end + count;
						modifytemp = false;
					}
					M->MEDdata[end + count] = dummy;
				}
				else
				{
					if (n == 0)
						if (m > 4)
							M->MEDdata[end + count] = median_of_5(
									M->MEDdata[M->MEDmap[i]],
									M->MEDdata[M->MEDmap[i + 1]],
									M->MEDdata[M->MEDmap[i + 2]],
									M->MEDdata[M->MEDmap[i + 3]],
									M->MEDdata[M->MEDmap[i + 4]],
									dimensions);
						else
							M->MEDdata[end + count] = median_of_3(
									M->MEDdata[M->MEDmap[i]],
									M->MEDdata[M->MEDmap[i + 1]],
									M->MEDdata[M->MEDmap[i + 2]],
									dimensions);
					else
						if (m > 4)
							M->MEDdata[end + count] = median_of_5(
									M->MEDdata[i],
									M->MEDdata[i + 1],
									M->MEDdata[i + 2],
									M->MEDdata[i + 3],
									M->MEDdata[i + 4],
									dimensions);
						else
							M->MEDdata[end + count] = median_of_3(
									M->MEDdata[i],
									M->MEDdata[i + 1],
									M->MEDdata[i + 2],
									dimensions);
				}
			}
		}
		start = end + 1;
		end += count;
		rstart = temp;
	}
}

/*============================================================================*/
/***                   PAGE::p_find_median_right	  		    ***/
/*============================================================================*/
// called by the left-hand page
Hcode PAGE::p_find_median_right( PAGE& right )
{
	int i, n, m, temp, start, end, count, rstart, nobj, size = page_hdr->size;
	HU_int	*key = new HU_int[dimensions];
	Hcode dummy( dimensions );

	for (i = 0; i < dimensions; i++)
		dummy.hcode[i] = 0;

	/* place right page's data's hilbert codes in MEDdata array */
	for (i = 1; i <= right.page_hdr->size; i++)
		keycopy( M->MEDdata[size + i - 1],  ENCODE( key, data[i], dimensions ) );
	delete [] key;

	temp = size % MEDIAN;
	rstart = size - temp;
	/* place left page's key in MEDdata beyond the last multiple of MEDIAN
	   entries in left page */
	for (; temp > 0 ; temp--)
		M->MEDdata[size - temp] = dummy;

	/* the first MEDdata element which notionally holds the right page's
	   hilbert codes, excluding those used to pad out the end of left page */
	start = 0;
	end = size + right.page_hdr->size - 1;

/* see lengthy comment at beginning of the MEDIAN of 5 version of
	p_find_median ()*/

	if (end - start <= THRESHOLD)
		errorexit("ERROR in p_find_median_right(): data pages are too small!\n");

	/* find the median of the two pages */
	for (n = 0; start != end; n++)
	{
		if (end - start <= THRESHOLD)
		{
			nobj = end - start + 1;

			sort( &(M->MEDdata[start]), &(M->MEDdata[start + nobj]), median_sortfun );

			return M->MEDdata[start + nobj / 2];
		}

		for (i = start, count = 0, temp = 0; i <= end; i += MEDIAN)
		{
			count++;
			if (i < rstart)
			{
				temp = end + count;
				M->MEDdata[end + count] = dummy;
				continue;
			}

			if (n == 0 && i == rstart)
			{
				M->MEDdata[end + count] = median_of_5(
					M->MEDdata[i     < size ? i    : M->MEDmap[i - size]],
					M->MEDdata[i + 1 < size ? i + 1: M->MEDmap[i - size +1]],
					M->MEDdata[i + 2 < size ? i + 2: M->MEDmap[i - size +2]],
					M->MEDdata[i + 3 < size ? i + 3: M->MEDmap[i - size +3]],
					M->MEDdata[i + 4 < size ? i + 4: M->MEDmap[i - size +4]],
					dimensions);
				continue;
			}

			if ((m = end + 1 - i) > 2)
			{
				if (n == 0)
					if (m > 4)
						M->MEDdata[end + count] = median_of_5(
								M->MEDdata[M->MEDmap[i - size]],
								M->MEDdata[M->MEDmap[i - size + 1]],
								M->MEDdata[M->MEDmap[i - size + 2]],
								M->MEDdata[M->MEDmap[i - size + 3]],
								M->MEDdata[M->MEDmap[i - size + 4]],
								dimensions);
					else
						M->MEDdata[end + count] = median_of_3(
								M->MEDdata[M->MEDmap[i - size]],
								M->MEDdata[M->MEDmap[i - size + 1]],
								M->MEDdata[M->MEDmap[i - size + 2]],
								dimensions);
				else
					if (m > 4)
						M->MEDdata[end + count] = median_of_5(
								M->MEDdata[i],
								M->MEDdata[i + 1],
								M->MEDdata[i + 2],
								M->MEDdata[i + 3],
								M->MEDdata[i + 4],
								dimensions);
					else
						M->MEDdata[end + count] = median_of_3(
								M->MEDdata[i],
								M->MEDdata[i + 1],
								M->MEDdata[i + 2],
								dimensions);
			}

		}
		start = end + 1;
		end += count;
		rstart = temp;
		if (rstart - start < MEDIAN)
			rstart = 0;
	}
}
