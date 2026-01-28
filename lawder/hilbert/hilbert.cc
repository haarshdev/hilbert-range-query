// Copyright (C) Jonathan Lawder 2001-2011

/*----------------------------------------------------------------------------*/
/*	THIS FILE CONTAINS CODE FOR THE HILBERT CURVE WHICH		      */
/*                            						      */
/*	USES BUTZ' CALCULATED MAPPING METHOD				      */
/*                            						      */
/*	ie MAPPING == BUTZ		 				      */
/*----------------------------------------------------------------------------*/

//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>

#ifdef __MSDOS__
	#include "..\gendefs.h"
	#include "..\utils\utils.h"
#else
	#include "../gendefs.h"
	#include "../utils/utils.h"
#endif

using namespace std;


/*
#if DIM == 32
	const U_int g_all_ones = ~0;
#elif DIM < 32
	const U_int g_all_ones = (1 << DIM) -1;
#endif
*/

/*============================================================================*/
/*                            						      */
/*                            ENCODING & DECODING FUNCTIONS		      */
/*                            						      */
/*============================================================================*/

/*============================================================================*/
/*                            ENCODE					      */
/*============================================================================*/
// point is the point to be encoded, hcode receives the hilbert code of point.
// hcode is returned.
// point and hcode must have storage allocated to them before this function
// is called.
HU_int* ENCODE( HU_int* hcode, const PU_int* const point, int DIMS )
{
	U_int	mask = (U_int)1 << WORDBITS - 1, element, temp1, temp2,
		A, W = 0, S, tS, T, tT, J, P = 0, xJ;
	int	i = NUMBITS * DIMS - DIMS, j;
	
	// initialise hcode
	memset( hcode, 0, sizeof(U_int) * DIMS );

	for (j = A = 0; j < DIMS; j++)
		if (point[j] & mask)
			A |= (1 << (DIMS-1-j));

	S = tS = A;

	P |= S & (1 << (DIMS-1));
	for (j = 1; j < DIMS; j++)
		if( S & 1 << (DIMS-1-j) ^ (P >> 1) & (1 << (DIMS-1-j)))
			P |= (1 << (DIMS-1-j));

	/* add in DIMS bits to hcode */
	element = i / WORDBITS;
	if (i % WORDBITS > WORDBITS - DIMS)
	{
		hcode[element] |= P << i % WORDBITS;
		hcode[element + 1] |= P >> WORDBITS - i % WORDBITS;
	}
	else
		hcode[element] |= P << i - element * WORDBITS;

	J = DIMS;
	for (j = 1; j < DIMS; j++)
		if ((P >> j & 1) == (P & 1))
			continue;
		else
			break;
	if (j != DIMS)
		J -= j;
	xJ = J - 1;

	if (P < 3)
		T = 0;
	else
		if (P % 2)
			T = (P - 1) ^ (P - 1) / 2;
		else
			T = (P - 2) ^ (P - 2) / 2;
	tT = T;

/*
printf ("           A   W  tS   S   P   J  xJ   T  tT\n\n");
printf ("        %4u%4u%4u%4u%4u%4u%4u%4u%4u\n\n", A, W, tS, S, P, J, xJ, T, tT);
*/

	for (i -= DIMS, mask >>= 1; i >=0; i -= DIMS, mask >>= 1)
	{
		for (j = A = 0; j < DIMS; j++)
			if (point[j] & mask)
				A |= (1 << (DIMS-1-j));

		W ^= tT;
		tS = A ^ W;
		if (xJ % DIMS != 0)
		{
			temp1 = tS << xJ % DIMS;
			temp2 = tS >> DIMS - xJ % DIMS;
			S = temp1 | temp2;
			S &= ((U_int)1 << DIMS) - 1;
		}
		else
			S = tS;

		P = S & (1 << (DIMS-1));
		for (j = 1; j < DIMS; j++)
			if( S & (1 << (DIMS-1-j)) ^ (P >> 1) & (1 << (DIMS-1-j)))
				P |= (1 << (DIMS-1-j));

		/* add in DIMS bits to hcode */
		element = i / WORDBITS;
		if (i % WORDBITS > WORDBITS - DIMS)
		{
			hcode[element] |= P << i % WORDBITS;
			hcode[element + 1] |= P >> WORDBITS - i % WORDBITS;
		}
		else
			hcode[element] |= P << i - element * WORDBITS;

		if (i > 0)
		{
			if (P < 3)
				T = 0;
			else
				if (P % 2)
					T = (P - 1) ^ (P - 1) / 2;
				else
					T = (P - 2) ^ (P - 2) / 2;

			if (xJ % DIMS != 0)
			{
				temp1 = T >> xJ % DIMS;
				temp2 = T << DIMS - xJ % DIMS;
				tT = temp1 | temp2;
				tT &= ((U_int)1 << DIMS) - 1;
			}
			else
				tT = T;

			J = DIMS;
			for (j = 1; j < DIMS; j++)
				if ((P >> j & 1) == (P & 1))
					continue;
				else
					break;
			if (j != DIMS)
				J -= j;

			xJ += J - 1;
		/*	J %= DIMS;*/
		}
/*
printf ("        %4u%4u%4u%4u%4u%4u%4u%4u%4u\n\n", A, W, tS, S, P, J, xJ, T, tT);
*/
	}
/*
for (i=0; i < DIMS; i++)
	printf("     %-13u", point[i]); printf("\n");
for (i=0; i < DIMS; i++)
/-*	printf("BUTZ %-13u", hcode[i]); printf("\n");*-/
	printf("XXX  %-13u", hcode[i]); printf("\n");
*/
	return hcode;
}

/*============================================================================*/
/*                            DECODE					      */
/*============================================================================*/
PU_int* DECODE ( PU_int* point, HU_int* hcode, int DIMS )
{
	U_int	mask = (U_int)1 << WORDBITS - 1, element, temp1, temp2,
		A, W = 0, S, tS, T, tT, J, P = 0, xJ;
	int	i = NUMBITS * DIMS - DIMS, j;


	/*--- P ---*/
	element = i / WORDBITS;
	P = hcode[element];
	if (i % WORDBITS > WORDBITS - DIMS)
	{
		temp1 = hcode[element + 1];
		P >>= i % WORDBITS;
		temp1 <<= WORDBITS - i % WORDBITS;
		P |= temp1;
	}
	else
		P >>= i % WORDBITS;	/* P is a DIMS bit hcode */

	/* the & masks out spurious highbit values */
	#if DIMS < WORDBITS
		P &= (1 << DIMS) -1;
	#endif

	/*--- xJ ---*/
	J = DIMS;
	for (j = 1; j < DIMS; j++)
		if ((P >> j & 1) == (P & 1))
			continue;
		else
			break;
	if (j != DIMS)
		J -= j;
	xJ = J - 1;

	/*--- S, tS, A ---*/
	A = S = tS = P ^ P / 2;


	/*--- T ---*/
	if (P < 3)
		T = 0;
	else
		if (P % 2)
			T = (P - 1) ^ (P - 1) / 2;
		else
			T = (P - 2) ^ (P - 2) / 2;

	/*--- tT ---*/
	tT = T;

	/*--- distrib bits to coords ---*/
	for (j = DIMS - 1; P > 0; P >>=1, j--)
		if (P & 1)
			point[j] |= mask;


	for (i -= DIMS, mask >>= 1; i >=0; i -= DIMS, mask >>= 1)
	{
		/*--- P ---*/
		element = i / WORDBITS;
		P = hcode[element];
		if (i % WORDBITS > WORDBITS - DIMS)
		{
			temp1 = hcode[element + 1];
			P >>= i % WORDBITS;
			temp1 <<= WORDBITS - i % WORDBITS;
			P |= temp1;
		}
		else
			P >>= i % WORDBITS;	/* P is a DIMS bit hcode */

		/* the & masks out spurious highbit values */
		#if DIMS < WORDBITS
			P &= (1 << DIMS) -1;
		#endif

		/*--- S ---*/
		S = P ^ P / 2;

		/*--- tS ---*/
		if (xJ % DIMS != 0)
		{
			temp1 = S >> xJ % DIMS;
			temp2 = S << DIMS - xJ % DIMS;
			tS = temp1 | temp2;
			tS &= ((U_int)1 << DIMS) - 1;
		}
		else
			tS = S;

		/*--- W ---*/
		W ^= tT;

		/*--- A ---*/
		A = W ^ tS;

		/*--- distrib bits to coords ---*/
		for (j = DIMS - 1; A > 0; A >>=1, j--)
			if (A & 1)
				point[j] |= mask;

		if (i > 0)
		{
			/*--- T ---*/
			if (P < 3)
				T = 0;
			else
				if (P % 2)
					T = (P - 1) ^ (P - 1) / 2;
				else
					T = (P - 2) ^ (P - 2) / 2;

			/*--- tT ---*/
			if (xJ % DIMS != 0)
			{
				temp1 = T >> xJ % DIMS;
				temp2 = T << DIMS - xJ % DIMS;
				tT = temp1 | temp2;
				tT &= ((U_int)1 << DIMS) - 1;
			}
			else
				tT = T;

			/*--- xJ ---*/
			J = DIMS;
			for (j = 1; j < DIMS; j++)
				if ((P >> j & 1) == (P & 1))
					continue;
				else
					break;
			if (j != DIMS)
				J -= j;
			xJ += J - 1;
		}
	}
	return point;
}

/*============================================================================*/
/*                            						      */
/*                            PARTIAL MATCH QUERY FUNCTIONS		      */
/*                            						      */
/*============================================================================*/

/*============================================================================*/
/*                            HB_matches_PMQ				      */
/*============================================================================*/
/* based on H_decode()
   returns 1 if decode(P) is a match to query, otherwise returns 0
   the last 3 pointer parameters get modified ONLY WHEN THERE IS A MATCH
   very similar to HB_matches_RQ
 */
static bool HB_matches_PMQ( U_int P, U_int Qsaf, U_int Qp,
				U_int *xJ, U_int *tT, U_int *W, int DIMS )
{
	U_int	J, S, tS, T, local_W, A, temp1, temp2;
	int		j;

	local_W = *W;

	/*--- S ---*/
	S = P ^ P / 2;

	/*--- tS ---*/
	if (*xJ % DIMS != 0)
	{
		temp1 = S >> *xJ % DIMS;
		temp2 = S << DIMS - *xJ % DIMS;
		tS = temp1 | temp2;
		tS &= ((U_int)1 << DIMS) - 1;
	}
	else
		tS = S;

	/*--- W ---*/
	local_W ^= *tT;

	/*--- A ---*/
	A = local_W ^ tS;

	if ((A & Qsaf) == Qp)
	{
		/*--- W ---*/
		*W = local_W;

		/*--- T ---*/
		if (P < 3)
			T = 0;
		else
			if (P % 2)
				T = (P - 1) ^ (P - 1) / 2;
			else
				T = (P - 2) ^ (P - 2) / 2;

		/*--- tT ---*/
		if (*xJ % DIMS != 0)
		{
			temp1 = T >> *xJ % DIMS;
			temp2 = T << DIMS - *xJ % DIMS;
			*tT = temp1 | temp2;
			*tT &= ((U_int)1 << DIMS) - 1;
		}
		else
			*tT = T;

		/*--- xJ ---*/
		J = DIMS;
		for (j = 1; j < DIMS; j++)
			if ((P >> j & 1) == (P & 1))
				continue;
			else
				break;
		if (j != DIMS)
			J -= j;
		*xJ += J - 1;

		return true;
	}

	return false;
}


/*----------------------------------------------------------------------------*/
// this is HPMQ_VER == 2
/*----------------------------------------------------------------------------*/
/*============================================================================*/
/*                            HB_nextmatch_PM				      */
/*============================================================================*/
static bool HB_nextmatch_PM ( PU_int *query, HU_int *match, HU_int *key,
			U_int Qsaf, int counter, U_int xJ, U_int tT, U_int W, int DIMS )
{
	U_int	mask = 1 << counter / DIMS, Qp, H, element, temp,
			local_xJ, local_tT, local_W;
	U_int	Lo, Hi, LoHi_H, HiLo_H, tSL, tSH, LoHi_C, temp1, temp2, diffbit,
			LoBAK, HiBAK, NUMPOINTS = (U_int)1 << DIMS;
	bool	Backup = false;
	int		i, j;

/* static int jkl = 0; */

	if (counter >= 0)
	{
		for (i = 0, Qp = 0; i < DIMS; i++)
			if (query[i] & mask)
				Qp |= (1 << (DIMS-1-i));

		/* extract DIMS bits from key */
		element = counter / WORDBITS;
		H = key[element];
		if (counter % WORDBITS > WORDBITS - DIMS)
		{
			temp = key[element + 1];
			H >>= counter % WORDBITS;
			temp <<= WORDBITS - counter % WORDBITS;
			H |= temp;
		}
		else
			H >>= counter % WORDBITS;
		/* the & masks out spurious highbit values */
		if (DIMS < WORDBITS)
			H &= (1 << DIMS) -1;

		local_xJ = xJ;
		local_tT = tT;
		local_W  = W;

		for (Lo = 0, Hi = NUMPOINTS-1 ; ; )
		{
			for (; Lo != Hi ;)
			{
				LoHi_H = Lo + (Hi - Lo) / 2;
				HiLo_H = LoHi_H + 1;

				tSL = LoHi_H ^ LoHi_H / 2;
				tSH = HiLo_H ^ HiLo_H / 2;

				if (local_xJ % DIMS != 0)
				{
					temp1 = tSL >> local_xJ % DIMS;
					temp2 = tSL << DIMS - local_xJ % DIMS;
					tSL   = temp1 | temp2;
					tSL  &= ((U_int)1 << DIMS) - 1;
					temp1 = tSH >> local_xJ % DIMS;
					temp2 = tSH << DIMS - local_xJ % DIMS;
					tSH   = temp1 | temp2;
					tSH  &= ((U_int)1 << DIMS) - 1;
				}
				diffbit = tSL ^ tSH;

				LoHi_C = local_W ^ local_tT ^ tSL;
/*				HiC = local_W ^ local_tT ^ tSH; */

				if (!(diffbit & Qsaf) || ((diffbit & LoHi_C) == (diffbit & Qp)))
				{
					if (LoHi_H >= H)
					{
						if (!(diffbit & Qsaf))
						{
							LoBAK = HiLo_H;
							HiBAK = Hi;
							Backup = true;
						}
						Hi = LoHi_H;
					}
					else
						if (!(diffbit & Qsaf))
							Lo = HiLo_H;
						else
							if (true == Backup)
							{
								Lo = LoBAK;
								Hi = HiBAK;
								Backup = false;
							}
							else
								return false;
				}
				else
					Lo = HiLo_H;
/* jkl++; */
			} /* inner 'for' loop */

			if (/*Lo < H || */false == HB_matches_PMQ( Lo, Qsaf, Qp, &local_xJ, &local_tT, &local_W, DIMS ))
			{
				if (false == Backup)
					return false;
				Lo = LoBAK;
				Hi = HiBAK;
				Backup = false;
				continue;
			}
/* check value of Lo here */
			if (Lo == H)
				if (false == HB_nextmatch_PM( query, match, key, Qsaf,
					counter - DIMS, local_xJ, local_tT, local_W, DIMS ))
				{
					if (false == Backup)
						return false;
					local_xJ = xJ;
					local_tT = tT;
					local_W  = W;
					Lo = LoBAK;
					Hi = HiBAK;
					Backup = false;
					continue;
				}

			/* add in a bit of the hcode */
			element = counter / WORDBITS;
			if (counter % WORDBITS > WORDBITS - DIMS)
			{
				match[element] |= Lo << counter % WORDBITS;
				match[element + 1] |=
						Lo >> WORDBITS - counter % WORDBITS;
			}
			else
				match[element] |=
						Lo << counter - element * WORDBITS;

			if (Lo == H)
				return true;
			break; /* ..and move on to the fast-forward code */
		}
	}
	else
		return true;

	/* now we can fast-forward:  */
	for (j = counter - DIMS; j >= 0; j -=DIMS)
	{
		mask >>= 1;
		for (i = Qp = 0; i < DIMS; i++)
			if (query[i] & mask)
				Qp |= (1 << (DIMS-1-i));

		for (Lo = 0, Hi = NUMPOINTS-1; Lo != Hi ;)
		{
			LoHi_H = Lo + (Hi - Lo) / 2;
			HiLo_H = LoHi_H + 1;

			tSL = LoHi_H ^ LoHi_H / 2;
			tSH = HiLo_H ^ HiLo_H / 2;

			if (local_xJ % DIMS != 0)
			{
				temp1 = tSL >> local_xJ % DIMS;
				temp2 = tSL << DIMS - local_xJ % DIMS;
				tSL   = temp1 | temp2;
				tSL  &= ((U_int)1 << DIMS) - 1;
				temp1 = tSH >> local_xJ % DIMS;
				temp2 = tSH << DIMS - local_xJ % DIMS;
				tSH   = temp1 | temp2;
				tSH  &= ((U_int)1 << DIMS) - 1;
			}
			diffbit = tSL ^ tSH;

			LoHi_C = local_W ^ local_tT ^ tSL;
			if (!(diffbit & Qsaf) || ((diffbit & LoHi_C) == (diffbit & Qp)))
				Hi = LoHi_H;
			else
				Lo = HiLo_H;
		}

		if (false == HB_matches_PMQ( Lo, Qsaf, Qp, &local_xJ, &local_tT, &local_W, DIMS ))
			errorexit("Error in HB_nextmatch_PM()\n");

		element = j / WORDBITS;
		if (j % WORDBITS > WORDBITS - DIMS)
		{
			match[element] |= Lo << j % WORDBITS;
			match[element+1] |= Lo >> WORDBITS - j % WORDBITS;
		}
		else
			match[element] |= Lo << j - element * WORDBITS;
	}
	return true;
}

/*============================================================================*/
/*                            H_nextmatch_PM				      */
/*============================================================================*/
/* in call to HB_nextmatch_PM(); the last 3 parameters
   are xJ, tT and W */
bool H_nextmatch_PM( PU_int *query, HU_int *match, HU_int *key, U_int Qsaf, int dimensions )
{
    return HB_nextmatch_PM( query, match, key, Qsaf, dimensions * (NUMBITS-1), 0, 0, 0, dimensions );
}


/*============================================================================*/
/*                            						      */
/*                            RANGE QUERY FUNCTIONS			      */
/*                            						      */
/*============================================================================*/


static bool HB_matches_RQ (U_int P, U_int Lo, U_int L_xor_H,
                      U_int *A, U_int *xJ, U_int *tT, U_int *W, int DIMS)
{
    U_int    J, S, tS, T, local_W, temp1, temp2;
    int      j;

    local_W = *W;

    /*--- calculate S ---*/
    S = P ^ P / 2;

    /*--- calculate tS ---*/
    if (*xJ % DIMS != 0)
    {
        temp1 = S >> *xJ % DIMS;
        temp2 = S << DIMS - *xJ % DIMS;
        tS = temp1 | temp2;
        tS &= ((U_int)1 << DIMS) - 1;
    }
    else
        tS = S;

    /*--- calculate W ---*/
    local_W ^= *tT;

    /*--- calculate A ---*/
    *A = local_W ^ tS;

    if ((*A & L_xor_H) == (Lo & L_xor_H))
    /* the n-point whose derived-key is P lies within the
       query range whose derived-key is Lo */
    {
        /*--- calculate W ---*/
        *W = local_W;

        /*--- calculate T ---*/
        if (P < 3)
            T = 0;
        else
            if (P % 2)
                T = (P - 1) ^ (P - 1) / 2;
            else
                T = (P - 2) ^ (P - 2) / 2;

        /*--- calculate tT ---*/
        if (*xJ % DIMS != 0)
        {
            temp1 = T >> *xJ % DIMS;
            temp2 = T << DIMS - *xJ % DIMS;
            *tT = temp1 | temp2;
            *tT &= ((U_int)1 << DIMS) - 1;
        }
        else
            *tT = T;

        /*--- calculate xJ ---*/
        J = DIMS;
        for (j = 1; j < DIMS; j++)
            if ((P >> j & 1) == (P & 1))
                continue;
            else
                break;
        if (j != DIMS)
            J -= j;
        *xJ += J - 1;
//	printf("HB_matches_RQ returning 'true' .............\n");

        return true;
    }

//	printf("HB_matches_RQ returning 'false'.......\n");

    return false;
}

/*----------------------------------------------------------------------------*/
#if 0 // this is HRQ_VER == 7
/*----------------------------------------------------------------------------*/

/*=============================================================*/
/*                    HB_nextmatch_RQ                          */
/*=============================================================*/
static U_int HB_nextmatch_RQ (Point LoBound, Point HiBound,
            Hcode *next_match, Hcode *page_key, int K,
            U_int xJ, U_int tT, U_int W,
            U_int LoMask, U_int HiMask)
{
    U_int   mask = 1 << K, qLo, qHi,
            L_xor_H, R_xor_L, R_xor_H, region,
            H, element, temp, local_xJ, local_tT, local_W,
            partioning_dimension, lomask, himask,
            Lo, Hi, LoHi_H, HiLo_H, LoHi_C, tSL, tSH,
            temp1, temp2, LoBAK, HiBAK, Backup = FALSE, i;
    int     j, N, retval = 0;
    Point   lobound, hibound;

    if (K >= 0)
    {
/* step 1: find the n-points of the quadrants in which the
   query lower and upper bound points lie */

        for (i = qLo = qHi = 0; i < DIM; i++)
        {
            if (LoBound.hcode[i] & mask)
                qLo |= 1<< (DIM-1-i);
            if (HiBound.hcode[i] & mask)
                qHi |= 1<< (DIM-1-i);
        }

/* step 2: extract DIM consecutive bits from page_key */

        element = (DIM * K) / 32;
        H = page_key->hcode[element];
        if ((DIM * K) % 32 > 32 - DIM)
        {
            temp = page_key->hcode[element + 1];
            H >>= (DIM * K) % 32;
            temp <<= 32 - (DIM * K) % 32;
            H |= temp;
        }
        else
            H >>= (DIM * K) % 32;
        /* the & masks out spurious highbit values */
        if (DIM < 32)
            H &= (1 << DIM) -1;

/* step 3: find in which dimensions the query lower and upper
   bounds have the same values */

        L_xor_H = qLo ^ qHi;
        L_xor_H = ~L_xor_H;

/* step 4: do a binary search of the derived-keys of the
   quadrants in the current-state to find the quadrant with the
   lowest derived-key that is greater than or equal to the value
   of H and which is within the query range */

        local_xJ = xJ;
        local_tT = tT;
        local_W  = W;

        for (Lo = 0, Hi = NUMPOINTS - 1 ; ; )
        {
            /* do a binary search of the state */
            for (; Lo != Hi ;)
            {
/* the middle 2 ordered derived-keys of quadrants in a state or
   part of a state. They enable us to divide it in 2; one is the
   upper bound of the lower part and the other is the lower
   bound of the upper part */

                LoHi_H = Lo + (Hi - Lo) / 2;
                HiLo_H = LoHi_H + 1;

/* find the partitioning-dimension */

                tSL = LoHi_H ^ LoHi_H / 2;
                tSH = HiLo_H ^ HiLo_H / 2;

                if (local_xJ % DIM != 0)
                {
                    temp1 = tSL >> local_xJ % DIM;
                    temp2 = tSL << DIM - local_xJ % DIM;
                    tSL   = temp1 | temp2;
                    tSL  &= ((U_int)1 << DIM) - 1;
                    temp1 = tSH >> local_xJ % DIM;
                    temp2 = tSH << DIM - local_xJ % DIM;
                    tSH   = temp1 | temp2;
                    tSH  &= ((U_int)1 << DIM) - 1;
                }

/* the dimension in which all points in the lower half of the
   state (or part) have the same value and the points in the
   upper half have the opposite value */

                partioning_dimension = tSL ^ tSH;

/* The n-point of the quadrant with the highest derived-key in
   the 'lower' half of the state (or part) */

                LoHi_C = local_W ^ local_tT ^ tSL;

/* check to see if the query intersects with the lower half of
   the state (or part) i.e. if its lower bound is within it (the
   first test checks whether the query intersects both halves -
   both parts of the test are needed. NB 'partioning_dimension'
   and 'LoHi_C' change as the loop iterates) */

                if (!(partioning_dimension & L_xor_H) ||
                    ((partioning_dimension & LoHi_C) ==
                     (partioning_dimension & qLo)))
                {
                    /* there is a match in the lower half of
                       the state; either the range spans both
                       halves (so far as the partioning_dimension
                       dim is concerned) or LoHi_C is within the
                       range... */

/* check to see if H is in the lower half of the state
   (or part) */

                    if (LoHi_H >= H)
                    {

/* H is in the lower half of the state (or part). If the query
   spans both halves of the state (or part), record the limits
   of the top half of the state and note that search may continue
   there by assigning 'TRUE' to 'Backup'. If a next-match is not
   subsequently found in the lower half of the state the search
   may in the upper half (this is back-tracking) */

                        if (!(partioning_dimension & L_xor_H))
                        {
                            LoBAK = HiLo_H;
                            HiBAK = Hi;
                            Backup = TRUE;
                        }

/* set (reduce) the 'Hi' value for the next iteration of the
   binary search */

                        Hi = LoHi_H;
                    }
                    else

/* ...but H is in the upper half of the state (or part),
   therefore check if the range also intersects with the upper
   half...if so, set (increase) the 'Lo' value for the next
   iteration of the binary search */

                        if (!(partioning_dimension & L_xor_H))
                            Lo = HiLo_H;
                        else

/*...otherwise H is in the upper half of the state (or part) but
   the query doesn't intersect with it; therefore, if we
   previously found a (larger) half state (or part) with which
   the query intersects, backtrack to it, changing the values of
   'Lo' & 'Hi' accordingly */

                            if (Backup)
                            {
                                Lo = LoBAK;
                                Hi = HiBAK;
                                Backup = FALSE;
                            }
                            else
                                return 0;
                }
                else

/* the query region does not intersect with the upper half of
   the state (or part), therefore set (increase) the value of
   'Lo' for the next iteration of the binary search */

                    Lo = HiLo_H;
            } /* inner 'for' loop */

/* check to see if the quadrant whose derived-key is 'Lo' lies
   within the query region and, if so, determine the
   characteristics of the next-state (using variables
   local_xJ, local_tT and local_W). 'region' is set to be
   the n-point corresponding to 'Lo'. If this function call
   fails, back-track if possible or else return to an earlier
   call to HB_nextmatch_RQ to resume the search, i.e. backtrack
   there */

            if (!HB_matches_RQ(Lo, qLo, L_xor_H, &region,
                        &local_xJ, &local_tT, &local_W))
            {
                if (!Backup)
                    return 0;
                Lo = LoBAK;
                Hi = HiBAK;
                Backup = FALSE;
                continue;
            }

/* save existing variable values in case search continues
   in another quadrant later */

            lobound = LoBound; hibound = HiBound;
            lomask  = LoMask;  himask  = HiMask;

            if (K > 0)
            {

/* restrict the query region (lobound and hibound) to that part
   which intersects the region containing the next-match, if it
   exists, i.e. 'Lo'. 'R_xor_L' has bits set to 1 if 'region'
   and 'qlo' have different values in the correponding
   dimensions. Since 'region' is within the query range, its
   bits will be 1 and qLo's will be 0 therefore in restricting
   the query range, we increase the value of the current order's
   bit in the query to 1 and set the lower bits to 0. Since at
   higher orders we are only going to encounter 0 bits in qLo
   values, it's simplest to set the whole value of the query
   lower bound to 0. A similar logic applies to 'R_xor_H' if the
   query region was wholly within the subspace which is defined
   by 'region' then both 'R_xor_L' and 'R_xor_H' would be 0 in
   all bits and the query region would not get changed */

                for (R_xor_L = region ^ qLo, N = DIM - 1;
                        R_xor_L;
                        R_xor_L >>= 1, N--)
                    if (R_xor_L & 1)
                    {
                        lobound.hcode[N] = 0;
                        lomask |= 1<< (DIM-1-N);
                    }
                for (R_xor_H = region ^ qHi, N = DIM - 1;
                        R_xor_H;
                        R_xor_H >>= 1, N--)
                    if (R_xor_H & 1)
                    {
                        hibound.hcode[N] = U_int_MAX;
                        himask |= 1<< (DIM-1-N);
                    }
            }

            if (Lo == H)
                if (lomask == g_all_ones && himask == g_all_ones)
                {
                    /* the page_key is a match to the query */
                    *next_match = *page_key;
                    return 2; /* a next-match has been found */
                }
                else
                    if (!(retval = HB_nextmatch_RQ(lobound,
                        hibound, next_match, page_key, K - 1,
                        local_xJ, local_tT, local_W,
                        lomask, himask)))
                    {
                        if (!Backup)
                            return 0;
                        local_xJ = xJ;
                        local_tT = tT;
                        local_W  = W;
                        Lo = LoBAK;
                        Hi = HiBAK;
                        Backup = FALSE;
                        continue;
                    }

            if (retval == 2)
/* the whole of the match has been determined; therefore no need
   for step 6 */
                return 2;

/* step 6: add in DIM bits of the current-quadrant's
   derived-key into the next-match */

            element = (DIM * K) / 32;
            if ((DIM * K) % 32 > 32 - DIM)
            {
                next_match->hcode[element] |=
                            Lo << (DIM * K) % 32;
                next_match->hcode[element + 1] |=
                            Lo >> 32 - (DIM * K) % 32;
            }
            else
                next_match->hcode[element] |=
                            Lo << (DIM * K) - element * 32;

            if (Lo == H)
                return 1; /* next_match tentatively found */
            break; /* we can move into fast-forward mode */
        }
    }
    else return 1; /* next_match tentatively found */

    /* now we can fast-forward:  */
    for (j = K - 1; j >= 0; j --)
    {

/* at his stage, we are no longer interseted in the value of bits
   in the page-key; a next-match is guaranteed to exist - all
   that is needed now is to find the lowest derived-key in the
   current query region. */

        if (lomask == g_all_ones && himask == g_all_ones)
            return 1; /* next_match tentatively found */

        mask >>= 1;

        for (i = qLo = qHi = 0; i < DIM; i++)
        {
            if (lobound.hcode[i] & mask)
                qLo |= 1<< (DIM-1-i);
            if (hibound.hcode[i] & mask)
                qHi |= 1<< (DIM-1-i);
        }
        L_xor_H = qLo ^ qHi;
        L_xor_H = ~L_xor_H;

        /* do a binary search of the state */
        for (Lo = 0, Hi = NUMPOINTS - 1; Lo != Hi ; )
        {
            LoHi_H = Lo + (Hi - Lo) / 2;
            HiLo_H = LoHi_H + 1;

            tSL = LoHi_H ^ LoHi_H / 2;
            tSH = HiLo_H ^ HiLo_H / 2;

            if (local_xJ % DIM != 0)
            {
                temp1 = tSL >> local_xJ % DIM;
                temp2 = tSL << DIM - local_xJ % DIM;
                tSL   = temp1 | temp2;
                tSL  &= ((U_int)1 << DIM) - 1;
                temp1 = tSH >> local_xJ % DIM;
                temp2 = tSH << DIM - local_xJ % DIM;
                tSH   = temp1 | temp2;
                tSH  &= ((U_int)1 << DIM) - 1;
            }
            partioning_dimension = tSL ^ tSH;

            LoHi_C = local_W ^ local_tT ^ tSL;

            if (!(partioning_dimension & L_xor_H) ||
                ((partioning_dimension & LoHi_C) ==
                 (partioning_dimension & qLo)))
                Hi = LoHi_H;
            else
                Lo = HiLo_H;
        }

        if (!HB_matches_RQ (Lo, qLo, L_xor_H, &region,
                    &local_xJ, &local_tT, &local_W))
            exit(1); /* an error has occurred */

        /* calc new values for lobound & hibound: ie reduce
           range */

        if (j > 0)
        {
            for (R_xor_L = region ^ qLo, N = DIM - 1; R_xor_L;
                    R_xor_L >>= 1, N--)
                if (R_xor_L & 1)
                {
                    lobound.hcode[N] = 0;
                    lomask |= 1<< (DIM-1-N);
                }
            for (R_xor_H = region ^ qHi, N = DIM - 1; R_xor_H;
                    R_xor_H >>= 1, N--)
                if (R_xor_H & 1)
                {
                    hibound.hcode[N] = U_int_MAX;
                    himask |= 1<< (DIM-1-N);
                }
        }

        /* add in a bit of the hcode */
        element = (DIM * j) / 32;
        if ((DIM * j) % 32 > 32 - DIM)
        {
            next_match->hcode[element] |= Lo << (DIM * j) % 32;
            next_match->hcode[element+1] |=
                                Lo >> 32 - (DIM * j) % 32;
        }
        else
            next_match->hcode[element] |=
                                Lo << (DIM * j) - element * 32;
    }
    return 1; /* next-match found */

}

/*----------------------------------------------------------------------------*/
#endif	// #if 0 - this is HRQ_VER == 7
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
// this is HRQ_VER == 8

/*=============================================================*/
/*                    HB_nextmatch_RQ                          */
/*=============================================================*/

#define jORDER WORDBITS

static bool HB_nextmatch_RQ( PU_int *LoBound, PU_int *HiBound,
            HU_int *next_match, HU_int *page_key, int K,
            U_int xJ, U_int tT, U_int W,
            U_int LoMask, U_int HiMask, int DIMS )
{
    U_int   mask, qLo, qHi,
            L_xor_H, R_xor_L, R_xor_H, region,
            H, element, temp,
            partioning_dimension,
            Lo, Hi, LoHi_H, HiLo_H, LoHi_C, tSL, tSH,
            temp1, temp2,
            NUMPOINTS = (U_int)1 << DIMS,
            g_all_ones = ( (U_int)1 << DIMS ) -1;
    int     i, j, N;

/* variables used for back-tracking */
	U_int	qLoBAK, qHiBAK,
			L_xor_HBAK,
			HBAK,
			xJBAK, tTBAK, WBAK,
			LoMaskBAK, HiMaskBAK,
			LoBAK, HiBAK;
	bool	Backup = false;
	// allocate storage to these and delete them when finished
	HU_int	*next_matchBAK = new HU_int[DIMS];
	PU_int	*LoBoundBAK = new PU_int[DIMS];
	PU_int	*HiBoundBAK = new PU_int[DIMS];
	int 	KBAK;

    while (K >= 0)
    {
	mask = 1 << K;
/* step 1: find the n-points of the quadrants in which the
   query lower and upper bound points lie */

        for (i = 0, qLo = qHi = 0; i < DIMS; i++)
        {
            if (LoBound[i] & mask)
            {
                qLo |= 1 << (DIMS-1-i);
            }
            if (HiBound[i] & mask)
            {
                qHi |= 1 << (DIMS-1-i);
            }
        }

/* step 2: extract DIMS consecutive bits from page_key */

        element = (DIMS * K) / jORDER;
        H = page_key[element];
        if ((DIMS * K) % jORDER > jORDER - DIMS)
        {
            temp = page_key[element + 1];
            H >>= (DIMS * K) % jORDER;
            temp <<= jORDER - (DIMS * K) % jORDER;
            H |= temp;
        }
        else
        {
            H >>= (DIMS * K) % jORDER;
        }
        /* the & masks out spurious highbit values */
        if (DIMS < jORDER)
        {
            H &= (1 << DIMS) -1;
        }

/* step 3: find in which dimensions the query lower and upper
   bounds have the same values */

        L_xor_H = qLo ^ qHi;
        L_xor_H = ~L_xor_H;

/* step 4: do a binary search of the derived-keys of the
   quadrants in the current-state to find the quadrant with the
   lowest derived-key that is greater than or equal to the value
   of H and which is within the query range */

        for (Lo = 0, Hi = NUMPOINTS - 1 ; Lo != Hi ; )
        {

/* the middle 2 ordered derived-keys of quadrants in a state or
   part of a state. They enable us to divide it in 2; one is the
   upper bound of the lower part and the other is the lower
   bound of the upper part */

            LoHi_H = Lo + (Hi - Lo) / 2;
            HiLo_H = LoHi_H + 1;

/* find the partitioning-dimension */

            tSL = LoHi_H ^ LoHi_H / 2;
            tSH = HiLo_H ^ HiLo_H / 2;

            if (xJ % DIMS != 0)
            {
/* right circular shift tSL and tSH by local_xJ bits */
                temp1 = tSL >> xJ % DIMS;
                temp2 = tSL << DIMS - xJ % DIMS;
                tSL   = temp1 | temp2;
                tSL  &= ((U_int)1 << DIMS) - 1;
                temp1 = tSH >> xJ % DIMS;
                temp2 = tSH << DIMS - xJ % DIMS;
                tSH   = temp1 | temp2;
                tSH  &= ((U_int)1 << DIMS) - 1;
            }

/* the dimension in which all points in the lower half of the
   state (or part) have the same value and the points in the
   upper half have the opposite value */

            partioning_dimension = tSL ^ tSH;

/* The n-point of the quadrant with the highest derived-key in
   the 'lower' half of the state (or part) */

            LoHi_C = W ^ tT ^ tSL;

/* check to see if the query intersects with the lower half of
   the state (or part) i.e. if its lower bound is within it (the
   first test checks whether the query intersects both halves -
   both parts of the test are needed. NB 'partioning_dimension'
   and 'LoHi_C' change as the loop iterates) */

            if (!(partioning_dimension & L_xor_H) ||
                ((partioning_dimension & LoHi_C) ==
                    (partioning_dimension & qLo)))
            {
                    /* there is a match in the lower half of
                       the state; either the range spans both
                       halves (so far as the partioning_dimension
                       dim is concerned) or LoHi_C is within the
                       range... */

/* check to see if H is in the lower half of the state
   (or part) */

                if (LoHi_H >= H)
                {

/* H is in the lower half of the state (or part). If the query
   spans both halves of the state (or part), record the limits
   of the top half of the state and note that search may continue
   there by assigning 'TRUE' to 'Backup'. If a next-match is not
   subsequently found in the lower half of the state the search
   may in the upper half (this is back-tracking) */

                    if (!(partioning_dimension & L_xor_H))
                    {
						qLoBAK = qLo; qHiBAK = qHi;
						L_xor_HBAK = L_xor_H;
						HBAK = H;
						xJBAK = xJ; tTBAK = tT; WBAK = W;
						KBAK = K;
						keycopy( next_matchBAK, next_match, DIMS );
						keycopy( LoBoundBAK, LoBound, DIMS );
						keycopy( HiBoundBAK, HiBound, DIMS );
						LoMaskBAK = LoMask; HiMaskBAK = HiMask;
						LoBAK = HiLo_H; HiBAK = Hi;
						Backup = true;
                    }

/* set (reduce) the 'Hi' value for the next iteration of the
   binary search */

                    Hi = LoHi_H;
                }
                else
                {

/* ...but H is in the upper half of the state (or part),
   therefore check if the range also intersects with the upper
   half...if so, set (increase) the 'Lo' value for the next
   iteration of the binary search */

                    if (!(partioning_dimension & L_xor_H))
                    {
                        Lo = HiLo_H;
                    }
                    else
                    {

/*...otherwise H is in the upper half of the state (or part) but
   the query doesn't intersect with it; therefore, if we
   previously found a (larger) half state (or part) with which
   the query intersects, backtrack to it, changing the values of
   'Lo' & 'Hi' accordingly */

                        if (true == Backup)
                        {
							qLo = qLoBAK; qHi = qHiBAK;
							L_xor_H = L_xor_HBAK;
							H = HBAK;
							xJ = xJBAK; tT = tTBAK; W = WBAK;
							K = KBAK;
							keycopy( next_match, next_matchBAK, DIMS );
							keycopy( LoBound, LoBoundBAK, DIMS );
							keycopy( HiBound, HiBoundBAK, DIMS );
							LoMask = LoMaskBAK; HiMask = HiMaskBAK;
							Lo = LoBAK; Hi = HiBAK;
							Backup = false;
                        }
                        else
                        {
                        	delete [] next_matchBAK;
                        	delete [] LoBoundBAK;
                        	delete [] HiBoundBAK;
                         	return false;
                        }
       		    }
       		}
            }
            else
            {

/* the query region does not intersect with the upper half of
   the state (or part), therefore set (increase) the value of
   'Lo' for the next iteration of the binary search */

                Lo = HiLo_H;
   	    }

	} /* end of binary search */

/* check to see if the quadrant whose derived-key is 'Lo' lies
   within the query region and, if so, determine the
   characteristics of the next-state (using variables
   local_xJ, local_tT and local_W). 'region' is set to be
   the n-point corresponding to 'Lo'. If this function call
   fails, back-track if possible or else return to an earlier
   call to HB_nextmatch_RQ to resume the search, i.e. backtrack
   there */

        if (false == HB_matches_RQ( Lo, qLo, L_xor_H, &region,
                    &xJ, &tT, &W, DIMS))
        {
            if (Backup == false)
            {
            	cerr << "ERROR in ......\n";
             	delete [] next_matchBAK;
             	delete [] LoBoundBAK;
             	delete [] HiBoundBAK;
                return false;
            }

			qLo = qLoBAK; qHi = qHiBAK;
			L_xor_H = L_xor_HBAK;
			H = HBAK;
			xJ = xJBAK; tT = tTBAK; W = WBAK;
			K = KBAK;
			keycopy( next_match, next_matchBAK, DIMS );
			keycopy( LoBound, LoBoundBAK, DIMS );
			keycopy( HiBound, HiBoundBAK, DIMS );
			LoMask = LoMaskBAK; HiMask = HiMaskBAK;
			Lo = LoBAK; Hi = HiBAK;
			Backup = false;
            continue;
        }

        if (K > 0)
        {

/* restrict the query region (lobound and hibound) to that part
   which intersects the region containing the next-match, if it
   exists, i.e. 'Lo'. 'R_xor_L' has bits set to 1 if 'region'
   and 'qlo' have different values in the correponding
   dimensions. Since 'region' is within the query range, its
   bits will be 1 and qLo's will be 0 therefore in restricting
   the query range, we increase the value of the current order's
   bit in the query to 1 and set the lower bits to 0. Since at
   higher orders we are only going to encounter 0 bits in qLo
   values, it's simplest to set the whole value of the query
   lower bound to 0. A similar logic applies to 'R_xor_H' if the
   query region was wholly within the subspace which is defined
   by 'region' then both 'R_xor_L' and 'R_xor_H' would be 0 in
   all bits and the query region would not get changed */

            for (R_xor_L = region ^ qLo, N = DIMS - 1;
                    R_xor_L;
                    R_xor_L >>= 1, N--)
            {
                if (R_xor_L & 1)
                {
                    LoBound[N] = 0;
                    LoMask |= 1 << (DIMS-1-N);
                }
            }
            for (R_xor_H = region ^ qHi, N = DIMS - 1;
                    R_xor_H;
                    R_xor_H >>= 1, N--)
            {
                if (R_xor_H & 1)
                {
                    HiBound[N] = ~0;	// UINT_MAX;
                    HiMask |= 1 << (DIMS-1-N);
                }
            }
        }

        if (Lo == H)
        {
            if (LoMask == g_all_ones && HiMask == g_all_ones)
            {
                    /* the page_key is a match to the query */
                keycopy( next_match, page_key, DIMS );
               	delete [] next_matchBAK;
               	delete [] LoBoundBAK;
               	delete [] HiBoundBAK;
                return true; // 2;  a next-match has been found
            }
        }

/* step 6: add in DIMS bits of the current-quadrant's
   derived-key into the next-match */

        element = (DIMS * K) / jORDER;
        if ((DIMS * K) % jORDER > jORDER - DIMS)
        {
            next_match[element] |=
                        Lo << (DIMS * K) % jORDER;
            next_match[element + 1] |=
                        Lo >> jORDER - (DIMS * K) % jORDER;
        }
        else
        {
            next_match[element] |=
                        Lo << (DIMS * K) - element * jORDER;
        }

	K--;

	if (Lo > H)
	{
            break; /* we can move into fast-forward mode */
 	}
    } /* end of while loop */

	if ( K < 0 )
	{
               	delete [] next_matchBAK;
               	delete [] LoBoundBAK;
               	delete [] HiBoundBAK;
		return true;
	}

	/* now we can fast-forward:  */
    for (j = K; j >= 0; j--)
    {

/* at his stage, we are no longer interseted in the value of bits
   in the page-key; a next-match is guaranteed to exist - all
   that is needed now is to find the lowest derived-key in the
   current query region. */

        if (LoMask == g_all_ones && HiMask == g_all_ones)
        {
		delete [] next_matchBAK;
		delete [] LoBoundBAK;
		delete [] HiBoundBAK;
		return true; /* next_match found */
        }

        mask = 1 << j;

        for (i = 0, qLo = qHi = 0; i < DIMS; i++)
        {
            if (LoBound[i] & mask)
                qLo |= 1 << (DIMS-1-i);
            if (HiBound[i] & mask)
                qHi |= 1 << (DIMS-1-i);
        }
        L_xor_H = qLo ^ qHi;
        L_xor_H = ~L_xor_H;

        /* do a binary search of the state */
        for (Lo = 0, Hi = NUMPOINTS - 1; Lo != Hi ; )
        {
            LoHi_H = Lo + (Hi - Lo) / 2;
            HiLo_H = LoHi_H + 1;

            tSL = LoHi_H ^ LoHi_H / 2;
            tSH = HiLo_H ^ HiLo_H / 2;

            if (xJ % DIMS != 0)
            {
                temp1 = tSL >> xJ % DIMS;
                temp2 = tSL << DIMS - xJ % DIMS;
                tSL   = temp1 | temp2;
                tSL  &= ((U_int)1 << DIMS) - 1;
                temp1 = tSH >> xJ % DIMS;
                temp2 = tSH << DIMS - xJ % DIMS;
                tSH   = temp1 | temp2;
                tSH  &= ((U_int)1 << DIMS) - 1;
            }
            partioning_dimension = tSL ^ tSH;

            LoHi_C = W ^ tT ^ tSL;

            if (!(partioning_dimension & L_xor_H) ||
                ((partioning_dimension & LoHi_C) ==
                 (partioning_dimension & qLo)))
            {
                Hi = LoHi_H;
            }
            else
            {
                Lo = HiLo_H;
            }
        }

        if (false == HB_matches_RQ(Lo, qLo, L_xor_H, &region,
                    &xJ, &tT, &W, DIMS ))
        {
            exit(1); /* an error has occurred */
        }

        /* calc new values for LoBound & HiBound: ie reduce
           range */

        if (j > 0)
        {
            for (R_xor_L = region ^ qLo, N = DIMS - 1; R_xor_L;
                    R_xor_L >>= 1, N--)
            {
                if (R_xor_L & 1)
                {
                    LoBound[N] = 0;
                    LoMask |= 1 << (DIMS-1-N);
                }
            }
            for (R_xor_H = region ^ qHi, N = DIMS - 1; R_xor_H;
                    R_xor_H >>= 1, N--)
            {
                if (R_xor_H & 1)
                {
                    HiBound[N] = ~0;	// U_int_MAX;
                    HiMask |= 1 << (DIMS-1-N);
                }
            }
        }

        /* add in a bit of the hcode */
        element = (DIMS * j) / jORDER;
        if ((DIMS * j) % jORDER > jORDER - DIMS)
        {
            next_match[element] |= Lo << (DIMS * j) % jORDER;
            next_match[element+1] |=
                                Lo >> jORDER - (DIMS * j) % jORDER;
        }
        else
        {
            next_match[element] |=
                                Lo << (DIMS * j) - element * jORDER;
        }
    }
	delete [] next_matchBAK;
	delete [] LoBoundBAK;
	delete [] HiBoundBAK;

	return true; /* next-match found */

}

/*============================================================================*/
/*                            H_nextmatch_RQ				      */
/*============================================================================*/
/* in call to HB_nextmatch_RQ_v5d(), state is not used; the last 3 parameters
   are xJ, tT and W */

bool H_nextmatch_RQ( PU_int *LB, PU_int *UB,
		HU_int *match, HU_int *key, int dimensions )
{
    bool retval;
//    int i;
    // these get changed by HB_nextmatch_RQ, therefore use copies
    PU_int *LoBound = new PU_int[dimensions];
    PU_int *HiBound = new PU_int[dimensions];

    keycopy( LoBound, LB, dimensions );
    keycopy( HiBound, UB, dimensions );

    retval = HB_nextmatch_RQ( LoBound, HiBound, match, key, NUMBITS-1, 0, 0, 0, 0, 0, dimensions );

    delete [] LoBound;
    delete [] HiBound;

	/*
	printf("key: ");
    for ( i = 0; i < DIM; i++ )
    {
        printf("%u ", key[i]);
    }
    printf("\n");

    if ( retval )
    {
        printf("match: ");
        for ( i = 0; i < DIM; i++ )
        {
            printf("%u ", match[i]);
        }
        printf("\n");
    }
    else printf("no match\n");
	 */

    return retval;
}


