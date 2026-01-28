// Copyright (C) Jonathan Lawder 2001-2011

#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include "btree.h"
#ifdef __MSDOS__
	#include "..\utils\utils.h"
#else
	#include "../utils/utils.h"
#endif


/*============================================================================*/
/*                            		                          	      */
/*                            BTnode	                          	      */
/*                            		                          	      */
/*============================================================================*/

/*============================================================================*/
/*                            BTnode::BTnode	                      	      */
/*============================================================================*/
// constructor
BTnode::BTnode( int dims, int n_entries ) {

	dimensions = dims;
	node_entries = n_entries;
	node_entry_size = sizeof(U_int) * dimensions + sizeof(X);
	int raw_data_bytes = node_entries * node_entry_size;

	// the data block that makes up a node (excluding parent pointer)
/*	raw_data = (unsigned char*)malloc( node_entries * node_entry_size );*/
	raw_data = new unsigned char[raw_data_bytes];
	// initialise the raw_data
	memset( raw_data, '\0', raw_data_bytes );

	// the array of pointers to keys in a node
/*	Hkey = (U_int**)malloc(sizeof(U_int*) * node_entries);*/
//	Hkey = new (U_int*)[sizeof(U_int*) * node_entries];
	Hkey = new U_int*[node_entries];

	// the array of pointers to lpage (downptr) associated with a leaf (inner)
/*	XX = (X**)malloc(sizeof(X*) * node_entries);*/
//	XX = new (X*)[sizeof(X*) * node_entries];
	XX = new X*[node_entries];

	// setup the pointers into raw_data
	// the number of U_ints in a node entry is dim + 1
	// NB Hkey[0] and XX[0] are not used; they overlap with the BTnodehdr
	// NB the BTnodehdr 'occupies' the same amount of space as a <Hkey,XX>
	btnodehdr = (BTnodehdr*)raw_data;
	for ( int i = 0; i < node_entries; i++ )
	{
		// Hkeys come before XX's in a node entry
		Hkey[i] = (U_int*)raw_data + i * (dimensions + 1);
		XX[i] = (X*)raw_data + i * (dimensions + 1) + dimensions;
	}
};

/*============================================================================*/
/*                            BTnode::~BTnode	                          	      */
/*============================================================================*/
// destructor
BTnode::~BTnode() {
//	delete btnodehdr; not needed; this just points into 'data'
	delete [] Hkey;
	delete [] XX;
	delete [] raw_data;
}

/*============================================================================*/
/*                            idxi_find_slot				      */
/*============================================================================*/
/* Find the SLOT (element) in a node that corresponds to a 'key' less
	 than or equal to the parameter 'key'. IN THIS RESPECT IT IS DIFFERENT
	 FROM db_pageinbuffer().
   This function can be used both for leaf and inner node searching.
   Note that since the first array element in a node is used for
	 info, 'parent' and 'next' pointers (if it's a leaf) and info,
	 'parent' and 'first' pointers (if it's an inner node), the 'low'
	 parameter is never < 1.
   Therefore, if return value is 0, this can mean one of 2 things:
	 leaf: the key is lower than the lowest key in the leaf.
	 inner: the 'first' pointer should be followed down to the next level.
   A non-zero positive number denotes that an EXACT match was found.
*/
int BTnode::idxi_find_slot( HU_int *key )
{
	int		mid, i, j, low = 1, high = btnodehdr->size;
	HU_int	*entry, *helement;	// no memory allocation needed
#ifdef JKLDEBUGxxx
	cout << "idxi_find_slot : high = " << high << " btnodehdr->size = " << btnodehdr->size ;
	if (in_HDR->flags & isROOT) cout << " - root\n"; else cout << " - not root\n";
#endif
	if (high < low)
		errorexit("ERROR in BTnode::idxi_find_slot()\n");

	while (high >= low)
	{
		mid = (high + low) / 2;
// CHECK CHECK CHECK !!! - just say 'helement = Hkey[mid]' ????
		entry = Hkey[mid];	// point at an entry in the node
		helement = entry;	// point at the first U_int in an entry
#ifdef JKLDEBUGxxx
	cout << "..." << *helement << "...\n";
#endif
		for (j = 0, i = dimensions - 1; i >= 0; i--)
		{
			if (key[i] == helement[i])
				continue;
			if (key[i] < helement[i])
			{	j = -1; break;} /* a < b */
			else
			{	j = 1; break;}  /* a > b */
		}
		if (j == 0)        /* they are equal */
			return mid;
		if (j == -1)         /* key is lower */
			high = mid - 1;
		else
			low = mid + 1;
	}
	/* key not found - return the next lower key's element as a negative no */
	if (j == -1)
		mid--;
	if (mid > 0)
		mid *= -1;
	return mid;
}

/*============================================================================*/
/*                            idxi_find_leaf				      */
/*============================================================================*/
// finds the LEAF which will contain the page which may contain a key
BTnode* BTnode::idxi_find_leaf( HU_int *key )
{
	int i;

	if (btnodehdr->flags & isLEAF)
		return this;
	i = idxi_find_slot( key );
	if (i == 0)
		return in_HDR->firstptr->idxi_find_leaf( key );
	if (i < 0)
		i *= -1;
	return in_ENTRY[i]->downptr->idxi_find_leaf( key );
}

/*============================================================================*/
/*                            idxi_find_prev				      */
/*============================================================================*/
// used by idx_get_prev()
int BTnode::idxi_find_prev( HU_int *key, int lpage, BTnode *left )
{
	int slot;
	BTnode *newnode, *myleft;

	slot = idxi_find_slot( key );

	/* if thisnode not a LEAF: set up parameters for next recursive call */
	if (!(lf_HDR->flags & isLEAF))
	{
		if (slot < 0)
			slot *= -1;

		/* setup myleft */
		if (slot == 0)
		{
			newnode = lf_HDR->firstptr;
			if (left == NULL)
				myleft = NULL;
			else
				myleft = left->in_ENTRY[left->lf_HDR->size]->downptr;
		}
		else
		{
			newnode = in_ENTRY[slot]->downptr;
			if (slot == 1)
				myleft = in_HDR->firstptr;
			else
				myleft = in_ENTRY[slot - 1]->downptr;
		}
		return newnode->idxi_find_prev( key, lpage, myleft );
	}
	else
	{
		if (slot < 1 || slot > lf_HDR->size)
			errorexit("ERROR 1 in idxi_find_prev(): key absent\n");
		if (lpage != lf_ENTRY[slot]->lpage)
			errorexit("ERROR 2 in idxi_find_prev(): key-page mismatch\n");
		if (slot > 1)
			return lf_ENTRY[slot - 1]->lpage;
		if (left == NULL)
			return -1;
		return left->lf_ENTRY[left->lf_HDR->size]->lpage;
	}
}

/*============================================================================*/
/*                            idxi_merge_nodes				      */
/*============================================================================*/
/* Move all data from the RIGHT node to the end of the LEFT node and delete
	 entry for right node in parent.
   NB if leaf, left grows by 'size' of right, else size+1. */

// called by the left node
int BTnode::idxi_merge_nodes( BTnode *right )
{
	int i, nbytes, nobj, slot, isleaf;
	BTnode *Parent;
	HU_int *key;	// no memory allocation needed

	if (!(this && right))
		errorexit("ERROR 1 in idxi_merge_nodes(): merge with NULL node(s)\n");
	isleaf = right->lf_HDR->flags & isLEAF;
	Parent = right->lf_HDR->parent;
	if (!Parent)
		errorexit("ERROR 2 in idxi_merge_nodes(): merging orphan node(s)\n");
	if (Parent != lf_HDR->parent)
		errorexit("ERROR 3 in idxi_merge_nodes(): merging nodes "
				  "with different parents\n");

	// remove entry for right node from parent
	key = right->Hkey[1];
	// find slot for right node in parent
	slot = Parent->idxi_find_slot( key );
	if (!isleaf)
		slot *= -1;
	if (slot < 1 || slot > Parent->in_HDR->size)
		errorexit("ERROR 4 in idxi_merge_nodes(): lowest key in (right "
				  "hand)\nnode to be merged not compatible with parent\n");

	if (!isleaf) /* fill in key in first free slot in left */
		keycopy ( Hkey[in_HDR->size + 1], Parent->Hkey[slot], dimensions );
//		left->X.ientry[left->X.in.size + 1].Hkey =
//			parent->X.ientry[slot].Hkey;
/*		BT->keycopy(&left->X.ientry[left->X.in.size + 1].Hkey,
			&parent->X.ientry[slot].Hkey);*/
	/* do the deletion in parent */
	nbytes = (Parent->in_HDR->size - slot) * node_entry_size;
	memmove( Parent->Hkey[slot], Parent->Hkey[slot + 1], nbytes );
/*	memset(&parent->X.lentry[parent->X.lf.size], NULL, sizeof (innerENTRY));*/
	Parent->lf_HDR->size--;

	/* move data from right node to left node */
	nobj = right->in_HDR->size;
	if (isleaf)
	{
		nbytes = nobj * node_entry_size;
		memcpy( Hkey[lf_HDR->size + 1], right->Hkey[1], nbytes);
		lf_HDR->nextptr = right->lf_HDR->nextptr;
		lf_HDR->size += nobj;
	}
	else
	{
		in_ENTRY[in_HDR->size + 1]->downptr = right->in_HDR->firstptr;
		/* adjust parents of children of right node */
		right->in_HDR->firstptr->in_HDR->parent = this;
		for ( i = 1; i <= right->in_HDR->size; i++)
			right->in_ENTRY[i]->downptr->in_HDR->parent = this;

		nbytes = nobj * node_entry_size;
		memcpy( Hkey[in_HDR->size + 2], right->Hkey[1], nbytes );
		in_HDR->size += nobj + 1;
	}

	delete right;
	return 1;
}

/*============================================================================*/
/*                            idxi_shift_from_left			      */
/*============================================================================*/
/* move 'excess' data from left to right */
// right is the implied parameter
int BTnode::idxi_shift_from_left( BTnode *left, BTnode *anchor )
{
	int i, slot, isleaf, lsize, rsize, numtomove, nbytes;
	HU_int *key;	// no memory allocation needed

	if (!(left && this && anchor))
		errorexit("ERROR 1 in idxi_shift_from_left(): at least one required "
				  "node is NULL\n");

	lsize = left->in_HDR->size;
	rsize = in_HDR->size;

	isleaf = lf_HDR->flags & isLEAF;
	if (isleaf)
	{
		if (!(lsize + rsize >= lf_minFANOUT))
			errorexit("ERROR 2 in idx_shift_from_left(): combined node size "
					  "insufficient\nfor rebalancing\n");

		numtomove = (lsize - lf_minFANOUT) / 2;
		numtomove = numtomove == 0 ? 1 : numtomove;

		key = Hkey[1];
		slot = anchor->idxi_find_slot( key );
		if (slot < 1 || slot > anchor->in_HDR->size)
			errorexit("ERROR 3 in idxi_shift_from_left(): lowest key in "
					  "right hand node not\ncompatible with anchor node\n");
		/* move all entries in right up to make room for left's values */
		nbytes = node_entry_size * rsize;
		memmove( Hkey[numtomove + 1],  Hkey[1], nbytes );
		/* move elements from left node to right */
		nbytes = node_entry_size * numtomove;
		memcpy( Hkey[1], left->Hkey[lsize - numtomove + 1], nbytes );
/*		memset(&left->X.lentry[lsize - numtomove + 1], NULL, nbytes);*/
		/* adjust anchor key value */
		keycopy ( anchor->Hkey[slot], Hkey[1], dimensions );
/*		BT->keycopy(&anchor->X.ientry[slot].Hkey, &right->X.lentry[1].Hkey);*/
		/* adjust sizes */
		left->lf_HDR->size -= numtomove;
		lf_HDR->size += numtomove;
	}
	else
	{
/*		numtomove = (lsize - BT->inminFANOUT) / 2; */
		if (!(lsize + rsize >= in_minFANOUT))
			errorexit("ERROR 4 in idx_shift_from_leftt(): combined node size "
					  "insufficient\nfor rebalancing\n");

		numtomove = (lsize - in_minFANOUT - 1) / 2;
		key = Hkey[1];
		slot = anchor->idxi_find_slot( key );
		slot *= -1;
		if (slot < 1 || slot > anchor->in_HDR->size)
			errorexit("ERROR 5 in idxi_shift_from_left(): lowest key in "
					  "right hand node not\ncompatible with anchor node\n");
		/* 1 move right entries to position numtomove + 2 */
		nbytes = node_entry_size * rsize;
		memmove( Hkey[numtomove + 2], Hkey[1], nbytes);
		/* 2 move anchor value to right node */
		keycopy ( Hkey[numtomove + 1], anchor->Hkey[slot], dimensions );
/*		BT->keycopy(&right->X.ientry[numtomove + 1].Hkey,
				&anchor->X.ientry[slot].Hkey);*/
		/* 3 move right's first pointer to numtomove+1 element */
		in_ENTRY[numtomove + 1]->downptr =in_HDR->firstptr;
		/* 4 move value from left to anchor */
		keycopy ( anchor->Hkey[slot], left->Hkey[lsize - numtomove], dimensions );
/*		BT->keycopy(&anchor->X.ientry[slot].Hkey,
				&left->X.ientry[lsize - numtomove].Hkey);*/
		/* 5 move pointer from left to right's first */
		in_HDR->firstptr = left->in_ENTRY[lsize - numtomove]->downptr;
		// 6 move high entries in left to low end of right
		nbytes = node_entry_size * numtomove;
		memcpy( Hkey[1], left->Hkey[lsize - numtomove + 1], nbytes );
		nbytes = node_entry_size * (numtomove + 1);		// ???? redundant????
/*		memset(&left->X.ientry[lsize - numtomove], NULL, nbytes);*/
		// adjust parents of imported elements' children
		in_HDR->firstptr->in_HDR->parent = this;
		for (i = 1; i <= numtomove; i++)
			in_ENTRY[i]->downptr->in_HDR->parent = this;
		// adjust sizes
		in_HDR->size += numtomove + 1;
		left->in_HDR->size -= numtomove + 1;
	}
	if (left->in_HDR->size < in_minFANOUT || in_HDR->size < in_minFANOUT)
		errorexit("ERROR 6 in idxi_shift_from_left(): undersized node\n");
	return 1;
}

/*============================================================================*/
/*                            idxi_shift_from_right			      */
/*============================================================================*/
/* Move 'excess' data from right to left.
   NB if leaf, left grows by numtomove, else by numtomove +1 */
// left is the implied parameter
int BTnode::idxi_shift_from_right( BTnode *right, BTnode *anchor )
{
	int i, slot, isleaf, lsize, rsize, numtomove, nbytes;
	HU_int *key;	// no memory allocation needed

	if (!(this && right && anchor))
		errorexit("ERROR 1 in idxi_shift_from_right(): at least one "
				  "required node is NULL\n");
	lsize = in_HDR->size;
	rsize = right->in_HDR->size;

	isleaf = right->lf_HDR->flags & isLEAF;
	if (isleaf)
	{
		if (!(lsize + rsize >= lf_minFANOUT))
			errorexit("ERROR 2 in idxi_shift_from_right(): combined node "
					  "size insufficient\nfor rebalancing\n");

		numtomove = (rsize - lf_minFANOUT) / 2;
		numtomove = numtomove == 0 ? 1 : numtomove;

		key = right->Hkey[1];
		slot = anchor->idxi_find_slot( key );
		if (slot < 1 || slot > anchor->in_HDR->size)
			errorexit("ERROR 3 in idxi_shift_from_right(): lowest key in "
					  "right hand node not\ncompatible with anchor node\n");
		// move elements from right to left
		nbytes = node_entry_size * numtomove;
		memcpy( Hkey[lsize + 1], right->Hkey[1], nbytes );
		// move high entries in right to low end (and empty vacated slots)
		nbytes = node_entry_size * (rsize - numtomove);
		memmove( right->Hkey[1], right->Hkey[numtomove + 1], nbytes );
		nbytes = node_entry_size * numtomove;
/*		memset(&right->X.lentry[rsize - numtomove + 1], NULL, nbytes);*/
		// adjust anchor key value
		keycopy ( anchor->Hkey[slot], right->Hkey[1], dimensions );
/*		BT->keycopy(&anchor->X.ientry[slot].Hkey, &right->X.lentry[1].Hkey);*/
		// adjust sizes
		lf_HDR->size += numtomove;
		right->lf_HDR->size -= numtomove;
	}
	else
	{
/*		numtomove = (rsize - BT->inminFANOUT) / 2; */
		if (!(lsize + rsize >= in_minFANOUT))
			errorexit("ERROR 4 in idxi_shift_from_right(): combined node "
					  "size insufficient\nfor rebalancing\n");

		numtomove = (rsize - in_minFANOUT - 1) / 2;

		key = right->Hkey[1];
		slot = anchor->idxi_find_slot( key );
		slot *= -1;
		if (slot < 1 || slot > anchor->in_HDR->size)
			errorexit("ERROR 5 in idxi_shift_from_right(): lowest key in "
					  "right hand node not\ncompatible with anchor node\n");
		// move anchor value to left node
		keycopy ( Hkey[lsize + 1], anchor->Hkey[slot], dimensions );
/*		BT->keycopy(&left->X.ientry[lsize + 1].Hkey, &anchor->X.ientry[slot].Hkey);*/
		// move right's first pointer to left
		in_ENTRY[lsize + 1]->downptr = right->in_HDR->firstptr;
		// assign new pointer to right's first
		right->in_HDR->firstptr = right->in_ENTRY[numtomove + 1]->downptr;
		// move value from right to anchor
		keycopy( anchor->Hkey[slot], right->Hkey[numtomove + 1], dimensions );
/*		BT->keycopy(&anchor->X.ientry[slot].Hkey,
				&right->X.ientry[numtomove + 1].Hkey);*/
		// move elements from right to left
		nbytes = node_entry_size * numtomove;
		memcpy( Hkey[lsize + 2], right->Hkey[1], nbytes );
		// move high entries in right to low end (and empty vacated slots)
		nbytes = node_entry_size * (rsize - numtomove - 1);
		memmove( right->Hkey[1], right->Hkey[numtomove + 2], nbytes );
		nbytes = node_entry_size * (numtomove + 1);	// ???? redundant ????
/*		memset(&right->X.ientry[rsize - numtomove], NULL, nbytes);*/
		// adjust parents of imported elements' children: BEFORE size mods
		for (i = 1; i <= numtomove + 1; i++)
			in_ENTRY[lsize + i]->downptr->in_HDR->parent = this;
		// adjust sizes
		in_HDR->size += numtomove + 1;
		right->in_HDR->size -= numtomove + 1;
	}
	if (in_HDR->size < in_minFANOUT || right->in_HDR->size < in_minFANOUT)
		errorexit("ERROR 6 in idxi_shift_from_right(): undersized node\n");
	return 1;
}

/*============================================================================*/
/*                            idxi_process_underflow			      */
/*============================================================================*/
int BTnode::idxi_process_underflow( BTnode *left, BTnode *right, BTnode *LAnchor, BTnode *RAnchor )
{
	int bal_left = 0, bal_right = 0;

	/* find out if rebalancing is possible */
	if (left)
	{
		if (left->lf_HDR->flags & isLEAF)
			bal_left = (left->lf_HDR->size > lf_minFANOUT);
		else
			bal_left = (left->in_HDR->size > in_minFANOUT);
		bal_left <<= 1;
	}
	if (right)
		if (right->lf_HDR->flags & isLEAF)
			bal_right = (right->lf_HDR->size > lf_minFANOUT);
		else
			bal_right = (right->in_HDR->size > in_minFANOUT);

	/* invoke rebalancing or merging */
	switch (bal_left | bal_right)
	{
		case 0:
			/* Must merge: both left and right are half empty.
			   At least one of left & right has the same parent as thisnode:
				 if both do, it doesn't matter which one we merge with.
			   Anchors not needed since we always merge into the leftmost
				 node and any higher level index adjustment will take place
				 in the common parent */
			if (right && btnodehdr->parent == right->btnodehdr->parent)
				idxi_merge_nodes( right );
			else
				left->idxi_merge_nodes( this );
			break;

		/* Cases 1 - 3: we can do a re-balance
		   Index entries in the anchors will need adjusting */
		case 1: /* can only balance with right */
			idxi_shift_from_right( right, RAnchor );
			break;
		case 2: /* can only balance with left */
			idxi_shift_from_left( left, LAnchor );
			break;
		case 3: /* balance with either */
			if (left->in_HDR->size == right->in_HDR->size)
				/* at least 1 of l or r will have thisnode's parent */
				if(left->in_HDR->parent == in_HDR->parent)
					idxi_shift_from_left( left, LAnchor);
				else
					idxi_shift_from_right( right, RAnchor );
			else /* balance with largest sibling */
				if (left->in_HDR->size > right->in_HDR->size)
					idxi_shift_from_left( left, LAnchor );
				else
					idxi_shift_from_right( right, RAnchor );
			break;
		default:
			errorexit("ERROR in idxi_process_underflow(): cannot process\n");
	}
	return 1;
}

/*============================================================================*/
/*                            idxi_delete_from_node			      */
/*============================================================================*/
// a return value of 0 triggers collapse of root
// uses idxi_process_underflow() which uses in_minFANOUT
int BTnode::idxi_delete_from_node( HU_int *key, int lpage,
				BTnode *left, BTnode *right,
				BTnode *LAnchor, BTnode *RAnchor )
{
	int slot, tempslot, nbytes;
	BTnode *newnode, *myleft, *myright, *myLA, *myRA;

	slot = idxi_find_slot( key );

	/* if thisnode not a LEAF: set up parameters for next recursive call */
	if (!(in_HDR->flags & isLEAF))
	{
		if (slot < 0)
			slot *= -1;

		/* setup myleft and myLA */
		if (slot == 0)
		{
			newnode = in_HDR->firstptr;
			myLA = LAnchor;
			if (left == NULL)
				myleft = NULL;
			else
				myleft = left->in_ENTRY[left->in_HDR->size]->downptr;
		}
		else
		{
			newnode = in_ENTRY[slot]->downptr;
			myLA = this;
			if (slot == 1)
				myleft = in_HDR->firstptr;
			else
				myleft = in_ENTRY[slot - 1]->downptr;
		}

		/* setup myright and myRA */
		if (slot == in_HDR->size)
		{
			myRA = RAnchor;
			if (right == NULL)
				myright = NULL;
			else
				myright = right->in_HDR->firstptr;
		}
		else
		{
			myRA = this;
			myright = in_ENTRY[slot + 1]->downptr;
		}
		newnode->idxi_delete_from_node( key, lpage, myleft, myright, myLA, myRA );
	}

	/* delete the entry from the leaf */
	else
	{
		if (slot < 1 || slot > lf_HDR->size)
			errorexit("ERROR 1 in idxi_delete_from_node(): key absent\n");

		if (lf_ENTRY[slot]->lpage != lpage)
			errorexit("ERROR 2 in idxi_delete_from_node() - attempting "
					  "to delete a\nnon-existent page entry\n");

		/* Update index entry in LAnchor - if necessary
		   A copy of a key which is the first one in any leaf always exists
			 as a key in some higher level node (if any), specifically, the
			 leaf's LAnchor.
		   If the first key in a leaf is deleted the copy in the LAnchor
			 must be updated so it becomes the new first entry in the leaf.*/
		if (slot == 1 && LAnchor != NULL)
		{
			tempslot = LAnchor->idxi_find_slot( key );
			if (tempslot < 1 || tempslot > LAnchor->in_HDR->size)
				errorexit("ERROR 3 in idxi_delete_from_node(): "
						  "key not in anchor\n");

			keycopy ( LAnchor->Hkey[tempslot], Hkey[2], dimensions ) ;
/*			BT->keycopy(&LAnchor->X.ientry[tempslot].Hkey,
					  &thisnode->X.ientry[2].Hkey);*/
		/* As a deletion in a leaf may cause cascading of deletions upwards,
			 the search key may be searched for in the LAnchor at some point
			 but it won't be there - unless we also modify the search key,
			 which, since it is a pointer, will take effect in all other
			 function calls. */
			keycopy ( key, Hkey[2], dimensions );
/*			BT->keycopy(key, &thisnode->X.ientry[2].Hkey);*/
		}

		nbytes = node_entry_size * (lf_HDR->size - slot);
		memmove( Hkey[slot], Hkey[slot + 1], nbytes );
/*		memset(&thisnode->X.lentry[thisnode->X.lf.size], NULL,
			sizeof (thisnode->X.lentry[0]));*/
		lf_HDR->size--;
	}

	/* Deal with thisnode if it is the root */
	if (in_HDR->flags & isROOT)
	{
/*		if (thisnode != g_BTroot)
			errorexit("ERROR 4 in idxi_delete_from_node\n"); */
		if (lf_HDR->flags & isLEAF)
			return 1;
		if (in_HDR->size == 0)
		/* this will trigger collapse of root */
			return 0;
		else
			return 1;
	}

	/* Deal with underflow in nodes other than the root */
	if (lf_HDR->flags & isLEAF)
	{
		if (lf_HDR->size < lf_minFANOUT)
			idxi_process_underflow( left, right, LAnchor, RAnchor );
		return 1;
	}
	else
		if (in_HDR->size < in_minFANOUT)
			idxi_process_underflow( left, right, LAnchor, RAnchor );

/* DON'T do any error checking for underflow here! */

	return 1;
}

/*============================================================================*/
/*                            idxi_setup_parents			      */
/*============================================================================*/
/* For setting up 'parent' pointers between nodes in a btree which
	 has been read into memory from file. */
// must be called originally by a BTree's root node
int BTnode::idxi_setup_parents()
{
	int i;

	if (in_HDR->flags & isROOT)
		in_HDR->parent = NULL;

	if (!(lf_HDR->flags & isLEAF))
	{
		in_HDR->firstptr->in_HDR->parent = this;
		for (i = 1; i <= in_HDR->size; i++)
			in_ENTRY[i]->downptr->btnodehdr->parent = this;
		in_HDR->firstptr->idxi_setup_parents();
		for (i = 1; i <= in_HDR->size; i++)
			in_ENTRY[i]->downptr->idxi_setup_parents();
	}
	return 1;
}


/*============================================================================*/
/*                            idxi_append				      */
/*============================================================================*/
// used by idxi_setup_nextptrs(): in constructing the node list
void BTnode::idxi_append( ListNode *tail )
{
//	tail->next = static_cast<ListNode*>(getstorage(sizeof(ListNode)));
	tail->next = new ListNode;
	tail = tail->next;
	tail->node = this;
	tail->next = NULL;
}

/*============================================================================*/
/*                            idxi_setup_nextptrs			      */
/*============================================================================*/
/* traverses the index breadth first: returns the number of leaf
   node entries */
// must be called by root node
int BTnode::idxi_setup_nextptrs()
{
	int i, j;
	ListNode *head, *tail, *temp;

//	head = static_cast<ListNode*>(getstorage(sizeof(ListNode)));
	head = new ListNode;
	head->node = this;
	head->next = NULL;
	tail = temp = head;

	while (!(head->node->in_HDR->flags & isLEAF))
	{
		head->node->in_HDR->firstptr->idxi_append( tail );
		tail = tail->next;
		for (i = 1; i <= head->node->in_HDR->size; i++)
		{
			head->node->in_ENTRY[i]->downptr->idxi_append( tail );
			tail = tail->next;
		}
		/* remove inner nodes as list is traversed */
		head = head->next;
//		free(temp);
		delete temp;
		temp = head;
	}
	/* now head points at node that points at first leaf */
	j = 0;
	while (head != tail)
	{
		j += head->node->lf_HDR->size;
		head->node->lf_HDR->nextptr = head->next->node;
		head = head->next;
//		free(temp);
		delete temp;
		temp = head;
	}
	j += head->node->lf_HDR->size;
	head->node->lf_HDR->nextptr = NULL;
//	free(head);
	delete head;

	return j;
}

/*============================================================================*/
/*                            idxi_read_file				      */
/*============================================================================*/
// this is assuming that the file contains at least one node - what if it doesn't?!!
BTnode* BTnode::idxi_read_file(fstream& f, int dims, int n_entries, int n_entry_size )
{
	int i;
/*	BTnode *node = static_cast<BTnode*>(getstorage(sizeof(BTnode)));*/
	BTnode *node = new BTnode(dims, n_entries );

	f.read( reinterpret_cast<char*>(node->raw_data), n_entries * n_entry_size );
	if (! f)
		errorexit("ERROR in idxi_read_file()\n");

	if (!(node->in_HDR->flags & isLEAF))
	{
		node->in_HDR->firstptr = idxi_read_file( f, dims, n_entries, n_entry_size );
		for (i = 1; i <= node->in_HDR->size; i++)
			node->in_ENTRY[i]->downptr = idxi_read_file( f, dims, n_entries, n_entry_size );
	}
	return node;
}

/*============================================================================*/
/*                            idxi_free_btree				      */
/*============================================================================*/
// Release btree from memory
void BTnode::idxi_free_btree()
{
	int i;

	if (!this)
		return;

	if (!(btnodehdr->flags & isLEAF))
	{
		btnodehdr->firstptr->idxi_free_btree();
		for (i = 1; i <= btnodehdr->size; i++)
			XX[i]->downptr->idxi_free_btree();
	}
	delete this;
}


/*============================================================================*/
/*                            		                          	      */
/*                            BTree	                          	      */
/*                            		                          	      */
/*============================================================================*/

/*============================================================================*/
/*                            BTree::BTree	                          	      */
/*============================================================================*/
// constructor
BTree::BTree()
{
}

/*============================================================================*/
/*                            BTree::BTree	                          	      */
/*============================================================================*/
// constructor
BTree::BTree( string db_name, int dims, int n_entries )
{
	name = db_name;
	root = NULL;
	dimensions = dims;
	node_entries = n_entries;
	node_entry_size = sizeof(U_int) * dimensions + sizeof(X);

//	idxfile = NULL;
}

/*============================================================================*/
/*                            BTree::~BTree	                          	      */
/*============================================================================*/
// destructor
BTree::~BTree()
{
	if (root)
		root->idxi_free_btree();
	root = NULL;
}

/*============================================================================*/
/*                            BTree::free_root	                          	      */
/*============================================================================*/
// or freeing root created in db_create()
void BTree::free_root()
{
	if (root)
		root->idxi_free_btree();
	root = NULL;
}

/*============================================================================*/
/*                            idxi_make_new_root			      */
/*============================================================================*/
// creates a new root which points at 2 children
void BTree::idxi_make_new_root( BTnode *right, HU_int *newkey )
{
	BTnode	*left;

	left = root;

	root = new BTnode( dimensions, node_entries );
	
	root->in_HDR->parent = NULL;
	root->in_HDR->flags = isROOT;
	root->in_HDR->size = 1;
	root->in_HDR->firstptr = left;
	root->in_ENTRY[1]->downptr = right;
	left->lf_HDR->flags &= ~isROOT;
	if (left->lf_HDR->flags & isLEAF)
	{
		keycopy ( root->Hkey[1], right->Hkey[1], dimensions );
		left->lf_HDR->parent = right->lf_HDR->parent = root ;
/*		BT->keycopy(&BT->root->X.ientry[1].Hkey, &right->X.lentry[1].Hkey);*/
	}
	else
	{
		keycopy ( root->Hkey[1], newkey, dimensions );
/*		BT->keycopy(&BT->root->X.ientry[1].Hkey, newkey);*/
		left->in_HDR->parent = right->in_HDR->parent = root;
	}
}

/*============================================================================*/
/*                            idxi_split_leaf				      */
/*============================================================================*/
/* splits a full leaf into 2, inserts new key into parent
   (New node size same as p or one bigger) */
// calls idxi_make_new_root()  - OK
void BTree::idxi_split_leaf( BTnode *p )
{
	BTnode *New;
	int numtomove, nbytes, size;

	/* make new node */
	New = new BTnode( dimensions, node_entries );
	New->lf_HDR->flags = isLEAF;
	New->lf_HDR->nextptr = p->lf_HDR->nextptr;
	New->lf_HDR->parent = p->lf_HDR->parent;
	p->lf_HDR->nextptr = New;
	/* deal with 'data' */
/*	split = p->X.lf.size / 2 + 1;
	nbytes = sizeof (p->X.lentry[0]) * (p->X.lf.size - split + 1);
	memmove(&New->X.lentry[1], &p->X.lentry[split], nbytes); */
	size = p->lf_HDR->size;
	numtomove = size / 2;
	nbytes = node_entry_size * numtomove;
	memcpy( New->Hkey[1], p->Hkey[size - numtomove + 1], nbytes );
		/* erase data from p */
/*	memset(&p->X.lentry[split], NULL, nbytes); */
/*	memset(&p->X.lentry[size - numtomove + 1], NULL, nbytes);*/
		/* adjust sizes - MUST deal with New first! */
/*	New->X.lf.size = p->X.lf.size - split + 1;
	p->X.lf.size = split - 1; */
	New->lf_HDR->size = numtomove;
	p->lf_HDR->size -= numtomove;

	/* insert new key into parent: p & New parents get adjusted as reqd.*/
	if (p->lf_HDR->flags & isROOT)
	{
		idxi_make_new_root( New, NULL );
	}
	else
	{
		idxi_insert_in_node( p->lf_HDR->parent, New->Hkey[1], 0, New );
	}
}

/*============================================================================*/
/*                            idxi_split_inner				      */
/*============================================================================*/
/* Splits a full inner node into 2, inserts new key into parent.
   New node takes the higher values of p and is the same size as p or
	 one smaller. */
// calls idxi_make_new_root() - OK
void BTree::idxi_split_inner( BTnode *p )
{
	BTnode *New = new BTnode( dimensions, node_entries );
	int promotee, nbytes, i;
	HU_int *promkey = new U_int[dimensions];

	/* make new node */
	New->in_HDR->flags = 0;
	New->in_HDR->parent = p->in_HDR->parent;

	promotee = p->in_HDR->size / 2 + 1;
	keycopy ( promkey, p->Hkey[promotee], dimensions );
/*	BT->keycopy(&promkey, &p->X.ientry[promotee].Hkey);*/
	New->in_HDR->firstptr = p->in_ENTRY[promotee]->downptr;

	/* deal with 'data' */
	nbytes = node_entry_size * (p->in_HDR->size - promotee);
	memcpy (New->Hkey[1], p->Hkey[promotee + 1], nbytes );
		/* erase data from p */
/*	memset(&p->X.ientry[promotee], NULL, nbytes + sizeof (p->X.ientry[0]));*/
		/* adjust sizes - MUST deal with New first! */
	New->in_HDR->size = p->in_HDR->size - promotee;
	p->in_HDR->size = promotee - 1;

	if (p->in_HDR->size < in_minFANOUT || New->in_HDR->size < in_minFANOUT)
		errorexit("ERROR in idxi_split_inner(): undersized node\n");

	/* insert new key into parent: p & New parents get adjusted as reqd.*/
	if (p->in_HDR->flags & isROOT)
		idxi_make_new_root( New, promkey );
	else
		idxi_insert_in_node( p->lf_HDR->parent, promkey, 0, New );

	/* modify parent pointers in children of New */
	New->in_HDR->firstptr->in_HDR->parent = New;
	for ( i = 1; i <= New->in_HDR->size; i++ )
		New->in_ENTRY[i]->downptr->in_HDR->parent = New;
}

/*============================================================================*/
/*                            idxi_insert_in_node			      */
/*============================================================================*/
/*  Generic function to insert into a btree node a <key, pointer> into
	  an inner node or a <key, page> into a leaf node.
	If 'p' is a leaf, we are inserting a <key, page>: q is not used and is
	  therefore NULL.
	If 'p' is an inner node, we are inserting a <key, pointer>: lpage is not
	  used and is therefore 0.
	'p' is the node which will have something inserted in it.
	'key' is the key which will be inserted in p.
	'lpage' is the 'data' which will be inserted in p if it's a leaf.
	'q' is the pointer which will be inserted in p if it's an inner.
	If after insertion, the node is full, it splits itself.
*/
// calls idxi_split_inner() and idxi_split_leaf() - OK
void BTree::idxi_insert_in_node( BTnode *p, HU_int *key, int lpage, BTnode *q )
{
	int slot, nbytes;

/*
	void idx_dump(BTree *);

 	if (slot > 0)
	{
		int i;
		FILE *f;

		if (p->X.lf.flags & isLEAF)
			printf("Index node is a leaf; Page no.: %u\n",lpage);
		else
			printf("Index node is internal\n");

		printf("Key :");
		for (i = DIM - 1; i >= 0; i--)
			printf("%-15u", key->hcode[i]);
		printf("\n");

		idx_dump(BT);

		errorexit("ERROR in idxi_insert_in_node(): key to be "
				  "inserted already in index\n");
	}
 */

	slot = p->idxi_find_slot( key );
 	slot *= -1;
	slot++; /* this is where the new entry goes */

	/* else it's an inner node */
	nbytes = node_entry_size * (p->btnodehdr->size - slot + 1);
	memmove( p->Hkey[slot + 1], p->Hkey[slot], nbytes );
/*	memset(&p->X.ientry[slot], NULL, sizeof (p->X.ientry[0]));*/
	/* insert entry */
	keycopy ( p->Hkey[slot], key, dimensions );
/*	BT->keycopy(&p->X.ientry[slot].Hkey, key);*/
	p->btnodehdr->size++;

	if (p->lf_HDR->flags & isLEAF)
	{
		p->lf_ENTRY[slot]->lpage = lpage;
		if (p->lf_HDR->size > lf_FANOUT)
			idxi_split_leaf( p );
	}
	else
	{
		p->in_ENTRY[slot]->downptr = q;
		if (p->in_HDR->size > in_FANOUT)
			idxi_split_inner( p );
	}
}

/*============================================================================*/
/*                            idxi_write_file				      */
/*============================================================================*/
/* For writing a btree to file:
   traverses btree depth first, writing nodes to file as it goes */
int BTree::idxi_write_file( BTnode *node )
{
	int i;

	if (!node)
		return 0;
	idxfile.write(reinterpret_cast<char*>(node->raw_data), node_entries * node_entry_size);
	if (!(node->in_HDR->flags & isLEAF))
	{
		idxi_write_file(node->in_HDR->firstptr);
		for (i = 1; i <= node->in_HDR->size; i++)
			idxi_write_file(node->in_ENTRY[i]->downptr);
	}
	return 1;
}

/*============================================================================*/
/*                            idx_read					      */
/*============================================================================*/
// Don't call this function if the index is empty!
// Returns the number of leaf node entries
// (ie the number of pages in the database)
int BTree::idx_read( string fname )
{
 	string filename;

	if ( fname == "" )
		filename = name + ".idx";
	else
		filename = fname;

//	string filename = name + ".idx";

	idxfile.open( filename.c_str(), ios::in | ios::binary );
	if (! idxfile)
		errorexit( "ERROR 1 in idx_read(): writing index to file \n" );

	root = root->idxi_read_file( idxfile, dimensions, node_entries, node_entry_size );

	if (! idxfile)
		errorexit( "ERROR 1 in idx_read(): reading index file\n" );

	/* make sure there's nothing more to read */
	(void) idxfile.get();
	if (! idxfile.eof())
		errorexit( "ERROR 2 in idx_read(): end of index file not read\n" );
	idxfile.close();

	root->idxi_setup_parents();
	int i = root->idxi_setup_nextptrs();

	return i;
}

/*============================================================================*/
/*                            idx_write					      */
/*============================================================================*/
/* Writes a btree index to file.
   'extn' is the file name extension if required: MUST include the dot (.) */
int BTree::idx_write( string fname )
{
	string filename;

	if ( fname == "" )
		filename = name + ".idx";
	else
		filename = fname;

	idxfile.open( filename.c_str(), ios::out | ios::binary );

	idxi_write_file( root );
//	if (idxfile.ferror())  ????		????
//		errorexit("ERROR in idx_write(): writing index to file \n");

	idxfile.close();
	return 1;
}

/*============================================================================*/
/*                            idx_get_last_page				      */
/*============================================================================*/
int BTree::idx_get_last_page()
{
	BTnode *p = root;

	while(!(p->in_HDR->flags & isLEAF))
		p = p->in_ENTRY[p->in_HDR->size]->downptr;

	return p->lf_ENTRY[p->lf_HDR->size]->lpage;
}

/*============================================================================*/
/*                            idx_get_next				      */
/*============================================================================*/
/* Similar to idx_search() except that it returns the lpage no. of
	 the page with the next highest key to the parameter key value -
	 or -1 if there is no higher key in the index. */
int BTree::idx_get_next( HU_int *key, int lpage )
{
	BTnode *p;
	int slot;

	if (!root || root->lf_HDR->size == 0)
		return -1;

	p = root->idxi_find_leaf( key );
	/* search the LEAF page */
	slot = p->idxi_find_slot( key );
	if (slot < 1)
		errorexit("ERROR 1 in idx_get_next(): key not in index\n");
	if (lpage != p->lf_ENTRY[slot]->lpage)
			errorexit("ERROR 2 in idx_get_next(): key-page mismatch\n");
	if (slot == p->lf_HDR->size)
	{
		if (p->lf_HDR->nextptr == NULL)
			return -1;  /* this is the last page in the database */
		return p->lf_HDR->nextptr->lf_ENTRY[1]->lpage;
	}
	return p->lf_ENTRY[slot + 1]->lpage;
}

/*============================================================================*/
/*                            idx_get_next_key				      */
/*============================================================================*/
/* Similar to idx_get_next() except that it returns the key of
	 the page with the next highest key to the parameter key value -
	 or returns the parameter key if there is no higher key in the index. */
HU_int* BTree::idx_get_next_key( HU_int *key, int lpage )
{
	BTnode *p;
	int slot;

	if (!root || root->lf_HDR->size == 0)
		return key;

	p = root->idxi_find_leaf( key );
	/* search the LEAF page */
	slot = p->idxi_find_slot( key );
	if (slot < 1)
		errorexit("ERROR 1 in idx_get_next_key(): key not in index\n");
	if (lpage != p->lf_ENTRY[slot]->lpage)
	{
			errorexit("ERROR 2 in idx_get_next_key(): key-page mismatch\n");
	}
	if (slot == p->lf_HDR->size)
	{
		if (p->lf_HDR->nextptr == NULL)
			return key;  /* this is the last page in the database */
		return p->lf_HDR->nextptr->Hkey[1];
	}
	return p->Hkey[slot + 1];
}

/*============================================================================*/
/*                            idx_get_prev				      */
/*============================================================================*/
/* Similar to idx_search() except that it returns the lpage no. of
	 the page with the next lowest key to the parameter key value -
	 or -1 if there is no lower key in the index. */
int BTree::idx_get_prev( HU_int *key, int lpage )
{
	if (!root || root->lf_HDR->size == 0)
		return -1;
	return root->idxi_find_prev( key, lpage, NULL );
}

/*============================================================================*/
/*                            idx_insert_key				      */
/*============================================================================*/
/*  top level call to insert a <key, lpage> pair into a leaf in the btree -
	takes care of the empty index and the full leaf situations */
int BTree::idx_insert_key( HU_int *key, int lpage )
{
	BTnode *p;

	if (!root)
	{
/*		root = static_cast<BTnode*>(getstorage(sizeof(BTnode)));*/
		root = new BTnode( dimensions, node_entries );
		root->lf_HDR->flags = (isROOT | isLEAF);
		root->lf_HDR->size = 1;
		keycopy ( root->Hkey[1], key, dimensions );
/*		BT->keycopy(&BT->root->X.lentry[1].Hkey, key);*/
		root->lf_ENTRY[1]->lpage = lpage;
		root->lf_HDR->nextptr = root->lf_HDR->parent = NULL;
		return 1;
	}
	p = root->idxi_find_leaf( key );

	idxi_insert_in_node( p, key, lpage, NULL );
	return 1;
}

/*============================================================================*/
/*                            idx_delete_key				      */
/*============================================================================*/
int BTree::idx_delete_key( HU_int *key, int lpage )
{
	int i;

	if (!root || /* g_BTroot->X.lf.flags & isLEAF && */
			root->lf_HDR->size == 0)
		errorexit("ERROR in idx_delete_key(): database is empty\n");
	i = root->idxi_delete_from_node( key, lpage, NULL, NULL, NULL, NULL );
	if (i == 0) /* we need to collapse root */
	{
		root = root->in_HDR->firstptr;
		root->in_HDR->flags |= isROOT;
		delete root->in_HDR->parent;
		root->in_HDR->parent = NULL;
	}
	if (root->lf_HDR->flags & isLEAF)
		if (root->lf_HDR->size == 0)
		{   /* root is an empty leaf */
			delete root;
			root = NULL;
		}
	return 1;
}

/*============================================================================*/
/*                            idx_search				      */
/*============================================================================*/
/* finds the LOGICAL PAGE which may contain a 'key' value */
int BTree::idx_search( HU_int *key )
{
	BTnode *p;
	int slot;

	if (root == NULL || root->lf_HDR->size == 0)
	{
		printf("Database is empty\n");
		return 0;
	}
	p = root->idxi_find_leaf( key );
	/* search the LEAF page */
	slot = p->idxi_find_slot( key );
	if (slot == 0)
		errorexit( "ERROR in idx_search : key is lower than any key in the database" );
//		return -1; key is lower than any key in the database
	if (slot < 0)
		slot *= -1;
#ifdef JKLDEBUGxxx
	cout << "idx_search : lpage = " << p->lf_ENTRY[slot]->lpage << endl;
#endif		
	return p->lf_ENTRY[slot]->lpage;
}

/*============================================================================*/
/*                            idx_dump					      */
/*============================================================================*/
// outputs the contents of the leaves of the btree
// the index must aready be in memory
void BTree::idx_dump( string filename )
{
	int		i, j;
	BTnode		*p = root;
	U_int		*h;
	fstream		f;

	f.open( filename.c_str(), ios::out );
	
	f << "From the index\n\n";

	if (!p)
	{
		f << "Index empty\n";
		return;
	}
	while (!(p->in_HDR->flags & isLEAF))
		p = p->in_HDR->firstptr;

	do
	{
		for (i = 1; i <= p->lf_HDR->size; i++)
		{
			f << "PAGE NO.: " << p->lf_ENTRY[i]->lpage << endl;
			f << "KEY : ";
			for (j = dimensions - 1, h = p->Hkey[i]; j >= 0; j--)
				f << setw(15) << h[j];
			f << endl;
		}
		p = p->lf_HDR->nextptr;
	} while (p);

	f.close();
}
