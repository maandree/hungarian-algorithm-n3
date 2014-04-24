// -*- mode: c, coding: utf-8 -*-

/**
 * 𝓞(n³) implementation of the Hungarian algorithm
 * 
 * Copyright (C) 2011, 2014  Mattias Andrée
 * 
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it
 * and/or modify it under the terms of the Do What The Fuck You Want
 * To Public License, Version 2, as published by Sam Hocevar. See
 * http://sam.zoy.org/wtfpl/COPYING for more details.
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


#ifndef RANDOM_DEVICE
#define RANDOM_DEVICE "/dev/urandom"
#endif


#define cell      long
#define CELL_STR  "%li"

#define llong    int_fast64_t
#define byte     int_fast8_t
#define boolean  int_fast8_t
#define null     0
#define true     1
#define false    0



//new cell[X]
#define new_cells(X)  new_longs(X)

//new boolean[X]
#define new_booleans(X)  new_longs(X)

//new byte[X]
#define new_bytes(X)  malloc((size_t)(X))

//new llong[X]
#define new_llongs(X)  malloc((size_t)(X) << 3)

//new long[X]
#if !(defined __LP64__ || defined __LLP64__)
#  define new_longs(X)  malloc((size_t)(X) << 2) /*32-bit*/
#else
#  define new_longs(X)  malloc((size_t)(X) << 3) /*64-bit*/
#endif

//new float[X]
#define new_floats(X)  malloc((size_t)(X) << 2)

//new double[X]
#define new_doubles(X)  malloc((size_t)(X) << 3)

//new ?[][X]
#define new_arrays(X)  new_longs(X)


#ifdef DEBUG
#  define debug(X) fprintf(stderr, "\033[31m%s\033[m\n", X)
#else
#  define debug(X) 
#endif




/**
 * Cell marking:  none
 */
#define UNMARKED  0L

/**
 * Cell marking:  marked
 */
#define MARKED    1L

/**
 * Cell marking:  prime
 */
#define PRIME     2L



/**
 * Bit set, a set of fixed number of bits/booleans
 */
typedef struct
{
    /**
     * The set of all limbs, a limb consist of 64 bits
     */
    cell* limbs;
    
    /**
     * Singleton array with the index of the first non-zero limb
     */
    long* first;
    
    /**
     * Array the the index of the previous non-zero limb for each limb
     */
    long* prev;
    
    /**
     * Array the the index of the next non-zero limb for each limb
     */
    long* next;
    
} BitSet;



long** kuhn_match(cell** table, long n, long m);
void kuhn_reduceRows(cell** t, long n, long m);
byte** kuhn_mark(cell** t, long n, long m);
boolean kuhn_isDone(byte** marks, boolean* colCovered, long n, long m);
long* kuhn_findPrime(cell** t, byte** marks, boolean* rowCovered, boolean* colCovered, long n, long m);
void kuhn_altMarks(byte** marks, long* altRow, long* altCol, long* colMarks, long* rowPrimes, long* prime, long n, long m);
void kuhn_addAndSubtract(cell** t, boolean* rowCovered, boolean* colCovered, long n, long m);
long** kuhn_assign(byte** marks, long n, long m);

BitSet new_BitSet(long size);
void BitSet_set(BitSet this, long i);
void BitSet_unset(BitSet this, long i);
long BitSet_any(BitSet this) __attribute__((pure));

long lb(llong x) __attribute__((const));



void print(cell** t, long n, long m, long** assignment);

int main(int argc, char** argv __attribute__((unused)))
{
    FILE* urandom = fopen(RANDOM_DEVICE, "r");
    unsigned int seed;
    fread(&seed, sizeof(unsigned int), 1, urandom);
    srand(seed);
    fclose(urandom);
    
    
    long n = 10, m = 15;
    
    long i;
    cell** t = new_arrays(n);
    cell** table = new_arrays(n);
    long j;
    if (argc < 2)
        for (i = 0; i < n; i++)
	{
	    *(t + i) = new_llongs(m);
	    *(table + i) = new_llongs(m);
	    for (j = 0; j < m; j++)
	        *(*(table + i) + j) = *(*(t + i) + j) = (cell)(random() & 63);
	}
    else
    {
        long x;
	scanf("%li", &n);
	scanf("%li", &m);
        t = new_arrays(n);
        table = new_arrays(n);
        for (i = 0; i < n; i++)
	{
	    *(t + i) = new_llongs(m);
	    *(table + i) = new_llongs(m);
	    for (j = 0; j < m; j++)
	    {
	        scanf(CELL_STR, &x);
	        *(*(table + i) + j) = *(*(t + i) + j) = x;
	    }
	}
    }
    
    printf("\nInput:\n\n");
    print(t, n, m, null);
    
    long** assignment = kuhn_match(table, n, m);
    printf("\nOutput:\n\n");
    print(t, n, m, assignment);
    
    cell sum = 0;
    for (i = 0; i < n; i++)
    {
        sum += *(*(t + *(*(assignment + i) + 0)) + *(*(assignment + i) + 1));
	free(*(assignment + i));
	free(*(table + i));
	free(*(t + i));
    }
    free(assignment);
    free(table);
    free(t);
    printf("\n\nSum: %li\n\n", sum);
    
    return 0;
}

void print(cell** t, long n, long m, long** assignment)
{
    long i, j;
    
    long** assigned = new_arrays(n);
    for (i = 0; i < n; i++)
    {
        *(assigned + i) = new_longs(m);
	for (j = 0; j < m; j++)
	    *(*(assigned + i) + j) = 0;
    }
    if (assignment != null)
        for (i = 0; i < n; i++)
	    (*(*(assigned + **(assignment + i)) + *(*(assignment + i) + 1)))++;
    
    for (i = 0; i < n; i++)
    {
	printf("    ");
	for (j = 0; j < m; j++)
	{
	    if (*(*(assigned + i) + j))
		printf("\033[%lim", 30 + *(*(assigned + i) + j));
	    printf("%5li%s\033[m   ", (cell)(*(*(t + i) + j)), (*(*(assigned + i) + j) ? "^" : " "));
        }
	printf("\n\n");
	
	free(*(assigned + i));
    }
    
    free(assigned);
}



/**
 * Calculates an optimal bipartite minimum weight matching using an
 * O(n³)-time implementation of The Hungarian Algorithm, also known
 * as Kuhn's Algorithm.
 * 
 * @param   table  The table in which to perform the matching
 * @param   n      The height of the table
 * @param   m      The width of the table
 * @return         The optimal assignment, an array of row–coloumn pairs
 */
long** kuhn_match(cell** table, long n, long m)
{
    long i;
    
    /* not copying table since it will only be used once */
    
    kuhn_reduceRows(table, n, m);
    byte** marks = kuhn_mark(table, n, m);
    
    boolean* rowCovered = new_booleans(n);
    boolean* colCovered = new_booleans(m);
    for (i = 0; i < n; i++)
    {
        *(rowCovered + i) = false;
        *(colCovered + i) = false;
    }
    for (i = n; i < m; i++)
        *(colCovered + i) = false;
    
    long* altRow = new_longs(n * m);
    long* altCol = new_longs(n * m);
    
    long* rowPrimes = new_longs(n);
    long* colMarks  = new_longs(m);
    
    long* prime;
    
    for (;;)
    {
	if (kuhn_isDone(marks, colCovered, n, m))
	    break;
	
        for (;;)
	{
	    prime = kuhn_findPrime(table, marks, rowCovered, colCovered, n, m);
	    if (prime != null)
	    {
		kuhn_altMarks(marks, altRow, altCol, colMarks, rowPrimes, prime, n, m);
		for (i = 0; i < n; i++)
		{
		    *(rowCovered + i) = false;
		    *(colCovered + i) = false;
		}
		for (i = n; i < m; i++)
		    *(colCovered + i) = false;
		free(prime);
		break;
	    }
	    kuhn_addAndSubtract(table, rowCovered, colCovered, n, m);
	}
    }
    
    free(rowCovered);
    free(colCovered);
    free(altRow);
    free(altCol);
    free(rowPrimes);
    free(colMarks);
    
    long** rc = kuhn_assign(marks, n, m);
    
    for (i = 0; i < n; i++)
        free(*(marks + i));
    free(marks);
    
    return rc;
}


/**
 * Reduces the values on each rows so that, for each row, the
 * lowest cells value is zero, and all cells' values is decrease
 * with the same value [the minium value in the row].
 * 
 * @param  t  The table in which to perform the reduction
 * @param  n  The table's height
 * @param  m  The table's width
 */
void kuhn_reduceRows(cell** t, long n, long m)
{
    long i, j;
    cell min;
    cell* ti;
    for (i = 0; i < n; i++)
    {
        ti = *(t + i);
        min = *ti;
	for (j = 1; j < m; j++)
	    if (min > *(ti + j))
	        min = *(ti + j);
	
	for (j = 0; j < m; j++)
	    *(ti + j) -= min;
    }
}


/**
 * Create a matrix with marking of cells in the table whose
 * value is zero [minimal for the row]. Each marking will
 * be on an unique row and an unique column.
 * 
 * @param   t  The table in which to perform the reduction
 * @param   n  The table's height
 * @param   m  The table's width
 * @return     A matrix of markings as described in the summary
 */
byte** kuhn_mark(cell** t, long n, long m)
{
    long i, j;
    byte** marks = new_arrays(n);
    byte* marksi;
    for (i = 0; i < n; i++)
    {
	marksi = *(marks + i) = new_bytes(m);
        for (j = 0; j < m; j++)
 	    *(marksi + j) = UNMARKED;
    }
    
    boolean* rowCovered = new_booleans(n);
    boolean* colCovered = new_booleans(m);
    for (i = 0; i < n; i++)
    {
        *(rowCovered + i) = false;
        *(colCovered + i) = false;
    }
    for (i = 0; i < m; i++)
        *(colCovered + i) = false;
    
    for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
	    if ((!*(rowCovered + i)) && (!*(colCovered + j)) && (*(*(t + i) + j) == 0))
	    {
	        *(*(marks + i) + j) = MARKED;
		*(rowCovered + i) = true;
		*(colCovered + j) = true;
	    }
    
    free(rowCovered);
    free(colCovered);
    return marks;
}


/**
 * Determines whether the marking is complete, that is
 * if each row has a marking which is on a unique column.
 *
 * @param   marks       The marking matrix
 * @param   colCovered  An array which tells whether a column is covered
 * @param   n           The table's height
 * @param   m           The table's width
 * @return              Whether the marking is complete
 */
boolean kuhn_isDone(byte** marks, boolean* colCovered, long n, long m)
{
    long i, j;
    for (j = 0; j < m; j++)
        for (i = 0; i < n; i++)
	    if (*(*(marks + i) + j) == MARKED)
	    {
	        *(colCovered + j) = true;
		break;
	    }
    
    long count = 0;
    for (j = 0; j < m; j++)
        if (*(colCovered + j))
	    count++;
    
    return count == n;
}


/**
 * Finds a prime
 * 
 * @param   t           The table
 * @param   marks       The marking matrix
 * @param   rowCovered  Row cover array
 * @param   colCovered  Column cover array
 * @param   n           The table's height
 * @param   m           The table's width
 * @return              The row and column of the found print, <code>null</code> will be returned if none can be found
 */
long* kuhn_findPrime(cell** t, byte** marks, boolean* rowCovered, boolean* colCovered, long n, long m)
{
    long i, j;
    BitSet zeroes = new_BitSet(n * m);
    
    for (i = 0; i < n; i++)
        if (!*(rowCovered + i))
	    for (j = 0; j < m; j++)
	        if ((!*(colCovered + j)) && (*(*(t + i) + j) == 0))
		  BitSet_set(zeroes, i * m + j);
    
    long p, row, col;
    boolean markInRow;
    
    for (;;)
    {
        p = BitSet_any(zeroes);
	if (p < 0)
        {
	    free(zeroes.limbs);
	    free(zeroes.first);
	    free(zeroes.next);
	    free(zeroes.prev);
	    return null;
	}
	
	row = p / m;
	col = p % m;
	
	*(*(marks + row) + col) = PRIME;
	
	markInRow = false;
	for (j = 0; j < m; j++)
	    if (*(*(marks + row) + j) == MARKED)
	    {
		markInRow = true;
		col = j;
	    }
	
	if (markInRow)
	{
	    *(rowCovered + row) = true;
	    *(colCovered + col) = false;
	    
	    for (i = 0; i < n; i++)
	        if ((*(*(t + i) + col) == 0) && (row != i))
		{
		    if ((!*(rowCovered + i)) && (!*(colCovered + col)))
		        BitSet_set(zeroes, i * m + col);
		    else
		        BitSet_unset(zeroes, i * m + col);
		}
	    
	    for (j = 0; j < m; j++)
	        if ((*(*(t + row) + j) == 0) && (col != j))
		{
		    if ((!*(rowCovered + row)) && (!*(colCovered + j)))
		        BitSet_set(zeroes, row * m + j);
		    else
		        BitSet_unset(zeroes, row * m + j);
		}
	    
	    if ((!*(rowCovered + row)) && (!*(colCovered + col)))
	        BitSet_set(zeroes, row * m + col);
	    else
	        BitSet_unset(zeroes, row * m + col);
	}
	else
	{
	    long* rc = new_longs(2);
	    *rc = row;
	    *(rc + 1) = col;
	    free(zeroes.limbs);
	    free(zeroes.first);
	    free(zeroes.next);
	    free(zeroes.prev);
	    return rc;
	}
    }
}


/**
 * Removes all prime marks and modifies the marking
 *
 * @param  marks      The marking matrix
 * @param  altRow     Marking modification path rows
 * @param  altCol     Marking modification path columns
 * @param  colMarks   Markings in the columns
 * @param  rowPrimes  Primes in the rows
 * @param  prime      The last found prime
 * @param  n          The table's height
 * @param  m          The table's width
 */
void kuhn_altMarks(byte** marks, long* altRow, long* altCol, long* colMarks, long* rowPrimes, long* prime, long n, long m)
{
    long index = 0, i, j;
    *altRow = *prime;
    *altCol = *(prime + 1);
    
    for (i = 0; i < n; i++)
    {
        *(rowPrimes + i) = -1;
        *(colMarks + i) = -1;
    }
    for (i = n; i < m; i++)
        *(colMarks + i) = -1;
    
    for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
	    if (*(*(marks + i) + j) == MARKED)
	        *(colMarks + j) = i;
	    else if (*(*(marks + i) + j) == PRIME)
	        *(rowPrimes + i) = j;
    
    long row, col;
    for (;;)
    {
        row = *(colMarks + *(altCol + index));
	if (row < 0)
	    break;
	
	index++;
	*(altRow + index) = row;
	*(altCol + index) = *(altCol + index - 1);
	
	col = *(rowPrimes + *(altRow + index));
	
	index++;
	*(altRow + index) = *(altRow + index - 1);
	*(altCol + index) = col;
    }
    
    byte* markx;
    for (i = 0; i <= index; i++)
    {
        markx = *(marks + *(altRow + i)) + *(altCol + i);
        if (*markx == MARKED)
	    *markx = UNMARKED;
	else
	    *markx = MARKED;
    }
    
    byte* marksi;
    for (i = 0; i < n; i++)
    {
        marksi = *(marks + i);
        for (j = 0; j < m; j++)
	    if (*(marksi + j) == PRIME)
	        *(marksi + j) = UNMARKED;
    }
}


/**
 * Depending on whether the cells' rows and columns are covered,
 * the the minimum value in the table is added, subtracted or
 * neither from the cells.
 *
 * @param  t           The table to manipulate
 * @param  rowCovered  Array that tell whether the rows are covered
 * @param  colCovered  Array that tell whether the columns are covered
 * @param  n           The table's height
 * @param  m           The table's width
 */
void kuhn_addAndSubtract(cell** t, boolean* rowCovered, boolean* colCovered, long n, long m)
{
    long i, j;
    cell min = 0x7FFFffffL;
    for (i = 0; i < n; i++)
        if (!*(rowCovered + i))
	    for (j = 0; j < m; j++)
	        if ((!*(colCovered + j)) && (min > *(*(t + i) + j)))
		    min = *(*(t + i) + j);
    
    for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
	{
	    if (*(rowCovered + i))
	        *(*(t + i) + j) += min;
	    if (*(colCovered + j) == false)
	        *(*(t + i) + j) -= min;
	}
}


/**
 * Creates a list of the assignment cells
 * 
 * @param   marks  Matrix markings
 * @param   n      The table's height
 * @param   m      The table's width
 * @return         The assignment, an array of row–coloumn pairs
 */
long** kuhn_assign(byte** marks, long n, long m)
{
    long** assignment = new_arrays(n);
    
    long i, j;
    for (i = 0; i < n; i++)
    {
        *(assignment + i) = new_longs(2);
        for (j = 0; j < m; j++)
	    if (*(*(marks + i) + j) == MARKED)
	    {
	        **(assignment + i) = i;
		*(*(assignment + i) + 1) = j;
	    }
    }
    
    return assignment;
}


/**
 * Constructor for BitSet
 *
 * @param   size  The (fixed) number of bits to bit set should contain
 * @return        The a unique BitSet instance with the specified size
 */
BitSet new_BitSet(long size)
{
    BitSet this;
    
    long c = size >> 6L;
    if (size & 63L)
        c++;
    
    this.limbs = new_llongs(c);
    this.prev = new_longs(c + 1L);
    this.next = new_longs(c + 1L);
    *(this.first = new_longs(1)) = 0;
    
    long i;
    for (i = 0; i < c; i++)
    {
        *(this.limbs + i) = 0LL;
        *(this.prev + i) = *(this.next + i) = 0L;
    }
    *(this.prev + c) = *(this.next + c) = 0L;
    
    return this;
}

/**
 * Turns on a bit in a bit set
 * 
 * @param  this  The bit set
 * @param  i     The index of the bit to turn on
 */
void BitSet_set(BitSet this, long i)
{
    long j = i >> 6L;
    llong old = *(this.limbs + j);
    
    *(this.limbs + j) |= 1LL << (llong)(i & 63L);
    
    if ((!*(this.limbs + j)) ^ (!old))
    {
        j++;
	*(this.prev + *(this.first)) = j;
	*(this.prev + j) = 0;
	*(this.next + j) = *(this.first);
	*(this.first) = j;
    }
}

/**
 * Turns off a bit in a bit set
 * 
 * @param  this  The bit set
 * @param  i     The index of the bit to turn off
 */
void BitSet_unset(BitSet this, long i)
{
    long j = i >> 6L;
    llong old = *(this.limbs + j);
    
    *(this.limbs + j) &= ~(1LL << (llong)(i & 63L));
    
    if ((!*(this.limbs + j)) ^ (!old))
    {
	j++;
	long p = *(this.prev + j);
	long n = *(this.next + j);
	*(this.prev + n) = p;
	*(this.next + p) = n;
	if (*(this.first) == j)
	    *(this.first) = n;
    }
}

/**
 * Gets the index of any set bit in a bit set
 * 
 * @param   this  The bit set
 * @return        The index of any set bit
 */
long BitSet_any(BitSet this)
{
    if (*(this.first) == 0L)
        return -1;
    
    long i = *(this.first) - 1L;
    return lb(*(this.limbs + i) & -*(this.limbs + i)) + (i << 6L);
}


/**
 * Calculates the floored binary logarithm of a positive integer
 *
 * @param   value  The integer whose logarithm to calculate
 * @return         The floored binary logarithm of the integer
 */
long lb(llong value)
{
    long rc = 0L;
    llong v = value;
    
    if (v & (int_fast64_t)0xFFFFFFFF00000000LL)  {  rc |= 32L;  v >>= 32LL;  }
    if (v & (int_fast64_t)0x00000000FFFF0000LL)  {  rc |= 16L;  v >>= 16LL;  }
    if (v & (int_fast64_t)0x000000000000FF00LL)  {  rc |=  8L;  v >>=  8LL;  }
    if (v & (int_fast64_t)0x00000000000000F0LL)  {  rc |=  4L;  v >>=  4LL;  }
    if (v & (int_fast64_t)0x000000000000000CLL)  {  rc |=  2L;  v >>=  2LL;  }
    if (v & (int_fast64_t)0x0000000000000002LL)     rc |=  1L;
    
    return rc;
}

