#ifndef HUNGARIAN_H
#define HUNGARIAN_H

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#define cell      long
#define CELL_STR  "%li"

#define llong    int_fast64_t
#define byte     int_fast8_t
#define boolean  int_fast8_t
#define null     0
#define true     1
#define false    0


ssize_t** kuhn_match(cell**, size_t, size_t);
void kuhn_reduceRows(cell**, size_t, size_t);
byte** kuhn_mark(cell**, size_t, size_t);
boolean kuhn_isDone(byte**, boolean*, size_t, size_t);
size_t* kuhn_findPrime(cell**, byte**, boolean*, boolean*, size_t, size_t);
void kuhn_altMarks(byte**, size_t*, size_t*, ssize_t*, ssize_t*, size_t*, size_t, size_t);
void kuhn_addAndSubtract(cell**, boolean*, boolean*, size_t, size_t);
ssize_t** kuhn_assign(byte**, size_t, size_t);
cell score_assignment(cell**, size_t, size_t);


/**
 * Bit set, a set of fixed number of bits/booleans
 */
typedef struct
{
    /**
     * The set of all limbs, a limb consist of 64 bits
     */
    llong* limbs;

    /**
     * Singleton array with the index of the first non-zero limb
     */
    size_t* first;

    /**
     * Array the the index of the previous non-zero limb for each limb
     */
    size_t* prev;

    /**
     * Array the the index of the next non-zero limb for each limb
     */
    size_t* next;

} BitSet;


BitSet new_BitSet(size_t);
void BitSet_set(BitSet, size_t);
void BitSet_unset(BitSet, size_t);
ssize_t BitSet_any(BitSet) __attribute__((pure));

size_t lb(llong x) __attribute__((const));
void print(cell**, size_t, size_t, ssize_t**);

#endif  /*  HUNGARIAN_H */
