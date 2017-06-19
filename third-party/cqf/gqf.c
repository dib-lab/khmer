/*
 * =====================================================================================
 *
 *       Filename:  gqf.c
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  2017-02-04 03:40:58 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (ppandey@cs.stonybrook.edu)
 *                  Rob Johnson (rob@cs.stonybrook.edu)
 *   Organization:  Stony Brook University
 *
 * =====================================================================================
 */

# include <assert.h>
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/mman.h>

#include "gqf.h"

/* Can be
	 0 (choose size at run-time),
	 8, 16, 32, or 64 (for optimized versions),
	 or other integer <= 56 (for compile-time-optimized bit-shifting-based versions)
	 */

/******************************************************************
 * Code for managing the metadata bits and slots w/o interpreting *
 * the content of the slots.
 ******************************************************************/

#define METADATA_WORD(qf,field,slot_index) (get_block((qf), (slot_index) / SLOTS_PER_BLOCK)->field[((slot_index) % SLOTS_PER_BLOCK) / 64])
#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) - 1ULL)
#define MAX_VALUE(nbits) ((1ULL << (nbits)) - 1)


static inline int popcnt(uint64_t val)
{
	asm("popcnt %[val], %[val]"
			: [val] "+r" (val)
			:
			: "cc");
	return val;
}

static inline int64_t bitscanreverse(uint64_t val)
{
	if (val == 0) {
		return -1;
	} else {
		asm("bsr %[val], %[val]"
				: [val] "+r" (val)
				:
				: "cc");
		return val;
	}
}

static inline int popcntv(const uint64_t val, int ignore)
{
	if (ignore % 64)
		return popcnt (val & ~BITMASK(ignore % 64));
	else
		return popcnt(val);
}

// Returns the number of 1s up to (and including) the pos'th bit
// Bits are numbered from 0
static inline int bitrank(uint64_t val, int pos) {
	val = val & ((2ULL << pos) - 1);
	asm("popcnt %[val], %[val]"
			: [val] "+r" (val)
			:
			: "cc");
	return val;
}

/**
 * Returns the position of the k-th 1 in the 64-bit word x.
 * k is 0-based, so k=0 returns the position of the first 1.
 *
 * Uses the broadword selection algorithm by Vigna [1], improved by Gog
 * and Petri [2] and Vigna [3].
 *
 * [1] Sebastiano Vigna. Broadword Implementation of Rank/Select
 *    Queries. WEA, 2008
 *
 * [2] Simon Gog, Matthias Petri. Optimized succinct data
 * structures for massive data. Softw. Pract. Exper., 2014
 *
 * [3] Sebastiano Vigna. MG4J 5.2.1. http://mg4j.di.unimi.it/
 * The following code is taken from
 * https://github.com/facebook/folly/blob/b28186247104f8b90cfbe094d289c91f9e413317/folly/experimental/Select64.h
 */
const uint8_t kSelectInByte[2048] = {
	8, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0,
	1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0,
	2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0,
	1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0,
	3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 7, 0,
	1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0,
	2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0,
	1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
	4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0,
	1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 8, 8, 8, 1,
	8, 2, 2, 1, 8, 3, 3, 1, 3, 2, 2, 1, 8, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2,
	2, 1, 8, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1,
	4, 3, 3, 1, 3, 2, 2, 1, 8, 6, 6, 1, 6, 2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1, 6, 4,
	4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1,
	3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 8, 7, 7, 1, 7, 2,
	2, 1, 7, 3, 3, 1, 3, 2, 2, 1, 7, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
	7, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3,
	3, 1, 3, 2, 2, 1, 7, 6, 6, 1, 6, 2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1, 6, 4, 4, 1,
	4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2,
	2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 8, 8, 8, 8, 8, 8, 8, 2,
	8, 8, 8, 3, 8, 3, 3, 2, 8, 8, 8, 4, 8, 4, 4, 2, 8, 4, 4, 3, 4, 3, 3, 2, 8, 8,
	8, 5, 8, 5, 5, 2, 8, 5, 5, 3, 5, 3, 3, 2, 8, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3,
	4, 3, 3, 2, 8, 8, 8, 6, 8, 6, 6, 2, 8, 6, 6, 3, 6, 3, 3, 2, 8, 6, 6, 4, 6, 4,
	4, 2, 6, 4, 4, 3, 4, 3, 3, 2, 8, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3, 3, 2,
	6, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2, 8, 8, 8, 7, 8, 7, 7, 2, 8, 7,
	7, 3, 7, 3, 3, 2, 8, 7, 7, 4, 7, 4, 4, 2, 7, 4, 4, 3, 4, 3, 3, 2, 8, 7, 7, 5,
	7, 5, 5, 2, 7, 5, 5, 3, 5, 3, 3, 2, 7, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3,
	3, 2, 8, 7, 7, 6, 7, 6, 6, 2, 7, 6, 6, 3, 6, 3, 3, 2, 7, 6, 6, 4, 6, 4, 4, 2,
	6, 4, 4, 3, 4, 3, 3, 2, 7, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3, 3, 2, 6, 5,
	5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 3, 8, 8, 8, 8, 8, 8, 8, 4, 8, 8, 8, 4, 8, 4, 4, 3, 8, 8, 8, 8, 8, 8,
	8, 5, 8, 8, 8, 5, 8, 5, 5, 3, 8, 8, 8, 5, 8, 5, 5, 4, 8, 5, 5, 4, 5, 4, 4, 3,
	8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 3, 8, 8, 8, 6, 8, 6, 6, 4, 8, 6,
	6, 4, 6, 4, 4, 3, 8, 8, 8, 6, 8, 6, 6, 5, 8, 6, 6, 5, 6, 5, 5, 3, 8, 6, 6, 5,
	6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7,
	7, 3, 8, 8, 8, 7, 8, 7, 7, 4, 8, 7, 7, 4, 7, 4, 4, 3, 8, 8, 8, 7, 8, 7, 7, 5,
	8, 7, 7, 5, 7, 5, 5, 3, 8, 7, 7, 5, 7, 5, 5, 4, 7, 5, 5, 4, 5, 4, 4, 3, 8, 8,
	8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 3, 8, 7, 7, 6, 7, 6, 6, 4, 7, 6, 6, 4,
	6, 4, 4, 3, 8, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 3, 7, 6, 6, 5, 6, 5,
	5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 5, 8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8, 5, 8, 5, 5, 4, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6,
	6, 4, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 5, 8, 8, 8, 6, 8, 6, 6, 5,
	8, 6, 6, 5, 6, 5, 5, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8,
	8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 4, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7,
	8, 7, 7, 5, 8, 8, 8, 7, 8, 7, 7, 5, 8, 7, 7, 5, 7, 5, 5, 4, 8, 8, 8, 8, 8, 8,
	8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 4,
	8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 5, 8, 7, 7, 6, 7, 6, 6, 5, 7, 6,
	6, 5, 6, 5, 5, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6,
	8, 6, 6, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7,
	8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 8,
	8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6,
	6, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7
};

static inline uint64_t _select64(uint64_t x, int k)
{
	if (k >= popcnt(x)) { return 64; }

	const uint64_t kOnesStep4  = 0x1111111111111111ULL;
	const uint64_t kOnesStep8  = 0x0101010101010101ULL;
	const uint64_t kMSBsStep8  = 0x80ULL * kOnesStep8;

	uint64_t s = x;
	s = s - ((s & 0xA * kOnesStep4) >> 1);
	s = (s & 0x3 * kOnesStep4) + ((s >> 2) & 0x3 * kOnesStep4);
	s = (s + (s >> 4)) & 0xF * kOnesStep8;
	uint64_t byteSums = s * kOnesStep8;

	uint64_t kStep8 = k * kOnesStep8;
	uint64_t geqKStep8 = (((kStep8 | kMSBsStep8) - byteSums) & kMSBsStep8);
	uint64_t place = popcnt(geqKStep8) * 8;
	uint64_t byteRank = k - (((byteSums << 8) >> place) & (uint64_t)(0xFF));
	return place + kSelectInByte[((x >> place) & 0xFF) | (byteRank << 8)];
}

// Returns the position of the rank'th 1.  (rank = 0 returns the 1st 1)
// Returns 64 if there are fewer than rank+1 1s.
static inline uint64_t bitselect(uint64_t val, int rank) {
#ifdef __SSE4_2_
	uint64_t i = 1ULL << rank;
	asm("pdep %[val], %[mask], %[val]"
			: [val] "+r" (val)
			: [mask] "r" (i));
	asm("tzcnt %[bit], %[index]"
			: [index] "=r" (i)
			: [bit] "g" (val)
			: "cc");
	return i;
#endif
	return _select64(val, rank);
}

static inline uint64_t bitselectv(const uint64_t val, int ignore, int rank)
{
	return bitselect(val & ~BITMASK(ignore % 64), rank);
}

#if BITS_PER_SLOT > 0
static inline qfblock * get_block(const QF *qf, uint64_t block_index)
{
	return &qf->blocks[block_index];
}
#else
static inline qfblock * get_block(const QF *qf, uint64_t block_index)
{
	return (qfblock *)(((char *)qf->blocks) + block_index * (sizeof(qfblock) + SLOTS_PER_BLOCK * qf->bits_per_slot / 8));
}
#endif

static inline int is_runend(const QF *qf, uint64_t index)
{
	return (METADATA_WORD(qf, runends, index) >> ((index % SLOTS_PER_BLOCK) % 64)) & 1ULL;
}

static inline int is_occupied(const QF *qf, uint64_t index)
{
	return (METADATA_WORD(qf, occupieds, index) >> ((index % SLOTS_PER_BLOCK) % 64)) & 1ULL;
}

#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64

static inline uint64_t get_slot(const QF *qf, uint64_t index)
{
	assert(index < qf->xnslots);
	return get_block(qf, index / SLOTS_PER_BLOCK)->slots[index % SLOTS_PER_BLOCK];
}

static inline void set_slot(const QF *qf, uint64_t index, uint64_t value)
{
	assert(index < qf->xnslots);
	get_block(qf, index / SLOTS_PER_BLOCK)->slots[index % SLOTS_PER_BLOCK] = value & BITMASK(qf->bits_per_slot);
}

#elif BITS_PER_SLOT > 0

/* Little-endian code ....  Big-endian is TODO */

static inline uint64_t get_slot(const QF *qf, uint64_t index)
{
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems to generate buggy code.  :/  */
	assert(index < qf->xnslots);
	uint64_t *p = (uint64_t *)&get_block(qf, index / SLOTS_PER_BLOCK)->slots[(index % SLOTS_PER_BLOCK) * BITS_PER_SLOT / 8];
	return (uint64_t)(((*p) >> (((index % SLOTS_PER_BLOCK) * BITS_PER_SLOT) % 8)) & BITMASK(BITS_PER_SLOT));
}

static inline void set_slot(const QF *qf, uint64_t index, uint64_t value)
{
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems to generate buggy code.  :/  */
	assert(index < qf->xnslots);
	uint64_t *p = (uint64_t *)&get_block(qf, index / SLOTS_PER_BLOCK)->slots[(index % SLOTS_PER_BLOCK) * BITS_PER_SLOT / 8];
	uint64_t t = *p;
	uint64_t mask = BITMASK(BITS_PER_SLOT);
	uint64_t v = value;
	int shift = ((index % SLOTS_PER_BLOCK) * BITS_PER_SLOT) % 8;
	mask <<= shift;
	v <<= shift;
	t &= ~mask;
	t |= v;
	*p = t;
}

#else

/* Little-endian code ....  Big-endian is TODO */

static inline uint64_t get_slot(const QF *qf, uint64_t index)
{
	assert(index < qf->xnslots);
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems to generate buggy code.  :/  */
	uint64_t *p = (uint64_t *)&get_block(qf, index / SLOTS_PER_BLOCK)->slots[(index % SLOTS_PER_BLOCK) * qf->bits_per_slot / 8];
	return (uint64_t)(((*p) >> (((index % SLOTS_PER_BLOCK) * qf->bits_per_slot) % 8)) & BITMASK(qf->bits_per_slot));
}

static inline void set_slot(const QF *qf, uint64_t index, uint64_t value)
{
	assert(index < qf->xnslots);
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems to generate buggy code.  :/  */
	uint64_t *p = (uint64_t *)&get_block(qf, index / SLOTS_PER_BLOCK)->slots[(index % SLOTS_PER_BLOCK) * qf->bits_per_slot / 8];
	uint64_t t = *p;
	uint64_t mask = BITMASK(qf->bits_per_slot);
	uint64_t v = value;
	int shift = ((index % SLOTS_PER_BLOCK) * qf->bits_per_slot) % 8;
	mask <<= shift;
	v <<= shift;
	t &= ~mask;
	t |= v;
	*p = t;
}

#endif

static inline uint64_t run_end(const QF *qf, uint64_t hash_bucket_index);

static inline uint64_t block_offset(const QF *qf, uint64_t blockidx)
{
	/* If we have extended counters and a 16-bit (or larger) offset
		 field, then we can safely ignore the possibility of overflowing
		 that field. */
	if (sizeof(qf->blocks[0].offset > 1) ||
			get_block(qf, blockidx)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
		return get_block(qf, blockidx)->offset;

	return run_end(qf, SLOTS_PER_BLOCK * blockidx - 1) - SLOTS_PER_BLOCK * blockidx + 1;
}

static inline uint64_t run_end(const QF *qf, uint64_t hash_bucket_index)
{
	uint64_t bucket_block_index       = hash_bucket_index / SLOTS_PER_BLOCK;
	uint64_t bucket_intrablock_offset = hash_bucket_index % SLOTS_PER_BLOCK;
	uint64_t bucket_blocks_offset = block_offset(qf, bucket_block_index);

	// uint64_t bucket_intrablock_rank   = bitrankv(get_block(qf, bucket_block_index)->occupieds, METADATA_WORDS_PER_BLOCK, bucket_intrablock_offset);
	uint64_t bucket_intrablock_rank   = bitrank(get_block(qf, bucket_block_index)->occupieds[0], bucket_intrablock_offset);

	if (bucket_intrablock_rank == 0) {
		if (bucket_blocks_offset <= bucket_intrablock_offset)
			return hash_bucket_index;
		else
			return SLOTS_PER_BLOCK * bucket_block_index + bucket_blocks_offset - 1;
	}

	/* FIXME: Should handle offset overflow  --- Ah screw it, just use 16-bit offset.  That'll never overflow.  */
	uint64_t runend_block_index  = bucket_block_index + bucket_blocks_offset / SLOTS_PER_BLOCK;
	uint64_t runend_ignore_bits  = bucket_blocks_offset % SLOTS_PER_BLOCK;
	uint64_t runend_rank         = bucket_intrablock_rank - 1;
	// uint64_t runend_block_offset = bitselectv(get_block(qf, runend_block_index)->runends, METADATA_WORDS_PER_BLOCK, runend_ignore_bits, runend_rank);
	uint64_t runend_block_offset = bitselectv(get_block(qf, runend_block_index)->runends[0], runend_ignore_bits, runend_rank);
	if (runend_block_offset == SLOTS_PER_BLOCK) {
		if (bucket_blocks_offset == 0 && bucket_intrablock_rank == 0) {
			/* The block begins in empty space, and this bucket is in that region of empty space */
			return hash_bucket_index;
		} else {
			do {
				// runend_rank        -= popcntv(get_block(qf, runend_block_index)->runends, METADATA_WORDS_PER_BLOCK, runend_ignore_bits);
				runend_rank        -= popcntv(get_block(qf, runend_block_index)->runends[0], runend_ignore_bits);
				runend_block_index++;
				runend_ignore_bits  = 0;
				// runend_block_offset = bitselectv(get_block(qf, runend_block_index)->runends, METADATA_WORDS_PER_BLOCK, runend_ignore_bits, runend_rank);
				runend_block_offset = bitselectv(get_block(qf, runend_block_index)->runends[0], runend_ignore_bits, runend_rank);
			} while (runend_block_offset == SLOTS_PER_BLOCK);
		}
	}

	uint64_t runend_index = SLOTS_PER_BLOCK * runend_block_index + runend_block_offset;
	if (runend_index < hash_bucket_index)
		return hash_bucket_index;
	else
		return runend_index;
}

/* find if the given slot is empty.
 * If the slot_offset is smaller than the boffset then it's never empty. */
static inline int offset_lower_bound(const QF *qf, uint64_t slot_index)
{
	const qfblock * b = get_block(qf, slot_index / SLOTS_PER_BLOCK);
	const uint64_t slot_offset = slot_index % SLOTS_PER_BLOCK;
	const uint64_t boffset = b->offset;
	const uint64_t occupieds = b->occupieds[0] & BITMASK(slot_offset+1);
	assert(SLOTS_PER_BLOCK == 64);
	if (boffset <= slot_offset) {
		const uint64_t runends = (b->runends[0] & BITMASK(slot_offset)) >> boffset;
		return popcnt(occupieds) - popcnt(runends);
	}
	return boffset - slot_offset + popcnt(occupieds);
}

static inline int is_empty(const QF *qf, uint64_t slot_index)
{
	return offset_lower_bound(qf, slot_index) == 0;
}

/* static inline int is_empty(const QF *qf, uint64_t slot_index) */
/* { */
/* 	return get_slot(qf, slot_index) == 0 */
/* 		&& !is_occupied(qf, slot_index)  */
/* 		&& !is_runend(qf, slot_index) */
/* 		&& run_end(qf, slot_index) == slot_index; */
/* } */

static inline int might_be_empty(const QF *qf, uint64_t slot_index)
{
	return !is_occupied(qf, slot_index)
		&& !is_runend(qf, slot_index);
}

static inline int probably_is_empty(const QF *qf, uint64_t slot_index)
{
	return get_slot(qf, slot_index) == 0
		&& !is_occupied(qf, slot_index)
		&& !is_runend(qf, slot_index);
}

static inline uint64_t find_first_empty_slot(QF *qf, uint64_t from)
{
	do {
		/* while (!probably_is_empty(qf, from)) */
		/*   from++; */
		int t = offset_lower_bound(qf, from);
		assert(t>=0);
		if (t == 0)
			break;
		from = from + t;
	} while(1);
	return from;
}

static inline uint64_t shift_into_b(const uint64_t a, const uint64_t b,
																		const int bstart, const int bend,
																		const int amount)
{
	const uint64_t a_component = bstart == 0 ? (a >> (64 - amount)) : 0;
	const uint64_t b_shifted_mask = BITMASK(bend - bstart) << bstart;
	const uint64_t b_shifted = ((b_shifted_mask & b) << amount) & b_shifted_mask;
	const uint64_t b_mask = ~b_shifted_mask;
	return a_component | b_shifted | (b & b_mask);
}

#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64

static inline void shift_remainders(QF *qf, uint64_t start_index, uint64_t empty_index)
{
	uint64_t start_block  = start_index / SLOTS_PER_BLOCK;
	uint64_t start_offset = start_index % SLOTS_PER_BLOCK;
	uint64_t empty_block  = empty_index / SLOTS_PER_BLOCK;
	uint64_t empty_offset = empty_index % SLOTS_PER_BLOCK;

//std::cout << start_index <<" "<< empty_index<<" "<< empty_index <<" "<< qf->xnslots <<std::endl;
	assert (start_index <= empty_index && empty_index < qf->xnslots);

	while (start_block < empty_block) {
		memmove(&get_block(qf, empty_block)->slots[1],
						&get_block(qf, empty_block)->slots[0],
						empty_offset * sizeof(qf->blocks[0].slots[0]));
		get_block(qf, empty_block)->slots[0] = get_block(qf, empty_block-1)->slots[SLOTS_PER_BLOCK-1];
		empty_block--;
		empty_offset = SLOTS_PER_BLOCK-1;
	}

	memmove(&get_block(qf, empty_block)->slots[start_offset+1],
					&get_block(qf, empty_block)->slots[start_offset],
					(empty_offset - start_offset) * sizeof(qf->blocks[0].slots[0]));
}

#else

#define REMAINDER_WORD(qf, i) ((uint64_t *)&(get_block(qf, (i)/qf->bits_per_slot)->slots[8 * ((i) % qf->bits_per_slot)]))

static inline void shift_remainders(QF *qf, const uint64_t start_index, const uint64_t empty_index)
{
	uint64_t last_word = (empty_index + 1) * qf->bits_per_slot / 64;
	const uint64_t first_word = start_index * qf->bits_per_slot / 64;
	int bend = ((empty_index + 1) * qf->bits_per_slot) % 64;
	const int bstart = (start_index * qf->bits_per_slot) % 64;

	while (last_word != first_word) {
		*REMAINDER_WORD(qf, last_word) = shift_into_b(*REMAINDER_WORD(qf, last_word-1),
																									*REMAINDER_WORD(qf, last_word),
																									0, bend, qf->bits_per_slot);
		last_word--;
		bend = 64;
	}
	*REMAINDER_WORD(qf, last_word) = shift_into_b(0, *REMAINDER_WORD(qf, last_word),
																								bstart, bend, qf->bits_per_slot);
}

#endif

static inline void qf_dump_block(const QF *qf, uint64_t i)
{
	uint64_t j;

	printf("%-192d", get_block(qf, i)->offset);
	printf("\n");

	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf("%02lx ", j);
	printf("\n");

	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf(" %d ", (get_block(qf, i)->occupieds[j/64] & (1ULL << (j%64))) ? 1 : 0);
	printf("\n");

	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf(" %d ", (get_block(qf, i)->runends[j/64] & (1ULL << (j%64))) ? 1 : 0);
	printf("\n");

#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32
	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf("%02x ", get_block(qf, i)->slots[j]);
#elif BITS_PER_SLOT == 64
	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf("%02lx ", get_block(qf, i)->slots[j]);
#else
	for (j = 0; j < SLOTS_PER_BLOCK * qf->bits_per_slot / 8; j++)
		printf("%02x ", get_block(qf, i)->slots[j]);
#endif

	printf("\n");

	printf("\n");
}

void qf_dump(const QF *qf)
{
	uint64_t i;

	printf("%lu %lu %lu\n",
				 qf->nblocks,
				 qf->ndistinct_elts,
				 qf->nelts);

	for (i = 0; i < qf->nblocks; i++) {
		qf_dump_block(qf, i);
	}

}

void qf_serialize(const QF *qf, const char *filename)
{
	uint64_t tmp_range;
	FILE *fout;
	fout = fopen(filename, "wb+");
	if (fout == NULL) {
		perror("Error opening file for serializing\n");
		exit(EXIT_FAILURE);
	}

	/* just a hack to handle __uint128_t value. Don't know a better to handle it
	 * right now */
	tmp_range = qf->range;

	fprintf(fout, "%lu ", qf->nslots);
	fprintf(fout, "%lu ", qf->xnslots);
	fprintf(fout, "%lu ", qf->key_bits);
	fprintf(fout, "%lu ", qf->value_bits);
	fprintf(fout, "%lu ", qf->key_remainder_bits);
	fprintf(fout, "%lu ", qf->bits_per_slot);
	fprintf(fout, "%lu ", tmp_range);
	fprintf(fout, "%lu ", qf->nblocks);
	fprintf(fout, "%lu ", qf->nelts);
	fprintf(fout, "%lu ", qf->ndistinct_elts);
	fprintf(fout, "%lu ", qf->noccupied_slots);

#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64
	assert(qf->nblocks == fwrite(qf->blocks, sizeof(qfblock), qf->nblocks, fout));
#else
	assert(qf->nblocks == fwrite(qf->blocks, sizeof(qfblock) + SLOTS_PER_BLOCK * qf->bits_per_slot / 8, qf->nblocks, fout));
#endif

	fclose(fout);
}

void qf_deserialize(QF *qf, const char *filename)
{
	uint64_t tmp_range;
	FILE *fin;
	fin = fopen(filename, "rb");
	if (fin == NULL) {
		perror("Error opening file for deserializing\n");
		exit(EXIT_FAILURE);
	}

	fscanf(fin, "%lu ", &qf->nslots);
	fscanf(fin, "%lu ", &qf->xnslots);
	fscanf(fin, "%lu ", &qf->key_bits);
	fscanf(fin, "%lu ", &qf->value_bits);
	fscanf(fin, "%lu ", &qf->key_remainder_bits);
	fscanf(fin, "%lu ", &qf->bits_per_slot);
	fscanf(fin, "%lu ", &tmp_range);
	fscanf(fin, "%lu ", &qf->nblocks);
	fscanf(fin, "%lu ", &qf->nelts);
	fscanf(fin, "%lu ", &qf->ndistinct_elts);
	fscanf(fin, "%lu ", &qf->noccupied_slots);

	/* just a hack to handle __uint128_t value. Don't know a better to handle it
	 * right now */
	qf->range = tmp_range;

	/* allocate the space for the actual qf blocks */
#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64
	qf->blocks = (qfblock *)calloc(qf->nblocks, sizeof(qfblock));
#else
	qf->blocks = (qfblock *)calloc(qf->nblocks, sizeof(qfblock) + SLOTS_PER_BLOCK * qf->bits_per_slot / 8);
#endif

#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64
	assert(qf->nblocks == fread(qf->blocks, sizeof(qfblock), qf->nblocks, fin));
#else
	assert(qf->nblocks == fread(qf->blocks, sizeof(qfblock) + SLOTS_PER_BLOCK * qf->bits_per_slot / 8, qf->nblocks, fin));
#endif

	fclose(fin);
}

static inline void find_next_n_empty_slots(QF *qf, uint64_t from, uint64_t n, uint64_t *indices)
{
	while (n) {
		indices[--n] = find_first_empty_slot(qf, from);
		from = indices[n] + 1;
	}
}

static inline void shift_slots(QF *qf, int64_t first, uint64_t last, uint64_t distance)
{
	int64_t i;
	if (distance == 1)
		shift_remainders(qf, first, last+1);
	else
		for (i = last; i >= first; i--)
			set_slot(qf, i + distance, get_slot(qf, i));
}

static inline void shift_runends(QF *qf, int64_t first, uint64_t last, uint64_t distance)
{
	assert(last < qf->xnslots && distance < 64);
	uint64_t first_word = first / 64;
	uint64_t bstart = first % 64;
	uint64_t last_word = (last + distance + 1) / 64;
	uint64_t bend = (last + distance + 1) % 64;

	if (last_word != first_word) {
		METADATA_WORD(qf, runends, 64*last_word) = shift_into_b(METADATA_WORD(qf, runends, 64*(last_word-1)),
																														METADATA_WORD(qf, runends, 64*last_word),
																														0, bend, distance);
		bend = 64;
		last_word--;
		while (last_word != first_word) {
			METADATA_WORD(qf, runends, 64*last_word) = shift_into_b(METADATA_WORD(qf, runends, 64*(last_word-1)),
																															METADATA_WORD(qf, runends, 64*last_word),
																															0, bend, distance);
			last_word--;
		}
	}
	METADATA_WORD(qf, runends, 64*last_word) = shift_into_b(0, METADATA_WORD(qf, runends, 64*last_word),
																													bstart, bend, distance);

}

static inline void insert_replace_slots_and_shift_remainders_and_runends_and_offsets(QF		*qf,
																																										 int		 operation,
																																										 uint64_t		 bucket_index,
																																										 uint64_t		 overwrite_index,
																																										 const uint64_t	*remainders,
																																										 uint64_t		 total_remainders,
																																										 uint64_t		 noverwrites)
{
	uint64_t empties[67];
	uint64_t i, j;
	uint64_t ninserts = total_remainders - noverwrites;
	uint64_t insert_index = overwrite_index + noverwrites;

	if (ninserts > 0) {
		/* First, shift things to create n empty spaces where we need them. */
		find_next_n_empty_slots(qf, insert_index, ninserts, empties);

		for (i = 0; i < ninserts - 1; i++)
			shift_slots(qf, empties[i+1] + 1, empties[i] - 1, i + 1);
		if (ninserts > 0)
			shift_slots(qf, insert_index, empties[ninserts - 1] - 1, ninserts);

		for (i = 0; i < ninserts - 1; i++)
			shift_runends(qf, empties[i+1] + 1, empties[i] - 1, i + 1);
		if (ninserts > 0)
			shift_runends(qf, insert_index, empties[ninserts - 1] - 1, ninserts);


		for (i = noverwrites; i < total_remainders - 1; i++)
			METADATA_WORD(qf, runends, overwrite_index + i) &= ~(1ULL << (((overwrite_index + i) % SLOTS_PER_BLOCK) % 64));

		switch (operation) {
			case 0: /* insert into empty bucket */
				assert (noverwrites == 0);
				METADATA_WORD(qf, runends, overwrite_index + total_remainders - 1) |= 1ULL << (((overwrite_index + total_remainders - 1) % SLOTS_PER_BLOCK) % 64);
				break;
			case 1: /* append to bucket */
				METADATA_WORD(qf, runends, overwrite_index + noverwrites - 1)      &= ~(1ULL << (((overwrite_index + noverwrites - 1) % SLOTS_PER_BLOCK) % 64));
				METADATA_WORD(qf, runends, overwrite_index + total_remainders - 1) |= 1ULL << (((overwrite_index + total_remainders - 1) % SLOTS_PER_BLOCK) % 64);
				break;
			case 2: /* insert into bucket */
				METADATA_WORD(qf, runends, overwrite_index + total_remainders - 1) &= ~(1ULL << (((overwrite_index + total_remainders - 1) % SLOTS_PER_BLOCK) % 64));
				break;
			default:
				fprintf(stderr, "Invalid operation %d\n", operation);
				abort();
		}

		if (ninserts > 0) {
			for (i = bucket_index / SLOTS_PER_BLOCK + 1; i <= empties[ninserts - 1]/SLOTS_PER_BLOCK; i++) {
				if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
					get_block(qf, i)->offset += ninserts;
			}
			for (j = 0; j < ninserts -1 ; j++)
				for (i = empties[ninserts - j - 1] / SLOTS_PER_BLOCK + 1; i <= empties[ninserts - j - 2]/SLOTS_PER_BLOCK; i++) {
					if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
						get_block(qf, i)->offset += ninserts - j - 1;
				}
		}
	}
	for (i = 0; i < total_remainders; i++)
		set_slot(qf, overwrite_index + i, remainders[i]);

	qf->noccupied_slots += ninserts;
}

static inline void remove_replace_slots_and_shift_remainders_and_runends_and_offsets(QF		        *qf,
																																										 int		 operation,
																																										 uint64_t		 bucket_index,
																																										 uint64_t		 overwrite_index,
																																										 const uint64_t	*remainders,
																																										 uint64_t		 total_remainders,
																																										 uint64_t		 old_length)
{
	uint64_t i;

	// Update the slots
	for (i = 0; i < total_remainders; i++)
		set_slot(qf, overwrite_index + i, remainders[i]);

	// If this is the last thing in its run, then we may need to set a new runend bit
	if (is_runend(qf, overwrite_index + old_length - 1)) {
	  if (total_remainders > 0) {
	    // If we're not deleting this entry entirely, then it will still the last entry in this run
	    METADATA_WORD(qf, runends, overwrite_index + total_remainders - 1) |= 1ULL << ((overwrite_index + total_remainders - 1) % 64);
	  } else if (overwrite_index > bucket_index &&
		     !is_runend(qf, overwrite_index - 1)) {
	    // If we're deleting this entry entirely, but it is not the first entry in this run,
	    // then set the preceding entry to be the runend
	    METADATA_WORD(qf, runends, overwrite_index - 1) |= 1ULL << ((overwrite_index - 1) % 64);
	  }
	}

	// shift slots back one run at a time
	uint64_t original_bucket = bucket_index;
	uint64_t current_bucket = bucket_index;
	uint64_t current_slot = overwrite_index + total_remainders;
	uint64_t current_distance = old_length - total_remainders;

	while (current_distance > 0) {
		if (is_runend(qf, current_slot + current_distance - 1)) {
			do {
				current_bucket++;
			} while (current_bucket < current_slot + current_distance &&
							 !is_occupied(qf, current_bucket));
		}

		if (current_bucket <= current_slot) {
			set_slot(qf, current_slot, get_slot(qf, current_slot + current_distance));
			if (is_runend(qf, current_slot) !=
					is_runend(qf, current_slot + current_distance))
				METADATA_WORD(qf, runends, current_slot) ^= 1ULL << (current_slot % 64);
			current_slot++;

		} else if (current_bucket <= current_slot + current_distance) {
			uint64_t i;
			for (i = current_slot; i < current_slot + current_distance; i++) {
				set_slot(qf, i, 0);
				METADATA_WORD(qf, runends, i) &= ~(1ULL << (i % 64));
			}

			current_distance = current_slot + current_distance - current_bucket;
			current_slot = current_bucket;
		} else {
			current_distance = 0;
		}
	}

	// reset the occupied bit of the hash bucket index if the hash is the
	// only item in the run and is removed completely.
	if (operation && !total_remainders)
		METADATA_WORD(qf, occupieds, bucket_index) &= ~(1ULL << (bucket_index % 64));

	// update the offset bits.
	// find the number of occupied slots in the original_bucket block.
	// Then find the runend slot corresponding to the last run in the
	// original_bucket block.
	// Update the offset of the block to which it belongs.
	uint64_t original_block = original_bucket / SLOTS_PER_BLOCK;
	while (1 && old_length > total_remainders) {	// we only update offsets if we shift/delete anything
		int32_t last_occupieds_bit = bitscanreverse(get_block(qf, original_block)->occupieds[0]);
		// there is nothing in the block
		if (last_occupieds_bit == -1) {
			if (get_block(qf, original_block + 1)->offset == 0)
				break;
			get_block(qf, original_block + 1)->offset = 0;
		} else {
			uint64_t last_occupieds_hash_index = SLOTS_PER_BLOCK * original_block + last_occupieds_bit;
			uint64_t runend_index = run_end(qf, last_occupieds_hash_index);
			// runend spans across the block
			// update the offset of the next block
			if (runend_index / SLOTS_PER_BLOCK == original_block) { // if the run ends in the same block
				if (get_block(qf, original_block + 1)->offset == 0)
					break;
				get_block(qf, original_block + 1)->offset = 0;
			} else if (runend_index / SLOTS_PER_BLOCK == original_block + 1) { // if the last run spans across one block
				if (get_block(qf, original_block + 1)->offset == (runend_index % SLOTS_PER_BLOCK) + 1)
					break;
				get_block(qf, original_block + 1)->offset = (runend_index % SLOTS_PER_BLOCK) + 1;
			} else { // if the last run spans across multiple blocks
				uint64_t i;
				for (i = original_block + 1; i < runend_index / SLOTS_PER_BLOCK - 1; i++)
					get_block(qf, i)->offset = SLOTS_PER_BLOCK;
				if (get_block(qf, runend_index / SLOTS_PER_BLOCK)->offset == (runend_index % SLOTS_PER_BLOCK) + 1)
					break;
				get_block(qf, runend_index / SLOTS_PER_BLOCK)->offset = (runend_index % SLOTS_PER_BLOCK) + 1;
			}
		}
		original_block++;
	}

	qf->noccupied_slots -= (old_length - total_remainders);
	if (!total_remainders) {
		qf->ndistinct_elts--;
	}
}

/*****************************************************************************
 * Code that uses the above to implement a QF with keys and inline counters. *
 *****************************************************************************/

/*
	 Counter format:
	 0 xs:    <empty string>
	 1 x:     x
	 2 xs:    xx
	 3 0s:    000
	 >2 xs:   xbc...cx  for x != 0, b < x, c != 0, x
	 >3 0s:   0c...c00  for c != 0
	 */
static inline uint64_t *encode_counter(QF *qf, uint64_t remainder, uint64_t counter, uint64_t *slots)
{
	uint64_t digit = remainder;
	uint64_t base = (1ULL << qf->bits_per_slot) - 1;
	uint64_t *p = slots;

	if (counter == 0)
		return p;

	*--p = remainder;

	if (counter == 1)
		return p;

	if (counter == 2) {
		*--p = remainder;
		return p;
	}

	if (counter == 3 && remainder == 0) {
		*--p = remainder;
		*--p = remainder;
		return p;
	}

	if (counter == 3 && remainder > 0) {
		*--p = 0;
		*--p = remainder;
		return p;
	}

	if (remainder == 0)
		*--p = remainder;
	else
		base--;

	if (remainder)
		counter -= 3;
	else
		counter -= 4;
	do {
		digit = counter % base;
		digit++; /* Zero not allowed */
		if (remainder && digit >= remainder)
			digit++; /* Cannot overflow since digit is mod 2^r-2 */
		*--p = digit;
		counter /= base;
	} while (counter);

	if (remainder && digit >= remainder)
		*--p = 0;

	*--p = remainder;

	return p;
}

/* Returns the index to the last element in the counter endoding.
REQUIRES: index points to first slot of a counter. */
static inline uint64_t decode_counter(const QF *qf, uint64_t index, uint64_t *remainder, uint64_t *count)
{
	uint64_t base;
	uint64_t rem;
	uint64_t cnt;
	uint64_t digit;
	uint64_t end;

	*remainder = rem = get_slot(qf, index);

	if (is_runend(qf, index)) { /* Entire run is "0" */
		*count = 1;
		return index;
	}

	digit = get_slot(qf, index + 1);

	if (is_runend(qf, index + 1)) {
		*count = digit == rem ? 2 : 1;
		return index + (digit == rem ? 1 : 0);
	}

	if (rem > 0 && digit >= rem) {
		*count = digit == rem ? 2 : 1;
		return index + (digit == rem ? 1 : 0);
	}

	if (rem > 0 && digit == 0 && get_slot(qf, index + 2) == rem) {
		*count = 3;
		return index + 2;
	}

	if (rem == 0 && digit == 0) {
		if (get_slot(qf, index + 2) == 0) {
			*count = 3;
			return index + 2;
		} else {
			*count = 2;
			return index + 1;
		}
	}

	cnt = 0;
	base = (1ULL << qf->bits_per_slot) - (rem ? 2 : 1);

	end = index + 1;
	while (digit != rem && !is_runend(qf, end)) {
		if (digit > rem)
			digit--;
		if (digit && rem)
			digit--;
		cnt = cnt * base + digit;

		end++;
		digit = get_slot(qf, end);
	}

	if (rem) {
		*count = cnt + 3;
		return end;
	}

	if (is_runend(qf, end) || get_slot(qf, end + 1) != 0) {
		*count = 1;
		return index;
	}

	*count = cnt + 4;
	return end + 1;
}

/* return the next slot which corresponds to a
 * different element
 * */
static inline uint64_t next_slot(QF *qf, uint64_t current)
{
	uint64_t rem = get_slot(qf, current);
	current++;

	while (get_slot(qf, current) == rem && current <= qf->nslots) {
		current++;
	}
	return current;
}

static inline void insert1(QF *qf, __uint128_t hash)
{
	uint64_t hash_remainder           = hash & BITMASK(qf->bits_per_slot);
	uint64_t hash_bucket_index        = hash >> qf->bits_per_slot;
	uint64_t hash_bucket_block_offset = hash_bucket_index % SLOTS_PER_BLOCK;

	if (is_empty(qf, hash_bucket_index)) {
		METADATA_WORD(qf, runends, hash_bucket_index) |= 1ULL << (hash_bucket_block_offset % 64);
		set_slot(qf, hash_bucket_index, hash_remainder);
		qf->noccupied_slots++;
		qf->ndistinct_elts++;
#ifdef LOG_NUM_SHIFTS
		shift_count[0]++;
#endif
	} else {
		uint64_t runend_index              = run_end(qf, hash_bucket_index);
		int operation = 0; /* Insert into empty bucket */
		uint64_t insert_index = runend_index + 1;
		uint64_t new_value = hash_remainder;

		/* printf("RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

		uint64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf, hash_bucket_index - 1) + 1;

		if (is_occupied(qf, hash_bucket_index)) {

			/* Find the counter for this remainder if it exists. */
			uint64_t current_remainder = get_slot(qf, runstart_index);
			uint64_t zero_terminator = runstart_index;

			/* The counter for 0 is special. */
			if (current_remainder == 0) {
				uint64_t t = runstart_index + 1;
				while (t < runend_index && get_slot(qf, t) != 0)
					t++;
				if (t < runend_index && get_slot(qf, t+1) == 0)
					zero_terminator = t+1; /* Three or more 0s */
				else if (runstart_index < runend_index && get_slot(qf, runstart_index + 1) == 0)
					zero_terminator = runstart_index + 1; /* Exactly two 0s */
				/* Otherwise, exactly one 0 (i.e. zero_terminator == runstart_index) */

				/* May read past end of run, but that's OK because loop below
					 can handle that */
				if (hash_remainder != 0) {
					runstart_index = zero_terminator + 1;
					current_remainder = get_slot(qf, runstart_index);
				}
			}

			/* Skip over counters for other remainders. */
			while (current_remainder < hash_remainder && runstart_index <= runend_index) {
				/* If this remainder has an extended counter, skip over it. */
				if (runstart_index < runend_index &&
						get_slot(qf, runstart_index + 1) < current_remainder) {
					runstart_index = runstart_index + 2;
					while (get_slot(qf, runstart_index) != current_remainder)
						runstart_index++;
					runstart_index++;

					/* This remainder has a simple counter. */
				} else {
					runstart_index++;
				}

				/* This may read past the end of the run, but the while loop
					 condition will prevent us from using the invalid result in
					 that case. */
				current_remainder = get_slot(qf, runstart_index);
			}

			/* If this is the first time we've inserted the new remainder,
				 and it is larger than any remainder in the run. */
			if (runstart_index > runend_index) {
				operation = 1;
				insert_index = runstart_index;
				new_value = hash_remainder;
				qf->ndistinct_elts++;

				/* This is the first time we're inserting this remainder, but
					 there are larger remainders already in the run. */
			} else if (current_remainder != hash_remainder) {
				operation = 2; /* Inserting */
				insert_index = runstart_index;
				new_value = hash_remainder;
				qf->ndistinct_elts++;

				/* Cases below here: we're incrementing the (simple or
					 extended) counter for this remainder. */

				/* If there's exactly one instance of this remainder. */
			} else if (runstart_index == runend_index ||
								 (hash_remainder > 0 && get_slot(qf, runstart_index + 1) > hash_remainder) ||
								 (hash_remainder == 0 && zero_terminator == runstart_index)) {
				operation = 2; /* Insert */
				insert_index = runstart_index;
				new_value = hash_remainder;

				/* If there are exactly two instances of this remainder. */
			} else if ((hash_remainder > 0 && get_slot(qf, runstart_index + 1) == hash_remainder) ||
								 (hash_remainder == 0 && zero_terminator == runstart_index + 1)) {
				operation = 2; /* Insert */
				insert_index = runstart_index + 1;
				new_value = 0;

				/* Special case for three 0s */
			} else if (hash_remainder == 0 && zero_terminator == runstart_index + 2) {
				operation = 2; /* Insert */
				insert_index = runstart_index + 1;
				new_value = 1;

				/* There is an extended counter for this remainder. */
			} else {

				/* Move to the LSD of the counter. */
				insert_index = runstart_index + 1;
				while (get_slot(qf, insert_index+1) != hash_remainder)
					insert_index++;

				/* Increment the counter. */
				uint64_t digit, carry;
				do {
					carry = 0;
					digit = get_slot(qf, insert_index);
					// Convert a leading 0 (which is special) to a normal encoded digit
					if (digit == 0) {
						digit++;
						if (digit == current_remainder)
							digit++;
					}

					// Increment the digit
					digit = (digit + 1) & BITMASK(qf->bits_per_slot);

					// Ensure digit meets our encoding requirements
					if (digit == 0) {
						digit++;
						carry = 1;
					}
					if (digit == current_remainder)
						digit = (digit + 1) & BITMASK(qf->bits_per_slot);
					if (digit == 0) {
						digit++;
						carry = 1;
					}

					set_slot(qf, insert_index, digit);
					insert_index--;
				} while(insert_index > runstart_index && carry);

				/* If the counter needs to be expanded. */
				if (insert_index == runstart_index && (carry > 0 || (current_remainder != 0 && digit >= current_remainder))) {
					operation = 2; /* insert */
					insert_index = runstart_index + 1;
					if (!carry)						/* To prepend a 0 before the counter if the MSD is greater than the rem */
						new_value = 0;
					else if (carry) {			/* Increment the new value because we don't use 0 to encode counters */
						new_value = 2;
						/* If the rem is greater than or equal to the new_value then fail*/
						if (current_remainder)
							assert(new_value < current_remainder);
					}
				} else {
					operation = -1;
				}
			}
		}

		if (operation >= 0) {
			uint64_t empty_slot_index = find_first_empty_slot(qf, runend_index+1);

#ifdef LOG_NUM_SHIFTS
			if (empty_slot_index-insert_index < len-1)
				shift_count[empty_slot_index-insert_index]++;
			else
				shift_count[len-1]++;
#endif
			shift_remainders(qf, insert_index, empty_slot_index);
			set_slot(qf, insert_index, new_value);
			shift_runends(qf, insert_index, empty_slot_index-1, 1);
			switch (operation) {
				case 0:
					METADATA_WORD(qf, runends, insert_index)   |= 1ULL << ((insert_index % SLOTS_PER_BLOCK) % 64);
					break;
				case 1:
					METADATA_WORD(qf, runends, insert_index-1) &= ~(1ULL << (((insert_index-1) % SLOTS_PER_BLOCK) % 64));
					METADATA_WORD(qf, runends, insert_index)   |= 1ULL << ((insert_index % SLOTS_PER_BLOCK) % 64);
					break;
				case 2:
					METADATA_WORD(qf, runends, insert_index)   &= ~(1ULL << ((insert_index % SLOTS_PER_BLOCK) % 64));
					break;
				default:
					fprintf(stderr, "Invalid operation %d\n", operation);
					abort();
			}

			/*
			 * Increment the offset for each block between the hash bucket index
			 * and block of the empty slot
			 * */
			uint64_t i;
			for (i = hash_bucket_index / SLOTS_PER_BLOCK + 1; i <= empty_slot_index/SLOTS_PER_BLOCK; i++) {
				if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
					get_block(qf, i)->offset++;
				//assert(get_block(qf, i)->offset != 0);
			}

			qf->noccupied_slots++;
		}
	}

	METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL << (hash_bucket_block_offset % 64);
	qf->nelts++;
}

static inline void insert(QF *qf, __uint128_t hash, uint64_t count)
{
	uint64_t hash_remainder           = hash & BITMASK(qf->bits_per_slot);
	uint64_t hash_bucket_index        = hash >> qf->bits_per_slot;
	uint64_t hash_bucket_block_offset = hash_bucket_index % SLOTS_PER_BLOCK;

	/* Empty slot */
	if (is_empty(qf, hash_bucket_index)) {
		METADATA_WORD(qf, runends, hash_bucket_index) |= 1ULL << (hash_bucket_block_offset % 64);
		set_slot(qf, hash_bucket_index, hash_remainder);
		qf->noccupied_slots++;

		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL << (hash_bucket_block_offset % 64);
		qf->nelts += 1;
		qf->ndistinct_elts++;

		/* This trick will, I hope, keep the fast case fast. */
		if (count > 1) {
			insert(qf, hash, count - 1);
		}

	} else { /* Non-empty slot */
		uint64_t new_values[67];
		int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf, hash_bucket_index - 1) + 1;

		if (!is_occupied(qf, hash_bucket_index)) { /* Empty bucket, but its slot is occupied. */
			uint64_t *p = encode_counter(qf, hash_remainder, count, &new_values[67]);
			insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf,
																																				0,
																																				hash_bucket_index,
																																				runstart_index,
																																				p,
																																				&new_values[67] - p,
																																				0);
			qf->ndistinct_elts++;

		} else { /* Non-empty bucket */

			uint64_t current_remainder, current_count, current_end;

			/* Find the counter for this remainder, if one exists. */
			current_end = decode_counter(qf, runstart_index, &current_remainder, &current_count);
			while (current_remainder < hash_remainder && !is_runend(qf, current_end)) {
				runstart_index = current_end + 1;
				current_end = decode_counter(qf, runstart_index, &current_remainder, &current_count);
			}

			/* If we reached the end of the run w/o finding a counter for this remainder,
				 then append a counter for this remainder to the run. */
			if (current_remainder < hash_remainder) {
				uint64_t *p = encode_counter(qf, hash_remainder, count, &new_values[67]);
				insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf,
																																					1, /* Append to bucket */
																																					hash_bucket_index,
																																					current_end + 1,
																																					p,
																																					&new_values[67] - p,
																																					0);
				qf->ndistinct_elts++;

				/* Found a counter for this remainder.  Add in the new count. */
			} else if (current_remainder == hash_remainder) {
				uint64_t *p = encode_counter(qf, hash_remainder, current_count + count, &new_values[67]);
				insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf,
																																					is_runend(qf, current_end) ? 1 : 2,
																																					hash_bucket_index,
																																					runstart_index,
																																					p,
																																					&new_values[67] - p,
																																					current_end - runstart_index + 1);

				/* No counter for this remainder, but there are larger
					 remainders, so we're not appending to the bucket. */
			} else {
				uint64_t *p = encode_counter(qf, hash_remainder, count, &new_values[67]);
				insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf,
																																					2, /* Insert to bucket */
																																					hash_bucket_index,
																																					runstart_index,
																																					p,
																																					&new_values[67] - p,
																																					0);
				qf->ndistinct_elts++;
			}
		}

		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL << (hash_bucket_block_offset % 64);
		qf->nelts += count;
	}
}

inline static void _remove(QF *qf, __uint128_t hash, uint64_t count)
{
	uint64_t hash_remainder           = hash & BITMASK(qf->bits_per_slot);
	uint64_t hash_bucket_index        = hash >> qf->bits_per_slot;
	uint64_t current_remainder, current_count, current_end;
	uint64_t new_values[67];

	/* Empty bucket */
	if (!is_occupied(qf, hash_bucket_index))
		return;

	uint64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf, hash_bucket_index - 1) + 1;
	uint64_t original_runstart_index = runstart_index;
	int only_item_in_the_run = 0;

	/*Find the counter for this remainder, if one exists.*/
	current_end = decode_counter(qf, runstart_index, &current_remainder, &current_count);
	while (current_remainder < hash_remainder && !is_runend(qf, current_end)) {
		runstart_index = current_end + 1;
		current_end = decode_counter(qf, runstart_index, &current_remainder, &current_count);
	}
	/* remainder not found in the given run */
	if (current_remainder != hash_remainder)
		return;

	if (original_runstart_index == runstart_index && is_runend(qf, current_end))
		only_item_in_the_run = 1;

	/* endode the new counter */
	uint64_t *p = encode_counter(qf, hash_remainder,
															 count > current_count ? 0 : current_count - count,
															 &new_values[67]);
	remove_replace_slots_and_shift_remainders_and_runends_and_offsets(qf,
																																		only_item_in_the_run,
																																		hash_bucket_index,
																																		runstart_index,
																																		p,
																																		&new_values[67] - p,
																																		current_end - runstart_index + 1);

	// update the nelements.
	qf->nelts -= count;
}

void qf_destroy(QF *qf)
{
	assert(qf->blocks != NULL);
	free(qf->blocks);
}

/***********************************************************************
 * Code that uses the above to implement key-value-counter operations. *
 ***********************************************************************/

void qf_init(QF *qf, uint64_t nslots, uint64_t key_bits, uint64_t value_bits)
{

	assert(popcnt(nslots) == 1); /* nslots must be a power of 2 */
	qf->nslots = nslots;
	qf->xnslots = nslots + 10*sqrt((double)nslots);
	qf->key_bits = key_bits;
	qf->value_bits = value_bits;

	qf->key_remainder_bits = key_bits;
	while (nslots > 1) {
		assert(qf->key_remainder_bits > 0);
		qf->key_remainder_bits--;
		nslots >>= 1;
	}

	qf->bits_per_slot = qf->key_remainder_bits + qf->value_bits;
	assert (BITS_PER_SLOT == 0 || BITS_PER_SLOT == qf->bits_per_slot);
	assert(qf->bits_per_slot > 1);

	qf->range = qf->nslots;
	qf->range <<= qf->bits_per_slot;
	qf->nblocks = (qf->xnslots + SLOTS_PER_BLOCK - 1) / SLOTS_PER_BLOCK;
	qf->nelts = 0;
	qf->ndistinct_elts = 0;
	qf->noccupied_slots = 0;
#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64
	qf->blocks = (qfblock *)calloc(qf->nblocks, sizeof(qfblock));
#else
	qf->blocks = (qfblock *)calloc(qf->nblocks, (sizeof(qfblock) + SLOTS_PER_BLOCK * qf->bits_per_slot / 8));
#endif
}

uint64_t qf_count_key_value(const QF *qf, uint64_t key, uint64_t value)
{
	__uint128_t hash = key;
	uint64_t hash_remainder   = hash & BITMASK(qf->bits_per_slot);
	int64_t hash_bucket_index = hash >> qf->bits_per_slot;

	if (!is_occupied(qf, hash_bucket_index))
		return 0;

	int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf, hash_bucket_index-1) + 1;
	if (runstart_index < hash_bucket_index)
		runstart_index = hash_bucket_index;

	/* printf("MC RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

	uint64_t current_remainder, current_count, current_end;
	do {
		current_end = decode_counter(qf, runstart_index, &current_remainder, &current_count);
		if (current_remainder == hash_remainder)
			return current_count;
		runstart_index = current_end + 1;
	} while (!is_runend(qf, current_end));

	return 0;
}

void qf_insert(QF *qf, uint64_t key, uint64_t value, uint64_t count)
{
	/*uint64_t hash = (key << qf->value_bits) | (value & BITMASK(qf->value_bits));*/
	if (count == 1)
		insert1(qf, key);
	else
		insert(qf, key, count);
}

/* count = 0 would remove the key completely. */
void qf_remove(QF *qf, uint64_t key, uint64_t value, uint64_t count)
{
	if (count == 0)
		count = qf_count_key_value(qf, key, value);

	_remove(qf, key, count);
}

/* initialize the iterator at the run corresponding
 * to the position index
 */
void qf_iterator(const QF *qf, QFi *qfi, uint64_t position)
{
	assert(position < qf->nslots);
	if (!is_occupied(qf, position)) {
		uint64_t block_index = position;
		uint64_t idx = bitselect(get_block(qf, block_index)->occupieds[0], 0);
		if (idx == 64) {
			while(idx == 64 && block_index < qf->nblocks) {
				block_index++;
				idx = bitselect(get_block(qf, block_index)->occupieds[0], 0);
			}
		}
		position = block_index * SLOTS_PER_BLOCK + idx;
	}

	qfi->qf = qf;
	qfi->run = position;
	qfi->current = position == 0 ? 0 : run_end(qfi->qf, position-1) + 1;
	if (qfi->current < position)
		qfi->current = position;
}

int qfi_get(QFi *qfi, uint64_t *key, uint64_t *value, uint64_t *count)
{
	assert(qfi->run <= qfi->qf->nslots);
	assert(qfi->current <= qfi->qf->xnslots);

	uint64_t current_remainder, current_count;
	decode_counter(qfi->qf, qfi->current, &current_remainder, &current_count);

	*key = (qfi->run << qfi->qf->bits_per_slot) | current_remainder;
	*value = 0;   // for now we are not using value
	*count = current_count;

	/*qfi->current = end_index;*/ 		//get should not change the current index
																		//of the iterator
	return 0;
}

int qfi_next(QFi *qfi)
{
	if (qfi_end(qfi))
		return 1;
	else {
		/* move to the end of the current counter*/
		uint64_t current_remainder, current_count;
		qfi->current = decode_counter(qfi->qf, qfi->current, &current_remainder, &current_count);

		if (!is_runend(qfi->qf, qfi->current)) {
			qfi->current++;
			if (qfi->current > qfi->qf->nslots)
				return 1;
			return 0;
		}
		else {
			uint64_t block_index = qfi->run / SLOTS_PER_BLOCK;
			uint64_t rank = bitrank(get_block(qfi->qf, block_index)->occupieds[0], qfi->run % SLOTS_PER_BLOCK);
			uint64_t next_run = bitselect(get_block(qfi->qf, block_index)->occupieds[0], rank);
			if (next_run == 64) {
				rank = 0;
				while (next_run == 64 && block_index < qfi->qf->nblocks) {
					block_index++;
					next_run = bitselect(get_block(qfi->qf, block_index)->occupieds[0], rank);
				}
			}
			if (block_index == qfi->qf->nblocks) {
				/* set the index values to max*/
				qfi->run = qfi->current = qfi->qf->xnslots;
				return 1;
			}
			qfi->run = block_index * SLOTS_PER_BLOCK + next_run;
			qfi->current++;
			if (qfi->current < qfi->run)
				qfi->current = qfi->run;
			return 0;
		}
	}
}

inline int qfi_end(QFi *qfi)
{
	if (qfi->current >= qfi->qf->xnslots /*&& is_runend(qfi->qf, qfi->current)*/)
		return 1;
	else
		return 0;
}

/*
 * Merge qfa and qfb into qfc
 */
/*
 * iterate over both qf (qfa and qfb)
 * simultaneously
 * for each index i
 * min(get_value(qfa, ia) < get_value(qfb, ib))
 * insert(min, ic)
 * increment either ia or ib, whichever is minimum.
 */
void qf_merge(const QF *qfa, const QF *qfb, QF *qfc)
{
	QFi qfia, qfib;
	qf_iterator(qfa, &qfia, 0);
	qf_iterator(qfb, &qfib, 0);

	uint64_t keya, valuea, counta, keyb, valueb, countb;
	qfi_get(&qfia, &keya, &valuea, &counta);
	qfi_get(&qfib, &keyb, &valueb, &countb);
	do {
		if (keya < keyb) {
			qf_insert(qfc, keya, valuea, counta);
			if (!qfi_next(&qfia))
				qfi_get(&qfia, &keya, &valuea, &counta);
		}
		else {
			qf_insert(qfc, keyb, valueb, countb);
			if (!qfi_next(&qfib))
				qfi_get(&qfib, &keyb, &valueb, &countb);
		}
	} while(!qfi_end(&qfia) && !qfi_end(&qfib));

	if (!qfi_end(&qfia)) {
		do {
			qfi_get(&qfia, &keya, &valuea, &counta);
			qf_insert(qfc, keya, valuea, counta);
		} while(!qfi_next(&qfia));
	}
	if (!qfi_end(&qfib)) {
		do {
			qfi_get(&qfib, &keyb, &valueb, &countb);
			qf_insert(qfc, keyb, valueb, countb);
		} while(!qfi_next(&qfib));
	}

	return;
}

/*
 * Merge an array of qfs into the resultant QF
 */
void qf_multi_merge(QF *qf_arr[], int nqf, QF *qfr)
{
	int i;
	QFi qfi_arr[nqf];
	int flag = 0;
	int smallest_i = 0;
	uint64_t smallest_key = UINT64_MAX;
	for (i=0; i<nqf; i++) {
		qf_iterator(qf_arr[i], &qfi_arr[i], 0);
	}

	while (!flag) {
		uint64_t keys[nqf];
		uint64_t values[nqf];
		uint64_t counts[nqf];
		for (i=0; i<nqf; i++)
			qfi_get(&qfi_arr[i], &keys[i], &values[i], &counts[i]);

		do {
			smallest_key = UINT64_MAX;
			for (i=0; i<nqf; i++) {
				if (keys[i] < smallest_key) {
					smallest_key = keys[i]; smallest_i = i;
				}
			}
			qf_insert(qfr, keys[smallest_i], values[smallest_i], counts[smallest_i]);
			if (!qfi_next(&qfi_arr[smallest_i]))
				qfi_get(&qfi_arr[smallest_i], &keys[smallest_i], &values[smallest_i], &counts[smallest_i]);
		} while(!qfi_end(&qfi_arr[smallest_i]));

		/* remove the qf that is exhausted from the array */
		if (smallest_i < nqf-1)
			memmove(&qfi_arr[smallest_i], &qfi_arr[smallest_i+1], (nqf-smallest_i-1)*sizeof(qfi_arr[0]));
		nqf--;
		if (nqf == 1)
			flag = 1;
	}
	if (!qfi_end(&qfi_arr[0])) {
		do {
			uint64_t key, value, count;
			qfi_get(&qfi_arr[0], &key, &value, &count);
			qf_insert(qfr, key, value, count);
		} while(!qfi_next(&qfi_arr[0]));
	}

	return;
}
