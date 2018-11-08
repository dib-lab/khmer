 #include <stdlib.h>
#if 0
# include <assert.h>
#else
# define assert(x)
#endif
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include "gqf.h"
#include "onDiskMQF.h"
#include <iostream>
#include <map>
#include "utils.h"
#include <type_traits>

namespace onDiskMQF_Namespace{
/******************************************************************
 * Code for managing the metadata bits and slots w/o interpreting *
 * the content of the slots.
 ******************************************************************/
#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))


/* Must be >= 6.  6 seems fastest. */
#define BLOCK_OFFSET_BITS (6)

//#define SLOTS_PER_BLOCK (1ULL << BLOCK_OFFSET_BITS)
//#define METADATA_WORDS_PER_BLOCK ((SLOTS_PER_BLOCK + 63) / 64)

#define NUM_SLOTS_TO_LOCK (1ULL<<16)
#define CLUSTER_SIZE (1ULL<<14)

#define METADATA_WORD(field,slot_index) (get_block((slot_index) / \
					SLOTS_PER_BLOCK)->field[((slot_index)  % SLOTS_PER_BLOCK) / 64])


uint64_t bitmaskLookup[]={0,1, 3, 7, 15, 31, 63, 127, 255, 511, 1023,
2047, 4095, 8191, 16383, 32767, 65535, 131071, 262143, 524287, 1048575,
2097151, 4194303, 8388607, 16777215, 33554431, 67108863, 134217727, 268435455, 536870911, 1073741823,
2147483647, 4294967295, 8589934591, 17179869183, 34359738367, 68719476735, 137438953471, 274877906943, 549755813887, 1099511627775,
2199023255551, 4398046511103, 8796093022207, 17592186044415, 35184372088831, 70368744177663, 140737488355327, 281474976710655, 562949953421311, 1125899906842623,
2251799813685247, 4503599627370495, 9007199254740991, 18014398509481983, 36028797018963967, 72057594037927935, 144115188075855871, 288230376151711743, 576460752303423487, 1152921504606846975,
2305843009213693951, 4611686018427387903, 9223372036854775807, 18446744073709551615};

//#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) \
												- 1ULL)
#define BITMASK(nbits) bitmaskLookup[nbits]

#define MAX_VALUE(nbits) ((1ULL << (nbits)) - 1)
#define BILLION 1000000000L

typedef struct __attribute__ ((__packed__)) qfblock {
	/* Code works with uint16_t, uint32_t, etc, but uint8_t seems just as fast as
   * anything else */
	uint8_t offset;
	uint64_t occupieds[METADATA_WORDS_PER_BLOCK];
	uint64_t runends[METADATA_WORDS_PER_BLOCK];
#if BITS_PER_SLOT == 8
	uint8_t  slots[SLOTS_PER_BLOCK];
#elif BITS_PER_SLOT == 16
	uint16_t  slots[SLOTS_PER_BLOCK];
#elif BITS_PER_SLOT == 32
	uint32_t  slots[SLOTS_PER_BLOCK];
#elif BITS_PER_SLOT == 64
	uint64_t  slots[SLOTS_PER_BLOCK];
#elif BITS_PER_SLOT != 0
	uint8_t   slots[SLOTS_PER_BLOCK * BITS_PER_SLOT / 8];
#else
	uint8_t   slots[];
#endif
} qfblock;

static __inline__ unsigned long long rdtsc(void)
{
	unsigned hi, lo;
	__asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
	return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

/**
 * Try to acquire a lock once and return even if the lock is busy.
 * If spin flag is set, then spin until the lock is available.
 */

static inline bool onDiskMQF_spin_lock(volatile int *lock, bool flag_spin)
{
	if (!flag_spin) {
		return !__sync_lock_test_and_set(lock, 1);
	} else {
		while (__sync_lock_test_and_set(lock, 1))
			while (*lock);
		return true;
	}

	return false;
}


static inline void onDiskMQF_spin_unlock(volatile int *lock)
{
	__sync_lock_release(lock);
	return;
}

template<uint64_t bitsPerSlot>
  bool _onDiskMQF<bitsPerSlot>::spin_lock(uint64_t hash_bucket_index, bool spin, bool flag)
{
	_onDiskMQF<bitsPerSlot> *cf=this;
	uint64_t hash_bucket_lock_offset  = hash_bucket_index % NUM_SLOTS_TO_LOCK;
	if (flag) {
		if (!onDiskMQF_spin_lock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK], spin))
			return false;
		if (NUM_SLOTS_TO_LOCK - hash_bucket_lock_offset <= CLUSTER_SIZE) {
			if (!onDiskMQF_spin_lock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1],
												spin)) {
				onDiskMQF_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
				return false;
			}
		}
	} else {
		/* take the lock for two lock-blocks; the lock-block in which the
		 * hash_bucket_index falls and the next lock-block */

		if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
				CLUSTER_SIZE) {
			if (!onDiskMQF_spin_lock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1], spin))
				return false;
		}
		if (!onDiskMQF_spin_lock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK], spin)) {
			if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
					CLUSTER_SIZE)
				onDiskMQF_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
			return false;
		}
		if (!onDiskMQF_spin_lock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1],
											spin)) {
			onDiskMQF_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
			if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
					CLUSTER_SIZE)
				onDiskMQF_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
			return false;
		}

	}
	return true;
}

template<uint64_t bitsPerSlot>
void _onDiskMQF<bitsPerSlot>::unlock(uint64_t hash_bucket_index, bool flag)
{
	_onDiskMQF<bitsPerSlot>* cf=this;
	uint64_t hash_bucket_lock_offset  = hash_bucket_index % NUM_SLOTS_TO_LOCK;
	if (flag) {
		if (NUM_SLOTS_TO_LOCK - hash_bucket_lock_offset <= CLUSTER_SIZE) {
			onDiskMQF_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1]);
		}
		onDiskMQF_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
	} else {
		onDiskMQF_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1]);
		onDiskMQF_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
		if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
				CLUSTER_SIZE)
			onDiskMQF_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
	}
}

template<uint64_t bitsPerSlot>
  void _onDiskMQF<bitsPerSlot>::modify_metadata(uint64_t *metadata, int cnt)
{
	_onDiskMQF<bitsPerSlot> *cf=this;
	onDiskMQF_spin_lock(&cf->mem->metadata_lock, true);
	*metadata = *metadata + cnt;
	onDiskMQF_spin_unlock(&cf->mem->metadata_lock);
	return;
}

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
inline void writeAndFreeBlocks(uint64_t memoryBufferIndex)
{
	// if(qf->blocksFilePos[memoryBufferIndex]!=0)
	// 	return;
  //
	// diskMQFStream.seek(qf->blocksFilePos[memoryBufferIndex]);
	// diskMQFStream.write((char*)blocks[memoryBufferIndex],(sizeof(qfblock) +
	// diskParams->nBlocksPerPointer*(8 * bits_per_slot )));
	// for(uint64_t i=0;i<diskParams->nBlocksPerPointer;i++){
	// 	qf->blocksFilePos[memoryBufferIndex]=0;
	// }
	// qf->blocksPointers[qf->reverseBlocksPointer[memoryBufferIndex]]=NULL;

}

// template<uint64_t bitsPerSlot>
// inline typename stxxl::vector<onDisk_qfblock<bitsPerSlot> >::iterator  get_block(onDiskMQF *qf, uint64_t block_index)
// {
// 	 // printf("block = %p\n",(void*)(((char *)qf->blocks) + block_index * (sizeof(qfblock) +
// 	 // 					 qf->metadata->bits_per_slot * 8 +
// 	 // 					8*qf->metadata->fixed_counter_size +
// 	 // 					8*qf->metadata->tag_bits
// 	 // 				 )) );
// 	 //printf("blocks start=%p\n",qf->blocks );
// 	// uint64_t BasepointerIndex= block_index / qf->diskParams->nBlocksPerPointer;
// 	// char* basePointer=(char*)qf->blocksPointers[BasepointerIndex];
// 	// if(basePointer==NULL)
// 	// {
// 	// 		//load
// 	// 		for(uint64_t i=0;i<qf->diskParams->nBlocksPerIOBatch;i+=qf->diskParams->nBlocksPerPointer)
// 	// 		{
// 	// 			writeAndFreeBlocks(qf->diskParams->memoryBufferPos+i);
// 	// 		}
// 	// }
// 	// uint64_t pointerShift= (block_index % qf->diskParams->nBlocksPerPointer)
// 	// 													*(sizeof(qfblock) + (8 * qf->metadata->bits_per_slot ));
// 	// return (qfblock*)(basePointer+pointerShift);
// 	//return (qfblock *)(((char *)qf->blocks) + block_index * (sizeof(qfblock) +
// 		//				 qf->metadata->bits_per_slot * 8
// 			//		 ));
//
// 			typename stxxl::vector<onDisk_qfblock<bitsPerSlot> >::iterator it = begin(qf->blocks);
//
// 			//return NULL;
//
// 		it+=block_index;
// 		//return (qfblock*)&(*it);
// 		return it;
// }
// template<uint64_t bitsPerSlot>
// inline typename stxxl::vector<onDisk_qfblock<bitsPerSlot> >::const_iterator  get_block_const(onDiskMQF *qf, uint64_t block_index)
// {
// 	 	typename stxxl::vector<onDisk_qfblock<bitsPerSlot> >::const_iterator it = begin(qf->blocks);
// 		it+=block_index;
// 		return it;
// }

template<uint64_t bitsPerSlot>
  inline int _onDiskMQF<bitsPerSlot>::is_runend(uint64_t index)
{
	return (METADATA_WORD( runends, index) >> ((index % SLOTS_PER_BLOCK) %
																								64)) & 1ULL;
}

template<uint64_t bitsPerSlot>
 inline int _onDiskMQF<bitsPerSlot>::is_occupied(uint64_t index)
{
	return (METADATA_WORD( occupieds, index) >> ((index % SLOTS_PER_BLOCK) %
																									64)) & 1ULL;
}







/* Little-endian code ....  Big-endian is TODO */
template<uint64_t bitsPerSlot>
  inline uint64_t _onDiskMQF<bitsPerSlot>::_get_slot(uint64_t index)
{
	_onDiskMQF<bitsPerSlot> *qf=this;
	assert(index < qf->metadata->xnslots);
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems
	 * to generate buggy code.  :/  */

	uint64_t *p = (uint64_t *)&qf->get_block_const( index /
																			 SLOTS_PER_BLOCK)->slots[(index %
																																SLOTS_PER_BLOCK)
																			 * qf->metadata->bits_per_slot / 8];
	return (uint64_t)(((*p) >> (((index % SLOTS_PER_BLOCK) *
															 qf->metadata->bits_per_slot) % 8)) &
										BITMASK(qf->metadata->bits_per_slot));
}

template<uint64_t bitsPerSlot>
  inline void _onDiskMQF<bitsPerSlot>::_set_slot(uint64_t index, uint64_t value)
{
	_onDiskMQF<bitsPerSlot> *qf=this;
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems
	 * to generate buggy code.  :/  */
	 //printf("ss %d\n",(index %SLOTS_PER_BLOCK)* qf->metadata->bits_per_slot / 8 );
	uint64_t *p = (uint64_t *)&qf->get_block( index /
																			 SLOTS_PER_BLOCK)->slots[(index %
																																SLOTS_PER_BLOCK)
																			 * qf->metadata->bits_per_slot / 8];

	int shift = ((index % SLOTS_PER_BLOCK) * qf->metadata->bits_per_slot) % 8;
	uint64_t mask = BITMASK(qf->metadata->bits_per_slot)<< shift;
	value <<= shift;
	*p=*p ^ ((*p ^ value) & mask);


}

template<uint64_t bitsPerSlot>
  inline void _onDiskMQF<bitsPerSlot>::set_slot(uint64_t index, uint64_t value){
	_onDiskMQF<bitsPerSlot> *qf=this;
	uint64_t original_value=_get_slot(index);
	value<<=qf->metadata->fixed_counter_size;
	uint64_t mask=BITMASK(qf->metadata->fixed_counter_size);
	original_value&=mask;
	value|=original_value;
	_set_slot(index,value);
}

template<uint64_t bitsPerSlot>
  inline uint64_t _onDiskMQF<bitsPerSlot>::get_slot(uint64_t index)
{
	_onDiskMQF<bitsPerSlot> *qf=this;
	uint64_t mask=BITMASK(qf->metadata->key_remainder_bits);
	uint64_t t=_get_slot(index);
	t >>= qf->metadata->fixed_counter_size;

	return t&mask;
}



template<uint64_t bitsPerSlot>
  inline uint64_t _onDiskMQF<bitsPerSlot>::get_fixed_counter(uint64_t index)
{
	_onDiskMQF<bitsPerSlot> *qf=this;
	uint64_t t=_get_slot(index);
	uint64_t mask=BITMASK(qf->metadata->fixed_counter_size);
	return t&mask;
	// uint64_t res=0;
	// uint64_t base=1;
	// uint64_t* p=(uint64_t*)((uint8_t*)qf->get_block( index /SLOTS_PER_BLOCK)->slots+
	// 															(SLOTS_PER_BLOCK * qf->metadata->bits_per_slot )/ 8);
  //
	// for(int i=qf->metadata->fixed_counter_size-1;i>=0;i--){
	// 	res+= base*(( p[i] >> ((index % SLOTS_PER_BLOCK) %64)) & 1ULL);
	// 	base*=2;
	// }
  //
	// return res;
}

template<uint64_t bitsPerSlot>
 inline  void _onDiskMQF<bitsPerSlot>::set_fixed_counter(uint64_t index,uint64_t value)
{
	_onDiskMQF<bitsPerSlot> *qf=this;
	uint64_t mask=BITMASK(qf->metadata->fixed_counter_size);
	uint64_t tvalue=_get_slot(index);
	tvalue=tvalue ^ ((tvalue ^ value) & mask);
	_set_slot(index,tvalue);
}

template<uint64_t bitsPerSlot>
  inline uint64_t _onDiskMQF<bitsPerSlot>::_get_tag(uint64_t index)
{
	_onDiskMQF<bitsPerSlot> *qf=this;
	uint64_t mask=BITMASK(qf->metadata->tag_bits);
	uint64_t t=_get_slot(index);
	t >>= (qf->metadata->fixed_counter_size+qf->metadata->key_remainder_bits);
	return t&mask;
}

template<uint64_t bitsPerSlot>
  inline void _onDiskMQF<bitsPerSlot>::set_tag(uint64_t index,uint64_t value)
{
	_onDiskMQF<bitsPerSlot> *qf=this;
	uint64_t original_value=_get_slot(index);
	value<<=(qf->metadata->fixed_counter_size+qf->metadata->key_remainder_bits);
	uint64_t mask=BITMASK(qf->metadata->fixed_counter_size+qf->metadata->key_remainder_bits);
	original_value&=mask;
	value|=original_value;
	_set_slot(index,value);

}

template<uint64_t bitsPerSlot>
  inline void _onDiskMQF<bitsPerSlot>::super_get(uint64_t index,uint64_t* slot,uint64_t *fcounter)
{
	_onDiskMQF<bitsPerSlot> *qf=this;
	uint64_t t=_get_slot(index);
	*fcounter=t & BITMASK(qf->metadata->fixed_counter_size);
	*slot= t >> qf->metadata->fixed_counter_size & BITMASK(qf->metadata->key_remainder_bits) ;
}

template<uint64_t bitsPerSlot>
  inline void _onDiskMQF<bitsPerSlot>::super_set(uint64_t index,uint64_t slot,uint64_t fcounter)
{
	_onDiskMQF<bitsPerSlot> *qf=this;
	uint64_t mask=BITMASK(qf->metadata->fixed_counter_size);
	slot<<=qf->metadata->fixed_counter_size;
	uint64_t newValue= slot ^ ((slot ^ fcounter) & mask);
	//cout<<"slot= "<<slot<<" fcounter= "<<fcounter<<" newValue= "<<newValue<<endl;
	uint64_t original_value=_get_slot(index);

	mask=BITMASK(qf->metadata->fixed_counter_size+
										qf->metadata->key_remainder_bits);
	newValue= original_value ^ ((original_value ^ newValue) & mask);

	_set_slot(index,newValue);


}


template<uint64_t bitsPerSlot>
  inline uint64_t _onDiskMQF<bitsPerSlot>::block_offset(uint64_t blockidx)
{
	_onDiskMQF<bitsPerSlot> *qf=this;
	/* If we have extended counters and a 16-bit (or larger) offset
		 field, then we can safely ignore the possibility of overflowing
		 that field. */
	if (sizeof(qf->get_block(0)->offset) > 1 ||
			qf->get_block( blockidx)->offset < BITMASK(8*sizeof(qf->get_block(0)->offset)))
		return qf->get_block( blockidx)->offset;

	return run_end( SLOTS_PER_BLOCK * blockidx - 1) - SLOTS_PER_BLOCK *
		blockidx + 1;
}

template<uint64_t bitsPerSlot>
  inline uint64_t _onDiskMQF<bitsPerSlot>::run_end(uint64_t hash_bucket_index)
{
	_onDiskMQF<bitsPerSlot> *qf=this;
	uint64_t bucket_block_index       = hash_bucket_index / SLOTS_PER_BLOCK;
	uint64_t bucket_intrablock_offset = hash_bucket_index % SLOTS_PER_BLOCK;
	uint64_t bucket_blocks_offset = block_offset( bucket_block_index);

	uint64_t bucket_intrablock_rank   = bitrank(qf->get_block(
																				bucket_block_index)->occupieds[0],
																				bucket_intrablock_offset);

	if (bucket_intrablock_rank == 0) {
		if (bucket_blocks_offset <= bucket_intrablock_offset)
			return hash_bucket_index;
		else
			return SLOTS_PER_BLOCK * bucket_block_index + bucket_blocks_offset - 1;
	}

	uint64_t runend_block_index  = bucket_block_index + bucket_blocks_offset /
		SLOTS_PER_BLOCK;
	uint64_t runend_ignore_bits  = bucket_blocks_offset % SLOTS_PER_BLOCK;
	uint64_t runend_rank         = bucket_intrablock_rank - 1;
	uint64_t runend_block_offset = bitselectv(qf->get_block(
																						runend_block_index)->runends[0],
																						runend_ignore_bits, runend_rank);
	if (runend_block_offset == SLOTS_PER_BLOCK) {
		if (bucket_blocks_offset == 0 && bucket_intrablock_rank == 0) {
			/* The block begins in empty space, and this bucket is in that region of
			 * empty space */
			return hash_bucket_index;
		} else {
			do {
				runend_rank        -= popcntv(qf->get_block(
																								runend_block_index)->runends[0],
																			runend_ignore_bits);
				runend_block_index++;
				runend_ignore_bits  = 0;
				runend_block_offset = bitselectv(qf->get_block(
																									 runend_block_index)->runends[0],
																				 runend_ignore_bits, runend_rank);
			} while (runend_block_offset == SLOTS_PER_BLOCK);
		}
	}

	uint64_t runend_index = SLOTS_PER_BLOCK * runend_block_index +
		runend_block_offset;
	if (runend_index < hash_bucket_index)
		return hash_bucket_index;
	else
		return runend_index;
}

//   inline uint64_t run_end2(onDiskMQF *qf, uint64_t hash_bucket_index)
// {
// 	uint8_t offset=hash_bucket_index % SLOTS_PER_BLOCK;
// 	uint64_t t= METADATA_WORD( runends, hash_bucket_index)
// 	&BITMASK(offset);
// 	uint64_t tt=hash_bucket_index / SLOTS_PER_BLOCK;
// 	hash_bucket_index+=SLOTS_PER_BLOCK-offset;
// 	while(t==0){
// 		t= METADATA_WORD( runends, hash_bucket_index);
// 		hash_bucket_index+=64;
// 		tt+=1;
// 	}
// 	uint64_t res=(tt*64)+SLOTS_PER_BLOCK-LOG2(t)-1;
// 	return res;
// }

template<uint64_t bitsPerSlot>
  inline int _onDiskMQF<bitsPerSlot>::offset_lower_bound(uint64_t slot_index)
{
	_onDiskMQF<bitsPerSlot> *qf=this;
	typename stxxl::vector<onDisk_qfblock<bitsPerSlot> >::const_iterator  b = qf->get_block_const( slot_index / SLOTS_PER_BLOCK);
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

template<uint64_t bitsPerSlot>
  inline int _onDiskMQF<bitsPerSlot>::is_empty2(uint64_t slot_index)
{
	return offset_lower_bound( slot_index) == 0;
}

template<uint64_t bitsPerSlot>
  inline int _onDiskMQF<bitsPerSlot>::might_be_empty(uint64_t slot_index)
{
	return !is_occupied( slot_index)
		&& !is_runend( slot_index);
}

template<uint64_t bitsPerSlot>
  inline int _onDiskMQF<bitsPerSlot>::probably_is_empty(uint64_t slot_index)
{
	_onDiskMQF<bitsPerSlot> *qf=this;
	return get_slot( slot_index) == 0
		&& !is_occupied( slot_index)
		&& !is_runend( slot_index);
}

template<uint64_t bitsPerSlot>
  inline uint64_t _onDiskMQF<bitsPerSlot>::find_first_empty_slot(uint64_t from)
{
	do {
		int t = offset_lower_bound( from);
		if(t<0)
		{
			cout<<"here"<<endl;
		}
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


#define REMAINDER_WORD(i) ((uint64_t *)&(get_block( (i)/metadata->bits_per_slot)->slots[8 * ((i) % metadata->bits_per_slot)]))

template<uint64_t bitsPerSlot>
  inline void _onDiskMQF<bitsPerSlot>::shift_remainders(const uint64_t start_index, const
																		uint64_t empty_index)
{
	_onDiskMQF<bitsPerSlot> *qf=this;
	uint64_t last_word = (empty_index + 1) * qf->metadata->bits_per_slot / 64;
	const uint64_t first_word = start_index * qf->metadata->bits_per_slot / 64;
	int bend = ((empty_index + 1) * qf->metadata->bits_per_slot) % 64;
	const int bstart = (start_index * qf->metadata->bits_per_slot) % 64;

	while (last_word != first_word) {
		*REMAINDER_WORD( last_word) = shift_into_b(*REMAINDER_WORD( last_word-1),
																									*REMAINDER_WORD( last_word),
																									0, bend, qf->metadata->bits_per_slot);
		last_word--;
		bend = 64;
	}
	*REMAINDER_WORD( last_word) = shift_into_b(0, *REMAINDER_WORD(
																																	 last_word),
																								bstart, bend,
																								qf->metadata->bits_per_slot);
}

template<uint64_t bitsPerSlot>
  inline void _onDiskMQF<bitsPerSlot>::dump_block(uint64_t i)
{
	_onDiskMQF<bitsPerSlot>* qf=this;
	uint64_t j;
 	printf("#Block %lu \n",i );
	printf("%-192d", qf->get_block( i)->offset);
	printf("\n");

	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf("%02lx ", j);
	printf("\n");

	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf(" %d ", (qf->get_block( i)->occupieds[j/64] & (1ULL << (j%64))) ? 1 : 0);
	printf("\n");

	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf(" %d ", (qf->get_block( i)->runends[j/64] & (1ULL << (j%64))) ? 1 : 0);
	printf("\n");

#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32
	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf("%02x ", qf->get_block( i)->slots[j]);
#elif BITS_PER_SLOT == 64
	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf("%02lx ", qf->get_block( i)->slots[j]);
#else
	//for (j = 0; j < SLOTS_PER_BLOCK * qf->metadata->bits_per_slot / 8; j++)
	//	printf("%02x ", qf->get_block( i)->slots[j]);

#endif

for(int j=0;j<64 && (j+64*i)<qf->metadata->xnslots;j++)
{
	printf("%lu ", get_slot(j+64*i));
}

	printf("\n fixed counter \n");

	for(int j=0;j<64&& (j+64*i)<qf->metadata->xnslots;j++)
	{
		printf("%lu ", get_fixed_counter(j+64*i));
	}

	printf("\n tags \n");

	for(int j=0;j<64&& (j+64*i)<qf->metadata->xnslots;j++)
	{
		printf("%lu ", _get_tag(j+64*i));
	}


	printf("\n");
	printf("\n");
}

template<uint64_t bitsPerSlot>
void _onDiskMQF<bitsPerSlot>::dump()
{
	_onDiskMQF<bitsPerSlot>* qf=this;
	uint64_t i;

	printf("%lu %lu %lu\n",
				 qf->metadata->nblocks,
				 qf->metadata->ndistinct_elts,
				 qf->metadata->nelts);

	for (i = 0; i < qf->metadata->nblocks; i++) {
		if(i==540){
			cout<<"last bloc"<<endl;
		}
		dump_block(i);
	}
	printf("End\n");




}

template<uint64_t bitsPerSlot>
  inline void _onDiskMQF<bitsPerSlot>::find_next_n_empty_slots(uint64_t from, uint64_t n,
																					 uint64_t *indices)
{
	while (n) {
		indices[--n] = find_first_empty_slot( from);
		from = indices[n] + 1;
	}
}

template<uint64_t bitsPerSlot>
  inline void _onDiskMQF<bitsPerSlot>::shift_slots(int64_t first, uint64_t last, uint64_t
															 distance)
{
	int64_t i;
	if (distance == 1)
		shift_remainders( first, last+1);
	else
		for (i = last; i >= first; i--)
			_set_slot( i + distance, _get_slot( i));
}

template<uint64_t bitsPerSlot>
  inline void _onDiskMQF<bitsPerSlot>::shift_runends(int64_t first, uint64_t last,
																 uint64_t distance)
{
	assert(last < metadata->xnslots && distance < 64);
	uint64_t first_word = first / 64;
	uint64_t bstart = first % 64;
	uint64_t last_word = (last + distance + 1) / 64;
	uint64_t bend = (last + distance + 1) % 64;

	if (last_word != first_word) {
		METADATA_WORD( runends, 64*last_word) = shift_into_b(METADATA_WORD( runends, 64*(last_word-1)),
																														METADATA_WORD( runends, 64*last_word),
																														0, bend, distance);
		bend = 64;
		last_word--;
		while (last_word != first_word) {
			METADATA_WORD( runends, 64*last_word) = shift_into_b(METADATA_WORD( runends, 64*(last_word-1)),
																															METADATA_WORD( runends, 64*last_word),
																															0, bend, distance);
			last_word--;
		}
	}
	METADATA_WORD( runends, 64*last_word) = shift_into_b(0, METADATA_WORD(
																																					 runends,
																																					 64*last_word),
																													bstart, bend, distance);

}

//   inline void shift_fixed_counters(onDiskMQF *qf, int64_t first, uint64_t last,
// 																 uint64_t distance)
// {
// 	assert(last < qf->metadata->xnslots && distance < 64);
// 	uint64_t first_word = first / 64;
// 	uint64_t bstart = first % 64;
// 	uint64_t last_word = (last + distance + 1) / 64;
// 	uint64_t bend = (last + distance + 1) % 64;
// 	uint64_t* curr, *prev;
// 	uint64_t tmp =last_word, tmp_bend=bend;
// 	for(int i=0;i<qf->metadata->fixed_counter_size;i++){
// 		last_word=tmp;
// 		bend=tmp_bend;
// 		if (last_word != first_word) {
// 			curr=(uint64_t*)((uint8_t*)qf->get_block( last_word)->slots+
// 			(SLOTS_PER_BLOCK * qf->metadata->bits_per_slot) / 8);
// 			prev=(uint64_t*)((uint8_t*)qf->get_block( last_word-1)->slots+
// 			(SLOTS_PER_BLOCK * qf->metadata->bits_per_slot) / 8);
// 			curr[i] = shift_into_b(prev[i],curr[i],0, bend, distance);
// 			bend = 64;
// 			last_word--;
// 			while (last_word != first_word) {
// 				curr=(uint64_t*)((uint8_t*)qf->get_block( last_word)->slots+(SLOTS_PER_BLOCK * qf->metadata->bits_per_slot / 8));
// 				prev=(uint64_t*)((uint8_t*)qf->get_block( last_word-1)->slots+(SLOTS_PER_BLOCK * qf->metadata->bits_per_slot / 8));
// 				curr[i] = shift_into_b(prev[i],curr[i],0, bend, distance);
// 				last_word--;
// 			}
// 		}
// 		curr=(uint64_t*)((uint8_t*)qf->get_block( last_word)->slots+(SLOTS_PER_BLOCK * qf->metadata->bits_per_slot / 8));
// 		curr[i] = shift_into_b(0,curr[i], bstart, bend, distance);
// 	}
//
// }

//   inline void shift_tags(QF *qf, int64_t first, uint64_t last,
// 																 uint64_t distance)
// {
// 	assert(last < qf->metadata->xnslots && distance < 64);
// 	uint64_t first_word = first / 64;
// 	uint64_t bstart = first % 64;
// 	uint64_t last_word = (last + distance + 1) / 64;
// 	uint64_t bend = (last + distance + 1) % 64;
// 	uint64_t* curr, *prev;
// 	uint64_t tmp =last_word, tmp_bend=bend;
// 	for(int i=0;i<qf->metadata->tag_bits;i++){
// 		last_word=tmp;
// 		bend=tmp_bend;
// 		if (last_word != first_word) {
// 			curr=(uint64_t*)((uint8_t*)qf->get_block( last_word)->slots+
// 											(8 *(qf->metadata->bits_per_slot + qf->metadata->fixed_counter_size )));
// 			prev=(uint64_t*)((uint8_t*)qf->get_block( last_word-1)->slots+
// 											(8 *(qf->metadata->bits_per_slot + qf->metadata->fixed_counter_size )));
// 			curr[i] = shift_into_b(prev[i],curr[i],0, bend, distance);
// 			bend = 64;
// 			last_word--;
// 			while (last_word != first_word) {
// 				curr=(uint64_t*)((uint8_t*)qf->get_block( last_word)->slots+
// 												(8 *(qf->metadata->bits_per_slot + qf->metadata->fixed_counter_size )));
// 				prev=(uint64_t*)((uint8_t*)qf->get_block( last_word-1)->slots+
// 												(8 *(qf->metadata->bits_per_slot + qf->metadata->fixed_counter_size )));
// 				curr[i] = shift_into_b(prev[i],curr[i],0, bend, distance);
// 				last_word--;
// 			}
// 		}
// 		curr=(uint64_t*)((uint8_t*)qf->get_block( last_word)->slots+
// 										(8 *(qf->metadata->bits_per_slot + qf->metadata->fixed_counter_size )));
// 		curr[i] = shift_into_b(0,curr[i], bstart, bend, distance);
// 	}
//
// }

template<uint64_t bitsPerSlot>
  inline void _onDiskMQF<bitsPerSlot>::insert_replace_slots_and_shift_remainders_and_runends_and_offsets(int		 operation,
																																										 uint64_t		 bucket_index,
																																										 uint64_t		 overwrite_index,
																																										 const uint64_t	*remainders,
																																										 const uint64_t	*fixed_size_counters,
																																										 uint64_t		 total_remainders,
																																										 uint64_t		 noverwrites)
{
	_onDiskMQF<bitsPerSlot> *qf=this;
	uint64_t empties[67];
	uint64_t i;
	int64_t ninserts = total_remainders - noverwrites;
	uint64_t insert_index = overwrite_index + noverwrites;
	if(qf->metadata->noccupied_slots+ninserts > qf->metadata->maximum_occupied_slots )
	{
		throw std::overflow_error("QF is 95% full, cannot insert more items.");
	}
	//printf("remainder =%lu ,overwrite_index = %lu , insert_index=%lu , operation=%d, noverwites=%lu total_remainders=%lu nnserts=%lu \n", remainders[0],overwrite_index,insert_index,operation,noverwrites,total_remainders,ninserts);
	if (ninserts > 0) {
		/* First, shift things to create n empty spaces where we need them. */
		//printf("shift %lu, ninserts=%lu\n",insert_index,ninserts );

		find_next_n_empty_slots( insert_index, ninserts, empties);
		for (i = 0; i < ninserts - 1; i++){
			shift_slots( empties[i+1] + 1, empties[i] - 1, i + 1);
		}
		shift_slots( insert_index, empties[ninserts - 1] - 1, ninserts);



		for (i = 0; i < ninserts - 1; i++)
			shift_runends( empties[i+1] + 1, empties[i] - 1, i + 1);
		shift_runends( insert_index, empties[ninserts - 1] - 1, ninserts);



		// for (i = 0; i < ninserts - 1; i++)
		// 	shift_fixed_counters( empties[i+1] + 1, empties[i] - 1, i + 1);
		// shift_fixed_counters( insert_index, empties[ninserts - 1] - 1, ninserts);
    //
    //
		// for (i = 0; i < ninserts - 1; i++)
		// 	shift_tags( empties[i+1] + 1, empties[i] - 1, i + 1);
		// shift_tags( insert_index, empties[ninserts - 1] - 1, ninserts);
    //



		for (i = noverwrites; i < total_remainders - 1; i++)
			METADATA_WORD( runends, overwrite_index + i) &= ~(1ULL <<
																													 (((overwrite_index
																															+ i) %
																														 SLOTS_PER_BLOCK)
																														% 64));
		// for (i = noverwrites; i < total_remainders - 1; i++)
		// 	set_fixed_counter(overwrite_index+i,0);


		switch (operation) {
			case 0: /* insert into empty bucket */
				assert (noverwrites == 0);
				METADATA_WORD( runends, overwrite_index + total_remainders - 1) |=
					1ULL << (((overwrite_index + total_remainders - 1) %
										SLOTS_PER_BLOCK) % 64);
				break;
			case 1: /* append to bucket */
				METADATA_WORD( runends, overwrite_index + noverwrites - 1)      &=
					~(1ULL << (((overwrite_index + noverwrites - 1) % SLOTS_PER_BLOCK) %
										 64));
				METADATA_WORD( runends, overwrite_index + total_remainders - 1) |=
					1ULL << (((overwrite_index + total_remainders - 1) %
										SLOTS_PER_BLOCK) % 64);
				break;
			case 2: /* insert into bucket */
				METADATA_WORD( runends, overwrite_index + total_remainders - 1) &=
					~(1ULL << (((overwrite_index + total_remainders - 1) %
											SLOTS_PER_BLOCK) % 64));
				break;
			default:
				fprintf(stderr, "Invalid operation %d\n", operation);
				abort();
		}

		uint64_t npreceding_empties = 0;
		for (i = bucket_index / SLOTS_PER_BLOCK + 1; i <= empties[0]/SLOTS_PER_BLOCK; i++) {
			while (npreceding_empties < ninserts &&
						 empties[ninserts - 1 - npreceding_empties]  / SLOTS_PER_BLOCK < i)
				npreceding_empties++;

			if (qf->get_block( i)->offset + ninserts - npreceding_empties < BITMASK(8*sizeof(qf->get_block(0)->offset)))
				qf->get_block( i)->offset += ninserts - npreceding_empties;
			else
				qf->get_block( i)->offset = (uint8_t) BITMASK(8*sizeof(qf->get_block(0)->offset));
		}
	}
	for (i = 0; i < total_remainders; i++){
		//set_slot( overwrite_index + i, remainders[i]);
		//set_fixed_counter(overwrite_index +i,fixed_size_counters[i]);
		super_set(overwrite_index+i,remainders[i],fixed_size_counters[i]);
//		printf("fixed counter = %lu\n",fixed_size_counters[i] );
	}

	modify_metadata( &qf->metadata->noccupied_slots, ninserts);
}


template<uint64_t bitsPerSlot>
  inline void _onDiskMQF<bitsPerSlot>::remove_replace_slots_and_shift_remainders_and_runends_and_offsets(int		 operation,
																																										 uint64_t		 bucket_index,
																																										 uint64_t		 overwrite_index,
																																										 const uint64_t	*remainders,
																																										 const uint64_t	*fcounters,
																																										 uint64_t		 total_remainders,
																																										 uint64_t		 old_length)
{
	_onDiskMQF<bitsPerSlot> *qf=this;
	uint64_t i;

	// Update the slots
	for (i = 0; i < total_remainders; i++){
		set_slot( overwrite_index + i, remainders[i]);
		set_fixed_counter( overwrite_index + i, fcounters[i]);
		if(qf->metadata->tag_bits>0)
			set_tag( overwrite_index + i, 0);
	}


	// If this is the last thing in its run, then we may need to set a new runend bit
	if (is_runend( overwrite_index + old_length - 1)) {
	  if (total_remainders > 0) {
	    // If we're not deleting this entry entirely, then it will still the last entry in this run
	    METADATA_WORD( runends, overwrite_index + total_remainders - 1) |= 1ULL << ((overwrite_index + total_remainders - 1) % 64);
	  } else if (overwrite_index > bucket_index &&
		     !is_runend( overwrite_index - 1)) {
	    // If we're deleting this entry entirely, but it is not the first entry in this run,
	    // then set the preceding entry to be the runend
	    METADATA_WORD( runends, overwrite_index - 1) |= 1ULL << ((overwrite_index - 1) % 64);
	  }
	}

	// shift slots back one run at a time
	uint64_t original_bucket = bucket_index;
	uint64_t current_bucket = bucket_index;
	uint64_t current_slot = overwrite_index + total_remainders;
	uint64_t current_distance = old_length - total_remainders;

	while (current_distance > 0) {
		if (is_runend( current_slot + current_distance - 1)) {
			do {
				current_bucket++;
			} while (current_bucket < current_slot + current_distance &&
							 !is_occupied( current_bucket));
		}

		if (current_bucket <= current_slot) {
			set_slot( current_slot, get_slot( current_slot + current_distance));
			set_fixed_counter( current_slot, get_fixed_counter( current_slot + current_distance));
			if(qf->metadata->tag_bits>0)
				set_tag( current_slot, _get_tag( current_slot + current_distance));

			if (is_runend( current_slot) !=
					is_runend( current_slot + current_distance))
				METADATA_WORD( runends, current_slot) ^= 1ULL << (current_slot % 64);
			current_slot++;

		} else if (current_bucket <= current_slot + current_distance) {
			uint64_t i;
			for (i = current_slot; i < current_slot + current_distance; i++) {
				set_slot( i, 0);
				set_fixed_counter(i,0);
				if(qf->metadata->tag_bits>0)
					set_tag(i,0);
				METADATA_WORD( runends, i) &= ~(1ULL << (i % 64));
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
		METADATA_WORD( occupieds, bucket_index) &= ~(1ULL << (bucket_index % 64));

	// update the offset bits.
	// find the number of occupied slots in the original_bucket block.
	// Then find the runend slot corresponding to the last run in the
	// original_bucket block.
	// Update the offset of the block to which it belongs.
	uint64_t original_block = original_bucket / SLOTS_PER_BLOCK;
	while (1 && old_length > total_remainders) {	// we only update offsets if we shift/delete anything
		int32_t last_occupieds_bit = bitscanreverse(qf->get_block( original_block)->occupieds[0]);
		// there is nothing in the block
		if (last_occupieds_bit == -1) {
			if (qf->get_block( original_block + 1)->offset == 0)
				break;
			qf->get_block( original_block + 1)->offset = 0;
		} else {
			uint64_t last_occupieds_hash_index = SLOTS_PER_BLOCK * original_block + last_occupieds_bit;
			uint64_t runend_index = run_end( last_occupieds_hash_index);
			// runend spans across the block
			// update the offset of the next block
			if (runend_index / SLOTS_PER_BLOCK == original_block) { // if the run ends in the same block
				if (qf->get_block( original_block + 1)->offset == 0)
					break;
				qf->get_block( original_block + 1)->offset = 0;
			} else if (runend_index / SLOTS_PER_BLOCK == original_block + 1) { // if the last run spans across one block
				if (qf->get_block( original_block + 1)->offset == (runend_index % SLOTS_PER_BLOCK) + 1)
					break;
				qf->get_block( original_block + 1)->offset = (runend_index % SLOTS_PER_BLOCK) + 1;
			} else { // if the last run spans across multiple blocks
				uint64_t i;
				for (i = original_block + 1; i < runend_index / SLOTS_PER_BLOCK - 1; i++)
					qf->get_block( i)->offset = SLOTS_PER_BLOCK;
				if (qf->get_block( runend_index / SLOTS_PER_BLOCK)->offset == (runend_index % SLOTS_PER_BLOCK) + 1)
					break;
				qf->get_block( runend_index / SLOTS_PER_BLOCK)->offset = (runend_index % SLOTS_PER_BLOCK) + 1;
			}
		}
		original_block++;
	}

	int num_slots_freed = old_length - total_remainders;
	modify_metadata( &qf->metadata->noccupied_slots, -num_slots_freed);
	/*qf->metadata->noccupied_slots -= (old_length - total_remainders);*/
	if (!total_remainders) {
		modify_metadata( &qf->metadata->ndistinct_elts, -1);
		/*qf->metadata->ndistinct_elts--;*/
	}
}

/*****************************************************************************
 * Code that uses the above to implement a QF with keys and inline counters. *
 *****************************************************************************/

/*
	 Counter format:
	 1- count <= 2^(fixed counter size)
	 	Slots :           [Remaining]
		Fixed counters :  [count-1]

	 2- count > 2^(fixed counter size)
	 Slots :           [Remaining] [first digit] [second digit] ... [last digit]
	 Fixed counters :  [Maximum]     [Maximum]   [Maximum]      ... [up to Maximum -1]

 */
template<uint64_t bitsPerSlot>
  inline uint64_t *_onDiskMQF<bitsPerSlot>::encode_counter(uint64_t remainder, uint64_t
																			 counter, uint64_t *slots, uint64_t *fixed_size_counters)
{
	_onDiskMQF<bitsPerSlot> *qf=this;
	//printf("inserting %lu repeated %lu\n",remainder,counter);
	const uint64_t slots_base = (1ULL << qf->metadata->key_remainder_bits) ;
	uint64_t *p  = slots;
	uint64_t *pf = fixed_size_counters;
	const uint64_t fixed_counter_max=(1ULL<<qf->metadata->fixed_counter_size)-1;

	if (counter == 0)
	 	return p;

	counter--;
	uint64_t fcounter_first=std::min(counter,fixed_counter_max);
	counter-=(fcounter_first);

	//printf("first fixed counter =%lu\n", fcounter_first);
	if(fcounter_first==fixed_counter_max){
		uint64_t max_count_in_fixed_counter=fixed_counter_max-1;// the fixed size count in the end of the counter should'nt be full
		do{
			*--p=counter%slots_base;
			*--pf=fixed_counter_max;
			//printf("vcount = %lu\n",counter%slots_base );
			counter >>= qf->metadata->key_remainder_bits;
		}	while(counter>max_count_in_fixed_counter);
		*(fixed_size_counters-1)=counter;// set the last counter
		//printf("last fixed counter = %lu\n",counter);
	}

	*--p = remainder;
	*--pf=fcounter_first;


	return p;
}

/* Returns the length of the encoding.
REQUIRES: index points to first slot of a counter. */
template<uint64_t bitsPerSlot>
  inline uint64_t _onDiskMQF<bitsPerSlot>::decode_counter(uint64_t index, uint64_t
																			*remainder, uint64_t *count)
{

	_onDiskMQF<bitsPerSlot> *qf=this;
	// *remainder  = get_slot( index);
	// uint64_t fcount=get_fixed_counter(index);
	uint64_t fcount;
	super_get(index,remainder,&fcount);
	uint64_t tmp_count= fcount+1;
	*count=0;
	uint64_t tmp_slot;
	const uint64_t fixed_count_max=(1ULL << qf->metadata->fixed_counter_size)-1;
	//printf("tmp count = %lu\n",tmp_count);
	if(fcount == fixed_count_max){
		uint64_t no_digits=0;
		do{
				index++;
				no_digits++;
				*count <<= qf->metadata->key_remainder_bits;

				//fcount= get_fixed_counter(index);
				//*count += get_slot( index);
				super_get(index,&tmp_slot,&fcount);
				*count+=tmp_slot;


		}while(fcount == fixed_count_max);
		*count += fcount<<(no_digits*qf->metadata->key_remainder_bits);
		//printf("fixed vcount= %lu\n", fcount<<(no_digits*qf->metadata->bits_per_slot + qf->metadata->fixed_counter_size));
	}
	*count += tmp_count;
	return index;

}

/* return the next slot which corresponds to a
 * different element
 * */
//   inline uint64_t next_slot(onDiskMQF *qf, uint64_t current)
// {
// 	uint64_t rem = get_slot( current);
// 	current++;
//
// 	while (get_slot( current) == rem && current <= qf->metadata->nslots) {
// 		current++;
// 	}
// 	return current;
// }
template<uint64_t bitsPerSlot>
  inline bool _onDiskMQF<bitsPerSlot>::insert1(__uint128_t hash, bool lock, bool spin)
{
	_onDiskMQF<bitsPerSlot> *qf=this;
	uint64_t hash_remainder           = hash & BITMASK(qf->metadata->key_remainder_bits);
	uint64_t hash_bucket_index        = hash >> qf->metadata->key_remainder_bits;
	uint64_t hash_bucket_block_offset = hash_bucket_index % SLOTS_PER_BLOCK;

	uint64_t fixed_count_max=((1ULL)<<qf->metadata->fixed_counter_size)-1;

//	uint64_t blockId=hash_bucket_index/64;
  // cout<<"before "<<(uint64_t)hash<<endl;
	// dump_block(blockId);


	if (lock) {
		if(qf->mem->general_lock)
			return false;
		if (!spin_lock( hash_bucket_index, spin, true))
			return false;
	}
	 /* might_be_empty( hash_bucket_index) && runend_index == hash_bucket_index */
	if (is_empty2( hash_bucket_index)) {
		METADATA_WORD( runends, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);
		// set_slot( hash_bucket_index, hash_remainder);
		// set_fixed_counter( hash_bucket_index, 0);

		super_set(hash_bucket_index,hash_remainder,0);
		METADATA_WORD( occupieds, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);

		modify_metadata( &qf->metadata->ndistinct_elts, 1);
		modify_metadata( &qf->metadata->noccupied_slots, 1);
		/*modify_metadata( &qf->metadata->nelts, 1);*/
	} else {
		uint64_t runend_index= run_end( hash_bucket_index);
		int operation = 0; /* Insert into empty bucket */
		uint64_t insert_index = runend_index + 1;
		uint64_t new_value = hash_remainder;
		uint64_t new_fixedCounter=0;
		/* printf("RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

		uint64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(
																																	 hash_bucket_index
																																	 - 1) + 1;

		if (is_occupied( hash_bucket_index)) {

			/* Find the counter for this remainder if it exists. */
			uint64_t current_remainder,current_fixed_counter;
			super_get( runstart_index,&current_remainder,&current_fixed_counter);
			//uint64_t zero_terminator = runstart_index;


			/* Skip over counters for other remainders. */
			while (current_remainder < hash_remainder && runstart_index <=
						 runend_index) {

			 while (runstart_index < runend_index &&
			 							 current_fixed_counter == fixed_count_max){
											 current_fixed_counter=get_fixed_counter( ++runstart_index);
										 }
        runstart_index++;

				/* This may read past the end of the run, but the while loop
					 condition will prevent us from using the invalid result in
					 that case. */
					 super_get(runstart_index,&current_remainder,&current_fixed_counter);
			}

			/* If this is the first time we've inserted the new remainder,

				 and it is larger than any remainder in the run. */
			if (runstart_index > runend_index) {
				operation = 1;
				insert_index = runstart_index;
				new_value = hash_remainder;
				new_fixedCounter=0;
				/*modify_metadata( &qf->metadata->ndistinct_elts, 1);*/

				/* This is the first time we're inserting this remainder, but
					 there are larger remainders already in the run. */
			} else if (current_remainder != hash_remainder) {
				operation = 2; /* Inserting */
				insert_index = runstart_index;
				new_value = hash_remainder;
				new_fixedCounter=0;
				/*modify_metadata( &qf->metadata->ndistinct_elts, 1);*/


				/* If there's exactly one instance of this remainder. */
			} else if (current_fixed_counter<fixed_count_max) {
				set_fixed_counter(runstart_index,current_fixed_counter+1);
				if(current_fixed_counter+1<fixed_count_max){
					operation=-1;
				}
				else{

					operation =runstart_index==runend_index? 1:2 ; /* Insert */

					insert_index = runstart_index+1;
					new_value = 0;
					new_fixedCounter=0;
				}
				/* There is an extended counter for this remainder. */
			} else
			 {
				 operation=-1;
				/* Move to the LSD of the counter. */
				insert_index = runstart_index + 1;
				super_get(insert_index,&current_remainder,&current_fixed_counter);
				while (current_fixed_counter==fixed_count_max)
						super_get(++insert_index,&current_remainder,&current_fixed_counter);


				uint64_t endCounterIndex=insert_index;

				uint64_t digit=current_remainder, carry;

				/* Increment the counter. */
				do {
					carry = 0;
					// Increment the digit
					digit = (digit + 1) & BITMASK(qf->metadata->key_remainder_bits);
					carry=(digit==0);
					set_slot( insert_index, digit);
					insert_index--;
					super_get(insert_index,&current_remainder,&current_fixed_counter);
					digit=current_remainder;
				} while(insert_index > runstart_index && carry);

				if(insert_index==runstart_index &&carry>0){
					current_fixed_counter=get_fixed_counter(endCounterIndex);
					set_fixed_counter(endCounterIndex,current_fixed_counter+1);
					if(current_fixed_counter+1<fixed_count_max){
					for(uint32_t i=runstart_index+1;i<=endCounterIndex;i++)
						set_slot(i,0);
					}
					else{
						operation=2;
						insert_index=runstart_index+1;
						new_value=fixed_count_max;
						new_fixedCounter=fixed_count_max;
						set_fixed_counter(endCounterIndex,0);
						//set_slot(runstart_index+1,1);
					}
				}
			}
		}

		if (operation >= 0) {
			uint64_t empty_slot_index = find_first_empty_slot( runend_index+1);

			shift_slots( insert_index, empty_slot_index-1,1);


			set_slot( insert_index, new_value);
			set_fixed_counter(insert_index,new_fixedCounter);

			shift_runends( insert_index, empty_slot_index-1, 1);

			switch (operation) {
				case 0:
					METADATA_WORD( runends, insert_index)   |= 1ULL << ((insert_index
																																	%
																																	SLOTS_PER_BLOCK)
																																 % 64);
					break;
				case 1:
					METADATA_WORD( runends, insert_index-1) &= ~(1ULL <<
																													(((insert_index-1) %
																														SLOTS_PER_BLOCK) %
																													 64));

					METADATA_WORD( runends, insert_index)   |= 1ULL << ((insert_index
																																	%
																																	SLOTS_PER_BLOCK)
																																 % 64);
					break;
				case 2:
					METADATA_WORD( runends, insert_index)   &= ~(1ULL <<
																													((insert_index %
																														SLOTS_PER_BLOCK) %
																													 64));
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
			for (i = hash_bucket_index / SLOTS_PER_BLOCK + 1; i <=
					 empty_slot_index/SLOTS_PER_BLOCK; i++) {
				if (qf->get_block( i)->offset < BITMASK(8*sizeof(qf->get_block(0)->offset)))
					qf->get_block( i)->offset++;
				assert(qf->get_block( i)->offset != 0);
			}
			modify_metadata( &qf->metadata->noccupied_slots, 1);
		}
		/*modify_metadata( &qf->metadata->nelts, 1);*/
		METADATA_WORD( occupieds, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);
	}



	if (lock) {
		unlock( hash_bucket_index, true);
	}

	return true;
}

template<uint64_t bitsPerSlot>
  inline bool _onDiskMQF<bitsPerSlot>::_insert( __uint128_t hash, uint64_t count, bool lock=false,bool spin=false)
{
	_onDiskMQF<bitsPerSlot> *qf=this;
	if(qf->metadata->maximum_count!=0){
		count=std::min(count,qf->metadata->maximum_count);
	}
	uint64_t hash_remainder           = hash & BITMASK(qf->metadata->key_remainder_bits);
	uint64_t hash_bucket_index        = hash >> qf->metadata->key_remainder_bits;
	uint64_t hash_bucket_block_offset = hash_bucket_index % SLOTS_PER_BLOCK;
	/*uint64_t hash_bucket_lock_offset  = hash_bucket_index % NUM_SLOTS_TO_LOCK;*/
	//printf("index= %lu remainder= %lu count=%lu\n",hash_bucket_index,hash_remainder,count);
	if(hash_bucket_index > qf->metadata->xnslots){
		throw std::out_of_range("Insert is called with hash index out of range");
	}
	if (lock) {
		if(qf->mem->general_lock)
			return false;
		if (!spin_lock( hash_bucket_index, spin, false))
			return false;
	}

	uint64_t runend_index             = run_end( hash_bucket_index);

	/* Empty slot */
	if (might_be_empty( hash_bucket_index) && runend_index ==
			hash_bucket_index) {
		METADATA_WORD( runends, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);
		//set_slot( hash_bucket_index, hash_remainder);
		//set_fixed_counter( hash_bucket_index, 0);

		super_set(hash_bucket_index,hash_remainder,0);
		METADATA_WORD( occupieds, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);

		modify_metadata( &qf->metadata->ndistinct_elts, 1);
		modify_metadata( &qf->metadata->noccupied_slots, 1);
		/*modify_metadata( &qf->metadata->nelts, 1);*/
		/* This trick will, I hope, keep the fast case fast. */
		if (count > 1) {
			_insert( hash, count - 1, false, false);
		}
	} else { /* Non-empty slot */
		uint64_t new_values[67];
		uint64_t new_fcounters[67];
		uint64_t total_remainders;
		int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(
																																	hash_bucket_index
																																	- 1) + 1;

		if (!is_occupied( hash_bucket_index)) { /* Empty bucket, but its slot is occupied. */
			uint64_t *p = encode_counter( hash_remainder, count, &new_values[67],&new_fcounters[67]);
			total_remainders=&new_values[67] - p;
			insert_replace_slots_and_shift_remainders_and_runends_and_offsets(
																																				0,
																																				hash_bucket_index,
																																				runstart_index,
																																				p,
																																				&new_fcounters[67]-total_remainders,
																																				&new_values[67] - p,
																																				0);
			modify_metadata( &qf->metadata->ndistinct_elts, 1);
		} else { /* Non-empty bucket */

			uint64_t current_remainder, current_count, current_end;

			/* Find the counter for this remainder, if one exists. */
			current_end = decode_counter( runstart_index, &current_remainder,
																	 &current_count);
			while (current_remainder < hash_remainder && current_end!=runend_index) {
				runstart_index = current_end + 1;
				current_end = decode_counter( runstart_index, &current_remainder,
																		 &current_count);
			}

			/* If we reached the end of the run w/o finding a counter for this remainder,
				 then append a counter for this remainder to the run. */
			if (current_remainder < hash_remainder) {
				uint64_t *p = encode_counter( hash_remainder, count, &new_values[67],&new_fcounters[67]);
				total_remainders=&new_values[67] - p;
				insert_replace_slots_and_shift_remainders_and_runends_and_offsets(
																																					1, /* Append to bucket */
																																					hash_bucket_index,
																																					current_end + 1,
																																					p,
																																					&new_fcounters[67]-total_remainders,
																																					&new_values[67] - p,
																																					0);
				modify_metadata( &qf->metadata->ndistinct_elts, 1);
				/* Found a counter for this remainder.  Add in the new count. */
			} else if (current_remainder == hash_remainder) {
				uint64_t tmp= current_count + count;
				if(qf->metadata->maximum_count!=0){
					tmp=std::min(tmp,qf->metadata->maximum_count);
				}

				uint64_t *p = encode_counter( hash_remainder, tmp, &new_values[67],&new_fcounters[67]);
				total_remainders=&new_values[67] - p;
				insert_replace_slots_and_shift_remainders_and_runends_and_offsets(
																																					current_end==runend_index ? 1 : 2,
																																					hash_bucket_index,
																																					runstart_index,
																																					p,
																																					&new_fcounters[67]-total_remainders,
																																					&new_values[67] - p,
																																					current_end - runstart_index + 1);
				/* No counter for this remainder, but there are larger
					 remainders, so we're not appending to the bucket. */
			} else {
				uint64_t *p = encode_counter( hash_remainder, count, &new_values[67],&new_fcounters[67]);
				total_remainders=&new_values[67] - p;
				insert_replace_slots_and_shift_remainders_and_runends_and_offsets(
																																					2, /* Insert to bucket */
																																					hash_bucket_index,
																																					runstart_index,
																																					p,
																																					&new_fcounters[67]-total_remainders,
																																					&new_values[67] - p,
																																					0);
				modify_metadata( &qf->metadata->ndistinct_elts, 1);
			}
		}
		METADATA_WORD( occupieds, hash_bucket_index) |= 1ULL << (hash_bucket_block_offset % 64);

		/*modify_metadata( &qf->metadata->nelts, count);*/
	}

	if (lock) {
		unlock( hash_bucket_index, false);
	}

	return true;
}

template<uint64_t bitsPerSlot>
bool _onDiskMQF<bitsPerSlot>::remove(uint64_t hash, uint64_t count , bool lock, bool spin)
{
	_onDiskMQF<bitsPerSlot>* qf=this;
	uint64_t hash_remainder           = hash & BITMASK(qf->metadata->key_remainder_bits);
	uint64_t hash_bucket_index        = hash >> qf->metadata->key_remainder_bits;
	uint64_t current_remainder, current_count, current_end;
	uint64_t new_values[67];
	uint64_t new_fcounters[67];

	if(hash_bucket_index > qf->metadata->xnslots){
		throw std::out_of_range("Remove function is called with hash index out of range");
	}


	/* Empty bucket */
	if (!is_occupied( hash_bucket_index)){
		return true;
	}

	uint64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end( hash_bucket_index - 1) + 1;
	uint64_t original_runstart_index = runstart_index;
	int only_item_in_the_run = 0;

	/*Find the counter for this remainder, if one exists.*/
	current_end = decode_counter( runstart_index, &current_remainder, &current_count);
	while (current_remainder < hash_remainder && !is_runend( current_end)) {
		runstart_index = current_end + 1;
		current_end = decode_counter( runstart_index, &current_remainder, &current_count);
	}
	/* remainder not found in the given run */
	if (current_remainder != hash_remainder){
		return true;
	}

	if (original_runstart_index == runstart_index && is_runend( current_end))
		only_item_in_the_run = 1;


	/* endode the new counter */
	uint64_t *p = encode_counter( hash_remainder,
															 count>current_count? 0 : current_count-count,
															 &new_values[67],&new_fcounters[67]);

	uint64_t total_reminders=&new_values[67] - p;
	// if(fcounter==0 && newcount==0){
	// 	total_reminders=0;
	// 	p=&new_values[67];
	// }
	if (lock) {
		if(qf->mem->general_lock)
			return false;
		if (!spin_lock( hash_bucket_index, spin, false))
		return false;
	}

	remove_replace_slots_and_shift_remainders_and_runends_and_offsets(
																																		only_item_in_the_run,
																																		hash_bucket_index,
																																		runstart_index,
																																		p,
																																		&new_fcounters[67]-total_reminders,
																																		total_reminders,
																																		current_end - runstart_index + 1);


  if (lock) {
  unlock( hash_bucket_index, true);
  }

	return true;
	// update the nelements.
	/*modify_metadata( &qf->metadata->nelts, -count);*/
	/*qf->metadata->nelts -= count;*/
}

/***********************************************************************
 * Code that uses the above to implement key-value-counter operations. *
 ***********************************************************************/
 void onDiskMQF::init( onDiskMQF *&qf, uint64_t nslots, uint64_t key_bits, uint64_t tag_bits,uint64_t fixed_counter_size ,const char * path){
	 uint64_t qbits=(uint64_t)log2((double)nslots);
	 uint64_t tmpslotsSize=key_bits-qbits+fixed_counter_size;
	 switch (tmpslotsSize) {
		 case 1:
		 		qf=new  onDiskMQF_Namespace::_onDiskMQF<1>(nslots,key_bits,tag_bits,fixed_counter_size,path);
	 			break;
		 case 2:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<2>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 3:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<3>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 4:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<4>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 5:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<5>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 6:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<6>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 7:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<7>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 8:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<8>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 9:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<9>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 10:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<10>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 11:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<11>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 12:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<12>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 13:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<13>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 14:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<14>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 15:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<15>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 16:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<16>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 17:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<17>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 18:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<18>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 19:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<19>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 20:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<20>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 21:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<21>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 22:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<22>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 23:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<23>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 24:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<24>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 25:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<25>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 26:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<26>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 27:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<27>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 28:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<28>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 29:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<29>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 30:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<30>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 31:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<31>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 32:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<32>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 33:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<33>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 34:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<34>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 35:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<35>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 36:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<36>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 37:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<37>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 38:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<38>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 39:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<39>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 40:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<40>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 41:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<41>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 42:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<42>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 43:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<43>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 44:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<44>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 45:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<45>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 46:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<46>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 47:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<47>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 48:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<48>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 49:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<49>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 50:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<50>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 51:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<51>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 52:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<52>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 53:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<53>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 54:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<54>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 55:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<55>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 56:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<56>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 57:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<57>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 58:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<58>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 59:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<59>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 60:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<60>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 61:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<61>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 62:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<62>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 63:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<63>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;
		 case 64:
			 qf=new  onDiskMQF_Namespace::_onDiskMQF<64>(nslots,key_bits,tag_bits,fixed_counter_size,path);
			 break;

	 }
	 //qf->insert(100,1,false,false);
	 //cout<<qf->count_key(100)<<endl;
 }

template<uint64_t bitsPerSlot>
_onDiskMQF<bitsPerSlot>::_onDiskMQF( uint64_t nslots, uint64_t key_bits, uint64_t tag_bits,uint64_t fixed_counter_size ,const char * path)
{
	//qf=(QF*)calloc(sizeof(QF),1);
	uint64_t num_slots, xnslots, nblocks;
	uint64_t key_remainder_bits, bits_per_slot;
	uint64_t size;


	if(popcnt(nslots) != 1){
		throw std::domain_error("nslots must be a power of 2");

	}
	num_slots = nslots;

	xnslots = nslots + 10*sqrt((double)nslots);
	nblocks = (xnslots + SLOTS_PER_BLOCK - 1) / SLOTS_PER_BLOCK;
	key_remainder_bits = key_bits;
	while (nslots > 1) {
		//assert(key_remainder_bits > 0);
		key_remainder_bits--;
		nslots >>= 1;
	}

	bits_per_slot = key_remainder_bits+fixed_counter_size+tag_bits ;
  if(bits_per_slot!=bitsPerSlot){
    cout<<"OnDisk MQF created wrong.bits per slot calculated form the parameters is not the same as the template parameter."<<endl;
  }
	size = nblocks * (sizeof(qfblock) + (8 * bits_per_slot )) ;
	stxxlBufferSize= (uint64_t)((double)(size)*0.2/(1024.0*1024.0));
	stxxlBufferSize=max(stxxlBufferSize,(uint64_t)16);
	mem = (qfmem *)calloc(sizeof(qfmem), 1);
	metadata = (qfmetadata *)calloc(sizeof(qfmetadata), 1);
	metadata->mem=false;
	metadata->size = size;
	metadata->seed = 2038074761;
	metadata->nslots = num_slots;
	metadata->xnslots = metadata->nslots +
		10*sqrt((double)metadata->nslots);
	metadata->key_bits = key_bits;
	metadata->tag_bits = tag_bits;
	metadata->fixed_counter_size = fixed_counter_size;
	metadata->key_remainder_bits = key_remainder_bits;
	metadata->bits_per_slot = bits_per_slot;

	metadata->range = metadata->nslots;
	metadata->range <<= metadata->key_remainder_bits;
	metadata->nblocks = (metadata->xnslots + SLOTS_PER_BLOCK - 1) /
		SLOTS_PER_BLOCK;
	metadata->nelts = 0;
	metadata->ndistinct_elts = 0;
	metadata->noccupied_slots = 0;
	metadata->maximum_occupied_slots=(uint64_t)((double)metadata->xnslots *0.95);
	metadata->num_locks = (metadata->xnslots/NUM_SLOTS_TO_LOCK)+2;
	metadata->maximum_count = 0;
	metadata->tags_map=NULL;




	/* initialize all the locks to 0 */
	mem->metadata_lock = 0;
	mem->general_lock = 0;
	mem->locks = (volatile int *)calloc(metadata->num_locks,
																					sizeof(volatile int));


	blocks=stxxl::vector<onDisk_qfblock<bitsPerSlot> >(metadata->nblocks,stxxlBufferSize/16);
	// for(uint64_t i=0;i<metadata->nblocks;i++)
	// 	blocks[i]=onDisk_qfblock<bitsPerSlot>();
  //
	// qf->diskParams=new diskParameters();
	// diskParameters* diskParams=qf->diskParams;
  //
	// diskParams->nBlocksPerPointer=8;
	// diskParams->nBlocksPerIOBatch=1024;
	// nslotsInMemory=qf->metadata->xnslots;
	// diskParams->nBlocksOnMemory=(nslotsInMemory/64)+1;
  //
  //
	// //make nBlocksPerIOBatch divisable by nBlocksPerPointer for code simplifications
	// diskParams->nBlocksPerIOBatch=(diskParams->nBlocksPerIOBatch / diskParams->nBlocksPerPointer)
	// 															* diskParams->nBlocksPerPointer;
  //
	// if(diskParams->nBlocksOnMemory<diskParams->nBlocksPerIOBatch)
	// {
	// 	diskParams->nBlocksOnMemory=diskParams->nBlocksPerIOBatch;
	// }
  //
  //
	// diskParams->nBlocksOnMemory=(diskParams->nBlocksOnMemory / diskParams->nBlocksPerIOBatch)
	// 															*diskParams->nBlocksPerIOBatch;
  //
	// diskParams->sizeOnMemory = diskParams->nBlocksOnMemory * (sizeof(qfblock) + (8 * bits_per_slot )) ;
	// //qf->blocks = (qfblock *)calloc(diskParams->sizeOnMemory, 1);
	// qf->blocksFilePos=(uint64_t*)calloc(diskParams->nBlocksOnMemory,sizeof(uint64_t));
  //
	// //qf->reverseBlocksPointer=(uint64_t*)calloc(diskParams->nBlocksOnMemory,sizeof(uint64_t));
  //
	// //make nBlocksonMemoy divisable by nBlocksPerIOBatch for code simplifications
	// diskParams->blocksPointersLen=diskParams->nBlocksOnMemory/diskParams->nBlocksPerPointer;
	// qf->blocksPointers=new qfblock*[diskParams->blocksPointersLen];
	// diskParams->memoryBufferPos=0;
	// for(uint64_t i=0;i<diskParams->blocksPointersLen;i++)
	// {
	// 	//qf->blocksPointers[i]=(qfblock*)((char*)qf->blocks
	// 	//+(i*diskParams->nBlocksPerPointer*(sizeof(qfblock) + (8 * bits_per_slot ))));
  //
	// 	qf->blocksPointers[i]=NULL;
	// }
  //
	// qf->diskMQFStream.open(path, ios::in | ios::out | ios::binary);
	// qf->diskMQFStream.seekp(size+sizeof(qfmetadata));
	// qf->diskMQFStream<<0;
  //
  //




}

/* The caller should call onDiskMQF_init on the dest QF before calling this function.
 */
template<uint64_t bitsPerSlot>
void _onDiskMQF<bitsPerSlot>::copy(onDiskMQF *dest)
{
	// _onDiskMQF<bitsPerSlot>* qf=this;
	throw std::logic_error("not implemented yet");
	// memcpy(dest->mem, src->mem, sizeof(qfmem));
	// memcpy(dest->metadata, src->metadata, sizeof(qfmetadata));
	// //memcpy(dest->blocks, src->blocks, src->metadata->size);
  //
	// if(src->metadata->tags_map!=NULL){
	// 	dest->metadata->tags_map=
	// 	new std::map<uint64_t, std::vector<int> >(*src->metadata->tags_map);
	// }
}

/* free up the memory if the QF is in memory.
 * else unmap the mapped memory from pagecache.
 *
 * It does not delete the file on disk for on-disk QF.
 */
template<uint64_t bitsPerSlot>
_onDiskMQF<bitsPerSlot>::~_onDiskMQF()
{
	_onDiskMQF<bitsPerSlot>* qf=this;
	//assert(qf->blocks != NULL);

	qf->metadata->noccupied_slots=0;
	if(qf->metadata->tags_map!=NULL){
		delete qf->metadata->tags_map;
		qf->metadata->tags_map=NULL;
	}
	if (qf->metadata->mem) {
		free(qf->mem);
		free(qf->metadata);
		qf->blocks.clear();
		//free(qf->blocks);
	} else {
	//msync(qf->metadata, qf->metadata->size + sizeof(qfmetadata),MS_SYNC);
//	munmap(qf->metadata, qf->metadata->size + sizeof(qfmetadata));
	//close(qf->mem->fd);
	}

}

// template<uint64_t bitsPerSlot>
// void _onDiskMQF<bitsPerSlot>::close()
// {
// 	_onDiskMQF<bitsPerSlot>* qf=this;
// 	//assert(qf->blocks != NULL);
// 	munmap(qf->metadata, qf->metadata->size + sizeof(qfmetadata));
// 	close(qf->mem->fd);
// }

/*
 * Will read the on-disk QF using mmap.
 * Data won't be copied in memory.
 *
 */
 template<uint64_t bitsPerSlot>
 void _onDiskMQF<bitsPerSlot>::read( const char *path)
 {
	 throw std::logic_error("not implemented yet");
	 _onDiskMQF<bitsPerSlot>* qf=this;
	 struct stat sb;
	 int ret;

	 qf->mem = (qfmem *)calloc(sizeof(qfmem), 1);
	 qf->mem->fd = open(path, O_RDWR, S_IRWXU);
	 if (qf->mem->fd < 0) {
		 perror("Couldn't open file:\n");
		 exit(EXIT_FAILURE);
	 }

	 ret = fstat (qf->mem->fd, &sb);
	 if ( ret < 0) {
		 perror ("fstat");
		 exit(EXIT_FAILURE);
	 }

	 if (!S_ISREG (sb.st_mode)) {
		 fprintf (stderr, "%s is not a file\n", path);
		 exit(EXIT_FAILURE);
	 }

	 qf->metadata = (qfmetadata *)mmap(NULL, sb.st_size, PROT_READ | PROT_WRITE, MAP_SHARED,
	 qf->mem->fd, 0);
	 qf->metadata->mem=false;
		 //qf->blocks = (qfblock *)(qf->metadata + 1);
		 qf->metadata->num_locks = (qf->metadata->xnslots/NUM_SLOTS_TO_LOCK)+2;
		 qf->mem->metadata_lock = 0;
		 qf->mem->locks = (volatile int *)calloc(qf->metadata->num_locks,
			 sizeof(volatile int));

	string tagsMapOutName=string(path)+".tags_map";
	if(file_exists(tagsMapOutName)){
		qf->metadata->tags_map=load_tags_map(tagsMapOutName.c_str());
	}

}

template<uint64_t bitsPerSlot>
void _onDiskMQF<bitsPerSlot>::reset()
{
	//assert(popcnt(nslots) == 1); /* nslots must be a power of 2 */
	_onDiskMQF<bitsPerSlot> *qf=this;
	throw std::logic_error("not implemented yet");
	qf->metadata->nelts = 0;
	qf->metadata->ndistinct_elts = 0;
	qf->metadata->noccupied_slots = 0;
	if(qf->metadata->tags_map!=NULL)
		qf->metadata->tags_map->clear();

	//memset(qf->blocks, 0, qf->metadata->nblocks*(sizeof(qfblock) + SLOTS_PER_BLOCK *
	//																	 qf->metadata->bits_per_slot / 8));

}

template<uint64_t bitsPerSlot>
void _onDiskMQF<bitsPerSlot>::serialize(const char *filename)
{
	throw std::logic_error("not implemented yet");
	_onDiskMQF<bitsPerSlot>* qf=this;
	FILE *fout;
	fout = fopen(filename, "wb+");
	if (fout == NULL) {
		perror("Error opening file for serializing\n");
		exit(EXIT_FAILURE);
	}
	fwrite(qf->metadata, sizeof(qfmetadata), 1, fout);
	/* we don't serialize the locks */
	//fwrite(qf->blocks, qf->metadata->size, 1, fout);
	fclose(fout);

	if(qf->metadata->tags_map!=NULL)
	{
		string tagsMapOutName=string(filename)+".tags_map";
		save_tags_map(qf->metadata->tags_map,tagsMapOutName.c_str());
	}
}


template<uint64_t bitsPerSlot>
void _onDiskMQF<bitsPerSlot>::deserialize(const char *filename)
{
	throw std::logic_error("not implemented yet");
	_onDiskMQF<bitsPerSlot>* qf=this;
	FILE *fin;
	fin = fopen(filename, "rb");
	if (fin == NULL) {
		perror("Error opening file for deserializing\n");
		exit(EXIT_FAILURE);
	}

	qf->mem = (qfmem *)calloc(sizeof(qfmem), 1);
	qf->metadata = (qfmetadata *)calloc(sizeof(qfmetadata), 1);

	fread(qf->metadata, sizeof(qfmetadata), 1, fin);

	/* initlialize the locks in the QF */
	qf->metadata->num_locks = (qf->metadata->xnslots/NUM_SLOTS_TO_LOCK)+2;
	qf->mem->metadata_lock = 0;
	/* initialize all the locks to 0 */
	qf->mem->locks = (volatile int *)calloc(qf->metadata->num_locks, sizeof(volatile int));

	//qf->blocks = (qfblock *)calloc(qf->metadata->size, 1);
	//fread(qf->blocks, qf->metadata->size, 1, fin);
	fclose(fin);

	string tagsMapOutName=string(filename)+".tags_map";
	if(file_exists(tagsMapOutName)){
		qf->metadata->tags_map=load_tags_map(tagsMapOutName.c_str());
	}


}
template<uint64_t bitsPerSlot>
uint64_t _onDiskMQF<bitsPerSlot>::add_tag(uint64_t key, uint64_t tag, bool lock, bool spin)
{
	_onDiskMQF<bitsPerSlot>* qf=this;
	if(qf->metadata->tag_bits==0){
		return 0;
	}
	__uint128_t hash = key;
	uint64_t hash_remainder   = hash & BITMASK(qf->metadata->key_remainder_bits);
	int64_t hash_bucket_index = hash >> qf->metadata->key_remainder_bits;


	if (!is_occupied( hash_bucket_index)){
		return 0;
	}

	int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(
																																hash_bucket_index-1)
		+ 1;

	if (runstart_index < hash_bucket_index)
		runstart_index = hash_bucket_index;

	/* printf("MC RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

	uint64_t current_remainder, current_count, current_end;
	do {
		current_end = decode_counter( runstart_index, &current_remainder,
																 &current_count);
		if (current_remainder == hash_remainder){
			if (lock) {
				if(qf->mem->general_lock)
					return false;
				if (!spin_lock( runstart_index, spin, false))
				return 0;
			}

			set_tag(runstart_index,tag);
			if (lock) {
				unlock( runstart_index, true);
			}


			return 1;
		}
		runstart_index = current_end + 1;
	} while (!is_runend( current_end));


	return 0;
}

template<uint64_t bitsPerSlot>
uint64_t _onDiskMQF<bitsPerSlot>::remove_tag(uint64_t key ,bool lock, bool spin)
{
_onDiskMQF<bitsPerSlot>* qf=this;
	if(qf->metadata->tag_bits==0){
		return 0;
	}

	__uint128_t hash = key;
	uint64_t hash_remainder   = hash & BITMASK(qf->metadata->key_remainder_bits);
	uint64_t hash_bucket_index = hash >> qf->metadata->key_remainder_bits;
	if(hash_bucket_index > qf->metadata->xnslots){
			throw std::out_of_range("onDiskMQF_remove_tag is called with hash index out of range");
		}

	if (!is_occupied( hash_bucket_index)){
		return 0;
	}
	uint64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(
																																hash_bucket_index-1)
		+ 1;
	if (runstart_index < hash_bucket_index)
		runstart_index = hash_bucket_index;

	/* printf("MC RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

	uint64_t current_remainder, current_count, current_end;
	do {
		current_end = decode_counter( runstart_index, &current_remainder,
																 &current_count);
		if (current_remainder == hash_remainder){
			if (lock) {
				if(qf->mem->general_lock)
					return false;
				if (!spin_lock( runstart_index, spin, false))
					return false;
				}
			set_tag(runstart_index,0);
			if (lock)
				unlock( runstart_index, true);
			return 1;
		}
		runstart_index = current_end + 1;
	} while (!is_runend( current_end));


	return 0;
}
template<uint64_t bitsPerSlot>
uint64_t _onDiskMQF<bitsPerSlot>::get_tag(uint64_t key)
{
	_onDiskMQF<bitsPerSlot>* qf=this;
	if(qf->metadata->tag_bits==0){
		return 0;
	}
	__uint128_t hash = key;
	uint64_t hash_remainder   = hash & BITMASK(qf->metadata->key_remainder_bits);
	uint64_t hash_bucket_index = hash >> qf->metadata->key_remainder_bits;
	if(hash_bucket_index > qf->metadata->xnslots){
			throw std::out_of_range("onDiskMQF_get_tag is called with hash index out of range");
		}
	if (!is_occupied( hash_bucket_index))
		return 0;

	uint64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(
																																hash_bucket_index-1)
		+ 1;
	if (runstart_index < hash_bucket_index)
		runstart_index = hash_bucket_index;

	/* printf("MC RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

	uint64_t current_remainder, current_count, current_end;
	do {
		current_end = decode_counter( runstart_index, &current_remainder,
																 &current_count);
		if (current_remainder == hash_remainder){
			return _get_tag(runstart_index);
		}
		runstart_index = current_end + 1;
	} while (!is_runend( current_end));

	return 0;
}

template<uint64_t bitsPerSlot>
bool _onDiskMQF<bitsPerSlot>::insert(uint64_t key, uint64_t count, bool
							 lock, bool spin)
{
	if(count==0)
	{
		return true;
	}
	/*uint64_t hash = (key << qf->metadata->tag_bits) | (value & BITMASK(qf->metadata->tag_bits));*/
	if (count == 1)
	 return insert1( key, lock, spin);
	else
	 return _insert( key, count, lock, spin);
}

template<uint64_t bitsPerSlot>
uint64_t _onDiskMQF<bitsPerSlot>::count_key(uint64_t key)
{
	_onDiskMQF<bitsPerSlot>* qf=this;
	__uint128_t hash = key;
	uint64_t hash_remainder   = hash & BITMASK(qf->metadata->key_remainder_bits);
	int64_t hash_bucket_index = hash >> qf->metadata->key_remainder_bits;

	if (!is_occupied( hash_bucket_index))
		return 0;

	int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(
																																hash_bucket_index-1)
		+ 1;
	if (runstart_index < hash_bucket_index)
		runstart_index = hash_bucket_index;

	/* printf("MC RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

	uint64_t current_remainder, current_count, current_end;
	do {
		current_end = decode_counter( runstart_index, &current_remainder,
																 &current_count);
		if (current_remainder == hash_remainder){
			return current_count;
		}

		runstart_index = current_end + 1;
	} while (!is_runend( current_end));
	return 0;
}




/* initialize the iterator at the run corresponding
 * to the position index
 */
template<uint64_t bitsPerSlot>
bool _onDiskMQF<bitsPerSlot>::getIterator(onDiskMQFIterator<bitsPerSlot> *qfi, uint64_t position)
{
	_onDiskMQF<bitsPerSlot>* qf=this;
	if(position > qf->metadata->xnslots){
		throw std::out_of_range("onDiskMQF_iterator is called with position out of range");
	}
	if (!is_occupied( position)) {
		uint64_t block_index = position;
		uint64_t idx = bitselect(qf->get_block( block_index)->occupieds[0], 0);
		if (idx == 64) {
			while(idx == 64 && block_index < qf->metadata->nblocks) {
				block_index++;
				idx = bitselect(qf->get_block( block_index)->occupieds[0], 0);
			}
		}
		position = block_index * SLOTS_PER_BLOCK + idx;
	}

	qfi->qf = qf;
	qfi->num_clusters = 0;
	qfi->run = position;
	qfi->current = position == 0 ? 0 : run_end(qfi->qf, position-1) + 1;
	if (qfi->current < position)
		qfi->current = position;



	if (qfi->current >= qf->metadata->nslots)
		return false;
	return true;
}

template<uint64_t bitsPerSlot>
int onDiskMQFIterator<bitsPerSlot>::get(uint64_t *key, uint64_t *value, uint64_t *count)
{
	onDiskMQFIterator<bitsPerSlot> *qfi=this;
	if(qfi->current > qfi->qf->metadata->xnslots){
		throw std::out_of_range("onDiskMQFIterator_get is called with hash index out of range");
	}
	uint64_t current_remainder, current_count;
	decode_counter(qfi->qf, qfi->current, &current_remainder, &current_count);
	*key = (qfi->run << qfi->qf->metadata->key_remainder_bits) | current_remainder;
	*value = _get_tag(qfi->qf,qfi->current);   // for now we are not using value
	*count = current_count;

	qfi->qf->metadata->ndistinct_elts++;
	qfi->qf->metadata->nelts += current_count;

	/*qfi->current = end_index;*/ 		//get should not change the current index
																		//of the iterator
	return 0;
}
template<uint64_t bitsPerSlot>
int onDiskMQFIterator<bitsPerSlot>::next()
{
	onDiskMQFIterator<bitsPerSlot> *qfi=this;
	if (onDiskMQFIterator_end(qfi))
		return 1;
	else {
		/* move to the end of the current counter*/
		uint64_t current_remainder, current_count;
		qfi->current = decode_counter(qfi->qf, qfi->current, &current_remainder,
																	&current_count);

		if (!is_runend(qfi->qf, qfi->current)) {
			qfi->current++;

			if (qfi->current > qfi->qf->metadata->xnslots)
				return 1;
			return 0;
		}
		else {

			uint64_t block_index = qfi->run / SLOTS_PER_BLOCK;
			uint64_t rank = bitrank(qfi->qf->get_block(block_index)->occupieds[0],
															qfi->run % SLOTS_PER_BLOCK);
			uint64_t next_run = bitselect(qfi->qf->get_block(block_index)->occupieds[0],
																		rank);
			if (next_run == 64) {
				rank = 0;
				while (next_run == 64 && block_index < qfi->qf->metadata->nblocks) {
					block_index++;
					next_run = bitselect(qfi->qf->get_block(block_index)->occupieds[0],
															 rank);
				}
			}
			if (block_index == qfi->qf->metadata->nblocks) {
				/* set the index values to max. */
				qfi->run = qfi->current = qfi->qf->metadata->xnslots;
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

template<uint64_t bitsPerSlot>
inline int onDiskMQFIterator<bitsPerSlot>::end()
{
	onDiskMQFIterator<bitsPerSlot> *qfi=this;
	if (qfi->current >= qfi->qf->metadata->xnslots /*&& is_runend(qfi->qf, qfi->current)*/)
		return 1;
	else
		return 0;
}


void unionFn(uint64_t  key_a, uint64_t  tag_a,uint64_t  count_a,
					   uint64_t  key_b, uint64_t  tag_b,uint64_t  count_b,
					   uint64_t *key_c, uint64_t *tag_c,uint64_t *count_c)
{
		if(count_a==0){
			*key_c=key_b;
			*tag_c=tag_a;
			*count_c=count_b;
		}
		else if(count_b==0){
			*key_c=key_a;
			*tag_c=tag_a;
			*count_c=count_a;
		}
		else{
			*key_c=key_a;
			*tag_c=tag_a;
			*count_c=count_a+count_b;
		}

}
void intersectFn(uint64_t  key_a, uint64_t  tag_a,uint64_t  count_a,
					   uint64_t  key_b, uint64_t  tag_b,uint64_t  count_b,
					   uint64_t *key_c, uint64_t *tag_c,uint64_t *count_c)
{
	*key_c=0;
	*tag_c=0;
	*count_c=0;
	if(count_a!=0 && count_b!=0){
			*key_c=key_a;
			*tag_c=tag_a;
			*count_c=std::min(count_a,count_b);
	}

}

void subtractFn(uint64_t  key_a, uint64_t  tag_a,uint64_t  count_a,
					   uint64_t  key_b, uint64_t  tag_b,uint64_t  count_b,
					   uint64_t *key_c, uint64_t *tag_c,uint64_t *count_c)
{
		if(count_b==0){
			*key_c=key_a;
			*tag_c=tag_a;
			*count_c=count_a;
		}
		else{
			*key_c=key_a;
			*tag_c=tag_a;
			*count_c=count_a<count_b ? 0:count_a-count_b;
		}

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
// void _onDiskMQF_merge(onDiskMQF *qfa, onDiskMQF *qfb, onDiskMQF *qfc,
// 	void(*mergeFn)(uint64_t   keya, uint64_t  tag_a,uint64_t  count_a,
// 						  	 uint64_t   keyb, uint64_t  tag_b,uint64_t  count_b,
// 							   uint64_t*  keyc, uint64_t* tag_c,uint64_t* count_c
// 							 ))
// {
// 	onDiskMQFIterator qfia, qfib;
// 	if(qfa->metadata->range != qfb->metadata->range ||
// 	qfb->metadata->range != qfc->metadata->range )
// 	{
// 		throw std::logic_error("Merging non compatible filters");
// 	}
// 	onDiskMQF_iterator(qfa, &qfia, 0);
// 	onDiskMQF_iterator(qfb, &qfib, 0);
//
// 	uint64_t keya, taga, counta, keyb, tagb, countb;
// 	uint64_t keyc,tagc, countc;
// 	onDiskMQFIterator_get(&qfia, &keya, &taga, &counta);
// 	onDiskMQFIterator_get(&qfib, &keyb, &tagb, &countb);
//
// 	do {
// 		if (keya < keyb) {
// 			mergeFn(keya,taga,counta,0,0,0,&keyc,&tagc,&countc);
// 			onDiskMQFIterator_next(&qfia);
// 			onDiskMQFIterator_get(&qfia, &keya, &taga, &counta);
// 		}
// 		else if(keya > keyb) {
// 			mergeFn(0,0,0,keyb,tagb,countb,&keyc,&tagc,&countc);
// 			onDiskMQFIterator_next(&qfib);
// 			onDiskMQFIterator_get(&qfib, &keyb, &tagb, &countb);
// 		}
// 		else{
// 			mergeFn(keya,taga,counta,keyb,tagb,countb,&keyc,&tagc,&countc);
// 			onDiskMQFIterator_next(&qfia);
// 			onDiskMQFIterator_next(&qfib);
// 			onDiskMQFIterator_get(&qfia, &keya, &taga, &counta);
// 			onDiskMQFIterator_get(&qfib, &keyb, &tagb, &countb);
// 		}
// 		if(countc!=0){
// 			onDiskMQF_insert(qfc, keyc, countc, true, true);
// 			onDiskMQF_add_tag(qfc,keya,tagc);
// 		}
//
// 	} while(!onDiskMQFIterator_end(&qfia) && !onDiskMQFIterator_end(&qfib));
//
// 	if (!onDiskMQFIterator_end(&qfia)) {
//
// 		do {
// 			onDiskMQFIterator_get(&qfia, &keya, &taga, &counta);
// 			mergeFn(keya,taga,counta,0,0,0,&keyc,&tagc,&countc);
// 			if(countc!=0){
// 				onDiskMQF_insert(qfc, keyc, countc, true, true);
// 				onDiskMQF_add_tag(qfc,keyc,tagc);
// 			}
// 		} while(!onDiskMQFIterator_next(&qfia));
// 	}
//
// 	if (!onDiskMQFIterator_end(&qfib)) {
// 		do {
// 			onDiskMQFIterator_get(&qfib, &keyb, &tagb, &countb);
// 			mergeFn(0,0,0,keyb,tagb,countb,&keyc,&tagc,&countc);
// 			if(countc!=0){
// 				onDiskMQF_insert(qfc, keyc, countc, true, true);
// 				onDiskMQF_add_tag(qfc,keyc,tagc);
// 			}
// 		} while(!onDiskMQFIterator_next(&qfib));
// 	}
//
// 	return;
// }
// void onDiskMQF_merge(onDiskMQF *qfa, onDiskMQF *qfb, onDiskMQF *qfc)
// {
// 	_onDiskMQF_merge(qfa,qfb,qfc,unionFn);
// }
//
// void onDiskMQF_intersect(onDiskMQF *qfa, onDiskMQF *qfb, onDiskMQF *qfc)
// {
// 	_onDiskMQF_merge(qfa,qfb,qfc,intersectFn);
// }
// void onDiskMQF_subtract(onDiskMQF *qfa, onDiskMQF *qfb, onDiskMQF *qfc)
// {
// 	_onDiskMQF_merge(qfa,qfb,qfc,subtractFn);
// }
//
// bool onDiskMQF_equals(onDiskMQF *qfa, onDiskMQF *qfb)
// {
// 	onDiskMQFIterator qfia, qfib;
// 	if(qfa->metadata->range != qfb->metadata->range  )
// 	{
// 		throw std::logic_error("comparing non compatible filters");
// 	}
// 	onDiskMQF_iterator(qfa, &qfia, 0);
// 	onDiskMQF_iterator(qfb, &qfib, 0);
//
// 	uint64_t keya, valuea, counta, keyb, valueb, countb;
// 	onDiskMQFIterator_get(&qfia, &keya, &valuea, &counta);
// 	onDiskMQFIterator_get(&qfib, &keyb, &valueb, &countb);
// 	do {
// 		if (keya != keyb) {
// 			return false;
// 		}
// 		else {
// 			if(counta!=countb || valuea != valueb)
// 			{
// 				return false;
// 			}
// 			onDiskMQFIterator_next(&qfib);
// 			onDiskMQFIterator_next(&qfia);
// 			onDiskMQFIterator_get(&qfia, &keya, &valuea, &counta);
// 			onDiskMQFIterator_get(&qfib, &keyb, &valueb, &countb);
// 		}
// 	} while(!onDiskMQFIterator_end(&qfia) && !onDiskMQFIterator_end(&qfib));
//
// 	if (!onDiskMQFIterator_end(&qfia) || !onDiskMQFIterator_end(&qfib)) {
// 		return false;
// 	}
//
// 	return true;
// }
//
//
//
// std::map<std::string, uint64_t> Tags_map;
// uint64_t last_index=0;
// void union_multi_Fn(uint64_t   key_arr[], uint64_t  tag_arr[],uint64_t  count_arr[]
// 	,std::map<uint64_t, std::vector<int> > ** inverted_indexes,int nqf,
// 							 uint64_t*  key_c, uint64_t* tag_c,uint64_t* count_c)
// {
//
// 	*count_c=0;
// 	for(int i=0;i<nqf;i++)
// 	{
// 		//printf("key =%lu, count=%lu\n", key_arr[i],count_arr[i]);
// 		if(count_arr[i]!=0)
// 		{
// 			*key_c=key_arr[i];
// 			*tag_c=tag_arr[i];
// 			*count_c+=count_arr[i];
// 		}
// 	}
//
// }
//
// void _onDiskMQF_multi_merge(onDiskMQF *onDiskMQF_arr[],int nqf, onDiskMQF *qfr,
// 	void(*mergeFn)(uint64_t   key_arr[], uint64_t  tag_arr[],uint64_t  count_arr[],
// 								std::map<uint64_t, std::vector<int> > ** inverted_indexes,int nqf,
// 							   uint64_t*  keyc, uint64_t* tag_c,uint64_t* count_c
// 							 ))
// {
// 	int i;
// 	uint64_t range=onDiskMQF_arr[0]->metadata->range;
// 	for (i=1; i<nqf; i++) {
// 		if(onDiskMQF_arr[i]->metadata->range!=range)
// 		{
// 			throw std::logic_error("Merging non compatible filters");
// 		}
// 	}
//
// 	onDiskMQFIterator *onDiskMQFIterator_arr[nqf];
//
// 	uint64_t smallest_key=UINT64_MAX,second_smallest_key;
// 	uint64_t keys[nqf];
// 	uint64_t tags[nqf];
// 	uint64_t counts[nqf];
// 	std::map<uint64_t, std::vector<int> > ** inverted_indexes=new std::map<uint64_t, std::vector<int> >*[nqf];
//
// 	for (i=0; i<nqf; i++) {
// 		onDiskMQFIterator_arr[i]=new onDiskMQFIterator();
// 		onDiskMQF_iterator(onDiskMQF_arr[i], onDiskMQFIterator_arr[i], 0);
// 		onDiskMQFIterator_get(onDiskMQFIterator_arr[i], &keys[i], &tags[i], &counts[i]);
// 		smallest_key=std::min(keys[i],smallest_key);
// 		inverted_indexes[i]=onDiskMQF_arr[i]->metadata->tags_map;
// 	}
//
// 	uint64_t keys_m[nqf];
// 	uint64_t tags_m[nqf];
// 	uint64_t counts_m[nqf];
//
//
// 	bool finish=false;
// 	while(!finish)
// 	{
// 		finish=true;
// 		second_smallest_key=UINT64_MAX;
// 		//printf("smallest_key = %llu\n",smallest_key );
// 		for(i=0;i<nqf;i++)
// 		{
// 			keys_m[i]=0;
// 			counts_m[i]=0;
// 			tags_m[i]=0;
//
// 			//printf(" key = %llu\n",keys[i]);
// 			if(keys[i]==smallest_key){
// 				keys_m[i]=keys[i];
// 				counts_m[i]=counts[i];
// 				tags_m[i]=tags[i];
// 				onDiskMQFIterator_next(onDiskMQFIterator_arr[i]);
// 				if(!onDiskMQFIterator_end(onDiskMQFIterator_arr[i]))
// 				{
// 					finish=false;
// 					onDiskMQFIterator_get(onDiskMQFIterator_arr[i], &keys[i], &tags[i], &counts[i]);
// 				}else{
// 					keys[i]=UINT64_MAX;
// 				}
// 			}
// 			second_smallest_key=std::min(second_smallest_key,keys[i]);
// 		}
// 		for (i = 0; i < nqf; i++) {
// 			if(keys[i]!=UINT64_MAX)
// 			{
// 				finish=false;
// 				break;
// 			}
// 		}
// 		//printf("second_smallest_key=%llu finish=%d\n",second_smallest_key,finish);
// 		uint64_t keyc,tagc, countc;
// 		mergeFn(keys_m,tags_m,counts_m,inverted_indexes,nqf,&keyc,&tagc,&countc);
//
// 		if(countc!=0){
// 			onDiskMQF_insert(qfr, keyc, countc, true, true);
// 			onDiskMQF_add_tag(qfr,keyc,tagc);
// 		}
// 		smallest_key=second_smallest_key;
// 	}
// 	// cout<<"before delete"<<endl;
// 	delete  inverted_indexes;
// 	for(i=0;i<nqf;i++)
// 	{
// 		delete onDiskMQFIterator_arr[i];
// 	}
//
// 	return;
// }
//
// /*
//  * Merge an array of qfs into the resultant QF
//  */
// void onDiskMQF_multi_merge(onDiskMQF *onDiskMQF_arr[], int nqf, onDiskMQF *qfr)
// {
// 	_onDiskMQF_multi_merge(onDiskMQF_arr,nqf,qfr,union_multi_Fn);
//
// }
//
// void inverted_union_multi_Fn(uint64_t   key_arr[], uint64_t  tag_arr[],uint64_t  count_arr[],
// 	std::map<uint64_t, std::vector<int> > ** inverted_indexes ,int nqf,
// 							 uint64_t*  key_c, uint64_t* tag_c,uint64_t* count_c)
// {
//
// 	std::string index_key="";
// 	*count_c=0;
// 	for(int i=0;i<nqf;i++)
// 	{
// 		if(count_arr[i]!=0)
// 		{
// 			*key_c=key_arr[i];
// 			*count_c+=count_arr[i];
// 			if(inverted_indexes==NULL){
// 				index_key+=std::to_string(i);
// 				index_key+=';';
// 			}
// 			else{
// 				auto it=inverted_indexes[i]->find(tag_arr[i]);
// 				for(auto k:it->second){
// 					index_key+=std::to_string(k);
// 					index_key+=';';
// 				}
// 			}
//
// 		}
// 	}
// 	index_key.pop_back();
// 	auto it=Tags_map.find(index_key);
// 	if(it==Tags_map.end())
// 	{
//
// 		Tags_map.insert(std::make_pair(index_key,Tags_map.size()));
// 		it=Tags_map.find(index_key);
// 	}
// 	*tag_c=it->second;
//
// }
//
// void inverted_union_multi_no_count_Fn(uint64_t   key_arr[], uint64_t  tag_arr[],uint64_t  count_arr[],
// 									std::map<uint64_t, std::vector<int> > ** inverted_indexes, int nqf,
// 							 uint64_t*  key_c, uint64_t* tag_c,uint64_t* count_c)
// {
// 	std::string index_key="";
// 	*count_c=0;
// 	for(int i=0;i<nqf;i++)
// 	{
// 		//printf("key =%lu, count=%lu\n", key_arr[i],count_arr[i]);
// 		if(count_arr[i]!=0)
// 		{
// 			*key_c=key_arr[i];
// 			if(inverted_indexes==NULL){
// 				index_key+=std::to_string(i);
// 				index_key+=';';
// 			}
// 			else{
// 				auto it=inverted_indexes[i]->find(tag_arr[i]);
// 				for(auto k:it->second){
// 					index_key+=std::to_string(k);
// 					index_key+=';';
// 				}
// 			}
//
// 		}
// 	}
// 	index_key.pop_back();
// 	auto it=Tags_map.find(index_key);
// 	if(it==Tags_map.end())
// 	{
//
// 		Tags_map.insert(std::make_pair(index_key,last_index));
// 		last_index++;
// 		it=Tags_map.find(index_key);
// 	}
// 	*count_c=it->second;
//
// }
//
//
// void onDiskMQF_invertable_merge(onDiskMQF *onDiskMQF_arr[], int nqf, onDiskMQF *qfr)
// {
// 	int i;
// 	int last_tag=0;
// 	Tags_map.clear();
// 	last_index=0;
// 	for(i=0;i<nqf;i++){
// 		if(onDiskMQF_arr[i]->metadata->tags_map==NULL){
// 			onDiskMQF_arr[i]->metadata->tags_map=new std::map<uint64_t, std::vector<int> >();
// 			vector<int> tmp(1);
// 			tmp[0]=last_index;
// 			Tags_map.insert(std::make_pair(std::to_string(i),last_index++));
// 			onDiskMQF_arr[i]->metadata->tags_map->insert(make_pair(0,tmp));
// 		}
// 		else{
// 			auto it=onDiskMQF_arr[i]->metadata->tags_map->begin();
// 			int updated_tags=0;
// 			while(it!=onDiskMQF_arr[i]->metadata->tags_map->end()){
// 				for(int j=0;j<it->second.size();j++){
// 					it->second[j]+=last_tag;
// 					auto it2=Tags_map.find(std::to_string(it->second[j]));
// 					if(it2==Tags_map.end()){
// 						Tags_map.insert(std::make_pair(std::to_string(it->second[j]),it->second[j]));
// 						updated_tags++;
// 					}
// 				}
// 				it++;
// 			}
// 		}
// 		last_tag+=Tags_map.size();
// 	}
//
//
//
// 	_onDiskMQF_multi_merge(onDiskMQF_arr,nqf,qfr,inverted_union_multi_Fn);
// 	qfr->metadata->tags_map=new std::map<uint64_t, std::vector<int> >();
// 	auto it=Tags_map.begin();
// 	while(it!=Tags_map.end()){
// 		std::vector<int> tmp=key_to_vector_int(it->first);
// 		qfr->metadata->tags_map->insert(std::make_pair(it->second,tmp));
// 		it++;
//
// 	}
//
//
// }
//
// void onDiskMQF_invertable_merge_no_count(onDiskMQF *onDiskMQF_arr[], int nqf, onDiskMQF *qfr)
// {
//
// 	int i;
// 	int last_tag=0;
// 	Tags_map.clear();
// 	last_index=0;
// 	for(i=0;i<nqf;i++){
// 		if(onDiskMQF_arr[i]->metadata->tags_map==NULL){
// 			onDiskMQF_arr[i]->metadata->tags_map=new std::map<uint64_t, std::vector<int> >();
// 			vector<int> tmp(1);
// 			tmp[0]=last_index;
// 			Tags_map.insert(std::make_pair(std::to_string(i),last_index++));
// 			onDiskMQF_arr[i]->metadata->tags_map->insert(make_pair(0,tmp));
// 		}
// 		else{
// 			auto it=onDiskMQF_arr[i]->metadata->tags_map->begin();
// 			int updated_tags=0;
// 			while(it!=onDiskMQF_arr[i]->metadata->tags_map->end()){
// 				for(int j=0;j<it->second.size();j++){
// 					it->second[j]+=last_tag;
// 					auto it2=Tags_map.find(std::to_string(it->second[j]));
// 					if(it2==Tags_map.end()){
// 						Tags_map.insert(std::make_pair(std::to_string(it->second[j]),it->second[j]));
// 						updated_tags++;
// 					}
// 				}
// 				it++;
// 			}
// 		}
// 		last_tag+=Tags_map.size();
// 	}
//
//
// 	_onDiskMQF_multi_merge(onDiskMQF_arr,nqf,qfr,inverted_union_multi_no_count_Fn);
//
// 	qfr->metadata->tags_map=new std::map<uint64_t, std::vector<int> >();
// 	auto it=Tags_map.begin();
// 	while(it!=Tags_map.end()){
// 		std::vector<int> tmp=key_to_vector_int(it->first);
// 		qfr->metadata->tags_map->insert(std::make_pair(it->second,tmp));
// 		it++;
// 	}
//
//
// }
//
// onDiskMQF* onDiskMQF_resize(onDiskMQF* qf, int newQ, const char * originalFilename, const char * newFilename)
// {
// 	throw std::logic_error("not implemented yet");
// 	if((int)qf->metadata->key_bits-newQ <2)
// 	{
// 		throw std::logic_error("Resize cannot be done. Slot size cannot be less than 2");
// 	}
//
// 	if(originalFilename)
// 	{
// 		onDiskMQF_serialize(originalFilename);
// 		onDiskMQF_destroy(qf);
// 		onDiskMQF_read(originalFilename);
// 	}
// 	onDiskMQF* newQF=(onDiskMQF *)calloc(sizeof(QF), 1);
// 	if(newFilename)
// 	{
// 		//onDiskMQF_init(newQF, (1ULL<<newQ),qf->metadata->key_bits, qf->metadata->tag_bits,qf->metadata->fixed_counter_size, false, newFilename, 2038074761);
// 	}
// 	else{
// 		//onDiskMQF_init(newQF, (1ULL<<newQ),qf->metadata->key_bits, qf->metadata->tag_bits,qf->metadata->fixed_counter_size, true, "" , 2038074761);
// 	}
// 	onDiskMQFIterator qfi;
// 	onDiskMQF_iterator( &qfi, 0);
//
//
// 	uint64_t keya, valuea, counta;
// 	onDiskMQFIterator_get(&qfi, &keya, &valuea, &counta);
//
// 	do {
// 			onDiskMQF_insert(newQF, keya, counta);
// 			onDiskMQF_add_tag(newQF,keya,valuea);
// 			onDiskMQFIterator_next(&qfi);
// 			onDiskMQFIterator_get(&qfi, &keya, &valuea, &counta);
// 	} while(!onDiskMQFIterator_end(&qfi));
// 	onDiskMQF_destroy(qf);
// 	return newQF;
//
//
// }
//
// /* find cosine similarity between two QFs. */
// uint64_t onDiskMQF_inner_product(onDiskMQF *qfa, onDiskMQF *qfb)
// {
// 	uint64_t acc = 0;
// 	onDiskMQFIterator qfi;
// 	onDiskMQF *onDiskMQF_mem, *onDiskMQF_disk;
//
// 	// create the iterator on the larger QF.
// 	if (qfa->metadata->size > qfb->metadata->size) {
// 		onDiskMQF_mem = qfb;
// 		onDiskMQF_disk = qfa;
// 	} else {
// 		onDiskMQF_mem = qfa;
// 		onDiskMQF_disk = qfb;
// 	}
//
// 	onDiskMQF_iterator(onDiskMQF_disk, &qfi, 0);
// 	do {
// 		uint64_t key = 0, value = 0, count = 0;
// 		uint64_t count_mem;
// 		onDiskMQFIterator_get(&qfi, &key, &value, &count);
// 		if ((count_mem = onDiskMQF_count_key(onDiskMQF_mem, key)) > 0) {
// 			acc += count*count_mem;
// 		}
// 	} while (!onDiskMQFIterator_next(&qfi));
//
// 	return acc;
// }
//
// /* find cosine similarity between two QFs. */
// // void onDiskMQF_intersect(QF *qfa, QF *qfb, QF *qfr)
// // {
// // 	QFi qfi;
// // 	QF *onDiskMQF_mem, *onDiskMQF_disk;
// //
// // 	// create the iterator on the larger QF.
// // 	if (qfa->metadata->size > qfb->metadata->size) {
// // 		onDiskMQF_mem = qfb;
// // 		onDiskMQF_disk = qfa;
// // 	} else {
// // 		onDiskMQF_mem = qfa;
// // 		onDiskMQF_disk = qfb;
// // 	}
// //
// // 	onDiskMQF_iterator(onDiskMQF_disk, &qfi, 0);
// // 	do {
// // 		uint64_t key = 0, value = 0, count = 0;
// // 		onDiskMQFIterator_get(&qfi, &key, &value, &count);
// // 		if (onDiskMQF_count_key(onDiskMQF_mem, key) > 0)
// // 			onDiskMQF_insert(qfr, key, count, false, false);
// // 	} while (!onDiskMQFIterator_next(&qfi));
// // }
//
// /* magnitude of a QF. */
// uint64_t onDiskMQF_magnitude(onDiskMQF *qf)
// {
// 	return sqrt(onDiskMQF_inner_product( qf));
// }
//
template<uint64_t bitsPerSlot>
int _onDiskMQF<bitsPerSlot>::space()
{
	_onDiskMQF<bitsPerSlot>* qf=this;
		return (int)(((double)qf->metadata->noccupied_slots/
								 (double)qf->metadata->xnslots
							 )* 100.0);
}

template<uint64_t bitsPerSlot>
bool onDiskMQF_Namespace::_onDiskMQF<bitsPerSlot>::general_lock(bool spin){
	_onDiskMQF<bitsPerSlot>* qf=this;
	if (!onDiskMQF_spin_lock(&qf->mem->general_lock, spin))
		return false;
	return true;
}
template<uint64_t bitsPerSlot>
void _onDiskMQF<bitsPerSlot>::general_unlock(){
	_onDiskMQF<bitsPerSlot>* qf=this;
	onDiskMQF_spin_unlock(&qf->mem->general_lock);
}
// void onDiskMQF_migrate(onDiskMQF* source, onDiskMQF* dest){
// 	onDiskMQFIterator source_i;
// 	if (onDiskMQF_iterator(source, &source_i, 0)) {
// 		do {
// 			uint64_t key = 0, value = 0, count = 0;
// 			onDiskMQFIterator_get(&source_i, &key, &value, &count);
// 			onDiskMQF_insert(dest, key, count, true, true);
// 			onDiskMQF_add_tag(dest,key,value);
// 		} while (!onDiskMQFIterator_next(&source_i));
// 	}
// }

template<uint64_t bitsPerSlot>
void onDiskMQF_Namespace::_onDiskMQF<bitsPerSlot>::migrateFromQF(QF* source){
	QFi source_i;
	if (qf_iterator(source, &source_i, 0)) {
		do {
			uint64_t key = 0, value = 0, count = 0;
			qfi_get(&source_i, &key, &value, &count);
			insert(key, count, true, true);
			add_tag(key,value);
		} while (!qfi_next(&source_i));
	}
}

};



#ifdef TEST
	#include "tests/lowLevelTests.hpp"
#endif
