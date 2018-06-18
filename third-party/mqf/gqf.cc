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
#include <iostream>
#include <map>
#include "utils.h"


/******************************************************************
 * Code for managing the metadata bits and slots w/o interpreting *
 * the content of the slots.
 ******************************************************************/

/* Must be >= 6.  6 seems fastest. */
#define BLOCK_OFFSET_BITS (6)

#define SLOTS_PER_BLOCK (1ULL << BLOCK_OFFSET_BITS)
#define METADATA_WORDS_PER_BLOCK ((SLOTS_PER_BLOCK + 63) / 64)

#define NUM_SLOTS_TO_LOCK (1ULL<<16)
#define CLUSTER_SIZE (1ULL<<14)

#define METADATA_WORD(qf,field,slot_index) (get_block((qf), (slot_index) / \
					SLOTS_PER_BLOCK)->field[((slot_index)  % SLOTS_PER_BLOCK) / 64])

#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) \
												- 1ULL)
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

#ifdef LOG_WAIT_TIME
static inline bool qf_spin_lock(QF *cf, volatile int *lock, uint64_t idx, bool
																flag_spin)
{
	struct timespec start, end;
	bool ret;

	clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
	if (!flag_spin) {
		ret = !__sync_lock_test_and_set(lock, 1);
		clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);
		cf->mem->wait_times[idx].locks_acquired_single_attempt++;
		cf->mem->wait_times[idx].total_time_single += BILLION * (end.tv_sec -
																												start.tv_sec) +
			end.tv_nsec - start.tv_nsec;
	} else {
		if (!__sync_lock_test_and_set(lock, 1)) {
			clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);
			cf->mem->wait_times[idx].locks_acquired_single_attempt++;
			cf->mem->wait_times[idx].total_time_single += BILLION * (end.tv_sec -
																													start.tv_sec) +
			end.tv_nsec - start.tv_nsec;
		} else {
			while (__sync_lock_test_and_set(lock, 1))
				while (*lock);
			clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);
			cf->mem->wait_times[idx].total_time_spinning += BILLION * (end.tv_sec -
																														start.tv_sec) +
				end.tv_nsec - start.tv_nsec;
		}
		ret = true;
	}
	cf->mem->wait_times[idx].locks_taken++;

	return ret;

	/*start = rdtsc();*/
	/*if (!__sync_lock_test_and_set(lock, 1)) {*/
		/*clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);*/
		/*cf->mem->wait_times[idx].locks_acquired_single_attempt++;*/
		/*cf->mem->wait_times[idx].total_time_single += BILLION * (end.tv_sec -
		 * start.tv_sec) + end.tv_nsec - start.tv_nsec;*/
	/*} else {*/
		/*while (__sync_lock_test_and_set(lock, 1))*/
			/*while (*lock);*/
		/*clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);*/
		/*cf->mem->wait_times[idx].total_time_spinning += BILLION * (end.tv_sec -
		 * start.tv_sec) + end.tv_nsec - start.tv_nsec;*/
	/*}*/

	/*end = rdtsc();*/
	/*cf->mem->wait_times[idx].locks_taken++;*/
	/*return;*/
}
#else
/**
 * Try to acquire a lock once and return even if the lock is busy.
 * If spin flag is set, then spin until the lock is available.
 */
static inline bool qf_spin_lock(volatile int *lock, bool flag_spin)
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
#endif

static inline void qf_spin_unlock(volatile int *lock)
{
	__sync_lock_release(lock);
	return;
}

static bool qf_lock(const QF *cf, uint64_t hash_bucket_index, bool spin, bool flag)
{
	uint64_t hash_bucket_lock_offset  = hash_bucket_index % NUM_SLOTS_TO_LOCK;
	if (flag) {
#ifdef LOG_WAIT_TIME
		if (!qf_spin_lock(cf, &cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK],
											hash_bucket_index/NUM_SLOTS_TO_LOCK, spin))
			return false;
		if (NUM_SLOTS_TO_LOCK - hash_bucket_lock_offset <= CLUSTER_SIZE) {
			if (!qf_spin_lock(cf, &cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1],
												hash_bucket_index/NUM_SLOTS_TO_LOCK+1, spin)) {
				qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
				return false;
			}
		}
#else
		if (!qf_spin_lock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK], spin))
			return false;
		if (NUM_SLOTS_TO_LOCK - hash_bucket_lock_offset <= CLUSTER_SIZE) {
			if (!qf_spin_lock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1],
												spin)) {
				qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
				return false;
			}
		}
#endif
	} else {
		/* take the lock for two lock-blocks; the lock-block in which the
		 * hash_bucket_index falls and the next lock-block */

#ifdef LOG_WAIT_TIME
		if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
				CLUSTER_SIZE) {
			if (!qf_spin_lock(cf, &cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1], spin))
				return false;
		}
		if (!qf_spin_lock(cf, &cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK], spin)) {
			if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
					CLUSTER_SIZE)
				qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
			return false;
		}
		if (!qf_spin_lock(cf, &cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1],
											spin)) {
			qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
			if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
					CLUSTER_SIZE)
				qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
			return false;
		}
#else
		if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
				CLUSTER_SIZE) {
			if (!qf_spin_lock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1], spin))
				return false;
		}
		if (!qf_spin_lock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK], spin)) {
			if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
					CLUSTER_SIZE)
				qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
			return false;
		}
		if (!qf_spin_lock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1],
											spin)) {
			qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
			if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
					CLUSTER_SIZE)
				qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
			return false;
		}
#endif
	}
	return true;
}

static void qf_unlock(const QF *cf, uint64_t hash_bucket_index, bool flag)
{
	uint64_t hash_bucket_lock_offset  = hash_bucket_index % NUM_SLOTS_TO_LOCK;
	if (flag) {
		if (NUM_SLOTS_TO_LOCK - hash_bucket_lock_offset <= CLUSTER_SIZE) {
			qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1]);
		}
		qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
	} else {
		qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1]);
		qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
		if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
				CLUSTER_SIZE)
			qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
	}
}

static void modify_metadata(QF *cf, uint64_t *metadata, int cnt)
{
#ifdef LOG_WAIT_TIME
	qf_spin_lock(cf, &cf->mem->metadata_lock,cf->num_locks, true);
#else
	qf_spin_lock(&cf->mem->metadata_lock, true);
#endif
	*metadata = *metadata + cnt;
	qf_spin_unlock(&cf->mem->metadata_lock);
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

#if BITS_PER_SLOT > 0
static inline qfblock * get_block(const QF *qf, uint64_t block_index)
{
	return &qf->blocks[block_index];
}
#else
static inline qfblock * get_block(const QF *qf, uint64_t block_index)
{
	 // printf("block = %p\n",(void*)(((char *)qf->blocks) + block_index * (sizeof(qfblock) +
	 // 					 qf->metadata->bits_per_slot * 8 +
	 // 					8*qf->metadata->fixed_counter_size +
	 // 					8*qf->metadata->tag_bits
	 // 				 )) );
	 //printf("blocks start=%p\n",qf->blocks );
	return (qfblock *)(((char *)qf->blocks) + block_index * (sizeof(qfblock) +
						 qf->metadata->bits_per_slot * 8
					 ));
}
#endif

static inline int is_runend(const QF *qf, uint64_t index)
{
	return (METADATA_WORD(qf, runends, index) >> ((index % SLOTS_PER_BLOCK) %
																								64)) & 1ULL;
}

static inline int is_occupied(const QF *qf, uint64_t index)
{
	return (METADATA_WORD(qf, occupieds, index) >> ((index % SLOTS_PER_BLOCK) %
																									64)) & 1ULL;
}

#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64

static inline uint64_t get_slot(const QF *qf, uint64_t index)
{
	assert(index < qf->metadata->xnslots);
	return get_block(qf, index / SLOTS_PER_BLOCK)->slots[index % SLOTS_PER_BLOCK];
}

static inline void set_slot(const QF *qf, uint64_t index, uint64_t value)
{
	assert(index < qf->metadata->xnslots);
	get_block(qf, index / SLOTS_PER_BLOCK)->slots[index % SLOTS_PER_BLOCK] =
		value & BITMASK(qf->metadata->bits_per_slot);
}

#elif BITS_PER_SLOT > 0

/* Little-endian code ....  Big-endian is TODO */

static inline uint64_t get_slot(const QF *qf, uint64_t index)
{
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems
	 * to generate buggy code.  :/  */
	assert(index < qf->metadata->xnslots);
	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 SLOTS_PER_BLOCK)->slots[(index %
																																SLOTS_PER_BLOCK)
																			 * BITS_PER_SLOT / 8];
	return (uint64_t)(((*p) >> (((index % SLOTS_PER_BLOCK) * BITS_PER_SLOT) %
															8)) & BITMASK(BITS_PER_SLOT));
}

static inline void set_slot(const QF *qf, uint64_t index, uint64_t value)
{
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems
	 * to generate buggy code.  :/  */
	assert(index < qf->metadata->xnslots);

	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 SLOTS_PER_BLOCK)->slots[(index %
																																SLOTS_PER_BLOCK)
																			 * BITS_PER_SLOT / 8];
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

static inline uint64_t _get_slot(const QF *qf, uint64_t index)
{
	assert(index < qf->metadata->xnslots);
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems
	 * to generate buggy code.  :/  */

	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 SLOTS_PER_BLOCK)->slots[(index %
																																SLOTS_PER_BLOCK)
																			 * qf->metadata->bits_per_slot / 8];
	return (uint64_t)(((*p) >> (((index % SLOTS_PER_BLOCK) *
															 qf->metadata->bits_per_slot) % 8)) &
										BITMASK(qf->metadata->bits_per_slot));
}

static inline void _set_slot(const QF *qf, uint64_t index, uint64_t value)
{
	assert(index < qf->metadata->xnslots);
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems
	 * to generate buggy code.  :/  */
	 //printf("ss %d\n",(index %SLOTS_PER_BLOCK)* qf->metadata->bits_per_slot / 8 );
	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 SLOTS_PER_BLOCK)->slots[(index %
																																SLOTS_PER_BLOCK)
																			 * qf->metadata->bits_per_slot / 8];

	uint64_t t = *p;
	uint64_t mask = BITMASK(qf->metadata->bits_per_slot);
	uint64_t v = value;
	int shift = ((index % SLOTS_PER_BLOCK) * qf->metadata->bits_per_slot) % 8;
	mask <<= shift;
	v <<= shift;
	t &= ~mask;
	t |= v;
	*p = t;
}

static inline void set_slot(const QF *qf, uint64_t index, uint64_t value){
	uint64_t original_value=_get_slot(qf,index);
	value<<=qf->metadata->fixed_counter_size;
	uint64_t mask=BITMASK(qf->metadata->fixed_counter_size);
	original_value&=mask;
	value|=original_value;
	_set_slot(qf,index,value);
}

static inline uint64_t get_slot(const QF *qf, uint64_t index)
{
	uint64_t mask=BITMASK(qf->metadata->key_remainder_bits);
	uint64_t t=_get_slot(qf,index);
	t >>= qf->metadata->fixed_counter_size;

	return t&mask;
}


#endif

static inline uint64_t get_fixed_counter(const QF *qf, uint64_t index)
{
	uint64_t t=_get_slot(qf,index);
	uint64_t mask=BITMASK(qf->metadata->fixed_counter_size);
	return t&mask;
	// uint64_t res=0;
	// uint64_t base=1;
	// uint64_t* p=(uint64_t*)((uint8_t*)get_block(qf, index /SLOTS_PER_BLOCK)->slots+
	// 															(SLOTS_PER_BLOCK * qf->metadata->bits_per_slot )/ 8);
  //
	// for(int i=qf->metadata->fixed_counter_size-1;i>=0;i--){
	// 	res+= base*(( p[i] >> ((index % SLOTS_PER_BLOCK) %64)) & 1ULL);
	// 	base*=2;
	// }
  //
	// return res;
}
 inline  void set_fixed_counter(const QF *qf, uint64_t index,uint64_t value)
{
	uint64_t mask=BITMASK(qf->metadata->fixed_counter_size);
	uint64_t t=value & mask;
	uint64_t original_value=_get_slot(qf,index);
	original_value&= ~mask;
	t|=original_value;
	_set_slot(qf,index,t);
	// //printf("bits per slot =%lu, fixed counter =%lu, tag_bits=%lu\n",
	// //qf->metadata->bits_per_slot,qf->metadata->fixed_counter_size,qf->metadata->tag_bits );
	// //printf("slots start=%p\n",(void*) get_block(qf, index /SLOTS_PER_BLOCK)->slots);
	// uint64_t* p=((uint64_t*)&get_block(qf, index /SLOTS_PER_BLOCK)->slots) +
	//  (qf->metadata->bits_per_slot) ;
	// printf("block start=%p\n", (void*)get_block(qf, index /SLOTS_PER_BLOCK));
	// printf("slots start=%p\n", (void*)&get_block(qf, index /SLOTS_PER_BLOCK)->slots);
	// //printf("fixed counter=%p index=%lu , value=%lu add=%d\n",(void*)p,index,value,(8 * qf->metadata->bits_per_slot) );
	// //printf("diff=%d\n",(char*)p-(char*) get_block(qf, index /SLOTS_PER_BLOCK) );
	// uint64_t bitmask=1ULL << ((index % SLOTS_PER_BLOCK) %64);
	//
	// int i= qf->metadata->fixed_counter_size-1 ;
	// //printf("fp=%p\n",&p[i]);
	// //printf("fp=%p\n",&p[0]);
	// for(;i>=0;i--){
	// 	//printf("fp=%p\n",&p[i]);
	// 	if(value%2){
	// 		p[i]|= bitmask;
	// 	}
	// 	else{
	// 		p[i]&= ~(bitmask);
	// 	}
	// 	value=value>>1;
	//
	// }
	// //printf("finish\n" );


}

static inline uint64_t get_tag(const QF *qf, uint64_t index)
{
	uint64_t mask=BITMASK(qf->metadata->tag_bits);
	uint64_t t=_get_slot(qf,index);
	t >>= (qf->metadata->fixed_counter_size+qf->metadata->key_remainder_bits);
	return t&mask;
}
static inline void set_tag(const QF *qf, uint64_t index,uint64_t value)
{
	uint64_t original_value=_get_slot(qf,index);
	value<<=(qf->metadata->fixed_counter_size+qf->metadata->key_remainder_bits);
	uint64_t mask=BITMASK(qf->metadata->fixed_counter_size+qf->metadata->key_remainder_bits);
	original_value&=mask;
	value|=original_value;
	_set_slot(qf,index,value);

}



static inline uint64_t run_end(const QF *qf, uint64_t hash_bucket_index);

static inline uint64_t block_offset(const QF *qf, uint64_t blockidx)
{
	/* If we have extended counters and a 16-bit (or larger) offset
		 field, then we can safely ignore the possibility of overflowing
		 that field. */
	if (sizeof(qf->blocks[0].offset) > 1 ||
			get_block(qf, blockidx)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
		return get_block(qf, blockidx)->offset;

	return run_end(qf, SLOTS_PER_BLOCK * blockidx - 1) - SLOTS_PER_BLOCK *
		blockidx + 1;
}

static inline uint64_t run_end(const QF *qf, uint64_t hash_bucket_index)
{
	uint64_t bucket_block_index       = hash_bucket_index / SLOTS_PER_BLOCK;
	uint64_t bucket_intrablock_offset = hash_bucket_index % SLOTS_PER_BLOCK;
	uint64_t bucket_blocks_offset = block_offset(qf, bucket_block_index);

	uint64_t bucket_intrablock_rank   = bitrank(get_block(qf,
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
	uint64_t runend_block_offset = bitselectv(get_block(qf,
																						runend_block_index)->runends[0],
																						runend_ignore_bits, runend_rank);
	if (runend_block_offset == SLOTS_PER_BLOCK) {
		if (bucket_blocks_offset == 0 && bucket_intrablock_rank == 0) {
			/* The block begins in empty space, and this bucket is in that region of
			 * empty space */
			return hash_bucket_index;
		} else {
			do {
				runend_rank        -= popcntv(get_block(qf,
																								runend_block_index)->runends[0],
																			runend_ignore_bits);
				runend_block_index++;
				runend_ignore_bits  = 0;
				runend_block_offset = bitselectv(get_block(qf,
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

static inline void shift_remainders(QF *qf, uint64_t start_index, uint64_t
																		empty_index)
{
	uint64_t start_block  = start_index / SLOTS_PER_BLOCK;
	uint64_t start_offset = start_index % SLOTS_PER_BLOCK;
	uint64_t empty_block  = empty_index / SLOTS_PER_BLOCK;
	uint64_t empty_offset = empty_index % SLOTS_PER_BLOCK;

	assert (start_index <= empty_index && empty_index < qf->metadata->xnslots);

	while (start_block < empty_block) {
		memmove(&get_block(qf, empty_block)->slots[1],
						&get_block(qf, empty_block)->slots[0],
						empty_offset * sizeof(qf->blocks[0].slots[0]));
		get_block(qf, empty_block)->slots[0] = get_block(qf,
																			empty_block-1)->slots[SLOTS_PER_BLOCK-1];
		empty_block--;
		empty_offset = SLOTS_PER_BLOCK-1;
	}

	memmove(&get_block(qf, empty_block)->slots[start_offset+1],
					&get_block(qf, empty_block)->slots[start_offset],
					(empty_offset - start_offset) * sizeof(qf->blocks[0].slots[0]));
}

#else

#define REMAINDER_WORD(qf, i) ((uint64_t *)&(get_block(qf, (i)/qf->metadata->bits_per_slot)->slots[8 * ((i) % qf->metadata->bits_per_slot)]))

static inline void shift_remainders(QF *qf, const uint64_t start_index, const
																		uint64_t empty_index)
{
	uint64_t last_word = (empty_index + 1) * qf->metadata->bits_per_slot / 64;
	const uint64_t first_word = start_index * qf->metadata->bits_per_slot / 64;
	int bend = ((empty_index + 1) * qf->metadata->bits_per_slot) % 64;
	const int bstart = (start_index * qf->metadata->bits_per_slot) % 64;

	while (last_word != first_word) {
		*REMAINDER_WORD(qf, last_word) = shift_into_b(*REMAINDER_WORD(qf, last_word-1),
																									*REMAINDER_WORD(qf, last_word),
																									0, bend, qf->metadata->bits_per_slot);
		last_word--;
		bend = 64;
	}
	*REMAINDER_WORD(qf, last_word) = shift_into_b(0, *REMAINDER_WORD(qf,
																																	 last_word),
																								bstart, bend,
																								qf->metadata->bits_per_slot);
}

#endif

static inline void qf_dump_block(const QF *qf, uint64_t i)
{
	uint64_t j;
 	printf("#Block %lu \n",i );
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
	//for (j = 0; j < SLOTS_PER_BLOCK * qf->metadata->bits_per_slot / 8; j++)
	//	printf("%02x ", get_block(qf, i)->slots[j]);

#endif

for(int j=0;j<64;j++)
{
	printf("%lu ", get_slot(qf,j+64*i));
}

	printf("\n fixed counter \n");

	for(int j=0;j<64;j++)
	{
		printf("%lu ", get_fixed_counter(qf,j+64*i));
	}

	printf("\n tags \n");

	for(int j=0;j<64;j++)
	{
		printf("%lu ", get_tag(qf,j+64*i));
	}


	printf("\n");
	printf("\n");
}

void qf_dump(const QF *qf)
{
	uint64_t i;

	printf("%lu %lu %lu\n",
				 qf->metadata->nblocks,
				 qf->metadata->ndistinct_elts,
				 qf->metadata->nelts);

	for (i = 0; i < qf->metadata->nblocks; i++) {
		qf_dump_block(qf, i);
	}
	printf("End\n");




}

static inline void find_next_n_empty_slots(QF *qf, uint64_t from, uint64_t n,
																					 uint64_t *indices)
{
	while (n) {
		indices[--n] = find_first_empty_slot(qf, from);
		from = indices[n] + 1;
	}
}

static inline void shift_slots(QF *qf, int64_t first, uint64_t last, uint64_t
															 distance)
{
	int64_t i;
	if (distance == 1)
		shift_remainders(qf, first, last+1);
	else
		for (i = last; i >= first; i--)
			_set_slot(qf, i + distance, _get_slot(qf, i));
}

static inline void shift_runends(QF *qf, int64_t first, uint64_t last,
																 uint64_t distance)
{
	assert(last < qf->metadata->xnslots && distance < 64);
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
	METADATA_WORD(qf, runends, 64*last_word) = shift_into_b(0, METADATA_WORD(qf,
																																					 runends,
																																					 64*last_word),
																													bstart, bend, distance);

}

// static inline void shift_fixed_counters(QF *qf, int64_t first, uint64_t last,
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
// 			curr=(uint64_t*)((uint8_t*)get_block(qf, last_word)->slots+
// 			(SLOTS_PER_BLOCK * qf->metadata->bits_per_slot) / 8);
// 			prev=(uint64_t*)((uint8_t*)get_block(qf, last_word-1)->slots+
// 			(SLOTS_PER_BLOCK * qf->metadata->bits_per_slot) / 8);
// 			curr[i] = shift_into_b(prev[i],curr[i],0, bend, distance);
// 			bend = 64;
// 			last_word--;
// 			while (last_word != first_word) {
// 				curr=(uint64_t*)((uint8_t*)get_block(qf, last_word)->slots+(SLOTS_PER_BLOCK * qf->metadata->bits_per_slot / 8));
// 				prev=(uint64_t*)((uint8_t*)get_block(qf, last_word-1)->slots+(SLOTS_PER_BLOCK * qf->metadata->bits_per_slot / 8));
// 				curr[i] = shift_into_b(prev[i],curr[i],0, bend, distance);
// 				last_word--;
// 			}
// 		}
// 		curr=(uint64_t*)((uint8_t*)get_block(qf, last_word)->slots+(SLOTS_PER_BLOCK * qf->metadata->bits_per_slot / 8));
// 		curr[i] = shift_into_b(0,curr[i], bstart, bend, distance);
// 	}
//
// }

// static inline void shift_tags(QF *qf, int64_t first, uint64_t last,
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
// 			curr=(uint64_t*)((uint8_t*)get_block(qf, last_word)->slots+
// 											(8 *(qf->metadata->bits_per_slot + qf->metadata->fixed_counter_size )));
// 			prev=(uint64_t*)((uint8_t*)get_block(qf, last_word-1)->slots+
// 											(8 *(qf->metadata->bits_per_slot + qf->metadata->fixed_counter_size )));
// 			curr[i] = shift_into_b(prev[i],curr[i],0, bend, distance);
// 			bend = 64;
// 			last_word--;
// 			while (last_word != first_word) {
// 				curr=(uint64_t*)((uint8_t*)get_block(qf, last_word)->slots+
// 												(8 *(qf->metadata->bits_per_slot + qf->metadata->fixed_counter_size )));
// 				prev=(uint64_t*)((uint8_t*)get_block(qf, last_word-1)->slots+
// 												(8 *(qf->metadata->bits_per_slot + qf->metadata->fixed_counter_size )));
// 				curr[i] = shift_into_b(prev[i],curr[i],0, bend, distance);
// 				last_word--;
// 			}
// 		}
// 		curr=(uint64_t*)((uint8_t*)get_block(qf, last_word)->slots+
// 										(8 *(qf->metadata->bits_per_slot + qf->metadata->fixed_counter_size )));
// 		curr[i] = shift_into_b(0,curr[i], bstart, bend, distance);
// 	}
//
// }

static inline void insert_replace_slots_and_shift_remainders_and_runends_and_offsets(QF		*qf,
																																										 int		 operation,
																																										 uint64_t		 bucket_index,
																																										 uint64_t		 overwrite_index,
																																										 const uint64_t	*remainders,
																																										 const uint64_t	*fixed_size_counters,
																																										 uint64_t		 total_remainders,
																																										 uint64_t		 noverwrites)
{
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

		find_next_n_empty_slots(qf, insert_index, ninserts, empties);
		for (i = 0; i < ninserts - 1; i++){
			shift_slots(qf, empties[i+1] + 1, empties[i] - 1, i + 1);
		}
		shift_slots(qf, insert_index, empties[ninserts - 1] - 1, ninserts);



		for (i = 0; i < ninserts - 1; i++)
			shift_runends(qf, empties[i+1] + 1, empties[i] - 1, i + 1);
		shift_runends(qf, insert_index, empties[ninserts - 1] - 1, ninserts);



		// for (i = 0; i < ninserts - 1; i++)
		// 	shift_fixed_counters(qf, empties[i+1] + 1, empties[i] - 1, i + 1);
		// shift_fixed_counters(qf, insert_index, empties[ninserts - 1] - 1, ninserts);
    //
    //
		// for (i = 0; i < ninserts - 1; i++)
		// 	shift_tags(qf, empties[i+1] + 1, empties[i] - 1, i + 1);
		// shift_tags(qf, insert_index, empties[ninserts - 1] - 1, ninserts);
    //



		for (i = noverwrites; i < total_remainders - 1; i++)
			METADATA_WORD(qf, runends, overwrite_index + i) &= ~(1ULL <<
																													 (((overwrite_index
																															+ i) %
																														 SLOTS_PER_BLOCK)
																														% 64));
		// for (i = noverwrites; i < total_remainders - 1; i++)
		// 	set_fixed_counter(qf,overwrite_index+i,0);


		switch (operation) {
			case 0: /* insert into empty bucket */
				assert (noverwrites == 0);
				METADATA_WORD(qf, runends, overwrite_index + total_remainders - 1) |=
					1ULL << (((overwrite_index + total_remainders - 1) %
										SLOTS_PER_BLOCK) % 64);
				break;
			case 1: /* append to bucket */
				METADATA_WORD(qf, runends, overwrite_index + noverwrites - 1)      &=
					~(1ULL << (((overwrite_index + noverwrites - 1) % SLOTS_PER_BLOCK) %
										 64));
				METADATA_WORD(qf, runends, overwrite_index + total_remainders - 1) |=
					1ULL << (((overwrite_index + total_remainders - 1) %
										SLOTS_PER_BLOCK) % 64);
				break;
			case 2: /* insert into bucket */
				METADATA_WORD(qf, runends, overwrite_index + total_remainders - 1) &=
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

			if (get_block(qf, i)->offset + ninserts - npreceding_empties < BITMASK(8*sizeof(qf->blocks[0].offset)))
				get_block(qf, i)->offset += ninserts - npreceding_empties;
			else
				get_block(qf, i)->offset = (uint8_t) BITMASK(8*sizeof(qf->blocks[0].offset));
		}
	}
	for (i = 0; i < total_remainders; i++){
		set_slot(qf, overwrite_index + i, remainders[i]);
		set_fixed_counter(qf,overwrite_index +i,fixed_size_counters[i]);
//		printf("fixed counter = %lu\n",fixed_size_counters[i] );
	}

	modify_metadata(qf, &qf->metadata->noccupied_slots, ninserts);
}

static inline void remove_replace_slots_and_shift_remainders_and_runends_and_offsets(QF		        *qf,
																																										 int		 operation,
																																										 uint64_t		 bucket_index,
																																										 uint64_t		 overwrite_index,
																																										 const uint64_t	*remainders,
																																										 const uint64_t	*fcounters,
																																										 uint64_t		 total_remainders,
																																										 uint64_t		 old_length)
{
	uint64_t i;

	// Update the slots
	for (i = 0; i < total_remainders; i++){
		set_slot(qf, overwrite_index + i, remainders[i]);
		set_fixed_counter(qf, overwrite_index + i, fcounters[i]);
		if(qf->metadata->tag_bits>0)
			set_tag(qf, overwrite_index + i, 0);
	}


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
			set_fixed_counter(qf, current_slot, get_fixed_counter(qf, current_slot + current_distance));
			if(qf->metadata->tag_bits>0)
				set_tag(qf, current_slot, get_tag(qf, current_slot + current_distance));

			if (is_runend(qf, current_slot) !=
					is_runend(qf, current_slot + current_distance))
				METADATA_WORD(qf, runends, current_slot) ^= 1ULL << (current_slot % 64);
			current_slot++;

		} else if (current_bucket <= current_slot + current_distance) {
			uint64_t i;
			for (i = current_slot; i < current_slot + current_distance; i++) {
				set_slot(qf, i, 0);
				set_fixed_counter(qf,i,0);
				if(qf->metadata->tag_bits>0)
					set_tag(qf,i,0);
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

	int num_slots_freed = old_length - total_remainders;
	modify_metadata(qf, &qf->metadata->noccupied_slots, -num_slots_freed);
	/*qf->metadata->noccupied_slots -= (old_length - total_remainders);*/
	if (!total_remainders) {
		modify_metadata(qf, &qf->metadata->ndistinct_elts, -1);
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
static inline uint64_t *encode_counter(QF *qf, uint64_t remainder, uint64_t
																			 counter, uint64_t *slots, uint64_t *fixed_size_counters)
{
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
static inline uint64_t decode_counter(const QF *qf, uint64_t index, uint64_t
																			*remainder, uint64_t *count)
{

	*remainder  = get_slot(qf, index);
	uint64_t fcount=get_fixed_counter(qf,index);
	uint64_t tmp_count= fcount+1;
	*count=0;
	const uint64_t fixed_count_max=(1ULL << qf->metadata->fixed_counter_size)-1;
	//printf("tmp count = %lu\n",tmp_count);
	if(fcount == fixed_count_max){
		uint64_t no_digits=0;
		do{
				index++;
				no_digits++;
				*count <<= qf->metadata->key_remainder_bits;

				fcount= get_fixed_counter(qf,index);
		//		printf("quer slot =%lu  fixed count= %lu\n", get_slot(qf, index),fcount);
				*count += get_slot(qf, index);

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
static inline uint64_t next_slot(QF *qf, uint64_t current)
{
	uint64_t rem = get_slot(qf, current);
	current++;

	while (get_slot(qf, current) == rem && current <= qf->metadata->nslots) {
		current++;
	}
	return current;
}

// static inline bool insert1(QF *qf, __uint128_t hash, bool lock, bool spin)
// {
// 	uint64_t hash_remainder           = hash & BITMASK(qf->metadata->bits_per_slot);
// 	uint64_t hash_bucket_index        = hash >> qf->metadata->bits_per_slot;
// 	uint64_t hash_bucket_block_offset = hash_bucket_index % SLOTS_PER_BLOCK;
//
// 	if (lock) {
// 		if (!qf_lock(qf, hash_bucket_index, spin, true))
// 			return false;
// 	}
// 	if (is_empty(qf, hash_bucket_index) /* might_be_empty(qf, hash_bucket_index) && runend_index == hash_bucket_index */) {
// 		METADATA_WORD(qf, runends, hash_bucket_index) |= 1ULL <<
// 			(hash_bucket_block_offset % 64);
// 		set_slot(qf, hash_bucket_index, hash_remainder);
// 		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL <<
// 			(hash_bucket_block_offset % 64);
//
// 		/*modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);*/
// 		modify_metadata(qf, &qf->metadata->noccupied_slots, 1);
// 		/*modify_metadata(qf, &qf->metadata->nelts, 1);*/
// 	} else {
// 		uint64_t runend_index              = run_end(qf, hash_bucket_index);
// 		int operation = 0; /* Insert into empty bucket */
// 		uint64_t insert_index = runend_index + 1;
// 		uint64_t new_value = hash_remainder;
//
// 		/* printf("RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */
//
// 		uint64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
// 																																	 hash_bucket_index
// 																																	 - 1) + 1;
//
// 		if (is_occupied(qf, hash_bucket_index)) {
//
// 			/* Find the counter for this remainder if it exists. */
// 			uint64_t current_remainder = get_slot(qf, runstart_index);
// 			uint64_t zero_terminator = runstart_index;
//
// 			/* The counter for 0 is special. */
// 			if (current_remainder == 0) {
// 				uint64_t t = runstart_index + 1;
// 				while (t < runend_index && get_slot(qf, t) != 0)
// 					t++;
// 				if (t < runend_index && get_slot(qf, t+1) == 0)
// 					zero_terminator = t+1; /* Three or more 0s */
// 				else if (runstart_index < runend_index && get_slot(qf, runstart_index
// 																													 + 1) == 0)
// 					zero_terminator = runstart_index + 1; /* Exactly two 0s */
// 				/* Otherwise, exactly one 0 (i.e. zero_terminator == runstart_index) */
//
// 				/* May read past end of run, but that's OK because loop below
// 					 can handle that */
// 				if (hash_remainder != 0) {
// 					runstart_index = zero_terminator + 1;
// 					current_remainder = get_slot(qf, runstart_index);
// 				}
// 			}
//
// 			/* Skip over counters for other remainders. */
// 			while (current_remainder < hash_remainder && runstart_index <=
// 						 runend_index) {
// 				/* If this remainder has an extended counter, skip over it. */
// 				if (runstart_index < runend_index &&
// 						get_slot(qf, runstart_index + 1) < current_remainder) {
// 					runstart_index = runstart_index + 2;
// 					while (runstart_index < runend_index &&
// 								 get_slot(qf, runstart_index) != current_remainder)
// 						runstart_index++;
// 					runstart_index++;
//
// 					/* This remainder has a simple counter. */
// 				} else {
// 					runstart_index++;
// 				}
//
// 				/* This may read past the end of the run, but the while loop
// 					 condition will prevent us from using the invalid result in
// 					 that case. */
// 				current_remainder = get_slot(qf, runstart_index);
// 			}
//
// 			/* If this is the first time we've inserted the new remainder,
// 				 and it is larger than any remainder in the run. */
// 			if (runstart_index > runend_index) {
// 				operation = 1;
// 				insert_index = runstart_index;
// 				new_value = hash_remainder;
// 				/*modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);*/
//
// 				/* This is the first time we're inserting this remainder, but
// 					 there are larger remainders already in the run. */
// 			} else if (current_remainder != hash_remainder) {
// 				operation = 2; /* Inserting */
// 				insert_index = runstart_index;
// 				new_value = hash_remainder;
// 				/*modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);*/
//
// 				/* Cases below here: we're incrementing the (simple or
// 					 extended) counter for this remainder. */
//
// 				/* If there's exactly one instance of this remainder. */
// 			} else if (runstart_index == runend_index ||
// 								 (hash_remainder > 0 && get_slot(qf, runstart_index + 1) >
// 									hash_remainder) ||
// 								 (hash_remainder == 0 && zero_terminator == runstart_index)) {
// 				operation = 2; /* Insert */
// 				insert_index = runstart_index;
// 				new_value = hash_remainder;
//
// 				/* If there are exactly two instances of this remainder. */
// 			} else if ((hash_remainder > 0 && get_slot(qf, runstart_index + 1) ==
// 									hash_remainder) ||
// 								 (hash_remainder == 0 && zero_terminator == runstart_index + 1)) {
// 				operation = 2; /* Insert */
// 				insert_index = runstart_index + 1;
// 				new_value = 0;
//
// 				/* Special case for three 0s */
// 			} else if (hash_remainder == 0 && zero_terminator == runstart_index + 2) {
// 				operation = 2; /* Insert */
// 				insert_index = runstart_index + 1;
// 				new_value = 1;
//
// 				/* There is an extended counter for this remainder. */
// 			} else {
//
// 				/* Move to the LSD of the counter. */
// 				insert_index = runstart_index + 1;
// 				while (get_slot(qf, insert_index+1) != hash_remainder)
// 					insert_index++;
//
// 				/* Increment the counter. */
// 				uint64_t digit, carry;
// 				do {
// 					carry = 0;
// 					digit = get_slot(qf, insert_index);
// 					// Convert a leading 0 (which is special) to a normal encoded digit
// 					if (digit == 0) {
// 						digit++;
// 						if (digit == current_remainder)
// 							digit++;
// 					}
//
// 					// Increment the digit
// 					digit = (digit + 1) & BITMASK(qf->metadata->bits_per_slot);
//
// 					// Ensure digit meets our encoding requirements
// 					if (digit == 0) {
// 						digit++;
// 						carry = 1;
// 					}
// 					if (digit == current_remainder)
// 						digit = (digit + 1) & BITMASK(qf->metadata->bits_per_slot);
// 					if (digit == 0) {
// 						digit++;
// 						carry = 1;
// 					}
//
// 					set_slot(qf, insert_index, digit);
// 					insert_index--;
// 				} while(insert_index > runstart_index && carry);
//
// 				/* If the counter needs to be expanded. */
// 				if (insert_index == runstart_index && (carry > 0 || (current_remainder
// 																														 != 0 && digit >=
// 																														 current_remainder)))
// 				{
// 					operation = 2; /* insert */
// 					insert_index = runstart_index + 1;
// 					if (!carry)						/* To prepend a 0 before the counter if the MSD is greater than the rem */
// 						new_value = 0;
// 					else if (carry) {			/* Increment the new value because we don't use 0 to encode counters */
// 						new_value = 2;
// 						/* If the rem is greater than or equal to the new_value then fail*/
// 						assert(new_value < current_remainder);
// 					}
// 				} else {
// 					operation = -1;
// 				}
// 			}
// 		}
//
// 		if (operation >= 0) {
// 			uint64_t empty_slot_index = find_first_empty_slot(qf, runend_index+1);
//
// 			shift_remainders(qf, insert_index, empty_slot_index);
//
// 			set_slot(qf, insert_index, new_value);
//
// 			shift_runends(qf, insert_index, empty_slot_index-1, 1);
// 			shift_fixed_counters(qf, insert_index, empty_slot_index-1, 1);
// 			switch (operation) {
// 				case 0:
// 					METADATA_WORD(qf, runends, insert_index)   |= 1ULL << ((insert_index
// 																																	%
// 																																	SLOTS_PER_BLOCK)
// 																																 % 64);
// 					break;
// 				case 1:
// 					METADATA_WORD(qf, runends, insert_index-1) &= ~(1ULL <<
// 																													(((insert_index-1) %
// 																														SLOTS_PER_BLOCK) %
// 																													 64));
// 					METADATA_WORD(qf, runends, insert_index)   |= 1ULL << ((insert_index
// 																																	%
// 																																	SLOTS_PER_BLOCK)
// 																																 % 64);
// 					break;
// 				case 2:
// 					METADATA_WORD(qf, runends, insert_index)   &= ~(1ULL <<
// 																													((insert_index %
// 																														SLOTS_PER_BLOCK) %
// 																													 64));
// 					break;
// 				default:
// 					fprintf(stderr, "Invalid operation %d\n", operation);
// 					abort();
// 			}
// 			/*
// 			 * Increment the offset for each block between the hash bucket index
// 			 * and block of the empty slot
// 			 * */
// 			uint64_t i;
// 			for (i = hash_bucket_index / SLOTS_PER_BLOCK + 1; i <=
// 					 empty_slot_index/SLOTS_PER_BLOCK; i++) {
// 				if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
// 					get_block(qf, i)->offset++;
// 				assert(get_block(qf, i)->offset != 0);
// 			}
// 			modify_metadata(qf, &qf->metadata->noccupied_slots, 1);
// 		}
// 		/*modify_metadata(qf, &qf->metadata->nelts, 1);*/
// 		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL <<
// 			(hash_bucket_block_offset % 64);
// 	}
//
// 	if (lock) {
// 		qf_unlock(qf, hash_bucket_index, true);
// 	}
//
// 	return true;
// }

static inline bool insert(QF *qf, __uint128_t hash, uint64_t count, bool lock=false,
													bool spin=false)
{
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
		if (!qf_lock(qf, hash_bucket_index, spin, false))
			return false;
	}

	uint64_t runend_index             = run_end(qf, hash_bucket_index);

	/* Empty slot */
	if (might_be_empty(qf, hash_bucket_index) && runend_index ==
			hash_bucket_index) {
		METADATA_WORD(qf, runends, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);
		set_slot(qf, hash_bucket_index, hash_remainder);
		set_fixed_counter(qf, hash_bucket_index, 0);
		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);

		modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);
		modify_metadata(qf, &qf->metadata->noccupied_slots, 1);
		/*modify_metadata(qf, &qf->metadata->nelts, 1);*/
		/* This trick will, I hope, keep the fast case fast. */
		if (count > 1) {
			insert(qf, hash, count - 1, false, false);
		}
	} else { /* Non-empty slot */
		uint64_t new_values[67];
		uint64_t new_fcounters[67];
		uint64_t total_remainders;
		int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
																																	hash_bucket_index
																																	- 1) + 1;

		if (!is_occupied(qf, hash_bucket_index)) { /* Empty bucket, but its slot is occupied. */
			uint64_t *p = encode_counter(qf, hash_remainder, count, &new_values[67],&new_fcounters[67]);
			total_remainders=&new_values[67] - p;
			insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf,
																																				0,
																																				hash_bucket_index,
																																				runstart_index,
																																				p,
																																				&new_fcounters[67]-total_remainders,
																																				&new_values[67] - p,
																																				0);
			modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);
		} else { /* Non-empty bucket */

			uint64_t current_remainder, current_count, current_end;

			/* Find the counter for this remainder, if one exists. */
			current_end = decode_counter(qf, runstart_index, &current_remainder,
																	 &current_count);
			while (current_remainder < hash_remainder && !is_runend(qf, current_end)) {
				runstart_index = current_end + 1;
				current_end = decode_counter(qf, runstart_index, &current_remainder,
																		 &current_count);
			}

			/* If we reached the end of the run w/o finding a counter for this remainder,
				 then append a counter for this remainder to the run. */
			if (current_remainder < hash_remainder) {
				uint64_t *p = encode_counter(qf, hash_remainder, count, &new_values[67],&new_fcounters[67]);
				total_remainders=&new_values[67] - p;
				insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf,
																																					1, /* Append to bucket */
																																					hash_bucket_index,
																																					current_end + 1,
																																					p,
																																					&new_fcounters[67]-total_remainders,
																																					&new_values[67] - p,
																																					0);
				modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);
				/* Found a counter for this remainder.  Add in the new count. */
			} else if (current_remainder == hash_remainder) {
				uint64_t tmp= current_count + count;
				if(qf->metadata->maximum_count!=0){
					tmp=std::min(tmp,qf->metadata->maximum_count);
				}
				uint64_t *p = encode_counter(qf, hash_remainder, tmp, &new_values[67],&new_fcounters[67]);
				total_remainders=&new_values[67] - p;
				insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf,
																																					is_runend(qf, current_end) ? 1 : 2,
																																					hash_bucket_index,
																																					runstart_index,
																																					p,
																																					&new_fcounters[67]-total_remainders,
																																					&new_values[67] - p,
																																					current_end - runstart_index + 1);
				/* No counter for this remainder, but there are larger
					 remainders, so we're not appending to the bucket. */
			} else {
				uint64_t *p = encode_counter(qf, hash_remainder, count, &new_values[67],&new_fcounters[67]);
				total_remainders=&new_values[67] - p;
				insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf,
																																					2, /* Insert to bucket */
																																					hash_bucket_index,
																																					runstart_index,
																																					p,
																																					&new_fcounters[67]-total_remainders,
																																					&new_values[67] - p,
																																					0);
				modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);
			}
		}
		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL << (hash_bucket_block_offset % 64);

		/*modify_metadata(qf, &qf->metadata->nelts, count);*/
	}

	if (lock) {
		qf_unlock(qf, hash_bucket_index, false);
	}

	return true;
}

 bool qf_remove(QF *qf, uint64_t hash, uint64_t count , bool lock, bool spin)
{
	uint64_t hash_remainder           = hash & BITMASK(qf->metadata->key_remainder_bits);
	uint64_t hash_bucket_index        = hash >> qf->metadata->key_remainder_bits;
	uint64_t current_remainder, current_count, current_end;
	uint64_t new_values[67];
	uint64_t new_fcounters[67];

	if(hash_bucket_index > qf->metadata->xnslots){
		throw std::out_of_range("Remove function is called with hash index out of range");
	}


	/* Empty bucket */
	if (!is_occupied(qf, hash_bucket_index)){
		return true;
	}

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
	if (current_remainder != hash_remainder){
		return true;
	}

	if (original_runstart_index == runstart_index && is_runend(qf, current_end))
		only_item_in_the_run = 1;


	/* endode the new counter */
	uint64_t *p = encode_counter(qf, hash_remainder,
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
		if (!qf_lock(qf, hash_bucket_index, spin, false))
		return false;
	}

	remove_replace_slots_and_shift_remainders_and_runends_and_offsets(qf,
																																		only_item_in_the_run,
																																		hash_bucket_index,
																																		runstart_index,
																																		p,
																																		&new_fcounters[67]-total_reminders,
																																		total_reminders,
																																		current_end - runstart_index + 1);


  if (lock) {
  qf_unlock(qf, hash_bucket_index, true);
  }

	return true;
	// update the nelements.
	/*modify_metadata(qf, &qf->metadata->nelts, -count);*/
	/*qf->metadata->nelts -= count;*/
}

/***********************************************************************
 * Code that uses the above to implement key-value-counter operations. *
 ***********************************************************************/

void qf_init(QF *qf, uint64_t nslots, uint64_t key_bits, uint64_t tag_bits,uint64_t fixed_counter_size,
						 bool mem, const char * path, uint32_t seed)
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
	//assert (BITS_PER_SLOT == 0 || BITS_PER_SLOT == qf->metadata->bits_per_slot);
	//assert(bits_per_slot > 1);
// #if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64
// 	size = nblocks * sizeof(qfblock) +  (64)*fixed_counter_size;
// #else
// 	size = nblocks * (sizeof(qfblock) + (SLOTS_PER_BLOCK * bits_per_slot / 8) +
// 	fixed_counter_size*8 + tag_bits*8 );
// #endif
//printf("bits per slot =%lu,key remainder bits =%lu, fixed counter =%lu, tag_bits=%lu\n",
//bits_per_slot,key_remainder_bits,fixed_counter_size,tag_bits );
size = nblocks * (sizeof(qfblock) + (8 * bits_per_slot )) ;

qf->mem = (qfmem *)calloc(sizeof(qfmem), 1);

	if (mem) {
		qf->metadata = (qfmetadata *)calloc(sizeof(qfmetadata), 1);
		qf->metadata->mem=mem;
		qf->metadata->size = size;
		qf->metadata->seed = seed;
		qf->metadata->nslots = num_slots;
		qf->metadata->xnslots = qf->metadata->nslots +
			10*sqrt((double)qf->metadata->nslots);
		qf->metadata->key_bits = key_bits;
		qf->metadata->tag_bits = tag_bits;
		qf->metadata->fixed_counter_size = fixed_counter_size;
		qf->metadata->key_remainder_bits = key_remainder_bits;
		qf->metadata->bits_per_slot = bits_per_slot;

		qf->metadata->range = qf->metadata->nslots;
		qf->metadata->range <<= qf->metadata->key_remainder_bits;
		qf->metadata->nblocks = (qf->metadata->xnslots + SLOTS_PER_BLOCK - 1) /
			SLOTS_PER_BLOCK;
		qf->metadata->nelts = 0;
		qf->metadata->ndistinct_elts = 0;
		qf->metadata->noccupied_slots = 0;
		qf->metadata->maximum_occupied_slots=(uint64_t)((double)qf->metadata->xnslots *0.95);
		qf->metadata->num_locks = (qf->metadata->xnslots/NUM_SLOTS_TO_LOCK)+2;
		qf->metadata->maximum_count = 0;
		qf->metadata->tags_map=NULL;
		qf->blocks = (qfblock *)calloc(size, 1);


	} else {

		qf->mem->fd = open(path, O_RDWR | O_CREAT | O_TRUNC, S_IRWXU);
		if (qf->mem->fd < 0) {
			perror("Couldn't open file:\n");
			exit(EXIT_FAILURE);
		}

		/* prashantpandey: Commenting out fallocate call to preallocate space for
		 * the file on disk because fallocate is not supported on MAC OS. Revisit
		 * it later. */
		// int ret;
		// ret = fallocate(qf->mem->fd, 0, 0, size+sizeof(qfmetadata));
		// if (ret < 0) {
		// 	perror("Couldn't fallocate file:\n");
		// 	exit(EXIT_FAILURE);
		// }

		// allocate space for mmaped file
		std::ofstream outputFile(path);
		outputFile.seekp(size+sizeof(qfmetadata));
		outputFile<<0;
		outputFile.close();



		qf->metadata = (qfmetadata *)mmap(NULL, size+sizeof(qfmetadata), PROT_READ |
																			PROT_WRITE, MAP_SHARED, qf->mem->fd, 0);
		qf->metadata->mem=mem;
		qf->metadata->size = size;
		qf->metadata->seed = seed;
		qf->metadata->nslots = num_slots;
		qf->metadata->xnslots = qf->metadata->nslots +
														10*sqrt((double)qf->metadata->nslots);
		qf->metadata->key_bits = key_bits;
		qf->metadata->tag_bits = tag_bits;
		qf->metadata->fixed_counter_size = fixed_counter_size;
		qf->metadata->key_remainder_bits = key_remainder_bits;
		qf->metadata->bits_per_slot = bits_per_slot;

		qf->metadata->range = qf->metadata->nslots;
		qf->metadata->range <<= qf->metadata->key_remainder_bits;
		qf->metadata->nblocks = (qf->metadata->xnslots + SLOTS_PER_BLOCK - 1) /
			SLOTS_PER_BLOCK;
		qf->metadata->nelts = 0;
		qf->metadata->ndistinct_elts = 0;
		qf->metadata->noccupied_slots = 0;
		qf->metadata->maximum_occupied_slots=(uint64_t)((double)qf->metadata->xnslots *0.95);
		qf->metadata->num_locks = (qf->metadata->xnslots/NUM_SLOTS_TO_LOCK)+2;
		qf->metadata->maximum_count = 0;
		qf->metadata->tags_map=NULL;
		qf->blocks = (qfblock *)(qf->metadata + 1);
	}

	/* initialize all the locks to 0 */
	qf->mem->metadata_lock = 0;
	qf->mem->general_lock = 0;
	qf->mem->locks = (volatile int *)calloc(qf->metadata->num_locks,
																					sizeof(volatile int));


#ifdef LOG_WAIT_TIME
	qf->mem->wait_times = (wait_time_data* )calloc(qf->metadata->num_locks+1,
																						sizeof(wait_time_data));
#endif
}

/* The caller should call qf_init on the dest QF before calling this function.
 */
void qf_copy(QF *dest, QF *src)
{
	memcpy(dest->mem, src->mem, sizeof(qfmem));
	memcpy(dest->metadata, src->metadata, sizeof(qfmetadata));
	memcpy(dest->blocks, src->blocks, src->metadata->size);

	if(src->metadata->tags_map!=NULL){
		dest->metadata->tags_map=
		new std::map<uint64_t, std::vector<int> >(*src->metadata->tags_map);
	}
}

/* free up the memory if the QF is in memory.
 * else unmap the mapped memory from pagecache.
 *
 * It does not delete the file on disk for on-disk QF.
 */
void qf_destroy(QF *qf)
{
	assert(qf->blocks != NULL);
	if(qf->metadata->tags_map!=NULL){
		delete qf->metadata->tags_map;
		qf->metadata->tags_map=NULL;
	}
	if (qf->metadata->mem) {
		free(qf->mem);
		free(qf->metadata);
		free(qf->blocks);
	} else {
	msync(qf->metadata, qf->metadata->size + sizeof(qfmetadata),MS_SYNC);
	munmap(qf->metadata, qf->metadata->size + sizeof(qfmetadata));
	close(qf->mem->fd);
	}

	//qf->metadata->noccupied_slots=0;
}

void qf_close(QF *qf)
{
	assert(qf->blocks != NULL);
	munmap(qf->metadata, qf->metadata->size + sizeof(qfmetadata));
	close(qf->mem->fd);
}

/*
 * Will read the on-disk QF using mmap.
 * Data won't be copied in memory.
 *
 */
 void qf_read(QF *qf, const char *path)
 {
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
		 qf->blocks = (qfblock *)(qf->metadata + 1);
		 qf->metadata->num_locks = (qf->metadata->xnslots/NUM_SLOTS_TO_LOCK)+2;
		 qf->mem->metadata_lock = 0;
		 qf->mem->locks = (volatile int *)calloc(qf->metadata->num_locks,
			 sizeof(volatile int));

	string tagsMapOutName=string(path)+".tags_map";
	if(file_exists(tagsMapOutName)){
		qf->metadata->tags_map=load_tags_map(tagsMapOutName.c_str());
	}

}

void qf_reset(QF *qf)
{
	assert(popcnt(nslots) == 1); /* nslots must be a power of 2 */

	qf->metadata->nelts = 0;
	qf->metadata->ndistinct_elts = 0;
	qf->metadata->noccupied_slots = 0;
	if(qf->metadata->tags_map!=NULL)
		qf->metadata->tags_map->clear();
#ifdef LOG_WAIT_TIME
	memset(qf->wait_times, 0, (qf->metadata->num_locks+1)*sizeof(wait_time_data));
#endif
#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64
	memset(qf->blocks, 0, qf->metadata->nblocks* sizeof(qfblock));
#else
	memset(qf->blocks, 0, qf->metadata->nblocks*(sizeof(qfblock) + SLOTS_PER_BLOCK *
																		 qf->metadata->bits_per_slot / 8));
#endif
}

void qf_serialize(const QF *qf, const char *filename)
{
	FILE *fout;
	fout = fopen(filename, "wb+");
	if (fout == NULL) {
		perror("Error opening file for serializing\n");
		exit(EXIT_FAILURE);
	}
	fwrite(qf->metadata, sizeof(qfmetadata), 1, fout);
	/* we don't serialize the locks */
	fwrite(qf->blocks, qf->metadata->size, 1, fout);
	fclose(fout);

	if(qf->metadata->tags_map!=NULL)
	{
		string tagsMapOutName=string(filename)+".tags_map";
		save_tags_map(qf->metadata->tags_map,tagsMapOutName.c_str());
	}
}



void qf_deserialize(QF *qf, const char *filename)
{
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

	qf->blocks = (qfblock *)calloc(qf->metadata->size, 1);
	fread(qf->blocks, qf->metadata->size, 1, fin);
	fclose(fin);

	string tagsMapOutName=string(filename)+".tags_map";
	if(file_exists(tagsMapOutName)){
		qf->metadata->tags_map=load_tags_map(tagsMapOutName.c_str());
	}


}
uint64_t qf_add_tag(const QF *qf, uint64_t key, uint64_t tag, bool lock, bool spin)
{
	if(qf->metadata->tag_bits==0){
		return 0;
	}
	__uint128_t hash = key;
	uint64_t hash_remainder   = hash & BITMASK(qf->metadata->key_remainder_bits);
	int64_t hash_bucket_index = hash >> qf->metadata->key_remainder_bits;


	if (!is_occupied(qf, hash_bucket_index)){
		return 0;
	}

	int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
																																hash_bucket_index-1)
		+ 1;

	if (runstart_index < hash_bucket_index)
		runstart_index = hash_bucket_index;

	/* printf("MC RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

	uint64_t current_remainder, current_count, current_end;
	do {
		current_end = decode_counter(qf, runstart_index, &current_remainder,
																 &current_count);
		if (current_remainder == hash_remainder){
			if (lock) {
				if(qf->mem->general_lock)
					return false;
				if (!qf_lock(qf, runstart_index, spin, false))
				return 0;
			}

			set_tag(qf,runstart_index,tag);
			if (lock) {
				qf_unlock(qf, runstart_index, true);
			}


			return 1;
		}
		runstart_index = current_end + 1;
	} while (!is_runend(qf, current_end));


	return 0;
}

uint64_t qf_remove_tag(const QF *qf, uint64_t key ,bool lock, bool spin)
{

	if(qf->metadata->tag_bits==0){
		return 0;
	}

	__uint128_t hash = key;
	uint64_t hash_remainder   = hash & BITMASK(qf->metadata->key_remainder_bits);
	int64_t hash_bucket_index = hash >> qf->metadata->key_remainder_bits;
	if(hash_bucket_index > qf->metadata->xnslots){
			throw std::out_of_range("qf_remove_tag is called with hash index out of range");
		}

	if (!is_occupied(qf, hash_bucket_index)){
		return 0;
	}
	int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
																																hash_bucket_index-1)
		+ 1;
	if (runstart_index < hash_bucket_index)
		runstart_index = hash_bucket_index;

	/* printf("MC RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

	uint64_t current_remainder, current_count, current_end;
	do {
		current_end = decode_counter(qf, runstart_index, &current_remainder,
																 &current_count);
		if (current_remainder == hash_remainder){
			if (lock) {
				if(qf->mem->general_lock)
					return false;
				if (!qf_lock(qf, runstart_index, spin, false))
					return false;
				}
			set_tag(qf,runstart_index,0);
			if (lock)
				qf_unlock(qf, runstart_index, true);
			return 1;
		}
		runstart_index = current_end + 1;
	} while (!is_runend(qf, current_end));


	return 0;
}

uint64_t qf_get_tag(const QF *qf, uint64_t key)
{
	if(qf->metadata->tag_bits==0){
		return 0;
	}
	__uint128_t hash = key;
	uint64_t hash_remainder   = hash & BITMASK(qf->metadata->key_remainder_bits);
	int64_t hash_bucket_index = hash >> qf->metadata->key_remainder_bits;
	if(hash_bucket_index > qf->metadata->xnslots){
			throw std::out_of_range("qf_get_tag is called with hash index out of range");
		}
	if (!is_occupied(qf, hash_bucket_index))
		return 0;

	int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
																																hash_bucket_index-1)
		+ 1;
	if (runstart_index < hash_bucket_index)
		runstart_index = hash_bucket_index;

	/* printf("MC RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

	uint64_t current_remainder, current_count, current_end;
	do {
		current_end = decode_counter(qf, runstart_index, &current_remainder,
																 &current_count);
		if (current_remainder == hash_remainder){
			return get_tag(qf,runstart_index);
		}
		runstart_index = current_end + 1;
	} while (!is_runend(qf, current_end));

	return 0;
}


bool qf_insert(QF *qf, uint64_t key, uint64_t count, bool
							 lock, bool spin)
{
	if(count==0)
	{
		return true;
	}
	/*uint64_t hash = (key << qf->metadata->tag_bits) | (value & BITMASK(qf->metadata->tag_bits));*/
	if (0 && count == 1)
	 return insert(qf, key,count, lock, spin);
	else
	 return insert(qf, key, count, lock, spin);
}

uint64_t qf_count_key(const QF *qf, uint64_t key)
{
	__uint128_t hash = key;
	uint64_t hash_remainder   = hash & BITMASK(qf->metadata->key_remainder_bits);
	int64_t hash_bucket_index = hash >> qf->metadata->key_remainder_bits;

	if (!is_occupied(qf, hash_bucket_index))
		return 0;

	int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
																																hash_bucket_index-1)
		+ 1;
	if (runstart_index < hash_bucket_index)
		runstart_index = hash_bucket_index;

	/* printf("MC RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

	uint64_t current_remainder, current_count, current_end;
	do {
		current_end = decode_counter(qf, runstart_index, &current_remainder,
																 &current_count);
		if (current_remainder == hash_remainder){
			return current_count;
		}

		runstart_index = current_end + 1;
	} while (!is_runend(qf, current_end));
	return 0;
}




/* initialize the iterator at the run corresponding
 * to the position index
 */
bool qf_iterator(QF *qf, QFi *qfi, uint64_t position)
{
	if(position > qf->metadata->xnslots){
		throw std::out_of_range("qf_iterator is called with position out of range");
	}
	if (!is_occupied(qf, position)) {
		uint64_t block_index = position;
		uint64_t idx = bitselect(get_block(qf, block_index)->occupieds[0], 0);
		if (idx == 64) {
			while(idx == 64 && block_index < qf->metadata->nblocks) {
				block_index++;
				idx = bitselect(get_block(qf, block_index)->occupieds[0], 0);
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

#ifdef LOG_CLUSTER_LENGTH
	qfi->c_info = (cluster_data* )calloc(qf->metadata->nslots/32, sizeof(cluster_data));
	qfi->cur_start_index = position;
	qfi->cur_length = 1;
#endif

	if (qfi->current >= qf->metadata->nslots)
		return false;
	return true;
}

int qfi_get(QFi *qfi, uint64_t *key, uint64_t *value, uint64_t *count)
{

	if(qfi->current > qfi->qf->metadata->xnslots){
		throw std::out_of_range("qfi_get is called with hash index out of range");
	}
	uint64_t current_remainder, current_count;
	decode_counter(qfi->qf, qfi->current, &current_remainder, &current_count);
	*key = (qfi->run << qfi->qf->metadata->key_remainder_bits) | current_remainder;
	*value = get_tag(qfi->qf,qfi->current);   // for now we are not using value
	*count = current_count;

	qfi->qf->metadata->ndistinct_elts++;
	qfi->qf->metadata->nelts += current_count;

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
		qfi->current = decode_counter(qfi->qf, qfi->current, &current_remainder,
																	&current_count);

		if (!is_runend(qfi->qf, qfi->current)) {
			qfi->current++;
#ifdef LOG_CLUSTER_LENGTH
			qfi->cur_length++;
#endif
			if (qfi->current > qfi->qf->metadata->xnslots)
				return 1;
			return 0;
		}
		else {
#ifdef LOG_CLUSTER_LENGTH
			/* save to check if the new current is the new cluster. */
			uint64_t old_current = qfi->current;
#endif
			uint64_t block_index = qfi->run / SLOTS_PER_BLOCK;
			uint64_t rank = bitrank(get_block(qfi->qf, block_index)->occupieds[0],
															qfi->run % SLOTS_PER_BLOCK);
			uint64_t next_run = bitselect(get_block(qfi->qf,
																							block_index)->occupieds[0],
																		rank);
			if (next_run == 64) {
				rank = 0;
				while (next_run == 64 && block_index < qfi->qf->metadata->nblocks) {
					block_index++;
					next_run = bitselect(get_block(qfi->qf, block_index)->occupieds[0],
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
#ifdef LOG_CLUSTER_LENGTH
			if (qfi->current > old_current + 1) { /* new cluster. */
				if (qfi->cur_length > 10) {
					qfi->c_info[qfi->num_clusters].start_index = qfi->cur_start_index;
					qfi->c_info[qfi->num_clusters].length = qfi->cur_length;
					qfi->num_clusters++;
				}
				qfi->cur_start_index = qfi->run;
				qfi->cur_length = 1;
			} else {
				qfi->cur_length++;
			}
#endif
			return 0;
		}
	}
}

inline int qfi_end(QFi *qfi)
{
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
void _qf_merge(QF *qfa, QF *qfb, QF *qfc,
	void(*mergeFn)(uint64_t   keya, uint64_t  tag_a,uint64_t  count_a,
						  	 uint64_t   keyb, uint64_t  tag_b,uint64_t  count_b,
							   uint64_t*  keyc, uint64_t* tag_c,uint64_t* count_c
							 ))
{
	QFi qfia, qfib;
	if(qfa->metadata->range != qfb->metadata->range ||
	qfb->metadata->range != qfc->metadata->range )
	{
		throw std::logic_error("Merging non compatible filters");
	}
	qf_iterator(qfa, &qfia, 0);
	qf_iterator(qfb, &qfib, 0);

	uint64_t keya, taga, counta, keyb, tagb, countb;
	uint64_t keyc,tagc, countc;
	qfi_get(&qfia, &keya, &taga, &counta);
	qfi_get(&qfib, &keyb, &tagb, &countb);

	do {
		if (keya < keyb) {
			mergeFn(keya,taga,counta,0,0,0,&keyc,&tagc,&countc);
			qfi_next(&qfia);
			qfi_get(&qfia, &keya, &taga, &counta);
		}
		else if(keya > keyb) {
			mergeFn(0,0,0,keyb,tagb,countb,&keyc,&tagc,&countc);
			qfi_next(&qfib);
			qfi_get(&qfib, &keyb, &tagb, &countb);
		}
		else{
			mergeFn(keya,taga,counta,keyb,tagb,countb,&keyc,&tagc,&countc);
			qfi_next(&qfia);
			qfi_next(&qfib);
			qfi_get(&qfia, &keya, &taga, &counta);
			qfi_get(&qfib, &keyb, &tagb, &countb);
		}
		if(countc!=0){
			qf_insert(qfc, keyc, countc, true, true);
			qf_add_tag(qfc,keya,tagc);
		}

	} while(!qfi_end(&qfia) && !qfi_end(&qfib));

	if (!qfi_end(&qfia)) {

		do {
			qfi_get(&qfia, &keya, &taga, &counta);
			mergeFn(keya,taga,counta,0,0,0,&keyc,&tagc,&countc);
			if(countc!=0){
				qf_insert(qfc, keyc, countc, true, true);
				qf_add_tag(qfc,keyc,tagc);
			}
		} while(!qfi_next(&qfia));
	}

	if (!qfi_end(&qfib)) {
		do {
			qfi_get(&qfib, &keyb, &tagb, &countb);
			mergeFn(0,0,0,keyb,tagb,countb,&keyc,&tagc,&countc);
			if(countc!=0){
				qf_insert(qfc, keyc, countc, true, true);
				qf_add_tag(qfc,keyc,tagc);
			}
		} while(!qfi_next(&qfib));
	}

	return;
}
void qf_merge(QF *qfa, QF *qfb, QF *qfc)
{
	_qf_merge(qfa,qfb,qfc,unionFn);
}

void qf_intersect(QF *qfa, QF *qfb, QF *qfc)
{
	_qf_merge(qfa,qfb,qfc,intersectFn);
}
void qf_subtract(QF *qfa, QF *qfb, QF *qfc)
{
	_qf_merge(qfa,qfb,qfc,subtractFn);
}

bool qf_equals(QF *qfa, QF *qfb)
{
	QFi qfia, qfib;
	if(qfa->metadata->range != qfb->metadata->range  )
	{
		throw std::logic_error("comparing non compatible filters");
	}
	qf_iterator(qfa, &qfia, 0);
	qf_iterator(qfb, &qfib, 0);

	uint64_t keya, valuea, counta, keyb, valueb, countb;
	qfi_get(&qfia, &keya, &valuea, &counta);
	qfi_get(&qfib, &keyb, &valueb, &countb);
	do {
		if (keya != keyb) {
			return false;
		}
		else {
			if(counta!=countb || valuea != valueb)
			{
				return false;
			}
			qfi_next(&qfib);
			qfi_next(&qfia);
			qfi_get(&qfia, &keya, &valuea, &counta);
			qfi_get(&qfib, &keyb, &valueb, &countb);
		}
	} while(!qfi_end(&qfia) && !qfi_end(&qfib));

	if (!qfi_end(&qfia) || !qfi_end(&qfib)) {
		return false;
	}

	return true;
}



std::map<std::string, uint64_t> Tags_map;
uint64_t last_index=0;
void union_multi_Fn(uint64_t   key_arr[], uint64_t  tag_arr[],uint64_t  count_arr[]
	,std::map<uint64_t, std::vector<int> > ** inverted_indexes,int nqf,
							 uint64_t*  key_c, uint64_t* tag_c,uint64_t* count_c)
{

	*count_c=0;
	for(int i=0;i<nqf;i++)
	{
		//printf("key =%lu, count=%lu\n", key_arr[i],count_arr[i]);
		if(count_arr[i]!=0)
		{
			*key_c=key_arr[i];
			*tag_c=tag_arr[i];
			*count_c+=count_arr[i];
		}
	}

}

void _qf_multi_merge(QF *qf_arr[],int nqf, QF *qfr,
	void(*mergeFn)(uint64_t   key_arr[], uint64_t  tag_arr[],uint64_t  count_arr[],
								std::map<uint64_t, std::vector<int> > ** inverted_indexes,int nqf,
							   uint64_t*  keyc, uint64_t* tag_c,uint64_t* count_c
							 ))
{
	int i;
	uint64_t range=qf_arr[0]->metadata->range;
	for (i=1; i<nqf; i++) {
		if(qf_arr[i]->metadata->range!=range)
		{
			throw std::logic_error("Merging non compatible filters");
		}
	}

	QFi *qfi_arr[nqf];

	uint64_t smallest_key=UINT64_MAX,second_smallest_key;
	uint64_t keys[nqf];
	uint64_t tags[nqf];
	uint64_t counts[nqf];
	std::map<uint64_t, std::vector<int> > ** inverted_indexes=new std::map<uint64_t, std::vector<int> >*[nqf];

	for (i=0; i<nqf; i++) {
		qfi_arr[i]=new QFi();
		qf_iterator(qf_arr[i], qfi_arr[i], 0);
		qfi_get(qfi_arr[i], &keys[i], &tags[i], &counts[i]);
		smallest_key=std::min(keys[i],smallest_key);
		inverted_indexes[i]=qf_arr[i]->metadata->tags_map;
	}

	uint64_t keys_m[nqf];
	uint64_t tags_m[nqf];
	uint64_t counts_m[nqf];


	bool finish=false;
	while(!finish)
	{
		finish=true;
		second_smallest_key=UINT64_MAX;
		//printf("smallest_key = %llu\n",smallest_key );
		for(i=0;i<nqf;i++)
		{
			keys_m[i]=0;
			counts_m[i]=0;
			tags_m[i]=0;

			//printf(" key = %llu\n",keys[i]);
			if(keys[i]==smallest_key){
				keys_m[i]=keys[i];
				counts_m[i]=counts[i];
				tags_m[i]=tags[i];
				qfi_next(qfi_arr[i]);
				if(!qfi_end(qfi_arr[i]))
				{
					finish=false;
					qfi_get(qfi_arr[i], &keys[i], &tags[i], &counts[i]);
				}else{
					keys[i]=UINT64_MAX;
				}
			}
			second_smallest_key=std::min(second_smallest_key,keys[i]);
		}
		for (i = 0; i < nqf; i++) {
			if(keys[i]!=UINT64_MAX)
			{
				finish=false;
				break;
			}
		}
		//printf("second_smallest_key=%llu finish=%d\n",second_smallest_key,finish);
		uint64_t keyc,tagc, countc;
		mergeFn(keys_m,tags_m,counts_m,inverted_indexes,nqf,&keyc,&tagc,&countc);

		if(countc!=0){
			qf_insert(qfr, keyc, countc, true, true);
			qf_add_tag(qfr,keyc,tagc);
		}
		smallest_key=second_smallest_key;
	}
	// cout<<"before delete"<<endl;
	delete  inverted_indexes;
	for(i=0;i<nqf;i++)
	{
		delete qfi_arr[i];
	}

	return;
}

/*
 * Merge an array of qfs into the resultant QF
 */
void qf_multi_merge(QF *qf_arr[], int nqf, QF *qfr)
{
	_qf_multi_merge(qf_arr,nqf,qfr,union_multi_Fn);

}

void inverted_union_multi_Fn(uint64_t   key_arr[], uint64_t  tag_arr[],uint64_t  count_arr[],
	std::map<uint64_t, std::vector<int> > ** inverted_indexes ,int nqf,
							 uint64_t*  key_c, uint64_t* tag_c,uint64_t* count_c)
{

	std::string index_key="";
	*count_c=0;
	for(int i=0;i<nqf;i++)
	{
		if(count_arr[i]!=0)
		{
			*key_c=key_arr[i];
			*count_c+=count_arr[i];
			if(inverted_indexes==NULL){
				index_key+=std::to_string(i);
				index_key+=';';
			}
			else{
				auto it=inverted_indexes[i]->find(tag_arr[i]);
				for(auto k:it->second){
					index_key+=std::to_string(k);
					index_key+=';';
				}
			}

		}
	}
	index_key.pop_back();
	auto it=Tags_map.find(index_key);
	if(it==Tags_map.end())
	{

		Tags_map.insert(std::make_pair(index_key,Tags_map.size()));
		it=Tags_map.find(index_key);
	}
	*tag_c=it->second;

}

void inverted_union_multi_no_count_Fn(uint64_t   key_arr[], uint64_t  tag_arr[],uint64_t  count_arr[],
									std::map<uint64_t, std::vector<int> > ** inverted_indexes, int nqf,
							 uint64_t*  key_c, uint64_t* tag_c,uint64_t* count_c)
{
	std::string index_key="";
	*count_c=0;
	for(int i=0;i<nqf;i++)
	{
		//printf("key =%lu, count=%lu\n", key_arr[i],count_arr[i]);
		if(count_arr[i]!=0)
		{
			*key_c=key_arr[i];
			if(inverted_indexes==NULL){
				index_key+=std::to_string(i);
				index_key+=';';
			}
			else{
				auto it=inverted_indexes[i]->find(tag_arr[i]);
				for(auto k:it->second){
					index_key+=std::to_string(k);
					index_key+=';';
				}
			}

		}
	}
	index_key.pop_back();
	auto it=Tags_map.find(index_key);
	if(it==Tags_map.end())
	{

		Tags_map.insert(std::make_pair(index_key,last_index));
		last_index++;
		it=Tags_map.find(index_key);
	}
	*count_c=it->second;

}


void qf_invertable_merge(QF *qf_arr[], int nqf, QF *qfr)
{
	int i;
	int last_tag=0;
	Tags_map.clear();
	last_index=0;
	for(i=0;i<nqf;i++){
		if(qf_arr[i]->metadata->tags_map==NULL){
			qf_arr[i]->metadata->tags_map=new std::map<uint64_t, std::vector<int> >();
			vector<int> tmp(1);
			tmp[0]=last_index;
			Tags_map.insert(std::make_pair(std::to_string(i),last_index++));
			qf_arr[i]->metadata->tags_map->insert(make_pair(0,tmp));
		}
		else{
			auto it=qf_arr[i]->metadata->tags_map->begin();
			int updated_tags=0;
			while(it!=qf_arr[i]->metadata->tags_map->end()){
				for(int j=0;j<it->second.size();j++){
					it->second[j]+=last_tag;
					auto it2=Tags_map.find(std::to_string(it->second[j]));
					if(it2==Tags_map.end()){
						Tags_map.insert(std::make_pair(std::to_string(it->second[j]),it->second[j]));
						updated_tags++;
					}
				}
				it++;
			}
		}
		last_tag+=Tags_map.size();
	}



	_qf_multi_merge(qf_arr,nqf,qfr,inverted_union_multi_Fn);
	qfr->metadata->tags_map=new std::map<uint64_t, std::vector<int> >();
	auto it=Tags_map.begin();
	while(it!=Tags_map.end()){
		std::vector<int> tmp=key_to_vector_int(it->first);
		qfr->metadata->tags_map->insert(std::make_pair(it->second,tmp));
		it++;

	}


}

void qf_invertable_merge_no_count(QF *qf_arr[], int nqf, QF *qfr)
{

	int i;
	int last_tag=0;
	Tags_map.clear();
	last_index=0;
	for(i=0;i<nqf;i++){
		if(qf_arr[i]->metadata->tags_map==NULL){
			qf_arr[i]->metadata->tags_map=new std::map<uint64_t, std::vector<int> >();
			vector<int> tmp(1);
			tmp[0]=last_index;
			Tags_map.insert(std::make_pair(std::to_string(i),last_index++));
			qf_arr[i]->metadata->tags_map->insert(make_pair(0,tmp));
		}
		else{
			auto it=qf_arr[i]->metadata->tags_map->begin();
			int updated_tags=0;
			while(it!=qf_arr[i]->metadata->tags_map->end()){
				for(int j=0;j<it->second.size();j++){
					it->second[j]+=last_tag;
					auto it2=Tags_map.find(std::to_string(it->second[j]));
					if(it2==Tags_map.end()){
						Tags_map.insert(std::make_pair(std::to_string(it->second[j]),it->second[j]));
						updated_tags++;
					}
				}
				it++;
			}
		}
		last_tag+=Tags_map.size();
	}


	_qf_multi_merge(qf_arr,nqf,qfr,inverted_union_multi_no_count_Fn);

	qfr->metadata->tags_map=new std::map<uint64_t, std::vector<int> >();
	auto it=Tags_map.begin();
	while(it!=Tags_map.end()){
		std::vector<int> tmp=key_to_vector_int(it->first);
		qfr->metadata->tags_map->insert(std::make_pair(it->second,tmp));
		it++;
	}


}

QF* qf_resize(QF* qf, int newQ, const char * originalFilename, const char * newFilename)
{
	if((int)qf->metadata->key_bits-newQ <2)
	{
		throw std::logic_error("Resize cannot be done. Slot size cannot be less than 2");
	}

	if(originalFilename)
	{
		qf_serialize(qf,originalFilename);
		qf_destroy(qf);
		qf_read(qf,originalFilename);
	}
	QF* newQF=(QF *)calloc(sizeof(QF), 1);
	if(newFilename)
	{
		qf_init(newQF, (1ULL<<newQ),qf->metadata->key_bits, qf->metadata->tag_bits,qf->metadata->fixed_counter_size, false, newFilename, 2038074761);
	}
	else{
		qf_init(newQF, (1ULL<<newQ),qf->metadata->key_bits, qf->metadata->tag_bits,qf->metadata->fixed_counter_size, true, "" , 2038074761);
	}
	QFi qfi;
	qf_iterator(qf, &qfi, 0);


	uint64_t keya, valuea, counta;
	qfi_get(&qfi, &keya, &valuea, &counta);

	do {
			qf_insert(newQF, keya, counta);
			qf_add_tag(newQF,keya,valuea);
			qfi_next(&qfi);
			qfi_get(&qfi, &keya, &valuea, &counta);
	} while(!qfi_end(&qfi));
	qf_destroy(qf);
	return newQF;


}

/* find cosine similarity between two QFs. */
uint64_t qf_inner_product(QF *qfa, QF *qfb)
{
	uint64_t acc = 0;
	QFi qfi;
	QF *qf_mem, *qf_disk;

	// create the iterator on the larger QF.
	if (qfa->metadata->size > qfb->metadata->size) {
		qf_mem = qfb;
		qf_disk = qfa;
	} else {
		qf_mem = qfa;
		qf_disk = qfb;
	}

	qf_iterator(qf_disk, &qfi, 0);
	do {
		uint64_t key = 0, value = 0, count = 0;
		uint64_t count_mem;
		qfi_get(&qfi, &key, &value, &count);
		if ((count_mem = qf_count_key(qf_mem, key)) > 0) {
			acc += count*count_mem;
		}
	} while (!qfi_next(&qfi));

	return acc;
}

/* find cosine similarity between two QFs. */
// void qf_intersect(QF *qfa, QF *qfb, QF *qfr)
// {
// 	QFi qfi;
// 	QF *qf_mem, *qf_disk;
//
// 	// create the iterator on the larger QF.
// 	if (qfa->metadata->size > qfb->metadata->size) {
// 		qf_mem = qfb;
// 		qf_disk = qfa;
// 	} else {
// 		qf_mem = qfa;
// 		qf_disk = qfb;
// 	}
//
// 	qf_iterator(qf_disk, &qfi, 0);
// 	do {
// 		uint64_t key = 0, value = 0, count = 0;
// 		qfi_get(&qfi, &key, &value, &count);
// 		if (qf_count_key(qf_mem, key) > 0)
// 			qf_insert(qfr, key, count, false, false);
// 	} while (!qfi_next(&qfi));
// }

/* magnitude of a QF. */
uint64_t qf_magnitude(QF *qf)
{
	return sqrt(qf_inner_product(qf, qf));
}


int qf_space(QF *qf)
{
		return (int)(((double)qf->metadata->noccupied_slots/
								 (double)qf->metadata->xnslots
							 )* 100.0);
}


bool qf_general_lock(QF* qf, bool spin){
	if (!qf_spin_lock(&qf->mem->general_lock, spin))
		return false;
	return true;
}
void qf_general_unlock(QF* qf){
	qf_spin_unlock(&qf->mem->general_lock);
}
void qf_migrate(QF* source, QF* dest){
	QFi source_i;
	if (qf_iterator(source, &source_i, 0)) {
		do {
			uint64_t key = 0, value = 0, count = 0;
			qfi_get(&source_i, &key, &value, &count);
			qf_insert(dest, key, count, true, true);
			qf_add_tag(dest,key,value);
		} while (!qfi_next(&source_i));
	}
}

#ifdef TEST
	#include "tests/lowLevelTests.hpp"
#endif
