/*
 * =====================================================================================
 *
 *       Filename:  gqf.h
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

#ifndef QF_H
#define QF_H

#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

#define BITS_PER_SLOT 8

/* Must be >= 6.  6 seems fastest. */
#define BLOCK_OFFSET_BITS (6)

#define SLOTS_PER_BLOCK (1ULL << BLOCK_OFFSET_BITS)
#define METADATA_WORDS_PER_BLOCK ((SLOTS_PER_BLOCK + 63) / 64)

	typedef struct __attribute__ ((__packed__)) qfblock {
		uint8_t offset; /* Code works with uint16_t, uint32_t, etc, but uint8_t seems just as fast as anything else */
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

	uint64_t shift_into_b2(uint64_t a, uint64_t b, int bstart, int bend, int amount);

#ifdef LOG_NUM_SHIFTS
#define len 5000
int shift_count[len];
#endif

	typedef struct quotient_filter {
		uint64_t nslots;
		uint64_t xnslots;
		uint64_t key_bits;
		uint64_t value_bits;
		uint64_t key_remainder_bits;
		uint64_t bits_per_slot;
		__uint128_t range;
		uint64_t nblocks;
		uint64_t nelts;
		uint64_t ndistinct_elts;
		uint64_t noccupied_slots;
		qfblock *blocks;
	} quotient_filter;


	typedef quotient_filter QF;

	typedef struct quotient_filter_iterator {
		const QF *qf;
		uint64_t run;
		uint64_t current;
	} quotient_filter_iterator;

	typedef quotient_filter_iterator QFi;

	void qf_init(QF *qf, uint64_t nslots, uint64_t key_bits, uint64_t value_bits);

	void qf_destroy(QF *qf);

	/* Increment the counter for this key/value pair by count. */
	void qf_insert(QF *qf, uint64_t key, uint64_t value, uint64_t count);

	/* Remove count instances of this key/value combination. */
	void qf_remove(QF *qf, uint64_t key, uint64_t value, uint64_t count);

	/* Remove all instances of this key/value pair. */
	void qf_delete_key_value(QF *qf, uint64_t key, uint64_t value);

	/* Remove all instances of this key. */
	void qf_delete_key(QF *qf, uint64_t key);

	/* Replace the association (key, oldvalue, count) with the association
		 (key, newvalue, count). If there is already an association (key,
		 newvalue, count'), then the two associations will be merged and
		 their counters will be summed, resulting in association (key,
		 newvalue, count' + count). */
	void qf_replace(QF *qf, uint64_t key, uint64_t oldvalue, uint64_t newvalue);

	/* Lookup the value associated with key.  Returns the count of that
		 key/value pair in the QF.  If it returns 0, then, the key is not
		 present in the QF. Only returns the first value associated with key
		 in the QF.  If you want to see others, use an iterator. */
	uint64_t qf_query(const QF *qf, uint64_t key, uint64_t *value);

	/* Return the number of times key has been inserted, with any value,
		 into qf. */
	uint64_t qf_count_key(const QF *qf, uint64_t key);

	/* Return the number of times key has been inserted, with the given
		 value, into qf. */
	uint64_t qf_count_key_value(const QF *qf, uint64_t key, uint64_t value);

	/* Initialize an iterator */
	void qf_iterator(const QF *qf, QFi *qfi, uint64_t position);

	/* Returns 0 if the iterator is still valid (i.e. has not reached the
		 end of the QF. */
	int qfi_get(QFi *qfi, uint64_t *key, uint64_t *value, uint64_t *count);

	/* Advance to next entry.  Returns whether or not another entry is
		 found.  */
	int qfi_next(QFi *qfi);

	/* Check to see if the if the end of the QF */
	int qfi_end(QFi *qfi);

	/* For debugging */
	void qf_dump(const QF *);

	/* write data structure of to the disk */
	void qf_serialize(const QF *qf, const char *filename);

	/* read data structure off the disk */
	void qf_deserialize(QF *qf, const char *filename);

	/* merge two QFs into the third one. */
	void qf_merge(const QF *qfa, const QF *qfb, QF *qfc);

	/* merge multiple QFs into the final QF one. */
	void qf_multi_merge(QF *qf_arr[], int nqf, QF *qfr);

#ifdef __cplusplus
}
#endif

#endif /* QF_H */
