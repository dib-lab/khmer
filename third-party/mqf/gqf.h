#ifndef QF_H
#define QF_H

#include <inttypes.h>
#include <stdbool.h>
#include <pthread.h>
#include <map>
#include <vector>

#ifdef __cplusplus
extern "C" {
#endif

/* Can be
	0 (choose size at run-time),
	8, 16, 32, or 64 (for optimized versions),
	or other integer <= 56 (for compile-time-optimized bit-shifting-based versions)
	*/
#define BITS_PER_SLOT 0

	struct __attribute__ ((__packed__)) qfblock;
	typedef struct qfblock qfblock;

	uint64_t shift_into_b2(uint64_t a, uint64_t b, int bstart, int bend, int amount);

	typedef struct {
		uint64_t total_time_single;
		uint64_t total_time_spinning;
		uint64_t locks_taken;
		uint64_t locks_acquired_single_attempt;
	} wait_time_data;

	typedef struct quotient_filter_mem {
		int fd;
		volatile int general_lock;
		volatile int metadata_lock;
		volatile int *locks;
		wait_time_data *wait_times;
	} quotient_filter_mem;

	typedef quotient_filter_mem qfmem;

	typedef struct quotient_filter_metadata {
		uint64_t size;
		uint32_t seed;
		uint64_t nslots;
		uint64_t xnslots;
		uint64_t key_bits;
		uint64_t tag_bits;
		uint64_t fixed_counter_size;
		uint64_t key_remainder_bits;
		uint64_t bits_per_slot;
		__uint128_t range;
		uint64_t nblocks;
		uint64_t nelts;
		uint64_t ndistinct_elts;
		uint64_t noccupied_slots;
		uint64_t maximum_occupied_slots;
		uint64_t num_locks;
		uint64_t maximum_count;
		bool mem;
		std::map<uint64_t, std::vector<int> > * tags_map;
	} quotient_filter_metadata;

	typedef quotient_filter_metadata qfmetadata;

	typedef struct quotient_filter {
		qfmem *mem;
		qfmetadata *metadata;
		qfblock *blocks;
	} quotient_filter;

	typedef quotient_filter QF;

	typedef struct {
		uint64_t start_index;
		uint16_t length;
	} cluster_data;

	typedef struct quotient_filter_iterator {
		QF *qf;
		uint64_t run;
		uint64_t current;
		uint64_t cur_start_index;
		uint16_t cur_length;
		uint32_t num_clusters;
		cluster_data *c_info;
	} quotient_filter_iterator;

	typedef quotient_filter_iterator QFi;

	/*!
	@breif initialize mqf .

	@param Qf* qf : pointer to the Filter.
	@param uint64_t nslots : Number of slots in the filter. Maximum number of items to be inserted depends on this number.
	@param uint64_t key_bits: Number of bits in the hash values. This number should equal log2(nslots) +r. Accuracy depends on r.
	@param uint64_t tag_bits: Number of bits in tag value.
	@param uint64_t fixed_counter_size: Fixed counter size. must be > 0.
	@param bool mem: Flag to create the filter on memeory. IF false, mmap is used.
	@param const char * path: In case of mmap. Path of the file used to pack the filter.
	@param uint32_t seed: useless value. To be removed
		  */
	void qf_init(QF *qf, uint64_t nslots, uint64_t key_bits, uint64_t tag_bits,uint64_t fixed_counter_size, bool mem, const char *path, uint32_t seed);

	void qf_reset(QF *qf);

	void qf_destroy(QF *qf);

	void qf_copy(QF *dest, QF *src);

	/*!
	 	@breif Increment the counter for this item by count.

		@param Qf* qf : pointer to the Filter
		@param uint64_t key : hash of the item to be insertedItems
		@param uint64_t count: Count to be added
		@param bool lock: For Multithreading, Lock the slot used by the current thread so that other threads can't change the value
		@param bool spin: For Multithreading, If there is a lock on the target slot. wait until the lock is freed and insert the count.

		@return bool: True if the item is inserted correctly.
	 */
	bool qf_insert(QF *qf, uint64_t key, uint64_t count,
								 bool lock=false, bool spin=false);


	/* Remove all instances of this key/value pair. */
	//void qf_delete_key_value(QF *qf, uint64_t key, uint64_t value);

	/* Remove all instances of this key. */
	//void qf_delete_key(QF *qf, uint64_t key);

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

	/*!
	@breif Return the number of times key has been inserted, with any value,
		 into qf.

	@param Qf* qf : pointer to the Filter.
	@param uint64_t key : hash of the item.

	@return uint64_t the count associated with the input key.
		  */
	uint64_t qf_count_key(const QF *qf, uint64_t key);

	/*!
	@breif Decrement the counter for this item by count.

	@param	Qf* qf : pointer to the Filter
	@param	uint64_t key : hash of the item to be removed
	@param	uint64_t count: Count to be removed
	@param	bool lock: For Multithreading, Lock the slot used by the current thread so that other threads can't change the value
	@param	bool spin: For Multithreading, If there is a lock on the target slot. wait until the lock is freed and insert the count.

	@return bool: Returns true if the item is removed successfully.
	 */
	bool qf_remove(QF *qf, uint64_t hash, uint64_t count,  bool lock=false, bool spin=false);


	/*!
		@breif Add Tag to item.

		@param Qf* qf : pointer to the Filter
		@param uint64_t key : hash of the item to be insertedItems
		@param uint64_t tag: tag to be added
		@param bool lock: For Multithreading, Lock the slot used by the current thread so that other threads can't change the value
		@param bool spin: For Multithreading, If there is a lock on the target slot. wait until the lock is freed and insert the count.

		@return bool: True if the item is inserted correctly.
	 */
	uint64_t qf_add_tag(const QF *qf, uint64_t key, uint64_t tag, bool lock=false, bool spin=false);
	/*!
	@breif Return the tag associated with a given item.

	@param Qf* qf : pointer to the Filter.
	@param uint64_t key : hash of the item.

	@return uint64_t the tag associated with the input key.
			*/
	uint64_t qf_get_tag(const QF *qf, uint64_t key);
	/*!
	@breif delete the tag associated with a given item.

	@param Qf* qf : pointer to the Filter.
	@param uint64_t key : hash of the item.

	@return bool: Returns true if the item is removed successfully.
			*/
	uint64_t qf_remove_tag(const QF *qf, uint64_t key, bool lock=false, bool spin=false);

	/* Initialize an iterator */
	bool qf_iterator(QF *qf, QFi *qfi, uint64_t position);

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

	/*! write data structure of to the disk */
	void qf_serialize(const QF *qf, const char *filename);

	/* read data structure off the disk */
	void qf_deserialize(QF *qf, const char *filename);

	/* mmap the QF from disk. */
	void qf_read(QF *qf, const char *path);

	/* merge two QFs into the third one. */
	void qf_merge(QF *qfa, QF *qfb, QF *qfc);

	void qf_intersect(QF *qfa, QF *qfb, QF *qfc);

	void qf_subtract(QF *qfa, QF *qfb, QF *qfc);
	/* merge multiple QFs into the final QF one. */
	void qf_multi_merge(QF *qf_arr[], int nqf, QF *qfr);


	/*! @breif Invertiable merge function adds tag for each key and creates index structure. The index is map of an integer and vector of integers where the integer is the value of the tags and vector on integers is the ids of the source filters.

	@param Qf* qf_arr : input array of filters
	@param int nqf: number of filters
	@param QF* qfr: pointer to the output filter.
	@param std::map<uint64_t, std::vector<int> > *inverted_index_ptr: Pointer to the output index.
	*/
	void qf_invertable_merge(QF *qf_arr[], int nqf, QF *qfr);
	void qf_invertable_merge_no_count(QF *qf_arr[], int nqf, QF *qfr);


	/*! @breif Resize the filter into a bigger or smaller one

	@param Qf* qf : pointer to the Filter
	@param uint64_t newQ: new number of slots(Q). the slot size will be recalculated to keep the range constant.
	@param string originalFilename(optional): dump the current filter to the disk to free space for the new filter. Filename is provided as the content of the string.
	@param string newFilename(optional): the new filter is created on disk. Filename is provided as the content of the string.

	@return QF: New Quotient Filter.
	*/
	QF* qf_resize(QF* qf, int newQ, const char * originalFilename=NULL, const char * newFilename=NULL);
	/* find cosine similarity between two QFs. */
	uint64_t qf_inner_product(QF *qfa, QF *qfb);

	/* magnitude of a QF. */
	uint64_t qf_magnitude(QF *qf);
	/* return the filled space(percent) */
	int qf_space(QF *qf);

	bool qf_equals(QF *qfa, QF *qfb);

	bool qf_general_lock(QF* qf, bool spin);
	void qf_general_unlock(QF* qf);

	void qf_migrate(QF* source, QF* destination);

#ifdef __cplusplus
}
#endif

#endif /* QF_H */
