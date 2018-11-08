#ifndef layeredMQF_H
#define layeredMQF_H

#include <inttypes.h>
#include <stdbool.h>
#include <pthread.h>
#include "gqf.h"
#ifdef __cplusplus
extern "C" {
#endif


	typedef class layeredMQF {
	public:
		QF* firstLayer_singletons;
		QF* secondLayer;
		layeredMQF(){
			firstLayer_singletons=new QF();
			secondLayer=new QF();
		}
		~layeredMQF()
		{
			delete firstLayer_singletons;
			delete secondLayer;
		}
	} layeredMQF;




	typedef struct layeredMQFIterator {
		QFi firstLayerIterator;
		QFi secondLayerIterator;
	} layeredMQF_iterator;


	void layeredMQF_init(layeredMQF *qf, uint64_t nslots_singletons ,uint64_t nslots, uint64_t key_bits, uint64_t value_bits,uint64_t fixed_counter_size, bool mem, const char *path, uint32_t seed);

	void layeredMQF_reset(layeredMQF *qf);

	void layeredMQF_destroy(layeredMQF *qf);

	void layeredMQF_copy(layeredMQF *dest, layeredMQF *src);

	/* Increment the counter for this key/value pair by count. */
	bool layeredMQF_insert(layeredMQF *qf, uint64_t key, uint64_t count,
								 bool lock, bool spin);

	/* Remove count instances of this key/value combination. */
	bool layeredMQF_remove(QF *qf, uint64_t hash, uint64_t count,  bool lock=false, bool spin=false);


		/*!
			@breif Add Tag to item.

			@param Qf* qf : pointer to the Filter
			@param uint64_t key : hash of the item to be insertedItems
			@param uint64_t tag: tag to be added
			@param bool lock: For Multithreading, Lock the slot used by the current thread so that other threads can't change the value
			@param bool spin: For Multithreading, If there is a lock on the target slot. wait until the lock is freed and insert the count.

			@return bool: True if the item is inserted correctly.
		 */
		uint64_t layeredMQF_add_tag(const QF *qf, uint64_t key, uint64_t tag, bool lock=false, bool spin=false);
		/*!
		@breif Return the tag associated with a given item.

		@param Qf* qf : pointer to the Filter.
		@param uint64_t key : hash of the item.

		@return uint64_t the tag associated with the input key.
				*/
		uint64_t layeredMQF_get_tag(const QF *qf, uint64_t key);
		/*!
		@breif delete the tag associated with a given item.

		@param Qf* qf : pointer to the Filter.
		@param uint64_t key : hash of the item.

		@return bool: Returns true if the item is removed successfully.
				*/
		uint64_t layeredMQF_remove_tag(const QF *qf, uint64_t key, bool lock=false, bool spin=false);



	/* Return the number of times key has been inserted, with any value,
		 into qf. */
	uint64_t layeredMQF_count_key(const layeredMQF *qf, uint64_t key);



	/* Initialize an iterator */
	bool layeredMQF_qf_iterator(layeredMQF *qf, layeredMQFIterator *qfi, uint64_t position);

	/* Returns 0 if the iterator is still valid (i.e. has not reached the
		 end of the QF. */
	int layeredMQF_qfi_get(layeredMQFIterator *qfi, uint64_t *key, uint64_t *value, uint64_t *count);

	/* Advance to next entry.  Returns whether or not another entry is
		 found.  */
	int layeredMQF_qfi_next(layeredMQFIterator *qfi);

	/* Check to see if the if the end of the QF */
	int layeredMQF_qfi_end(layeredMQFIterator *qfi);

	/* For debugging */
	void layeredMQF_dump(const layeredMQF *);

	/* write data structure of to the disk */
	void layeredMQF_serialize(const layeredMQF *qf, const char *filename);

	/* read data structure off the disk */
	void layeredMQF_deserialize(layeredMQF *qf, const char *filename);

	/* mmap the QF from disk. */
	void layeredMQF_read(layeredMQF *qf, const char *path);

	/* merge two QFs into the third one. */
	//void layeredMQF_merge(layeredMQF *layeredMQFa, layeredMQF *layeredMQFb, QF *qfc);

	/* merge multiple QFs into the final QF one. */
	//void qf_multi_merge(QF *qf_arr[], int nqf, QF *qfr);

	/* find cosine similarity between two QFs. */
//	uint64_t qf_inner_product(QF *qfa, QF *qfb);

	/* magnitude of a QF. */
	//uint64_t qf_magnitude(QF *qf);

	int layeredMQF_space(layeredMQF *qf);

#ifdef __cplusplus
}
#endif

#endif /* layeredMQF_H */
