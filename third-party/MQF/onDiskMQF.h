#ifndef onDiskMQF_H
#define onDiskMQF_H

#include <inttypes.h>
#include <stdbool.h>
#include <pthread.h>
#include <map>
#include <vector>
#include "gqf.h"
#include <fstream>
#include <stxxl/vector>


namespace onDiskMQF_Namespace{
/* Can be
	0 (choose size at run-time),
	8, 16, 32, or 64 (for optimized versions),
	or other integer <= 56 (for compile-time-optimized bit-shifting-based versions)
	*/
	#define SLOTS_PER_BLOCK (64)
	#define METADATA_WORDS_PER_BLOCK 1

	struct __attribute__ ((__packed__)) qfblock;
	typedef struct qfblock qfblock;

	template<uint64_t bitsPerSlot>
	class _onDiskMQF;


	typedef struct diskParameters {
		uint64_t nBlocksOnMemory;
		uint64_t sizeOnMemory;
		uint8_t nBlocksPerPointer;
		uint64_t memoryBufferPos;
		uint64_t nBlocksPerIOBatch;
		uint64_t blocksPointersLen;
	} diskParameters;


//	const uint8_t bitsPerSlot=15;
	template<uint64_t bitsPerSlot=64>
	class onDisk_qfblock{
	public:
		uint8_t offset;
		uint64_t occupieds[METADATA_WORDS_PER_BLOCK];
		uint64_t runends[METADATA_WORDS_PER_BLOCK];
		uint8_t slots[8*bitsPerSlot];
		onDisk_qfblock(){
			offset=0;
			memset(occupieds,0,METADATA_WORDS_PER_BLOCK*sizeof(uint64_t));
			memset(runends,0,METADATA_WORDS_PER_BLOCK*sizeof(uint64_t));
			memset(slots,0,8*bitsPerSlot);
		}
		// template<uint64_t bits2>
		// operator onDisk_qfblock<bits2>() {
    //
		// }
		// onDisk_qfblock(const onDisk_qfblock &obj) noexcept
		// {
		// 	offset=obj.offset;
		// 	memcpy(occupieds, obj.occupieds, METADATA_WORDS_PER_BLOCK*sizeof(uint64_t) );
		// 	memcpy(runends, obj.runends, METADATA_WORDS_PER_BLOCK*sizeof(uint64_t) );
		// 	memcpy(slots,obj.slots,8*bitsPerSlot);
		// }
		// onDisk_qfblock& operator=( onDisk_qfblock const& obj) noexcept {
		// 	this->offset=obj.offset;
		// 	memcpy(this->occupieds, obj.occupieds, METADATA_WORDS_PER_BLOCK*sizeof(uint64_t) );
		// 	memcpy(this->runends, obj.runends, METADATA_WORDS_PER_BLOCK*sizeof(uint64_t) );
		// 	memcpy(this->slots,obj.slots,8*bitsPerSlot);
	 	// 	return *this;
 		// }
	};


  template<uint64_t bitsPerSlot>
	class onDiskMQFIterator {
	public:
		_onDiskMQF<bitsPerSlot> *qf;
		uint64_t run;
		uint64_t current;
		uint64_t cur_start_index;
		uint16_t cur_length;
		uint32_t num_clusters;
		cluster_data *c_info;
		/* Returns 0 if the iterator is still valid (i.e. has not reached the
			 end of the QF. */
		int get(uint64_t *key, uint64_t *value, uint64_t *count);

		/* Advance to next entry.  Returns whether or not another entry is
			 found.  */
		int next();

		/* Check to see if the if the end of the QF */
		int end();

	};

	class onDiskMQF {
	public:
		qfmem *mem;
		qfmetadata *metadata;
		uint64_t stxxlBufferSize;
		uint64_t* blocksFilePos;
		uint64_t reverseBlocksPointer;
		std::fstream diskMQFStream;
		qfblock **blocksPointers;
		diskParameters* diskParams;
		static void init( onDiskMQF *&qf, uint64_t nslots, uint64_t key_bits, uint64_t tag_bits,uint64_t fixed_counter_size ,const char * path);

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


	virtual void reset()=0;



	virtual void copy(onDiskMQF*dest)=0;

	/*!
		@breif Increment the counter for this item by count.

		@param Qf* qf : pointer to the Filter
		@param uint64_t key : hash of the item to be insertedItems
		@param uint64_t count: Count to be added
		@param bool lock: For Multithreading, Lock the slot used by the current thread so that other threads can't change the value
		@param bool spin: For Multithreading, If there is a lock on the target slot. wait until the lock is freed and insert the count.

		@return bool: True if the item is inserted correctly.
	 */
	virtual bool insert( uint64_t key, uint64_t count,
								 bool lock=false, bool spin=false)=0;


	/* Remove all instances of this key/value pair. */
	//void onDiskMQF_delete_key_value(onDiskMQF*qf, uint64_t key, uint64_t value);

	/* Remove all instances of this key. */
	//void onDiskMQF_delete_key(onDiskMQF*qf, uint64_t key);

	/* Replace the association (key, oldvalue, count) with the association
		 (key, newvalue, count). If there is already an association (key,
		 newvalue, count'), then the two associations will be merged and
		 their counters will be summed, resulting in association (key,
		 newvalue, count' + count). */
	//void replace( uint64_t key, uint64_t oldvalue, uint64_t newvalue);

	/* Lookup the value associated with key.  Returns the count of that
		 key/value pair in the QF.  If it returns 0, then, the key is not
		 present in the QF. Only returns the first value associated with key
		 in the QF.  If you want to see others, use an iterator. */
	//uint64_t query( uint64_t key, uint64_t *value);

	/*!
	@breif Return the number of times key has been inserted, with any value,
		 into qf.

	@param Qf* qf : pointer to the Filter.
	@param uint64_t key : hash of the item.

	@return uint64_t the count associated with the input key.
			*/
	virtual uint64_t count_key( uint64_t key)=0;

	/*!
	@breif Decrement the counter for this item by count.

	@param	Qf* qf : pointer to the Filter
	@param	uint64_t key : hash of the item to be removed
	@param	uint64_t count: Count to be removed
	@param	bool lock: For Multithreading, Lock the slot used by the current thread so that other threads can't change the value
	@param	bool spin: For Multithreading, If there is a lock on the target slot. wait until the lock is freed and insert the count.

	@return bool: Returns true if the item is removed successfully.
	 */
	virtual bool remove( uint64_t hash, uint64_t count,  bool lock=false, bool spin=false)=0;


	/*!
		@breif Add Tag to item.

		@param Qf* qf : pointer to the Filter
		@param uint64_t key : hash of the item to be insertedItems
		@param uint64_t tag: tag to be added
		@param bool lock: For Multithreading, Lock the slot used by the current thread so that other threads can't change the value
		@param bool spin: For Multithreading, If there is a lock on the target slot. wait until the lock is freed and insert the count.

		@return bool: True if the item is inserted correctly.
	 */
	virtual uint64_t add_tag( uint64_t key, uint64_t tag, bool lock=false, bool spin=false)=0;
	/*!
	@breif Return the tag associated with a given item.

	@param Qf* qf : pointer to the Filter.
	@param uint64_t key : hash of the item.

	@return uint64_t the tag associated with the input key.
			*/
	virtual uint64_t get_tag(uint64_t key)=0;
	/*!
	@breif delete the tag associated with a given item.

	@param Qf* qf : pointer to the Filter.
	@param uint64_t key : hash of the item.

	@return bool: Returns true if the item is removed successfully.
			*/
	virtual uint64_t remove_tag(uint64_t key, bool lock=false, bool spin=false)=0;

	/* Initialize an iterator */
	//virtual bool getIterator(onDiskMQFIterator<bitsPerSlot> *qfi, uint64_t position)=0;



	/* For debugging */
	virtual void dump()=0;
	virtual void dump_block(uint64_t i)=0;
	//
	// /*! write data structure of to the disk */
	virtual void serialize(const char *filename)=0;
	//
	// /* read data structure off the disk */
	virtual void deserialize(const char *filename)=0;
	//
	// /* mmap the QF from disk. */
	virtual void read( const char *path)=0;
	//
	// /* merge two QFs into the third one. */
	// static void merge(onDiskMQF*qfa, onDiskMQF*qfb, onDiskMQF*qfc);
	//
	// void onDiskMQF_intersect(onDiskMQF *qfa, onDiskMQF *qfb, onDiskMQF *qfc);
	//
	// void onDiskMQF_subtract(onDiskMQF *qfa, onDiskMQF *qfb, onDiskMQF*qfc);
	// /* merge multiple QFs into the final QF one. */
	// void onDiskMQF_multi_merge(onDiskMQF *qf_arr[], int nqf, onDiskMQF *qfr);
	//

	/*! @breif Invertiable merge function adds tag for each key and creates index structure. The index is map of an integer and vector of integers where the integer is the value of the tags and vector on integers is the ids of the source filters.

	@param Qf* qf_arr : input array of filters
	@param int nqf: number of filters
	@param QF* qfr: pointer to the output filter.
	@param std::map<uint64_t, std::vector<int> > *inverted_index_ptr: Pointer to the output index.
	*/
	// void onDiskMQF_invertable_merge(onDiskMQF *qf_arr[], int nqf, onDiskMQF *qfr);
	// void onDiskMQF_invertable_merge_no_count(onDiskMQF *qf_arr[], int nqf, onDiskMQF *qfr);
	//

	/*! @breif Resize the filter into a bigger or smaller one

	@param Qf* qf : pointer to the Filter
	@param uint64_t newQ: new number of slots(Q). the slot size will be recalculated to keep the range constant.
	@param string originalFilename(optional): dump the current filter to the disk to free space for the new filter. Filename is provided as the content of the string.
	@param string newFilename(optional): the new filter is created on disk. Filename is provided as the content of the string.

	@return QF: New Quotient Filter.
	*/
	//onDiskMQF* onDiskMQF_resize(onDiskMQF* qf, int newQ, const char * originalFilename=NULL, const char * newFilename=NULL);
	/* find cosine similarity between two QFs. */
	//uint64_t onDiskMQF_inner_product(onDiskMQF *qfa, onDiskMQF *qfb);

	/* magnitude of a QF. */
	//uint64_t onDiskMQF_magnitude(onDiskMQF *qf);
	/* return the filled space(percent) */
	virtual int space()=0;

	//bool onDiskMQF_equals(onDiskMQF *qfa, onDiskMQF *qfb);

	virtual bool general_lock(bool spin)=0;
	virtual void general_unlock()=0;

	//void onDiskMQF_migrate(onDiskMQF* source, onDiskMQF* destination);
	virtual void migrateFromQF(QF* source)=0;
	};

	template<uint64_t bitsPerSlot=64>
	class _onDiskMQF: public onDiskMQF {
	private:
		bool spin_lock(uint64_t hash_bucket_index, bool spin, bool flag);
		void unlock(uint64_t hash_bucket_index, bool flag);
		void modify_metadata(uint64_t *metadata, int cnt);

		int is_runend(uint64_t index);
		int is_occupied(uint64_t index);
		uint64_t _get_slot(uint64_t index);
		void _set_slot(uint64_t index, uint64_t value);
		void set_slot(uint64_t index, uint64_t value);
		uint64_t get_slot(uint64_t index);
		uint64_t get_fixed_counter(uint64_t index);
		void set_fixed_counter(uint64_t index,uint64_t value);
		uint64_t _get_tag(uint64_t index);
		void set_tag(uint64_t index,uint64_t value);
		void super_get(uint64_t index,uint64_t* slot,uint64_t *fcounter);
		void super_set(uint64_t index,uint64_t slot,uint64_t fcounter);
		uint64_t run_end(uint64_t hash_bucket_index);
		uint64_t block_offset(uint64_t blockidx);
		int offset_lower_bound(uint64_t slot_index);
		int is_empty2(uint64_t slot_index);
		int might_be_empty(uint64_t slot_index);
		int probably_is_empty(uint64_t slot_index);
		uint64_t find_first_empty_slot(uint64_t from);
		void shift_remainders(const uint64_t start_index, const uint64_t empty_index);
		void find_next_n_empty_slots(uint64_t from, uint64_t n,uint64_t *indices);
		void shift_slots(int64_t first, uint64_t last, uint64_t distance);
		void shift_runends(int64_t first, uint64_t last, uint64_t distance);
		void insert_replace_slots_and_shift_remainders_and_runends_and_offsets(
							int operation,uint64_t bucket_index,uint64_t overwrite_index,
							const uint64_t *remainders,const uint64_t	*fixed_size_counters,
							uint64_t total_remainders,uint64_t noverwrites);
		void remove_replace_slots_and_shift_remainders_and_runends_and_offsets(
							int operation,uint64_t bucket_index,uint64_t overwrite_index,
							const uint64_t *remainders,const uint64_t	*fixed_size_counters,
							uint64_t total_remainders,uint64_t noverwrites);
		uint64_t *encode_counter(uint64_t remainder, uint64_t counter, uint64_t *slots, uint64_t *fixed_size_counters);
		uint64_t decode_counter(uint64_t index, uint64_t *remainder, uint64_t *count);
		bool insert1(__uint128_t hash, bool lock, bool spin);
		bool _insert(__uint128_t hash, uint64_t count, bool lock, bool spin);

	public:
		stxxl::vector<onDisk_qfblock<bitsPerSlot> > blocks;
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

	_onDiskMQF(uint64_t nslots, uint64_t key_bits, uint64_t tag_bits,uint64_t fixed_counter_size,const char * path);
	void reset() override;

	~_onDiskMQF();
	inline typename stxxl::vector<onDisk_qfblock<bitsPerSlot> >::iterator  get_block(uint64_t block_index)
	{
			typename stxxl::vector<onDisk_qfblock<bitsPerSlot> >::iterator it = this->blocks.begin();
			it+=block_index;
			return it;
	}
	inline typename stxxl::vector<onDisk_qfblock<bitsPerSlot> >::const_iterator get_block_const( uint64_t block_index)
	{
			typename stxxl::vector<onDisk_qfblock<bitsPerSlot> >::const_iterator it = this->blocks.cbegin();
			it+=block_index;
		  return it;
	}

	void copy(onDiskMQF *dest)override;

	/*!
	 	@breif Increment the counter for this item by count.

		@param Qf* qf : pointer to the Filter
		@param uint64_t key : hash of the item to be insertedItems
		@param uint64_t count: Count to be added
		@param bool lock: For Multithreading, Lock the slot used by the current thread so that other threads can't change the value
		@param bool spin: For Multithreading, If there is a lock on the target slot. wait until the lock is freed and insert the count.

		@return bool: True if the item is inserted correctly.
	 */
	bool insert( uint64_t key, uint64_t count,
								 bool lock=false, bool spin=false)override;


	/* Remove all instances of this key/value pair. */
	//void onDiskMQF_delete_key_value(onDiskMQF*qf, uint64_t key, uint64_t value);

	/* Remove all instances of this key. */
	//void onDiskMQF_delete_key(onDiskMQF*qf, uint64_t key);

	/* Replace the association (key, oldvalue, count) with the association
		 (key, newvalue, count). If there is already an association (key,
		 newvalue, count'), then the two associations will be merged and
		 their counters will be summed, resulting in association (key,
		 newvalue, count' + count). */
	//void replace( uint64_t key, uint64_t oldvalue, uint64_t newvalue);

	/* Lookup the value associated with key.  Returns the count of that
		 key/value pair in the QF.  If it returns 0, then, the key is not
		 present in the QF. Only returns the first value associated with key
		 in the QF.  If you want to see others, use an iterator. */
	//uint64_t query( uint64_t key, uint64_t *value);

	/*!
	@breif Return the number of times key has been inserted, with any value,
		 into qf.

	@param Qf* qf : pointer to the Filter.
	@param uint64_t key : hash of the item.

	@return uint64_t the count associated with the input key.
		  */
	uint64_t count_key( uint64_t key)override;

	/*!
	@breif Decrement the counter for this item by count.

	@param	Qf* qf : pointer to the Filter
	@param	uint64_t key : hash of the item to be removed
	@param	uint64_t count: Count to be removed
	@param	bool lock: For Multithreading, Lock the slot used by the current thread so that other threads can't change the value
	@param	bool spin: For Multithreading, If there is a lock on the target slot. wait until the lock is freed and insert the count.

	@return bool: Returns true if the item is removed successfully.
	 */
	bool remove( uint64_t hash, uint64_t count,  bool lock=false, bool spin=false)override;


	/*!
		@breif Add Tag to item.

		@param Qf* qf : pointer to the Filter
		@param uint64_t key : hash of the item to be insertedItems
		@param uint64_t tag: tag to be added
		@param bool lock: For Multithreading, Lock the slot used by the current thread so that other threads can't change the value
		@param bool spin: For Multithreading, If there is a lock on the target slot. wait until the lock is freed and insert the count.

		@return bool: True if the item is inserted correctly.
	 */
	uint64_t add_tag( uint64_t key, uint64_t tag, bool lock=false, bool spin=false)override;
	/*!
	@breif Return the tag associated with a given item.

	@param Qf* qf : pointer to the Filter.
	@param uint64_t key : hash of the item.

	@return uint64_t the tag associated with the input key.
			*/
	uint64_t get_tag(uint64_t key)override;
	/*!
	@breif delete the tag associated with a given item.

	@param Qf* qf : pointer to the Filter.
	@param uint64_t key : hash of the item.

	@return bool: Returns true if the item is removed successfully.
			*/
	uint64_t remove_tag(uint64_t key, bool lock=false, bool spin=false)override;

	/* Initialize an iterator */
	bool getIterator(onDiskMQFIterator<bitsPerSlot> *qfi, uint64_t position);



	/* For debugging */
	void dump()override;
	void dump_block(uint64_t i)override;
  //
	// /*! write data structure of to the disk */
	 void serialize(const char *filename)override;
  //
	// /* read data structure off the disk */
	void deserialize(const char *filename)override;
  //
	// /* mmap the QF from disk. */
	 void read( const char *path)override;
  //
	// /* merge two QFs into the third one. */
	// static void merge(onDiskMQF*qfa, onDiskMQF*qfb, onDiskMQF*qfc);
  //
	// void onDiskMQF_intersect(onDiskMQF *qfa, onDiskMQF *qfb, onDiskMQF *qfc);
  //
	// void onDiskMQF_subtract(onDiskMQF *qfa, onDiskMQF *qfb, onDiskMQF*qfc);
	// /* merge multiple QFs into the final QF one. */
	// void onDiskMQF_multi_merge(onDiskMQF *qf_arr[], int nqf, onDiskMQF *qfr);
  //

	/*! @breif Invertiable merge function adds tag for each key and creates index structure. The index is map of an integer and vector of integers where the integer is the value of the tags and vector on integers is the ids of the source filters.

	@param Qf* qf_arr : input array of filters
	@param int nqf: number of filters
	@param QF* qfr: pointer to the output filter.
	@param std::map<uint64_t, std::vector<int> > *inverted_index_ptr: Pointer to the output index.
	*/
	// void onDiskMQF_invertable_merge(onDiskMQF *qf_arr[], int nqf, onDiskMQF *qfr);
	// void onDiskMQF_invertable_merge_no_count(onDiskMQF *qf_arr[], int nqf, onDiskMQF *qfr);
  //

	/*! @breif Resize the filter into a bigger or smaller one

	@param Qf* qf : pointer to the Filter
	@param uint64_t newQ: new number of slots(Q). the slot size will be recalculated to keep the range constant.
	@param string originalFilename(optional): dump the current filter to the disk to free space for the new filter. Filename is provided as the content of the string.
	@param string newFilename(optional): the new filter is created on disk. Filename is provided as the content of the string.

	@return QF: New Quotient Filter.
	*/
	//onDiskMQF* onDiskMQF_resize(onDiskMQF* qf, int newQ, const char * originalFilename=NULL, const char * newFilename=NULL);
	/* find cosine similarity between two QFs. */
	//uint64_t onDiskMQF_inner_product(onDiskMQF *qfa, onDiskMQF *qfb);

	/* magnitude of a QF. */
	//uint64_t onDiskMQF_magnitude(onDiskMQF *qf);
	/* return the filled space(percent) */
	int space()override;

	//bool onDiskMQF_equals(onDiskMQF *qfa, onDiskMQF *qfb);

	bool general_lock(bool spin)override;
	void general_unlock()override;

	//void onDiskMQF_migrate(onDiskMQF* source, onDiskMQF* destination);
	void migrateFromQF(QF* source) override;
};

}
#endif /* onDiskMQF_H */
