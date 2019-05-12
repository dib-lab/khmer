#include <inttypes.h>
#include <stdbool.h>
#include <pthread.h>
#include "gqf.h"
#include "LayeredMQF.h"
#include <iostream>

using namespace std;
void layeredMQF_init(layeredMQF *qf, uint64_t nslots_singletons ,uint64_t nslots
	, uint64_t key_bits, uint64_t value_bits,uint64_t fixed_counter_size, bool mem
	, const char *path, uint32_t seed){

		if(qf==NULL)
		{
			qf= new layeredMQF();
		}

		qf_init(qf->firstLayer_singletons,nslots_singletons,key_bits,value_bits,1,0,mem,path,seed);
		qf_init(qf->secondLayer,nslots,key_bits,value_bits,fixed_counter_size,0,mem,path,seed);
}

void layeredMQF_reset(layeredMQF *qf){
	qf_reset(qf->firstLayer_singletons);
	qf_reset(qf->secondLayer);
}

void layeredMQF_destroy(layeredMQF *qf){
	qf_destroy(qf->firstLayer_singletons);
	qf_destroy(qf->secondLayer);
}

void layeredMQF_copy(layeredMQF *dest, layeredMQF *src){
	qf_copy(dest->firstLayer_singletons,src->firstLayer_singletons);
	qf_copy(dest->secondLayer,src->secondLayer);
}

/* Increment the counter for this key/value pair by count. */

bool layeredMQF_insert(layeredMQF *qf, uint64_t key, uint64_t count,
							 bool lock, bool spin){
	if(count==0)
		return true;
	bool inSecondLayer=qf_count_key(qf->secondLayer,key)>0;
	if(inSecondLayer)
	{
		return qf_insert(qf->secondLayer,key,count,lock,spin);
	}
	uint64_t CountinFirstLayer=qf_count_key(qf->firstLayer_singletons,key);
	if(CountinFirstLayer>1)
	{
		cerr<<"First Layer has items > 1"<<endl;
	}
	if(CountinFirstLayer>0)
	{
		if(qf_remove(qf->firstLayer_singletons,key,CountinFirstLayer,lock,spin))
			return qf_insert(qf->secondLayer,key,count+1,lock,spin);
		else{
			return false;
		}
	}
	if(count==1)
		return qf_insert(qf->firstLayer_singletons,key,count,lock,spin);
	else
		return qf_insert(qf->secondLayer,key,count,lock,spin);


}

/* Remove count instances of this key/value combination. */
bool LayeredMQF_remove(layeredMQF *qf, uint64_t hash, uint64_t count,  bool lock, bool spin){
	bool res=false;
	res|=qf_remove(qf->firstLayer_singletons,hash,count,lock,spin);
	res|=qf_remove(qf->secondLayer,hash,count,lock,spin);
	return res;
}




/* Return the number of times key has been inserted, with any value,
	 into qf. */
uint64_t layeredMQF_count_key(const layeredMQF *qf, uint64_t key){
	uint64_t res=qf_count_key(qf->secondLayer,key);
	if(res==0)
	{
		res=qf_count_key(qf->firstLayer_singletons,key);
	}
	return res;
}

int layeredMQF_space(layeredMQF *qf){
	return max(qf_space(qf->firstLayer_singletons),qf_space(qf->secondLayer));
}
//
// /* Initialize an iterator */
// bool layeredMQF_iterator(layeredMQF *qf, layeredMQFIterator* qfi, uint64_t position){
// 	if(qfi ==NULL)
// 		qfi=new layeredMQFIterator();
// 	bool res=false;
// 	res|=qf_iterator(qf->firstLayer_singletons,qfi->firstLayerIterator,position);
// 	res|=qf_iterator(qf->secondLayer,qfi->secondLayerIterator,position);
// 	return res;
// }
//
// /* Returns 0 if the iterator is still valid (i.e. has not reached the
// 	 end of the QF. */
// int layeredMQF_qfi_get(layeredMQFIterator*qfi, uint64_t *key, uint64_t *value, uint64_t *count){
//
// }
//
// /* Advance to next entry.  Returns whether or not another entry is
// 	 found.  */
// int layeredMQF_qfi_next(layeredMQFIterator*qfi);
//
// /* Check to see if the if the end of the QF */
// int layeredMQF_qfi_end(layeredMQFIterator*qfi);
//
// /* For debugging */
// void layeredMQF_dump(const layeredMQF *);
//
// /* write data structure of to the disk */
// void layeredMQF_serialize(const layeredMQF *qf, const char *filename);
//
// /* read data structure off the disk */
// void layeredMQF_deserialize(layeredMQF *qf, const char *filename);
//
// /* mmap the QF from disk. */
// void layeredMQF_read(layeredMQF *qf, const char *path);

/* merge two QFs into the third one. */
//void layeredMQF_merge(layeredMQF *layeredMQFa, layeredMQF *layeredMQFb, QF *qfc);

/* merge multiple QFs into the final QF one. */
//void qf_multi_merge(QF *qf_arr[], int nqf, QF *qfr);

/* find cosine similarity between two QFs. */
//	uint64_t qf_inner_product(QF *qfa, QF *qfb);

/* magnitude of a QF. */
//uint64_t qf_magnitude(QF *qf);
