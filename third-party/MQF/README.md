# MQF
[![Build Status](https://travis-ci.org/shokrof/MQF.svg?branch=mqfDevelopmenet)](https://travis-ci.org/shokrof/MQF)
[![codecov](https://codecov.io/gh/shokrof/MQF/branch/mqfDevelopmenet/graph/badge.svg)](https://codecov.io/gh/shokrof/MQF)

[Approximate membership query (AMQ)](http://www.cs.cmu.edu/~lblum/flac/Presentations/Szabo-Wexler_ApproximateSetMembership.pdf) data structures provide approximate representation for data using a smaller amount of memory compared to the real data size. As the name suggests, AMQ answers if a particular element exists or not in a given dataset but with possible false positive errors. AMQ has many examples such as [bloom filter](https://en.wikipedia.org/wiki/Bloom_filter), [count-min sketch](https://en.wikipedia.org/wiki/Count%E2%80%93min_sketch), [Quotient Filter](https://en.wikipedia.org/wiki/Quotient_filter), and Counting Quotient Filter([CQF](https://github.com/splatlab/cqf)). Here, we are proposing a new AMQ data structure called Mixed Counting Quotient Filter (MQF).

CQF splits the hash-bits of an item into two components: quotient and remaining parts. Quotient Part is used to determine the target slot. The remaining is inserted into the target slot. The insertion algorithm uses a variant of linear probing to resolve collisions. CQF allows counting the number of instances inserted by using slots following the item's reminder as a counter. CQF has relatively complex scheme to encode counters so that it can distinguish counters from item's reminder.

![alt text](https://raw.githubusercontent.com/shokrof/MQF/mqfDevelopmenet/QuotientFilter_MQF.png)

MQF, Mixed Quotient Filter, is a variant of CQF. MQF uses the same probing technique as CQF. MQF has more metadata called fixed size counters and different encoding for the counters. The improvement makes mqf more memory efficient for wider range of zipifan distribution.

When an item is inserted more than one time, MQF first use the fixed size counter to the number of insertions(count). After the fixed size counter is full, MQF store the count in the following slot and fixed-size counter. MQF uses the necessary number of slots to store the count.

In other words, Fixed-size counters is used in counting and marking the slots used for counting. Fixed-size counters for all slots related to the same item store the counter's maximum value except the last one should store value strictly less than the maximum. When the maximum is reached in the last fixed counter, a new slot is added with empty fixed-size counter.




## Advantages:
  - MQF supports counting and removing items. MQF uses variable size counters; therefore, It is memory efficient when count data that follows zipfian distribution where most the items occur once or twice but few items can happen in very high counts..
  - MQF has lower bits per element than Bloom filter and Count-min sketches ([Ref](https://www3.cs.stonybrook.edu/~ppandey/files/p775-pandey.pdf)).
  - MQF  has good data locality which makes it efficient when running on secondary storage.
  - MQF supports add/remove tags for the item.
  - MQF can be iterated to get the items and counts inserted in the filter.
  - MQF supports merging function. Two or more filters can be merged into one filter.
  - MQF can be resized to bigger/ smaller filter.


## Documentation
### Building
MQF only requires make and g++ to be installed.
```bash
apt-get install make g++
make NH=1
make test NH=1
./mqf_test
```
### Initialization
MQF Initialization requires the estimation of some parameters: number of slots, Key bits, fixed counter size, and tag size.

Fixed-size counter size is estimated from the shape of data distribution. If most of the items are singletons. The fixed-size counter should be limited to 1 bit. However, If a big portion of items is repeated more than one time, a bigger fixed-size counter will hold more counts and thus save slots.

The number of slots is estimated from the number of items expected to be inserted into the filter. Slots are used for inserting the remaining part of the hash of the items and the count. After calculating the required number of slots, multiply the result by 1.05 because MQF can't be filled by more than 95% of its capacity. Then, round the result to the nearest bigger power of two.

Key bits equal to log2(number of slots) + the remaining part bits. the remaining part bits is estimated from the desired accuracy using the formula below.

![eqn](https://raw.githubusercontent.com/shokrof/MQF/mqfDevelopmenet/r_eqn.gif)

Tag size is straightforward to be estimated. it can be set to zero if tags are not necessary.


1. qf_init
Initialize mqf .
```c++
void qf_init(QF *qf, uint64_t nslots, uint64_t key_bits, uint64_t tag_bits,uint64_t fixed_counter_size, bool mem, const char *path, uint32_t seed);
```

  * Qf* qf : pointer to the Filter.
  * uint64_t nslots : Number of slots in the filter. Should be of power of two. Maximum number of items to be inserted depends on this number.
  * uint64_t key_bits: Number of bits in the hash values.
  * uint64_t tag_bits: Number of bits in tag value.
  * uint64_t fixed_counter_size: Fixed counter size. must be > 0.
  * bool mem: Flag to create the filter on memory. IF false, mmap is used.
  * const char * path: In case of mmap. Path of the file used to pack the filter.
  * uint32_t seed: useless value. To be removed
2. qf_destroy
3. estimate

### Functions Supported
1. Insert :
Increment the counter for this item by count.
  ```c++
  bool qf_insert(QF *qf, uint64_t key,
     uint64_t count,bool lock, bool spin);
  ```

  * Qf* qf : pointer to the Filter
  * uint64_t key : hash of the item to be inserted.
  * uint64_t count: Count to be added
  * bool lock: For Multithreading, Lock the * slot used by the current thread so that other threads can't change the value
  * bool spin: For Multithreading, If there is a lock on the target slot. wait until the lock is freed and insert the count.
  * returns True if the insertion succeeded.

2. Count:
 Return the number of times key has been inserted, with any value, into qf.
 ```c++
 uint64_t qf_count_key(const QF *qf, uint64_t key);
 ```
 * Qf* qf : pointer to the Filter
 * uint64_t key : hash of the item to be counted.
 * returns the number of times the item is inserted.
3. Remove:
Decrement the counter for this item by count.
```c++
bool qf_remove(QF *qf, uint64_t hash, uint64_t count, bool lock, bool spin);
```
  * Qf* qf : pointer to the Filter
  * uint64_t key : hash of the item to be removed
  * uint64_t count: Count to be removed
  * bool lock: For Multithreading, Lock the slot used by the current thread so that other threads can't change the value
  * bool spin: For Multithreading, If there is a lock on the target slot. wait until the lock is freed and insert the count.

4. Add/Remove tag to elements
```c++
uint64_t qf_add_tag(const QF *qf, uint64_t key, uint64_t tag, bool lock, bool spin);
uint64_t qf_get_tag(const QF *qf, uint64_t key);
uint64_t qf_remove_tag(const QF *qf, uint64_t key, bool lock, bool spin);
```
  * Qf* qf : pointer to the Filter
  * uint64_t key : hash of the item.
  * uint64_t tag: tag for the item.
  * bool lock: For Multithreading, Lock the slot used by the current thread so that other threads can't change the value
  * bool spin: For Multithreading, If there is a lock on the target slot. wait until the lock is freed and insert the count.

5. Resize:
 resize the filter into a bigger or smaller one
 ```c++
 QF* qf_resize(QF* qf, int newQ, const char * originalFilename=NULL, const char * newFilename=NULL);
 ```
 * Qf* qf : pointer to the Filter
 * uint64_t newQ: new number of slots(Q). the slot size will be recalculated to keep the range constant.
 * string originalFilename(optional): dump the current filter to the disk to free space for the new filter. Filename is provided as the content of the string.
 * string newFilename(optional): the new filter is created on disk. Filename is provided as the content of the string.
 * returns a pointer to the new filter

6. Merge: merge more than one filter into a final one.
```c++
void qf_merge(QF *qfa, QF *qfb, QF *qfc);
void qf_multi_merge(QF *qf_arr[], int nqf, QF *qfr);
```
7. Invertible Merge: Invertible merge offers addiotinal functionality to normal merge. Original source filter can be queried for each key.
Invertiable merge function adds tag for each key and creates index structure. The index is map of an integer and vector of integers where the integer is the value of the tags and vector on integers is the ids of the source filters.
```c++
void qf_invertable_merge(QF *qf_arr[], int nqf, QF *qfr,std::map<uint64_t, std::vector<int> > *inverted_index_ptr);
```
  * Qf* qf_arr : input array of filters
  * int nqf: number of filters
  * QF* qfr: pointer to the output filter.
  * map (uint64_t,vector(int) )    inverted_index_ptr: Pointer to the output index.




7. Compare:
check if two filters have the same items, counts and tags.
```c++
bool qf_equals(QF *qfa, QF *qfb);
```
8. Intersect
calculate the intersection between two filters.
```c++
void qf_intersect(QF *qfa, QF *qfb, QF *qfc);
```
9. Subtract
subtract the second filter from the first.
```c++
void qf_subtract(QF *qfa, QF *qfb, QF *qfc);
```
10. Space:
returns the space percent occupied by the inserted items.
```c++
int qf_space(QF *qf);
```

### Miscellaneous Functions
1. Capacity
2. Copy
3. Serialize/ Deserialize
4. MMap read
