# MQF
[![Build Status](https://travis-ci.org/shokrof/MQF.svg?branch=mqfDevelopmenet)](https://travis-ci.org/shokrof/MQF)
[![codecov](https://codecov.io/gh/shokrof/MQF/branch/mqfDevelopmenet/graph/badge.svg)](https://codecov.io/gh/shokrof/MQF)

MQF, Mixed Quotient Filter, is approximate membership query data structure that supports many useful functions. MQF is a variant of [CQF](https://github.com/splatlab/cqf). MQF has lower bits per element than Bloom filter and Count-min sketches. MQF also has good data locality which makes it efficient when running on secondary storage. Moreover. It supports removing, iterating, merging ,and resizing.

## Documentation
### Building
MQF onlye requires make and g++ to be installed.
```bash
apt-get install make g++
make NH=1
make test NH=1
./mqf_test
```
### Initialization
1. qf_init
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

7. Compare:
check if two filters have the same items, counts and tags.
```c++
bool qf_equals(QF *qfa, QF *qfb);
```
8. Intersect
calculate the the intersection between two filters.
```c++
void qf_intersect(QF *qfa, QF *qfb, QF *qfc);
```
9. Subtract
subtract the second filter from the first.
```c++
void qf_subtract(QF *qfa, QF *qfb, QF *qfc);
```
10. Space:
returns the space  percent occupied by the inserted items.
```c++
int qf_space(QF *qf);
```

### Miscellaneous Functions
1. Capacity
2. Copy
3. Serialize/ Deserialize
4. MMap read
