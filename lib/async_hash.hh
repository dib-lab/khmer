#ifndef ASYNC_HASH_HH
#define ASYNC_HASH_HH

#include "hashtable.hh"
#include "cqueue.hh"
#include <iostream>
#include <thread>

namespace khmer {

class CountingHash;
class Hashtable;
class Hashbits;

class AsyncHashWriter {
    
    friend class Hashtable; 
    
    khmer::Hashtable * _ht;
    uint32_t _n_writer_threads;
    CQueue<HashIntoType> hash_q;
    std::thread * consume_hash_thread;
    HashIntoType _stop_sig = -1;

    public:

        AsyncHashWriter (khmer::Hashtable * ht, uint32_t n_writer_threads):
                        _ht(ht),
                        _n_writer_threads(n_writer_threads) {
        }

        ~AsyncHashWriter() {
        }

        void start();
        void stop();
        void enqueue(HashIntoType khash);
        void consume_hash(CQueue<HashIntoType>& q);

        /*
        void write_table(HashIntoType khash, unsigned int table, HashIntoType table_size) {
            const HashIntoType bin = khash % table_size;
            if (_ht->_max_count > _ht->_counts[table][bin]) {
                _ht->_counts[table][bin] += 1;
            } else {

            }
        }
        */
};
}; // namspace khmer

#endif
