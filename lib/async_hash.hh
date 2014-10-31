#ifndef ASYNC_HASH_HH
#define ASYNC_HASH_HH

#include "hashtable.hh"
//#include "cqueue.hh"
#include <iostream>
#include <thread>
#include <boost/lockfree/queue.hpp>

namespace khmer {

class CountingHash;
class Hashtable;
class Hashbits;

class AsyncHashWriter {
    
    friend class Hashtable; 
    
    khmer::Hashtable * _ht;

    //CQueue<HashIntoType> hash_q;
    //CQueue<const char*> kmer_q;
    boost::lockfree::queue<HashIntoType> hash_q;
    boost::lockfree::queue<const char *> kmer_q;

    std::thread * writer_thread;
    std::vector<std::thread> hasher_threads;
    //HashIntoType _stop_sig = -1;

    bool writers_running, hashers_running;

    unsigned int n_writers;
    unsigned int n_hashers;

    public:

        AsyncHashWriter (khmer::Hashtable * ht):
                        _ht(ht) {
            writers_running = false;
            hashers_running = false;
            n_writers = 1;
        }

        ~AsyncHashWriter() {
            if (writers_running || hashers_running) stop();
        }

        void start(unsigned int n_hasher_threads);
        void stop();
        void enqueue_kmer(const char * kmer);
        void consume_kmer(boost::lockfree::queue<const char *>& q);
        void enqueue_string(const std::string &s);
        void enqueue_hash(HashIntoType khash);
        void consume_hash(boost::lockfree::queue<HashIntoType>& q);

};
}; // namspace khmer

#endif
