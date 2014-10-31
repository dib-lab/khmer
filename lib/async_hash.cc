#include "hashtable.hh"
#include "async_hash.hh"
#include <boost/lockfree/queue.hpp>

namespace khmer {

void AsyncHashWriter::start(unsigned int n_hasher_threads) {
    n_hashers = n_hasher_threads;
    std::cout << "INIT AsyncHashWriter, n_writers=" << n_writers << " n_hashers=" << n_hashers << std::endl;

    //std::cout << "Thread " << std::this_thread::get_id() << " spawning hash consumer" << std::endl;
    writers_running = true;
    writer_thread = new std::thread(&AsyncHashWriter::consume_hash, this, std::ref(hash_q));

    hashers_running = true;
    for (unsigned int t=0; t<n_hashers; ++t) {
        std::cout << "AsyncHashWriter spawn consume_kmer" << std::endl;
        hasher_threads.push_back(std::thread(&AsyncHashWriter::consume_kmer, this, std::ref(kmer_q)));
        //hasher_threads.push_back(*hthread);
    }
}

void AsyncHashWriter::stop() {
    //hash_q.push(_stop_sig);
    writers_running = false;
    writer_thread->join();

    hashers_running = false;
    auto hthread = hasher_threads.begin();
    while (hthread != hasher_threads.end()) {
        hthread->join();
        hthread++;
    }
    //std::cout << "Thread " << std::this_thread::get_id() << " destroyed hash consumer" << std::endl;
}

void AsyncHashWriter::enqueue_string(const std::string &s) {
    char * sp = new char [s.length() + 1];
    std::strcpy (sp, s.c_str());
    for (unsigned int i=0; i<s.length()-_ht->_ksize+1; ++i) {
        enqueue_kmer(sp);
        ++sp;
    }
}

void AsyncHashWriter::enqueue_kmer(const char * kmer) {
    kmer_q.push(kmer);
}

void AsyncHashWriter::consume_kmer(boost::lockfree::queue<const char*>& q) {
    const char * kmer;
    unsigned long long total_consumed = 0;
    //std::cout << "Thread " << std::this_thread::get_id() << " waiting on hashes" << std::endl;
    while(1) {
        if (q.pop(kmer)) {
            enqueue_hash(_hash(kmer, _ht->_ksize));
            total_consumed++;
        } else {
            if (!hashers_running) return;
        }
        //if (total_consumed % 1000000 == 0) {
        //    std::cout << "Thread " << std::this_thread::get_id() << ": consumed " << 
        //        total_consumed << ", " << q.size()  << " in queue" << std::endl;
        //}
    }
}

void AsyncHashWriter::enqueue_hash(HashIntoType khash) {
    hash_q.push(khash);
}

void AsyncHashWriter::consume_hash(boost::lockfree::queue<HashIntoType>& q) {
    HashIntoType khash;
    unsigned long long total_consumed = 0;
    //std::cout << "Thread " << std::this_thread::get_id() << " waiting on hashes" << std::endl;
    while(1) {
        if (q.pop(khash)) {
            _ht->count(khash);
            total_consumed++;
        } else {
            if(!writers_running) return;
        }
        //if (total_consumed % 100000 == 0) {
        //    std::cout << "async consumed " << total_consumed << std::endl;
        //}
    }
}
}
