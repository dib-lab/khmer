#include "hashtable.hh"
#include "async_hash.hh"
#include <boost/lockfree/queue.hpp>

namespace khmer {

/////
//
// class Async methods
//
/////

void Async::start(unsigned int n_threads) {
    _n_hashers = n_threads;

    _workers_running = true;
    for (unsigned int t=0; t<_workers; ++t) {
        std::cout << "Async spawn worker" << std::endl;
        _worker_threads.push_back(std::thread(&Async::consume, this, std::ref(_in_queue)));
    }
}

void Async::stop() {

    _workers_running = false;
    auto ethread = _worker_threads.begin();
    while (ethread != _worker_threads.end()) {
        wthread->join();
        wthread++;
    }
}

void Async::push(T item) {
    _in_queue.push(item);
}

void Async::set_input(Queue<T>& new_q) {
    _in_queue = new_q;
}

/////
//
// class AsyncWriter methods
//
/////

void AsyncWriter::push(HashIntoType item) {
    _in_queue.push(item);
}

void AsyncWriter::consume(Queue<HashIntoType>& q) {
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
    }
}

/////
//
// class AsyncHasher methods
//
/////

bool AsyncHasher::pop(HashIntoType& khash) {
    return _out_queue.pop(khash);
}

void AsyncHasher::set_output(HashQueue& new_q) {
    _out_queue = new_q;
}

// Warning: this method will always consume the entire
// character array as decomposed k-mers! So, if you want to consume
// an individual k-mer, you have to push appropriately sized arrays
// onto the queue. You shouldn't do this though, because more queue
// accesses is more overhead.

void AsyncHashWriter::consume(CharQueue& q) {
    const char * kmer;
    khmer::KMerIterator kmers;
    unsigned long long total_consumed = 0;
    //std::cout << "Thread " << std::this_thread::get_id() << " waiting on hashes" << std::endl;
    while(1) {
        if (q.pop(kmer)) {
            kmers(kmer, _ksize);
            while(!kmers.done()) {
                _out_queue.push(kmers.next());
                total_consumed++;
            }
        } else {
            if (!hashers_running) return;
        }
    }
}


/////
//
// AsyncSequenceProcessor
//
/////

void AsyncSequenceProcessor::start(unsigned int n_processor_threads) {
    _n_processor_threads = n_processor_threads;
    std::cout << "INIT AsyncHashWriter, n_writers=" << n_writers << " n_hashers=" << n_hashers << std::endl;

    _hasher = new khmer::AsyncHashWriter(_ht);
    

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


};
