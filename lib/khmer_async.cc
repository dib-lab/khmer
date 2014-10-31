#include "hashtable.hh"
#include "khmer_async.hh"
#include <boost/lockfree/queue.hpp>

using namespace khmer;
using namespace khmer::read_parsers;
using namespace boost::lockfree;

/////
//
// class AsyncWriter methods
//
/////

void AsyncWriter::start() {
    Async<HashIntoType>::start(1);
}

void AsyncWriter::consume(HashQueue * q) {
    HashIntoType khash;
    unsigned long long total_consumed = 0;
    //std::cout << "Thread " << std::this_thread::get_id() << " waiting on hashes" << std::endl;
    while(1) {
        if (q->pop(khash)) {
            _ht->count(khash);
            total_consumed++;
        } else {
            if(!_workers_running) return;
        }
    }
}

/////
//
// class AsyncHasher methods
//
/////

bool AsyncHasher::pop(HashIntoType& khash) {
    return _out_queue->pop(khash);
}

void AsyncHasher::set_output(HashQueue* new_q) {
    _out_queue = new_q;
}

// Warning: this method will always consume the entire
// character array as decomposed k-mers! So, if you want to consume
// an individual k-mer, you have to push appropriately sized arrays
// onto the queue. You shouldn't do this though, because more queue
// accesses is more overhead.

void AsyncHasher::consume(CharQueue * q) {
    const char * kmer;
    unsigned long long total_consumed = 0;
    //std::cout << "Thread " << std::this_thread::get_id() << " waiting on hashes" << std::endl;
    while(1) {
        if (q->pop(kmer)) {
            khmer::KMerIterator kmers(kmer, _ksize);
            while(!kmers.done()) {
                _out_queue->push(kmers.next());
                total_consumed++;
            }
        } else {
            if (!_workers_running) return;
        }
    }
}


/////
//
// AsyncSequenceProcessor
//
/////

void AsyncSequenceProcessor::start(unsigned int n_threads) {
    _writer->start();
    Async<Read*>::start(n_threads);
}

void AsyncSequenceProcessor::stop() {
    // First stop the sequence processor, flushing its queue
    Async<Read*>::stop();
    // Then stop the writer, writing out all the flushed hashes
    _writer->stop();
}

bool AsyncSequenceProcessor::pop(Read* read) {
    return _out_queue->pop(read);
}

void AsyncSequenceProcessor::set_output(ReadQueue* new_q) {
    _out_queue = new_q;
}

/////
//
// AsyncDiginorm
//
/////

void AsyncDiginorm::start(unsigned int cutoff, unsigned int n_threads) {
    _cutoff = cutoff;
    AsyncSequenceProcessor::start(n_threads);
}

void AsyncDiginorm::consume(ReadQueue * q) {
    Read* read;
    BoundedCounterType median;
    float average, stddev;
    //std::cout << "Thread " << std::this_thread::get_id() << " waiting on hashes" << std::endl;
    while(1) {
        if (q->pop(read)) {
            _ht->get_median_count(read->sequence, median, average, stddev);
            if (median < _cutoff)
                _out_queue->push(read);
        } else {
            if (!_workers_running) return;
        }
    }
}
