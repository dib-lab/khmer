#include "async_writers.hh"

using namespace khmer;

/////
//
// class AsyncHashWriter methods
//
/////

void AsyncHashWriter::start() {
    _n_written = 0;
    AsyncConsumer<HashIntoType>::start(1);
}

unsigned int AsyncHashWriter::ksize() {
    return _ht->ksize();
}

void AsyncHashWriter::consume() {
    HashIntoType khash;
    //std::cout << "Thread " << std::this_thread::get_id() << " waiting on hashes" << std::endl;
    while(1) {
        if (_in_queue->pop(khash)) {
            _ht->count(khash);
            _n_written++;
        } else {
            if(!_workers_running && (_n_written >= _n_pushed)) return;
        }
    }
}

unsigned int AsyncHashWriter::n_written() {
    return _n_written;
}

/////
//
// AsyncSequenceWriter methods
//
/////


void AsyncSequenceWriter::start() {
    _n_written = 0;
    _ksize = _ht->ksize();
    AsyncConsumer<const char *>::start(1);
}

unsigned int AsyncSequenceWriter::ksize() {
    return _ht->ksize();
}

void AsyncSequenceWriter::consume() {
    const char * sequence;
    HashIntoType kmer;
    CountFunc kmer_counter(_ht, _ksize);
    //std::cout << "Thread " << std::this_thread::get_id() << " waiting on hashes" << std::endl;
    while(_workers_running) {
        if (_in_queue->pop(sequence)) {
            khmer::KMerIterator kmers(sequence, _ksize);
            try {
                while(!kmers.done()) {
                    kmer = kmers.next();
                    _ht->count(kmer);
                }
            } catch (khmer_exception &e) {
                _exc_handler->push(std::make_exception_ptr(e));
                std::cout << e.what() << " " << std::endl;
                std::cout << "ERROR in AsyncSequenceWriter: " << sequence << std::endl;
                exit(1);
                //flush_queue<const char *>(_in_queue);
                return;
            }
            _n_written++;
        }
    }
    #if(VERBOSITY)
    std::cout << "Exit AsyncSequenceWriter consume thread." << std::endl;
    #endif
    _in_queue->consume_all(kmer_counter);
}

unsigned int AsyncSequenceWriter::n_written() {
    return _n_written;
}

/////
//
// class AsyncHasher methods
//
/////


// Warning: this method will always consume the entire
// character array as decomposed k-mers! So, if you want to consume
// an individual k-mer, you have to push appropriately sized arrays
// onto the queue. You shouldn't do this though, because more queue
// accesses is more overhead.

void AsyncHasher::consume() {
    const char * kmer;
    unsigned long long total_consumed = 0;
    //std::cout << "Thread " << std::this_thread::get_id() << " waiting on hashes" << std::endl;
    while(_workers_running) {
        if (_in_queue->pop(kmer)) {
            khmer::KMerIterator kmers(kmer, _ksize);
            while(!kmers.done()) {
                while(!(_out_queue->bounded_push(kmers.next()))) if (!_workers_running) return;
                total_consumed++;
            }
        }
    }
}

/////
//
// class AsyncSequenceParser
//
/////


