#include "async_sequence_processor.hh"

using namespace khmer;
using namespace khmer::read_parsers;

/////
//
// AsyncSequenceProcessor
//
/////

void AsyncSequenceProcessor::write(const char * sequence) {
    khmer::KMerIterator kmers(sequence, _ksize);
    HashIntoType kmer;
    try {
        while(!kmers.done()) {
            kmer = kmers.next();
            _ht->count_ts(kmer);
        }
    } catch (khmer_exception &e) {
        _exc_handler->push(std::make_exception_ptr(e));
        std::cout << e.what() << " " << std::endl;
        std::cout << "ERROR in AsyncCounter: " << sequence << std::endl;
        return;
    }
}

void AsyncSequenceProcessor::start(const std::string &filename,
                                    bool paired,
                                    unsigned int n_threads) {
    _ksize = _ht->ksize();
    aparser = new AsyncParser();
    aparser->start(filename, paired);
    _n_processed = 0;
    AsyncConsumerProducer<ReadBatchPtr,ReadBatchPtr>::start(n_threads);
}

void AsyncSequenceProcessor::stop() {
    AsyncConsumerProducer<ReadBatchPtr,ReadBatchPtr>::stop();
}

unsigned int AsyncSequenceProcessor::n_processed() {
    return _n_processed;
}


/////
//
// AsyncSequenceProcessorTester
//
/////

bool AsyncSequenceProcessorTester::iter_stop() {
    if(!workers_running() && (n_popped() >= n_processed()))
        return true;
    return false;
}

void AsyncSequenceProcessorTester::consume() {
    ReadBatchPtr batch;
    const char * sp;
    while(_workers_running) {
        if(_in_queue->pop(batch)) {

            sp = copy_seq(batch->first());
            write(sp);            
            if (paired) {
                sp = copy_seq(batch->second());
                write(sp);
            }
            while(!(_out_queue->bounded_push(batch))) if (!_workers_running) return;
            __sync_fetch_and_add(&_n_processed, _batchsize); 
        } else {
            if (!is_parsing() && (_n_processed >= _n_parsed)) {
                _workers_running = false;
            }
        }
    }
}


