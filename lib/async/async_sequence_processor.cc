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
    aparser = new AsyncSequenceParser();
    aparser->start(filename, paired);
    _paired = paired;
    _n_processed = 0;
    AsyncConsumerProducer<ReadBatchPtr,ReadBatchPtr>::start(n_threads);
}

void AsyncSequenceProcessor::stop() {
    aparser->stop();
    AsyncConsumerProducer<ReadBatchPtr,ReadBatchPtr>::stop();
}

unsigned int AsyncSequenceProcessor::n_processed() {
    return _n_processed;
}

unsigned int AsyncSequenceProcessor::n_parsed() {
    return aparser->n_parsed();
}

bool AsyncSequenceProcessor::is_paired() {
    return _paired;
}

unsigned int AsyncSequenceProcessor::parser_queue_load() {
    return aparser->queue_load();
}

/////
//
// AsyncSequenceProcessorTester
//
/////

bool AsyncSequenceProcessorTester::iter_stop() {
    if((_STATE == STATE_WAIT) && !(aparser->has_output()))
        return true;
    return false;
}

void AsyncSequenceProcessorTester::consume() {
    ReadBatchPtr batch;
    const char * sp;
    while(_STATE == STATE_RUNNING) {
        if(aparser->pop(batch)) {
            sp = copy_seq(batch->first());
            write(sp);            
            if (_paired) {
                sp = copy_seq(batch->second());
                write(sp);
            }
            while(!(_out_queue->bounded_push(batch))) if (_STATE != STATE_RUNNING) return;
            __sync_fetch_and_add(&_n_processed, _batchsize); 
        } else {
            if ((aparser->get_state() == STATE_WAIT) && !(aparser->has_output())) {
                _STATE = STATE_WAIT;
            }
        }
    }
}


