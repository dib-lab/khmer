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

    if(!(_ht->is_threadsafe())) {
        _ht->init_threadstuff();
    }

    aparser = new AsyncSequenceParser();
    aparser->start(filename, paired);
    while(!aparser->check_running());

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

bool AsyncSequenceProcessor::iter_stop() {
    if(_STATE == STATE_WAIT ||
        aparser->get_state() == STATE_DORMANT) {
        
        #if(VERBOSITY)
        lock_stdout();
        std::cout << "AsyncSP iter_stop()" << std::endl;
        unlock_stdout();
        #endif
        
        return true;
    }
    return false;
}

/////
//
// AsyncSequenceProcessorTester
//
/////



void AsyncSequenceProcessorTester::consume() {
    ReadBatchPtr batch;
    const char * sp; 
    while(!check_running()); 
    while(check_running()) {
        if(aparser->pop(batch)) {
            sp = copy_seq(batch->first());
            write(sp);            
            if (_paired) {
                sp = copy_seq(batch->second());
                write(sp);
            }
            while(!(_out_queue->bounded_push(batch))) if (!check_running()) return;
            __sync_fetch_and_add(&_n_processed, _batchsize);
            if (_n_processed % 5000 == 0) std::cout << "processed " << _n_processed << std::endl;
        } else {
            if (((aparser->get_state() == STATE_WAIT) 
                && !(aparser->has_output()))
                || aparser->get_state() == STATE_DORMANT) {

                set_global_state(STATE_WAIT);
            }
        }
    }
}


