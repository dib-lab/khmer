#include "khmer_async.hh"
#include "async_sequence_processor.hh"

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
    this->paired = paired;
    _batchsize = paired ? 2 : 1;
    _reader_thread = new std::thread(&AsyncSequenceProcessor::read_iparser,
                                    this, filename);
    _n_processed = 0;
    AsyncConsumerProducer<ReadBatchPtr,ReadBatchPtr>::start(n_threads);
}

void AsyncSequenceProcessor::stop() {
    _parsing_reads = false;
    // First stop the sequence processor, flushing its queue
    AsyncConsumerProducer<ReadBatchPtr,ReadBatchPtr>::stop();
    //flush_queue<ReadBatchPtr>(_out_queue);
}

inline  imprint_paired(IParser * parser) {
    ReadPair read_pair;
    //ReadBatchPtr batch;
    //ReadPtr first = new Read();
    //ReadPtr second = new Read();
    parser->imprint_next_read_pair(read_pair);
    std::shared_ptr<Read> first( new Read(read_pair.first) );
    std::shared_ptr<Read> second( new Read(read_pair.second) );
    ReadBatchPtr batch = new ReadBatch(first, second);
    return batch;
}

inline ReadBatchPtr imprint_single(IParser * parser) {
    //ReadPtr read = new Read();
    std::shared_ptr<Read> read(new Read());
    //ReadBatchPtr batch;
    parser->imprint_next_read(*read);
    ReadBatchPtr batch = new ReadBatch(read);
    return batch;
}

void AsyncSequenceProcessor::read_iparser(const std::string &filename) {
    #if(VERBOSITY)
    std::cout << "Spawned iparser thread..." << std::endl;
    #endif
    _n_parsed = 0;
    _parsing_reads = true;

    IParser * parser = IParser::get_parser(filename);
    
    ReadBatchPtr (*imprint) (IParser *);
    if (paired) {
        imprint = &imprint_paired;
    } else {
        imprint = &imprint_single;
    }

    ReadBatchPtr batch;
    while(!parser->is_complete() && _parsing_reads) {
        try {
            batch = imprint(parser);
        } catch (...) {
            // Exception, probably from read pairing; by default, we
            // just transfer the exception and die
            _exc_handler->push( std::current_exception() );
            // Should probably flush queue here
            return;
        }
        __sync_fetch_and_add(&_n_parsed, _batchsize);

        while(!(_in_queue->bounded_push(batch))) {
            // Somebody turned us off; flush the queue and return gracefully
            if (!_parsing_reads) {
                #if(VERBOSITY)
                lock_stdout();
                std::cout << "Flush readparser queue" << std::endl;
                unlock_stdout();
                #endif
                delete parser;
                flush_queue<ReadBatchPtr>(_in_queue);
                return;
            }
        }

        #if(VERBOSITY)
        if (_n_parsed % 10000 == 0) std::cout << "...parsed " << _n_parsed << std::endl;
        #endif
    }
    // When we're done, set our status to false to notify the other threads
    _parsing_reads = false;
    #if(VERBOSITY)
    lock_stdout();
    std::cout << "Finished parsing " << _n_parsed << " reads." << std::endl;
    unlock_stdout();
    #endif
    delete parser;
}

bool AsyncSequenceProcessor::is_paired() {
    return paired;
}

bool AsyncSequenceProcessor::is_parsing() {
    return _parsing_reads;
}

unsigned int AsyncSequenceProcessor::n_parsed() {
    return _n_parsed;
}

unsigned int AsyncSequenceProcessor::n_processed() {
    return _n_processed;
}

unsigned int AsyncSequenceProcessor::reader_queue_load() {
    return _n_parsed - _n_processed;
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


