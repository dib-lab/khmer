#include "async_parser.hh"

using namespace khmer;
using namespace khmer::read_parsers;

inline ReadBatchPtr imprint_paired(IParser * parser) {
    ReadPair read_pair;
    parser->imprint_next_read_pair(read_pair);
    ReadPtr first( new Read(read_pair.first) );
    ReadPtr second( new Read(read_pair.second) );
    ReadBatchPtr batch = new ReadBatch(first, second);
    return batch;
}

inline ReadBatchPtr imprint_single(IParser * parser) {
    std::shared_ptr<Read> read(new Read());
    parser->imprint_next_read(*read);
    ReadBatchPtr batch = new ReadBatch(read);
    return batch;
}

void AsyncParser::start(const std::string &filename,
                        bool paired) {
    _batchsize = paired ? 2 : 1;
    _n_parsed = 0;
    _current_filename = filename;
    AsyncProducer<ReadBatchPtr>::start(1);
}

void AsyncParser::consume(const std::string &filename) {
    
    set_global_state(STATE_RUNNING);

    IParser * parser = IParser::get_parser(filename);
    // Use a function ptr to decide which imprint function to use
    // so as to avoid unnecssary branching
    ReadBatchPtr (*imprint) (IParser *);
    if (paired) {
        imprint = &imprint_paired;
    } else {
        imprint = &imprint_single;
    }

    ReadBatchPtr batch;
    while(!parser->is_complete() && (_STATE == STATE_RUNNING)) {
        try {
            batch = imprint(parser);
        } catch (...) {
            // Exception, probably from read pairing; by default, we
            // just transfer the exception and die
            _exc_handler->push( std::current_exception() );
            // TODO: flush queue here
            return;
        }
        __sync_fetch_and_add(&_n_parsed, _batchsize);

        while(!(_out_queue->bounded_push(batch))) {
            // Somebody turned us off; flush the queue and return gracefully
            if (_STATE != STATE_RUNNING) {
                #if(VERBOSITY)
                lock_stdout();
                std::cout << "Flush readparser queue" << std::endl;
                unlock_stdout();
                #endif
                delete parser;
                flush_queue<ReadBatchPtr>(_out_queue);
                return;
            }
        }

        #if(VERBOSITY)
        if (_n_parsed % 10000 == 0) std::cout << "...parsed " << _n_parsed << std::endl;
        #endif
    }
    // When we're done, set our status to false to notify the other threads
    set_global_state(STATE_WAIT);
    #if(VERBOSITY)
    lock_stdout();
    std::cout << "Finished parsing " << _n_parsed << " reads." << std::endl;
    unlock_stdout();
    #endif
    delete parser;
}


unsigned int AsyncParser::n_parsed() {
    return _n_parsed;
}
