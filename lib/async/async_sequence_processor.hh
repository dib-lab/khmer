#ifndef KHMER_ASYNC_SP_HH
#define KHMER_ASYNC_SP_HH

#include "khmer_async.hh"
#include "async_models.hh"
#include "async_parser.hh"

namespace khmer {

// General read processing model
// TODO: Split off the read_iparser functionality into its own AsyncParser class
class AsyncSequenceProcessor: public AsyncConsumerProducer<ReadBatchPtr, ReadBatchPtr> {

    protected:

        khmer::Hashtable * _ht;
        unsigned int _ksize;
        AsyncSequenceParser * aparser;
        unsigned int _n_processed;
        bool _paired;


    public:

        const char * name = "AsyncSequenceProcessor";

        AsyncSequenceProcessor (khmer::Hashtable * ht):
                                khmer::AsyncConsumerProducer<ReadBatchPtr,ReadBatchPtr>(),
                                _ht(ht) {
        }

        void start(const std::string &filename,
                    bool paired,
                    unsigned int n_threads);
        void stop();

        virtual void consume() = 0;
        unsigned int n_processed();
        unsigned int n_parsed();
        void write(const char * sequence);
        unsigned int parser_queue_load();

        bool iter_stop();

        void set_global_state(int state) {
            AsyncConsumerProducer<ReadBatchPtr, ReadBatchPtr>::set_global_state(state);
            #if(VERBOSITY)
            lock_stdout();
            std::cout << name << " STATE=" << _STATE << std::endl;
            unlock_stdout();
            #endif 
        }
        bool is_paired();
};

class AsyncSequenceProcessorTester: public AsyncSequenceProcessor {

    public:

        const char * name = "AsyncSequenceProcessorTester";

        AsyncSequenceProcessorTester (khmer::Hashtable * ht):
                        khmer::AsyncSequenceProcessor(ht) {
        }

        virtual void consume();
};
}; // namespace khmer
#endif
