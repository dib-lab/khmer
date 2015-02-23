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
        AsyncParser * aparser;

    public:

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
        void write(const char * sequence);

        virtual bool iter_stop() = 0;

        bool is_paired();
};

class AsyncSequenceProcessorTester: public AsyncSequenceProcessor {

    public:

        AsyncSequenceProcessorTester (khmer::Hashtable * ht):
                        khmer::AsyncSequenceProcessor(ht) {
        }

        virtual void consume();
        virtual bool iter_stop();
};
}; // namespace khmer
#endif
