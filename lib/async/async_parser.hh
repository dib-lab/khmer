#ifndef KHMER_ASYNC_PARSER_HH
#define KHMER_ASYNC_PARSER_HH

#include "khmer_async.hh"
#include "async_models.hh"

namespace khmer {

class AsyncSequenceParser: public AsyncProducer<ReadBatchPtr> {
    
    protected:

        bool paired = false;
        unsigned int _n_parsed;
        unsigned int _n_processed;

        AsyncSequenceParser ():
            khmer::AsyncProducer<ReadBatchPtr>() {
        }

        void start(const std::string &filename

}

} // namespace khmer
#endif
