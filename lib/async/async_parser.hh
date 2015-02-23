#ifndef KHMER_ASYNC_PARSER_HH
#define KHMER_ASYNC_PARSER_HH

#include "khmer_async.hh"
#include "async_models.hh"

namespace khmer {

class AsyncSequenceParser: public AsyncProducer<ReadBatchPtr> {
    
    protected:

        bool paired = false;
        unsigned int _n_parsed;
        unsigned int _batchsize;
        const std::string _current_filename;

        AsyncSequenceParser ():
            khmer::AsyncProducer<ReadBatchPtr>() {
        }

        void start(const std::string &filename, bool paired);
        void consume();
        unsigned int queue_load();
        unsigned int get_batchsize() { return _batchsize; };

}

} // namespace khmer
#endif
