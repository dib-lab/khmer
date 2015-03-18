#ifndef KHMER_ASYNC_PARSER_HH
#define KHMER_ASYNC_PARSER_HH

#include "khmer_async.hh"
#include "async_models.hh"

namespace khmer {

class AsyncSequenceParser: public AsyncProducer<ReadBatchPtr> {
    
    protected:

        bool _paired;
        unsigned int _n_parsed;
        unsigned int _batchsize;
        std::string _current_filename;

    public:

        const char * name = "AsyncSequenceParser";

        AsyncSequenceParser ():
            khmer::AsyncProducer<ReadBatchPtr>() {
        }

        void start(const std::string &filename, bool paired);
        void consume();
        unsigned int queue_load();
        unsigned int get_batchsize() { return _batchsize; };
        unsigned int n_parsed();
        
        void set_global_state(int state) {
            AsyncProducer<ReadBatchPtr>::set_global_state(state);
            #if(VERBOSITY)
            lock_stdout();
            std::cout << name << " STATE=" << _STATE << std::endl;
            unlock_stdout();
            #endif 
        }
};

} // namespace khmer
#endif
