#ifndef KHMER_ASYNC_DIGINORM_HH
#define KHMER_ASYNC_DIGINORM_HH

#include "khmer_async.hh"
#include "async_sequence_processor.hh"

namespace khmer {

class AsyncDiginorm: public AsyncSequenceProcessor {

    protected:

        unsigned int _cutoff;
        unsigned int _n_kept;
        #if(TIMING)
        double output_push_wait_global, reader_pop_wait_global, write_wait_global;
        uint32_t _timing_lock;
        #endif

    public:

        AsyncDiginorm (khmer::Hashtable * ht):
                        khmer::AsyncSequenceProcessor(ht) {
            #if(TIMING)
            _timing_lock = 0;
            output_push_wait_global = reader_pop_wait_global = write_wait_global = 0.0;
            #endif
        }

        void start(const std::string &filename,
                    unsigned int cutoff,
                    bool paired,
                    unsigned int n_threads);
        void stop();

        unsigned int n_kept();
        virtual void consume();

        bool filter_single(ReadPtr read);
        bool filter_paired(ReadBatchPtr read);
        unsigned int output_queue_load();
};

}; //namespace khmer

#endif
