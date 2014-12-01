#ifndef KHMER_ASYNC_DIGINORM_HH
#define KHMER_ASYNC_DIGINORM_HH

#include "khmer_async.hh"
#include "async_sequence_processor.hh"

namespace khmer {

class AsyncDiginorm: public AsyncSequenceProcessor {

    protected:

        unsigned int _cutoff;
        unsigned int _n_hashes_pushed;
        unsigned int _n_kept;

    public:

        AsyncDiginorm (khmer::Hashtable * ht):
                        khmer::AsyncSequenceProcessor(ht) {
        }

        void start(const std::string &filename,
                    unsigned int cutoff,
                    bool paired,
                    unsigned int n_threads);

        unsigned int n_kept();
        virtual void consume();
        bool iter_stop();

        bool filter_single(ReadPtr read);
        bool filter_paired(ReadBatchPtr read);
        unsigned int output_queue_load();
};

}; //namespace khmer

#endif
