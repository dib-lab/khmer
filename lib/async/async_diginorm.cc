#include "async_diginorm.hh"

/////
//
// AsyncDiginorm
//
/////

void AsyncDiginorm::start(const std::string &filename,
                            unsigned int cutoff,
                            bool paired,
                            unsigned int n_threads) {
    _cutoff = cutoff;
    _n_kept = 0;
    _n_hashes_pushed = 0;
    AsyncSequenceProcessor::start(filename, paired, n_threads);
}

unsigned int AsyncDiginorm::n_kept() {
    return _n_kept;
}

unsigned int AsyncDiginorm::output_queue_load() {
    return _n_kept - _n_popped;
}

bool AsyncDiginorm::iter_stop() {
    if (!workers_running() && (n_popped() >= n_kept()))
        return true;
    return false;
}

bool AsyncDiginorm::filter_single(ReadPtr read) {
    BoundedCounterType median = 0;
    float average = 0, stddev = 0;
    _ht->get_median_count(read->sequence, median, average, stddev);
    if (median < _cutoff) return false;
    return true;
}

bool AsyncDiginorm::filter_paired(ReadBatchPtr batch) {
    bool filter_first = false;
    filter_first = filter_single(batch->first());
    return filter_first && filter_single(batch->second());
}



void AsyncDiginorm::consume() {

    ReadBatchPtr batch;
    bool filter;

    #if(TIMING)
    timespec start_t, end_t;
    double output_push_wait = 0.0, reader_pop_wait = 0.0, write_wait = 0.0;
    #endif

    #if(VERBOSITY)
    lock_stdout();
    std::thread::id tid = std::this_thread::get_id();
    std::cout << "Thread " << tid << " : " << std::endl 
        << "\tCUTOFF: " << _cutoff << "\tK: " << _ht->ksize() << std::endl;
    unlock_stdout();
    #endif

    while(_workers_running) {
        TSTART()
        if (_in_queue->pop(batch)) {
            TEND(reader_pop_wait)
            if (paired) {
                filter = filter_paired(batch);
            } else {
                filter = filter_single(batch->first());
            }

            if (!filter) {
                __sync_fetch_and_add(&_n_kept, _batchsize);

                //sp = copy_seq(batch->first());
                TSTART()
                write(batch->first()->sequence.c_str());
                TEND(write_wait)
                if (paired) {
                    //sp = copy_seq(batch->second());
                    TSTART()
                    write(batch->second()->sequence.c_str());
                    TEND(write_wait)
                }

                TSTART()
                while(!(_out_queue->bounded_push(batch))) if (!_workers_running) return;
                TEND(output_push_wait)

            } else {
                delete batch;
            }
            __sync_fetch_and_add(&_n_processed, _batchsize);
        } else {
            if (!is_parsing() && (_n_processed >= _n_parsed)) {
                #if(VERBOSITY)
                lock_stdout();
                std::cout << "Done processing (" << _n_processed 
                    << " reads). Waiting for writeout" << std::endl;
                unlock_stdout();
                #endif
                #if(TIMING)
                lock_stdout();
                std::cout << reader_pop_wait << "\t" << output_push_wait << "\t" << write_wait << std::endl;
                unlock_stdout();
                #endif
                _workers_running = false;
            }
        }
    }

    #if(VERBOSITY)
    lock_stdout();
    std::cout << "\nReturning from AsyncDiginorm worker..." << std::endl;
    unlock_stdout();
    #endif
}
