#include "async_diginorm.hh"

/////
//
// AsyncDiginorm
//
/////

using namespace khmer;

void AsyncDiginorm::start(const std::string &filename,
                            unsigned int cutoff,
                            bool paired,
                            unsigned int n_threads) {
    _cutoff = cutoff;
    _n_kept = 0;

    AsyncSequenceProcessor::start(filename, paired, n_threads);
}

void AsyncDiginorm::stop() {
    #if(VERBOSITY)
    lock_stdout();
    std::cout << "Done processing (" << _n_processed 
        << " reads). Waiting for writeout" << std::endl;
    unlock_stdout();
    #endif
    #if(TIMING)
    lock_stdout();
    std::cout << "\n** TIMINGS **\n\tRead time: " << reader_pop_wait_global << "\n\tOutput time: " << output_push_wait_global << "\n\tWrite time: " << write_wait_global << std::endl << std::endl;
    unlock_stdout();
    #endif
    AsyncSequenceProcessor::stop();
}

unsigned int AsyncDiginorm::n_kept() {
    return _n_kept;
}

unsigned int AsyncDiginorm::output_queue_load() {
    return _n_kept - _n_popped;
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
        << "\tCUTOFF: " << _cutoff << "\tK: " << _ksize << std::endl;
    unlock_stdout();
    #endif

    //print("spin on checking_running")
    //while(!check_running());

    print("about to enter while loop in consume")
    std::cout << "aparser state is " << aparser->get_state() << std::endl;
    while(check_running()) {
        TSTART()
        if (aparser->pop(batch)) {
            TEND(reader_pop_wait)
            if (_paired) {
                filter = filter_paired(batch);
            } else {
                filter = filter_single(batch->first());
            }

            if (!filter) {
                __sync_fetch_and_add(&_n_kept, _batchsize);

                TSTART()
                write(batch->first()->sequence.c_str());
                TEND(write_wait)
                
		if (_paired) {
                    TSTART()
                    write(batch->second()->sequence.c_str());
                    TEND(write_wait)
                }

                TSTART()
                while(!(_out_queue->bounded_push(batch))) if (!check_running()) return;
                TEND(output_push_wait)

            } else {
                delete batch;
            }
            __sync_fetch_and_add(&_n_processed, _batchsize);
        } else {
            if (((aparser->get_state() == STATE_WAIT) && !(aparser->has_output()))
                || aparser->get_state() == STATE_DORMANT) {

                set_global_state(STATE_WAIT);
            }
        }
    }

    #if(VERBOSITY)
    lock_stdout();
    std::cout << "\nReturning from AsyncDiginorm worker. STATE: " << 
        _STATE << "; aparser STATE: " << aparser->get_state() << 
        "n_processed: " << _n_processed << std::endl;
    unlock_stdout();
    #endif

    #if(TIMING)
    while(!__sync_bool_compare_and_swap( &_timing_lock, 0, 1));
    output_push_wait_global += output_push_wait;
    reader_pop_wait_global += reader_pop_wait;
    write_wait_global += write_wait;
    __sync_bool_compare_and_swap( &_timing_lock, 1, 0);
    #endif
}
