#include "hashtable.hh"
#include "khmer_async.hh"
#include <atomic>
#include <chrono>
#include <time.h>
#include <boost/lockfree/queue.hpp>


using namespace khmer;
using namespace khmer::read_parsers;
using namespace boost::lockfree;

timespec timediff(timespec start, timespec end)
{
    timespec temp;
    if ((end.tv_nsec-start.tv_nsec)<0) {
	temp.tv_sec = end.tv_sec-start.tv_sec-1;
	temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
	temp.tv_sec = end.tv_sec-start.tv_sec;
	temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    return temp;
}

/////
//
// class AsyncWriter methods
//
/////

void AsyncWriter::start() {
    _n_pushed = 0;
    _n_written = 0;
    Async<HashIntoType>::start(1);
}

unsigned int AsyncWriter::ksize() {
    return _ht->ksize();
}

void AsyncWriter::consume(HashQueue * q) {
    HashIntoType khash;
    //std::cout << "Thread " << std::this_thread::get_id() << " waiting on hashes" << std::endl;
    while(1) {
        if (q->pop(khash)) {
            _ht->count(khash);
            _n_written++;
        } else {
            if(!_workers_running && (_n_written >= _n_pushed)) return;
        }
    }
}

bool AsyncWriter::push(HashIntoType &khash) {
    if(_in_queue->bounded_push(khash)) {
        __sync_fetch_and_add(&_n_pushed, 1);
        return true;
    }
    return false;
}

unsigned int AsyncWriter::n_pushed() {
    return _n_pushed;
}

unsigned int AsyncWriter::n_written() {
    return _n_written;
}

/////
//
// AsyncSequenceWriter methods
//
/////


void AsyncSequenceWriter::start() {
    _n_pushed = 0;
    _n_written = 0;
    _ksize = _ht->ksize();
    Async<const char *>::start(1);
}

unsigned int AsyncSequenceWriter::ksize() {
    return _ht->ksize();
}

void AsyncSequenceWriter::consume(CharQueue * q) {
    const char * sequence;
    HashIntoType kmer;
    CountFunc kmer_counter(_ht, _ksize);
    //std::cout << "Thread " << std::this_thread::get_id() << " waiting on hashes" << std::endl;
    while(_workers_running) {
        if (q->pop(sequence)) {
            khmer::KMerIterator kmers(sequence, _ksize);
            try {
                while(!kmers.done()) {
                    kmer = kmers.next();
                    _ht->count(kmer);
                }
            } catch (khmer_exception &e) {
                std::cout << e.what() << " " << std::endl;
                std::cout << "ERROR in AsyncSequenceWriter: " << sequence << std::endl;
                exit(1);
            }
            _n_written++;
        }
    }
    #if(VERBOSITY)
    std::cout << "Exit AsyncSequenceWriter consume thread." << std::endl;
    #endif
    q->consume_all(kmer_counter);
}


bool AsyncSequenceWriter::push(const char * &sequence) {
    if(_in_queue->bounded_push(sequence)) {
        __sync_fetch_and_add(&_n_pushed, 1);
        return true;
    }
    return false;
}

unsigned int AsyncSequenceWriter::n_pushed() {
    return _n_pushed;
}

unsigned int AsyncSequenceWriter::n_written() {
    return _n_written;
}

/////
//
// class AsyncHasher methods
//
/////

bool AsyncHasher::pop(HashIntoType& khash) {
    return _out_queue->pop(khash);
}

void AsyncHasher::set_output(HashQueue* new_q) {
    _out_queue = new_q;
}

// Warning: this method will always consume the entire
// character array as decomposed k-mers! So, if you want to consume
// an individual k-mer, you have to push appropriately sized arrays
// onto the queue. You shouldn't do this though, because more queue
// accesses is more overhead.

void AsyncHasher::consume(CharQueue * q) {
    const char * kmer;
    unsigned long long total_consumed = 0;
    //std::cout << "Thread " << std::this_thread::get_id() << " waiting on hashes" << std::endl;
    while(_workers_running) {
        if (q->pop(kmer)) {
            khmer::KMerIterator kmers(kmer, _ksize);
            while(!kmers.done()) {
                _out_queue->bounded_push(kmers.next());
                total_consumed++;
            }
        }
    }
}


/////
//
// AsyncSequenceProcessor
//
/////

void AsyncSequenceProcessor::start(const std::string &filename,
                                    unsigned int n_threads) {

    _reader_thread = new std::thread(&AsyncSequenceProcessor::read_iparser, 
                                    this, filename);

    _n_processed = 0;
    _n_popped = 0;

    _writer->start();
    Async<Read*>::start(n_threads);
}

void AsyncSequenceProcessor::stop() {
    _parsing_reads = false;
    // First stop the sequence processor, flushing its queue
    Async<Read*>::stop();
    // Then stop the writer, writing out all the flushed hashes
    _writer->stop();
    flush_queue<Read*>(_out_queue);
}

void AsyncSequenceProcessor::read_iparser(const std::string &filename) {
    #if(VERBOSITY)
    std::cout << "Spawned iparser thread..." << std::endl;
    #endif
    _n_parsed = 0;
    _parsing_reads = true;

    IParser * parser = IParser::get_parser(filename);

    while(!parser->is_complete() && _parsing_reads) {
        Read * read = new Read();
        parser->imprint_next_read(*read);
        __sync_fetch_and_add(&_n_parsed, 1);
        while(!(_in_queue->bounded_push(read))) {
            if (!_parsing_reads) {
                flush_queue<Read*>(_in_queue);
                return;
            }
        }
        #if(VERBOSITY)
        if (_n_parsed % 10000 == 0) std::cout << "...parsed " << _n_parsed << std::endl;
        #endif
    }
    // When we're done, set our status to false to notify the other threads
    _parsing_reads = false;
    #if(VERBOSITY)
    lock_stdout();
    std::cout << "Finished parsing " << _n_parsed << " reads." << std::endl;
    unlock_stdout();
    #endif
    delete parser;
}

void AsyncSequenceProcessor::set_output(ReadQueue* new_q) {
    _out_queue = new_q;
}

bool AsyncSequenceProcessor::is_parsing() {
    return _parsing_reads;
}

bool AsyncSequenceProcessor::has_output() {
    return !_out_queue->empty();
}

bool AsyncSequenceProcessor::pop(Read* &read) {
    if(_out_queue->pop(read)) {
        _n_popped++;
        return true;
    }
    return false;
}

unsigned int AsyncSequenceProcessor::n_popped() {
    return _n_popped;
}

unsigned int AsyncSequenceProcessor::n_parsed() {
    return _n_parsed;
}

unsigned int AsyncSequenceProcessor::n_processed() {
    return _n_processed;
}

unsigned int AsyncSequenceProcessor::n_written() {
    return _writer->n_written();
}

/////
//
// AsyncDiginorm
//
/////

void AsyncDiginorm::start(const std::string &filename,
                            unsigned int cutoff,
                            unsigned int n_threads) {
    _cutoff = cutoff;
    _n_kept = 0;
    _n_hashes_pushed = 0;
    AsyncSequenceProcessor::start(filename, n_threads);
}

unsigned int AsyncDiginorm::n_kept() {
    return _n_kept;
}

void AsyncDiginorm::consume(ReadQueue * q) {

    Read* read;
    BoundedCounterType median;
    float average, stddev;

    #if(VERBOSITY)
    timespec start_t, end_t;
    double output_push_wait = 0.0, hash_push_wait = 0.0;

    lock_stdout();
    std::thread::id tid = std::this_thread::get_id();
    std::cout << "Thread " << tid << " : " << std::endl 
        << "\tCUTOFF: " << _cutoff << "\tK: " << _ht->ksize() << std::endl;
    unlock_stdout();
    #endif

    std::string sequence;
    const char * sp;

    while(_workers_running) {
        if (q->pop(read)) {
            _ht->get_median_count(read->sequence, median, average, stddev);
            //std::cout << "Read " << read->sequence << "\n\tMED: " << median << 
            //    " (kept " << _n_kept << " thus far)" << std::endl;
            if (median < _cutoff) {
                __sync_fetch_and_add(&_n_kept, 1);

                std::string sequence(read->sequence);
                sp = sequence.c_str();

                #if(VERBOSITY)
                clock_gettime(CLOCK_MONOTONIC, &start_t);
                #endif
                while(!(_out_queue->bounded_push(read))) if (!_workers_running) return;
                #if(VERBOSITY)
                clock_gettime(CLOCK_MONOTONIC, &end_t);
                output_push_wait += (timediff(start_t,end_t).tv_nsec / 1000000000.0);
                clock_gettime(CLOCK_MONOTONIC, &start_t);
                #endif
                while(!(_writer->push(sp))) if (!_workers_running) return;;
                __sync_fetch_and_add(&_n_hashes_pushed, 1);
                #if(VERBOSITY)
                clock_gettime(CLOCK_MONOTONIC, &end_t);
                hash_push_wait += (timediff(start_t,end_t).tv_nsec / 1000000000.0);
                #endif
            } else {
                delete read;
            }
            __sync_fetch_and_add(&_n_processed, 1);
        } else {
            if (!is_parsing() && (_n_processed >= _n_parsed)) {
                #if(VERBOSITY)
                lock_stdout();
                std::cout << "Done processing (" << _n_processed 
                    << " reads). Waiting for writeout" << std::endl;
                unlock_stdout();
                #endif
                while(_n_hashes_pushed > _writer->n_written()) if (!_workers_running) return;
                _workers_running = false;
            }
        }
    }

    #if(VERBOSITY)
    lock_stdout();
    std::cout << "\nReturning from AsyncDiginorm worker..." << std::endl;
    std::cout << "\t(hashes pushed, hashes written): " << 
        _writer->n_pushed() << ", " << _writer->n_written() << std::endl;
    std::cout << "\tOutput push wait: " << 
        output_push_wait << 
        "s" << std::endl;
    std::cout << "\tHash push wait: " << 
        hash_push_wait << 
        "s" << std::endl;
    unlock_stdout();
    #endif
}
