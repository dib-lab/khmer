#include <Python.h>

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
// class AsyncHashWriter methods
//
/////

void AsyncHashWriter::start() {
    _n_written = 0;
    AsyncConsumer<HashIntoType>::start(1);
}

unsigned int AsyncHashWriter::ksize() {
    return _ht->ksize();
}

void AsyncHashWriter::consume() {
    HashIntoType khash;
    //std::cout << "Thread " << std::this_thread::get_id() << " waiting on hashes" << std::endl;
    while(1) {
        if (_in_queue->pop(khash)) {
            _ht->count(khash);
            _n_written++;
        } else {
            if(!_workers_running && (_n_written >= _n_pushed)) return;
        }
    }
}

unsigned int AsyncHashWriter::n_written() {
    return _n_written;
}

/////
//
// AsyncSequenceWriter methods
//
/////


void AsyncSequenceWriter::start() {
    _n_written = 0;
    _ksize = _ht->ksize();
    AsyncConsumer<const char *>::start(1);
}

unsigned int AsyncSequenceWriter::ksize() {
    return _ht->ksize();
}

void AsyncSequenceWriter::consume() {
    const char * sequence;
    HashIntoType kmer;
    CountFunc kmer_counter(_ht, _ksize);
    //std::cout << "Thread " << std::this_thread::get_id() << " waiting on hashes" << std::endl;
    while(_workers_running) {
        if (_in_queue->pop(sequence)) {
            khmer::KMerIterator kmers(sequence, _ksize);
            try {
                while(!kmers.done()) {
                    kmer = kmers.next();
                    _ht->count(kmer);
                }
            } catch (khmer_exception &e) {
                _exc_handler->push(std::make_exception_ptr(e));
                std::cout << e.what() << " " << std::endl;
                std::cout << "ERROR in AsyncSequenceWriter: " << sequence << std::endl;
                exit(1);
                //flush_queue<const char *>(_in_queue);
                return;
            }
            _n_written++;
        }
    }
    #if(VERBOSITY)
    std::cout << "Exit AsyncSequenceWriter consume thread." << std::endl;
    #endif
    _in_queue->consume_all(kmer_counter);
}

unsigned int AsyncSequenceWriter::n_written() {
    return _n_written;
}

/////
//
// class AynscBundledHashWriter
//
/////

void AsyncBundledHashWriter::start() {
    _n_written = 0;
    _ksize = _ht->ksize();
    AsyncConsumer<std::vector<HashIntoType>*>::start(1);
}


unsigned int AsyncBundledHashWriter::ksize() {
    return _ht->ksize();
}

unsigned int AsyncBundledHashWriter::n_written() {
    return _n_written;
}

void AsyncBundledHashWriter::count(std::vector<HashIntoType> * bundle) {
    try {
        for (std::vector<HashIntoType>::iterator it=bundle->begin(); it!=bundle->end(); ++it) {
            _ht->count(*it);
        }
    } catch (khmer_exception &e) {
        _exc_handler->push(std::make_exception_ptr(e));
    }
}

void AsyncBundledHashWriter::consume() {
    
    std::vector<HashIntoType> * bundle;
    while(_workers_running) {
        if(_in_queue->pop(bundle)) {
            count(bundle);
            _n_written++;
        }
    }
    _in_queue->consume_all( [this](std::vector<HashIntoType> * bundle){ count(bundle); });
}


/////
//
// class AsyncHasher methods
//
/////


// Warning: this method will always consume the entire
// character array as decomposed k-mers! So, if you want to consume
// an individual k-mer, you have to push appropriately sized arrays
// onto the queue. You shouldn't do this though, because more queue
// accesses is more overhead.

void AsyncHasher::consume() {
    const char * kmer;
    unsigned long long total_consumed = 0;
    //std::cout << "Thread " << std::this_thread::get_id() << " waiting on hashes" << std::endl;
    while(_workers_running) {
        if (_in_queue->pop(kmer)) {
            khmer::KMerIterator kmers(kmer, _ksize);
            while(!kmers.done()) {
                while(!(_out_queue->bounded_push(kmers.next()))) if (!_workers_running) return;
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
                                    bool paired,
                                    unsigned int n_threads) {

    this->paired = paired;
    _batchsize = paired ? 2 : 1;
    _reader_thread = new std::thread(&AsyncSequenceProcessor::read_iparser, 
                                    this, filename);
    _n_processed = 0;
    _writer->start();
    AsyncConsumerProducer<ReadBatch*,ReadBatch*>::start(n_threads);
}

void AsyncSequenceProcessor::stop() {
    _parsing_reads = false;
    // First stop the sequence processor, flushing its queue
    AsyncConsumerProducer<ReadBatch*,ReadBatch*>::stop();
    // Then stop the writer, writing out all the flushed hashes
    _writer->stop();
    flush_queue<ReadBatch*>(_out_queue);
}

inline ReadBatch * imprint_paired(IParser * parser) {
    ReadPair read_pair;
    ReadBatch * batch;
    parser->imprint_next_read_pair(read_pair);
    batch = new ReadBatch(new Read(read_pair.first), new Read(read_pair.second));
    return batch;
}

inline ReadBatch * imprint_single(IParser * parser) {
    Read * read = new Read();
    ReadBatch * batch;
    parser->imprint_next_read(*read);
    batch = new ReadBatch(read);
    return batch;
}

void AsyncSequenceProcessor::read_iparser(const std::string &filename) {
    #if(VERBOSITY)
    std::cout << "Spawned iparser thread..." << std::endl;
    #endif
    _n_parsed = 0;
    _parsing_reads = true;

    IParser * parser = IParser::get_parser(filename);
    
    ReadBatch * (*imprint) (IParser *);
    if (paired) {
        imprint = &imprint_paired;
    } else {
        imprint = &imprint_single;
    }

    ReadBatch * batch;
    while(!parser->is_complete() && _parsing_reads) {
        try {
            batch = imprint(parser);
        } catch (...) {
            // Exception, probably from read pairing; by default, we
            // just transfer the exception and die
            _exc_handler->push( std::current_exception() );
            // Should probably flush queue here
            return;
        }
        __sync_fetch_and_add(&_n_parsed, _batchsize);

        while(!(_in_queue->bounded_push(batch))) {
            // Somebody turned us off; flush the queue and return gracefully
            if (!_parsing_reads) {
                #if(VERBOSITY)
                lock_stdout();
                std::cout << "Flush readparser queue" << std::endl;
                unlock_stdout();
                #endif
                delete parser;
                flush_queue<ReadBatch*>(_in_queue);
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

bool AsyncSequenceProcessor::is_paired() {
    return paired;
}

bool AsyncSequenceProcessor::is_parsing() {
    return _parsing_reads;
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


std::vector<HashIntoType> * AsyncSequenceProcessor::get_bundle(Read * read) {
    std::vector<HashIntoType> * bundle = new std::vector<HashIntoType>();
    HashIntoType kmer;
    khmer::KMerIterator kmers(read->sequence.c_str(), _ht->ksize());
    while(!kmers.done()) {
        kmer = kmers.next();
        bundle->push_back(kmer);
    }
    return bundle;
}

/////
//
// AsyncSequenceProcessorTester
//
/////

bool AsyncSequenceProcessorTester::iter_stop() {
    if(!workers_running() && (n_popped() >= n_processed()))
        return true;
    return false;
}

void AsyncSequenceProcessorTester::consume() {
    ReadBatch * batch;
    std::vector<HashIntoType> * bundle;
    while(_workers_running) {
        if(_in_queue->pop(batch)) {

            bundle = get_bundle(batch->first());
            while(!(_writer->push(bundle))) if (!_workers_running) return;
            if (paired) {
                bundle = get_bundle(batch->second());
                while(!(_writer->push(bundle))) if (!_workers_running) return;
            }
            while(!(_out_queue->bounded_push(batch))) if (!_workers_running) return;
            __sync_fetch_and_add(&_n_processed, _batchsize); 
        } else {
            if (!is_parsing() && (_n_processed >= _n_parsed)) {
                while(_writer->n_pushed() > _writer->n_written()) if (!_workers_running) return;
                _workers_running = false;
            }
        }
    }
}
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

bool AsyncDiginorm::iter_stop() {
    if (!workers_running() && (n_popped() >= n_kept()))
        return true;
    return false;
}

bool AsyncDiginorm::filter_single(Read* read) {
    BoundedCounterType median = 0;
    float average = 0, stddev = 0;
    _ht->get_median_count(read->sequence, median, average, stddev);
    if (median < _cutoff) return false;
    return true;
}

bool AsyncDiginorm::filter_paired(ReadBatch * batch) {
    bool filter_first = false;
    filter_first = filter_single(batch->first());
    return filter_first && filter_single(batch->second());
}


void AsyncDiginorm::consume() {

    ReadBatch* batch;
    bool filter;

    #if(TIMING)
    timespec start_t, end_t;
    double output_push_wait = 0.0, writer_push_wait = 0.0, reader_pop_wait = 0.0;
    #endif

    #if(VERBOSITY)
    lock_stdout();
    std::thread::id tid = std::this_thread::get_id();
    std::cout << "Thread " << tid << " : " << std::endl 
        << "\tCUTOFF: " << _cutoff << "\tK: " << _ht->ksize() << std::endl;
    unlock_stdout();
    #endif

    //const char * sp;
    std::vector<HashIntoType> * bundle;

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

                bundle = get_bundle(batch->first());
                
                TSTART()
                while(!(_writer->push(bundle))) if (!_workers_running) return;
                TEND(writer_push_wait)
                if (paired) {
                    //sp = copy_seq(batch->second());
                    bundle = get_bundle(batch->second());
                    TSTART()
                    while(!(_writer->push(bundle))) if (!_workers_running) return;
                    TEND(writer_push_wait)
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
                std::cout << reader_pop_wait << "\t" << writer_push_wait << "\t" << output_push_wait << std::endl;
                unlock_stdout();
                #endif
                while(_writer->n_pushed() > _writer->n_written()) if (!_workers_running) return;
                _workers_running = false;
            }
        }
    }

    #if(VERBOSITY)
    lock_stdout();
    std::cout << "\nReturning from AsyncDiginorm worker..." << std::endl;
    std::cout << "\t(hashes pushed, hashes written): " << 
        _writer->n_pushed() << ", " << _writer->n_written() << std::endl;
    unlock_stdout();
    #endif
}
