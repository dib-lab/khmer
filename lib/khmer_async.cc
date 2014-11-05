#include "hashtable.hh"
#include "khmer_async.hh"
#include <atomic>
#include <boost/lockfree/queue.hpp>

using namespace khmer;
using namespace khmer::read_parsers;
using namespace boost::lockfree;

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
    if(_in_queue->push(khash)) {
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
    while(1) {
        if (q->pop(kmer)) {
            khmer::KMerIterator kmers(kmer, _ksize);
            while(!kmers.done()) {
                _out_queue->push(kmers.next());
                total_consumed++;
            }
        } else {
            if (!_workers_running) return;
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

    _parsed_count = 0;
    _parsing_reads = true;
    _reader_thread = new std::thread(&AsyncSequenceProcessor::read_iparser, 
                                    this, filename);

    _processed_count = 0;
    _n_popped = 0;

    _writer->start();
    Async<Read*>::start(n_threads);
}

void AsyncSequenceProcessor::stop() {
    // First stop the sequence processor, flushing its queue
    Async<Read*>::stop();
    // Then stop the writer, writing out all the flushed hashes
    _writer->stop();
}

void AsyncSequenceProcessor::read_iparser(const std::string &filename) {
    std::cout << "Spawned iparser thread..." << std::endl;

    khmer:: Config &the_config = khmer::get_active_config();
    IParser * parser = IParser::get_parser(filename, 1, 
            the_config.get_reads_input_buffer_size(), 
            the_config.get_reads_parser_trace_level());

    while(!parser->is_complete()) {
        Read * readp = new Read();
        parser->imprint_next_read(*readp);
        _parsed_count++;
        while(!(_in_queue->push(readp)));
        if (_parsed_count % 10000 == 0) std::cout << "...parsed " << _parsed_count << std::endl;
    }
    // When we're done, set our status to false to notify the other threads
    _parsing_reads = false;
    std::cout << "Finished parsing..." << std::endl;
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

    unsigned int ksize = _writer->ksize();
    HashIntoType khash;

    std::thread::id tid = std::this_thread::get_id();
    std::cout << "Thread " << tid << " : " << std::endl 
        << "\tCUTOFF: " << _cutoff << "\tK: " << ksize << std::endl;

    while(1) {
        if (q->pop(read)) {
            _ht->get_median_count(read->sequence, median, average, stddev);
            //std::cout << "Read " << read->sequence << "\n\tMED: " << median << 
            //    " (kept " << _n_kept << " thus far)" << std::endl;
            if (median < _cutoff) {
                __sync_fetch_and_add(&_n_kept, 1);
                while(!(_out_queue->push(read)));
                khmer::KMerIterator kmers(read->sequence.c_str(), ksize);
                while(!kmers.done()) {
                    khash = kmers.next();
                    while(!(_writer->push(khash)));
                    __sync_fetch_and_add(&_n_hashes_pushed, 1);
                }
            } else {
                delete read;
            }
            __sync_fetch_and_add(&_processed_count, 1);
        } else {
            if (!is_parsing() && _processed_count >= _parsed_count) {
                std::cout << "Finished processing..." << std::endl;
                std::cout << "(hashes pushed, hashes written): " << 
                    _writer->n_pushed() << ", " << _writer->n_written() << std::endl;
                while(_n_hashes_pushed > _writer->n_written());
                _workers_running = false;
            }
            if (!_workers_running) {
                std::cout << "Returning from AsyncDiginorm worker..." << std::endl;
                return;
            }
        }
    }
}
