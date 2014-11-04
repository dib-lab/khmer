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
    Async<HashIntoType>::start(1);
}

unsigned int AsyncWriter::ksize() {
    return _ht->ksize();
}

void AsyncWriter::consume(HashQueue * q) {
    HashIntoType khash;
    unsigned long long total_consumed = 0;
    //std::cout << "Thread " << std::this_thread::get_id() << " waiting on hashes" << std::endl;
    while(1) {
        if (q->pop(khash)) {
            _ht->count(khash);
            total_consumed++;
        } else {
            if(!_workers_running) return;
        }
    }
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

void AsyncSequenceProcessor::start(unsigned int n_threads) {
    _writer->start();
    Async<Read*>::start(n_threads);
}

void AsyncSequenceProcessor::stop() {
    // First stop the sequence processor, flushing its queue
    Async<Read*>::stop();
    // Then stop the writer, writing out all the flushed hashes
    _writer->stop();
}

bool AsyncSequenceProcessor::pop(Read* read) {
    return _out_queue->pop(read);
}

void AsyncSequenceProcessor::set_output(ReadQueue* new_q) {
    _out_queue = new_q;
}

/////
//
// AsyncDiginorm
//
/////

void AsyncDiginorm::read_iparser(const std::string &filename) {
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
}

void AsyncDiginorm::start(const std::string &filename,
                            unsigned int cutoff,
                            unsigned int n_threads) {
    _cutoff = cutoff;

    _parsed_count = 0;
    _parsing_reads = true;
    reader_thread = new std::thread(&AsyncDiginorm::read_iparser, this, filename);

    _processed_count = 0;
    _n_kept = 0;
    _processing_reads = true;
    AsyncSequenceProcessor::start(n_threads);
    _n_popped = 0;
}

bool AsyncDiginorm::is_parsing() {
    return _parsing_reads;
}

bool AsyncDiginorm::is_processing() {
    return _processing_reads;
}

unsigned int AsyncDiginorm::reads_kept() {
    return _n_kept;
}

bool AsyncDiginorm::has_output() {
    return !_out_queue->empty();
}

bool AsyncDiginorm::pop(Read* &read) {
    if(_out_queue->pop(read)) {
        //std::cout << "POP " << read->sequence << std::endl;
        _n_popped++;
        return true;
    }
    return false;
}

unsigned int AsyncDiginorm::reads_popped() {
    return _n_popped;
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
                //if (_n_kept % 1000 == 0) std::cout << tid << " kept " << read->sequence << std::endl;
                khmer::KMerIterator kmers(read->sequence.c_str(), ksize);
                while(!kmers.done()) {
                    khash = kmers.next();
                    _writer->push(khash);
                }
            }
            __sync_fetch_and_add(&_processed_count, 1);
        } else {
            if (!is_parsing() && _processed_count >= _parsed_count) {
                std::cout << "Finished processing..." << std::endl;
                _processing_reads = false;
                _workers_running = false;
                return;
            }
            if (!_workers_running) return;
        }
    }
}
