#ifndef KHMER_ASYNC_HH
#define KHMER_ASYNC_HH

#include "khmer.hh"
#include "hashtable.hh"
#include "read_parsers.hh"
#include <iostream>
#include <thread>
#include <exception>
#include <stdexcept>
#include <mutex>
#include <boost/lockfree/queue.hpp>

#define VERBOSITY 0

using namespace boost::lockfree;

typedef khmer::read_parsers::Read Read;

typedef boost::lockfree::capacity<50000> Cap;

typedef queue<khmer::HashIntoType, Cap> HashQueue;
typedef queue<const char *, Cap> CharQueue;
typedef queue<Read*, Cap> ReadQueue;

namespace khmer {

class CountingHash;
class Hashtable;
class Hashbits;


// Keeps a vector of exceptions which are transferred from other threads
// This is accessed from either a Python-facing function or a main C thread,
// using check_and_rethrow, which throws the most recent exception.
// TODO: Define better behavior for multiple exceptions
class AsyncExceptionHandler {
    
    protected:

    std::vector<std::exception_ptr> exceptions;
    std::mutex exceptions_mutex;

    public:

    AsyncExceptionHandler() {};

    ~AsyncExceptionHandler() {
        std::lock_guard<std::mutex> lock(exceptions_mutex);
        exceptions.clear();
    }

    void push(std::exception_ptr exc) {
        std::lock_guard<std::mutex> lock(exceptions_mutex);
        exceptions.push_back(exc);
    }

    void reset() {
        std::lock_guard<std::mutex> lock(exceptions_mutex);
        exceptions.clear();
    }

    void check_and_rethrow() {
        std::lock_guard<std::mutex> lock(exceptions_mutex);
        if (exceptions.size() > 0) {
            auto exc = exceptions.back();
            exceptions.pop_back();
            std::rethrow_exception(exc);
        }
    }
};


// General deletion functor
template <class T> class Delete {
    void operator()(T ptr) { delete ptr; }    
};

// Destroy all objects on a queue
template <typename T> void flush_queue(queue<T, Cap>* q) {
    q->consume_all([](T ptr){delete ptr;});
}

// avoids needing to define two processors for paired
// and unpaired reads by pushing both types to the same
// input queue
class ReadBatch {
    
    protected:

        Read ** batch;
        bool _size;

    public:

        ReadBatch (Read * first, Read * second) {
            batch = new Read*[2];
            batch[0] = first;
            batch[1] = second;
            _size = true;
        }
        ReadBatch (Read * read) {
            batch = new Read*[1];
            batch[0] = read;
            _size = false;
        }

        ~ReadBatch () {
            delete batch[0];
            if (_size) delete batch[1];
            delete batch;
        }

        Read * first() {
            return batch[0];
        }

        Read * second() {
            return batch[_size];
        }

        unsigned int size() {
            return (unsigned int) _size + 1;
        }
};

// Functor to call count() on the give sequence's k-mers
// Used to either consume all sequences on a queue (with consume_all) or
// to consume_one instead of pop
class CountFunc {
    public:

    khmer::Hashtable * _ht;
    unsigned int _ksize;

    CountFunc(khmer::Hashtable * ht, unsigned int ksize): _ht(ht), _ksize(ksize) {}

    void operator()(const char * sequence) {
        HashIntoType kmer;
        khmer::KMerIterator kmers(sequence, _ksize);
        try {
            while(!kmers.done()) {
                kmer = kmers.next();
                _ht->count(kmer);
            }
        } catch (khmer_exception &e) {
            std::cout << e.what() << " " << std::endl;
            std::cout << "ERROR in AsyncSequenceWriter: " << sequence << std::endl;
        }
    }
};

// Generic Async class
// Manage thread state such as exceptions
// This should be the base for all async operations
class Async {

    protected:
   
        uint32_t _stdout_spin_lock; 
        unsigned int _n_workers;
        std::vector<std::thread> _worker_threads;
        bool _workers_running;
        AsyncExceptionHandler * _exc_handler;
        unsigned int _batchsize;

    public:

        Async() {
            //_exc_handler();
            _stdout_spin_lock = 0;
            _workers_running = false;
            _batchsize = 1;
            _exc_handler = new AsyncExceptionHandler();
        }

        ~Async() {
            if(_workers_running) stop();
        }

        virtual void consume() = 0;

        void start(unsigned int n_threads) {
            _n_workers = n_threads;
            _workers_running = true;
            for (unsigned int t=0; t<_n_workers; ++t) {
                #if(VERBOSITY)
                std::cout << "Async spawn worker" << std::endl;
                #endif
                _worker_threads.push_back(std::thread(&Async::consume, this));
            }
        }

        void stop() {
            _workers_running = false;
            auto wthread = _worker_threads.begin();
            while (wthread != _worker_threads.end()) {
                if(wthread->joinable()) wthread->join();
                wthread++;
            }
        }

        bool workers_running() {
            return _workers_running;
        }

        // Define a lock for writing to standard out
        void lock_stdout() {
            while(!__sync_bool_compare_and_swap( &_stdout_spin_lock, 0, 1 ));
        }

        void unlock_stdout() {
            __sync_bool_compare_and_swap( &_stdout_spin_lock, 1, 0 );
        }

        void check_and_rethrow() {
            _exc_handler->check_and_rethrow();
        }

        // Register a new AsyncExceptionHandler
        void register_exception_handler(AsyncExceptionHandler * &exc_handler) {
            _exc_handler = exc_handler;
        }
};

// Consumer model: draws items off an input queue and does stuff
template <class T> class AsyncConsumer: public virtual Async {

    protected:

        queue<T, Cap> * _in_queue;
        unsigned int _n_pushed;

    public:
        AsyncConsumer(): 
            Async() {
            _in_queue = new queue<T, Cap>();
        }

        void start(unsigned int n_threads) {
            _n_pushed = 0;
            Async::start(n_threads);
        }

        bool push(T& item) {
            if( _in_queue->bounded_push(item)) {
                __sync_add_and_fetch(&_n_pushed, _batchsize);
                return true;
            }
            return false;
        }

        unsigned int n_pushed() {
            return _n_pushed;
        }   

        void set_input(queue<T, Cap>* new_q) {
            _in_queue = new_q;
        }
};

// Producer model: draw stuff from <something> and push it to an output queue
template <class T> class AsyncProducer: public virtual Async {

    protected:
    
        queue<T, Cap> * _out_queue;
        unsigned int _n_popped;

    public:
        AsyncProducer(): 
            Async() {
            _out_queue = new queue<T, Cap>();
        }

        void start(unsigned int n_threads) {
            _n_popped = 0;
            Async::start(n_threads);
        }

        bool pop(T& item) {
            if(_out_queue->pop(item)) {
                __sync_fetch_and_add(&_n_popped, _batchsize);
                return true;
            }
            return false;
        }

        unsigned int n_popped() {
            return _n_popped;
        }

        void set_output(queue<T, Cap>* new_q) {
            _out_queue = new_q;
        }

        bool has_output() {
            return !_out_queue->empty();
        }

};

// ConsumerProducer model: draw from an input queue, do stuff
// push to an output queue
template <class T, class V> class AsyncConsumerProducer:
    public AsyncConsumer<T>, public AsyncProducer<V> {

    public:
        AsyncConsumerProducer<T,V>() :
            Async(),
            AsyncConsumer<T>(),
            AsyncProducer<V>()
        {
        }

        void start(unsigned int n_threads) {
            AsyncProducer<V>::_n_popped = 0;
            AsyncConsumer<T>::_n_pushed = 0;
            Async::start(n_threads);
        }
};


class AsyncHashWriter: public AsyncConsumer<HashIntoType> {
    
    protected:

        khmer::Hashtable * _ht;
        unsigned int _n_written;

    public:

        AsyncHashWriter (khmer::Hashtable * ht):
                     khmer::AsyncConsumer<HashIntoType>(), 
                     _ht(ht) {
        }

        unsigned int ksize();
        void start();
        virtual void consume();
        unsigned int n_written();


};

class AsyncSequenceWriter: public AsyncConsumer<const char *> {
    
    protected:

        khmer::Hashtable * _ht;
        unsigned int _n_written;
        unsigned int _ksize;

    public:

        AsyncSequenceWriter (khmer::Hashtable * ht):
                     khmer::AsyncConsumer<const char *>(), 
                     _ht(ht) {
        }

        unsigned int ksize();
        void start();
        virtual void consume();
        unsigned int n_written();


};

/*
class AsyncSequenceParser: public Async<Read *> {
    
    protected:

        khmer::Hashtable * _ht;
        unsigned int _n_written;
        unsigned int _n_pushed;
        unsigned int _ksize;

    public:

        AsyncSequenceParser ():
                     khmer::Async<Read *>()
        {
        }

        void start();
        virtual void consume(ReadQueue * q);
        bool push(Read * &read);
        unsigned int n_pushed();
        unsigned int n_written();


};
*/

class AsyncHasher: public AsyncConsumerProducer<const char *, HashIntoType> {

    protected:
    
        unsigned int _ksize;

    public:

        AsyncHasher (unsigned int ksize):
                     khmer::AsyncConsumerProducer<const char *, HashIntoType>(),
                     _ksize(ksize) {
        }

        virtual void consume();
};

// General read processing model
// TODO: Split off the read_iparser functionality into its own AsyncParser class
class AsyncSequenceProcessor: public AsyncConsumerProducer<ReadBatch*, ReadBatch*> {

    protected:

        khmer::Hashtable * _ht;
        khmer::AsyncSequenceWriter * _writer;

        std::thread * _reader_thread;
        bool paired = false;

        unsigned int _n_parsed;
        unsigned int _n_processed;

        bool _parsing_reads;
        bool _processing_reads;

    public:

        AsyncSequenceProcessor (khmer::Hashtable * ht):
                                khmer::AsyncConsumerProducer<ReadBatch*,ReadBatch*>(),
                                _ht(ht) {
            _writer = new AsyncSequenceWriter(_ht);
        }

        void start(const std::string &filename,
                    bool paired,
                    unsigned int n_threads);
        void stop();

        virtual void consume() = 0;
        unsigned int n_processed();

        virtual bool iter_stop() = 0;

        bool is_parsing();
        void read_iparser(const std::string &filename);
        bool is_paired();
        unsigned int n_parsed();
        unsigned int n_written();
};

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

        bool filter_single(Read * read);
        bool filter_paired(ReadBatch * read);
};
};
#endif
