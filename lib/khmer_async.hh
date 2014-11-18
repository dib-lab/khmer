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

// Subclass read_parsers::Read for paired reads;
// avoids needing to define two processors for paired
// and unpaired reads by pushing both types to the same
// input queue
class PairedRead : public Read {
    public:

    Read * second;
    PairedRead (Read * first, Read * second) {
        name = first->name;
        annotations = first->annotations;
        sequence = first->sequence;
        accuracy = first->accuracy;
        this->second = second;
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

// General Async class
template <class T> class Async {

    protected:
   
        uint32_t _stdout_spin_lock; 
        unsigned int _n_workers;
        std::vector<std::thread> _worker_threads;
        bool _workers_running;
        AsyncExceptionHandler _exc_handler;

    public:

        queue<T, Cap> * _in_queue; 

        Async() {
            //_exc_handler();
            _stdout_spin_lock = 0;
            _workers_running = false;
            _in_queue = new queue<T, Cap>();
        }

        ~Async() {
            if(_workers_running) stop();
        }

        virtual void consume(queue<T, Cap>* q) = 0;

        void start(unsigned int n_threads) {
            _n_workers = n_threads;

            _workers_running = true;
            for (unsigned int t=0; t<_n_workers; ++t) {
                #if(VERBOSITY)
                std::cout << "Async spawn worker" << std::endl;
                #endif
                _worker_threads.push_back(std::thread(&Async<T>::consume, this, _in_queue));
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

        bool push(T& item) {
            return _in_queue->bounded_push(item);
        }

        void set_input(queue<T, Cap>* new_q) {
            _in_queue = new_q;
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
            _exc_handler.check_and_rethrow();
        }

        // Register a new AsyncExceptionHandler
        void register_exception_handler(AsyncExceptionHandler &exc_handler) {
            _exc_handler = exc_handler;
        }
};

class AsyncHashWriter: public Async<HashIntoType> {
    
    friend class Hashtable;
    friend class AsyncHasher;
    
    protected:

        khmer::Hashtable * _ht;
        unsigned int _n_written;
        unsigned int _n_pushed;

    public:

        AsyncHashWriter (khmer::Hashtable * ht):
                     khmer::Async<HashIntoType>(), 
                     _ht(ht) {
        }

        unsigned int ksize();
        void start();
        virtual void consume(HashQueue * q);
        bool push(HashIntoType &khash);
        unsigned int n_pushed();
        unsigned int n_written();


};

class AsyncSequenceWriter: public Async<const char *> {
    
    friend class AsyncHasher;
    
    protected:

        khmer::Hashtable * _ht;
        unsigned int _n_written;
        unsigned int _n_pushed;
        unsigned int _ksize;

    public:

        AsyncSequenceWriter (khmer::Hashtable * ht):
                     khmer::Async<const char *>(), 
                     _ht(ht) {
        }

        unsigned int ksize();
        void start();
        virtual void consume(CharQueue * q);
        bool push(const char * &sequence);
        unsigned int n_pushed();
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

class AsyncHasher: public Async<const char *> {

    friend class AsyncHashWriter;

    protected:
    
        unsigned int _ksize;

    public:

        HashQueue * _out_queue;

        AsyncHasher (unsigned int ksize):
                     khmer::Async<const char *>(),
                     _ksize(ksize) {
            _out_queue = new HashQueue();
        }

        virtual void consume(CharQueue * q);

        bool pop(HashIntoType& khash);
        void set_output(HashQueue* new_q);
};

class AsyncSequenceProcessor: public Async<Read*> {

    protected:

        khmer::Hashtable * _ht;
        khmer::AsyncSequenceWriter * _writer;

        std::thread * _reader_thread;
        bool paired = false;

        unsigned int _n_parsed;
        unsigned int _n_processed;
        // Number popped from output queue
        unsigned int _n_popped;

        bool _parsing_reads;
        bool _processing_reads;

    public:

        ReadQueue * _out_queue;

        AsyncSequenceProcessor (khmer::Hashtable * ht):
                                khmer::Async<Read*>(),
                                _ht(ht) {
            _writer = new AsyncSequenceWriter(_ht);
            _out_queue = new ReadQueue();
        }

        void start(const std::string &filename,
                    bool paired,
                    unsigned int n_threads);
        void stop();

        virtual void consume(ReadQueue* q) = 0;
        unsigned int n_processed();

        bool pop(Read * &read);
        unsigned int n_popped();
        bool has_output();
        void set_output(ReadQueue* new_q);
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
        virtual void consume(ReadQueue* q);
        bool iter_stop();
};
};
#endif
