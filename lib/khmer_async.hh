#ifndef KHMER_ASYNC_HH
#define KHMER_ASYNC_HH

#include "khmer.hh"
#include "hashtable.hh"
#include "read_parsers.hh"
#include <iostream>
#include <thread>
#include <boost/lockfree/queue.hpp>

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

template <class T> class Async {

    protected:
    
        unsigned int _n_workers;
        std::vector<std::thread> _worker_threads;
        bool _workers_running;

    public:

        queue<T, Cap> * _in_queue; 

        Async() {
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
                std::cout << "Async spawn worker" << std::endl;
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
            return _in_queue->push(item);
        }

        void set_input(queue<T, Cap>* new_q) {
            _in_queue = new_q;
        }
};

class AsyncWriter: public Async<HashIntoType> {
    
    friend class Hashtable;
    friend class AsyncHasher;
    
    protected:

        khmer::Hashtable * _ht;
        unsigned int _n_written;
        unsigned int _n_pushed;

    public:

        AsyncWriter (khmer::Hashtable * ht):
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

class AsyncHasher: public Async<const char *> {

    friend class AsyncWriter;

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
        khmer::AsyncWriter * _writer;


    public:

        ReadQueue * _out_queue;

        AsyncSequenceProcessor (khmer::Hashtable * ht):
                                khmer::Async<Read*>(),
                                _ht(ht) {
            _writer = new AsyncWriter(_ht);
            _out_queue = new ReadQueue();
        }

        void start(unsigned int n_threads);
        void stop();

        bool pop(Read* read);
        virtual void consume(ReadQueue* q) = 0;
        void set_output(ReadQueue* new_q);
};


class AsyncDiginorm: public AsyncSequenceProcessor {

    protected:

        bool _parsing_reads;
        bool _processing_reads;
        unsigned int _cutoff;
        std::thread * reader_thread;

        unsigned int _parsed_count;
        unsigned int _processed_count;
        unsigned int _n_kept;
        unsigned int _n_popped;

        unsigned int _n_hashes_pushed;

    public:

        AsyncDiginorm (khmer::Hashtable * ht):
                        khmer::AsyncSequenceProcessor(ht) {
        }
        void start(const std::string &filename,
                    unsigned int cutoff,
                    unsigned int n_threads);

        bool pop(Read * &read);
        bool has_output();

        bool is_processing();
        bool is_parsing();

        unsigned int reads_kept();
        unsigned int reads_popped();

        void read_iparser(const std::string &filename);
        virtual void consume(ReadQueue* q);
};
};
#endif
