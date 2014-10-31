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

typedef queue<khmer::HashIntoType> HashQueue;
typedef queue<const char *> CharQueue;
typedef queue<Read*> ReadQueue;

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

        queue<T> * _in_queue; 

        Async() {
            _workers_running = false;
            _in_queue = new queue<T>;
        }

        ~Async() {
            if(_workers_running) stop();
        }

        virtual void consume(queue<T>* q) = 0;

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
                wthread->join();
                wthread++;
            }
        }

        void push(T& item) {
            _in_queue->push(item);
        }

        void set_input(queue<T>* new_q) {
            _in_queue = new_q;
        }
};

class AsyncWriter: public Async<HashIntoType> {
    
    friend class Hashtable;
    friend class AsyncHasher;
    
    protected:

        khmer::Hashtable * _ht;

    public:

        AsyncWriter (khmer::Hashtable * ht):
                     khmer::Async<HashIntoType>(), 
                     _ht(ht) {
        }

        void start();
        virtual void consume(HashQueue * q);


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

        unsigned int _cutoff;

    public:

        void start(unsigned int cutoff, unsigned int n_threads);
        virtual void consume(ReadQueue* q);
};
};
#endif
