#ifndef KHMER_ASYNC_HH
#define KHMER_ASYNC_HH

#include "hashtable.hh"
#include <iostream>
#include <thread>
#include <boost/lockfree/queue.hpp>

typedef boost::lockfree::queue<HashIntoType> HashQueue;
typedef boost::lockfree::queue<const char *> CharQueue;
typedef boost::lockfree::queue<khmer::Read> ReadQueue;

typedef boost::lockfree::queue Queue;

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

        Queue<T> _in_queue; 

        Async() {
            _workers_running = false;
        }

        ~Async() {
            if(_workers_running) stop();
        }

        void start(unsigned int n_threads);
        void stop();

        void push(T item);
        virtual void consume(Queue<T>& q) = 0;
        
        void set_input(Queue<T>& new_q);
};

class AsyncWriter: public Async<HashIntoType> {
    
    friend class Hashtable;
    friend class AsyncHasher;
    
    protected:

        khmer::Hashtable * _ht;

    public:

        AsyncWriter (khmer::Hashtable * ht):
                     _ht(ht),
                     khmer::Async() {
        }

        void push(HashIntoType item);
        void consume(HashQueue& q);


};

class AsyncHasher: public Async<const char *> {

    friend class AsyncWriter;

    protected:
    
        unsigned int _ksize;

    public:

        HashQueue _out_queue;

        AsyncHasher (unsigned int ksize): 
                     _ksize(ksize),
                     khmer::Async() {
        }

        void push(const char * item);
        void consume(CharQueue& q);

        bool pop(HashIntoType& khash);
        void set_output(HashQueue& new_q);
}

class AsyncSequenceProcessor {

    private:

        khmer::Hashtable _ht;
        khmer::AsyncHashWriter _hasher;
        ReadQueue _in_queue;
        ReadQueue _out_queue;
        std::vector<std::thread> _processor_threads;
        bool _processors_running;

    public:

        AsyncSequenceProcessor (khmer::Hashtable * ht):
                                _ht(ht) {
            _processors_running = false;
        }

        ~AsyncSequenceProcessor() {
            if (processors_running) stop();
        }

        void start(unsigned int n_processor_threads);
        void stop();

        void push(khmer::Read read);
        bool pop(khmer::Read read);

        virtual void processor(ReadQueue q);
};

#endif
