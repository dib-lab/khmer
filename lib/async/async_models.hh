#ifndef KHMER_ASYNC_MODELS_HH
#define KHMER_ASYNC_MODELS_HH

#include "khmer.hh"
#include "khmer_async.hh"

namespace khmer {

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

/*
struct NoMoreQueuesAvailable : public khmer_exception {
};

template <class T> class AsyncRRProducer: public virtual Async {

    protected:
        std::vector<queue<T, Cap> * > _out_queues;
        bool _open_slots;
        uint32_t _acquire_q_spinlock;
    
    public:
        AsyncRRProducer():
            Async() {
            _open_slots = 0;
            _acquire_q_spinlock = 0;
        }

        void start(unsigned int n_threads) {
            _out_queues.clear();
            while(!__sync_bool_compare_and_swap( &_acquire_q_spinlock, 0, 1));
            for (unsigned int q=0; q<n_threads; ++g) {
                _out_queues.push_back(new std::vector<T, Cap>());
            }
            _open_slots = n_threads;
            __sync_bool_compare_and_swap( &_acquire_q_spinlock, 1, 0);
        }

        std::vector<T, Cap> * acquire_queue() {
            while(!__sync_bool_compare_and_swap ( &_acquire_q_spinlock, 0 1));
            if (_open_slots < 1) {
                throw NoMoreQueuesAvailable("No more queues to acquire");
            }
            _open_slots--;
            std::vector<T, Cap> * q = _out_queues[_open_slots];
            __sync_bool_compare_compare_and_swap( &_acquire_q_spinlock, 1, 0);
            return q;
        }
*/
        


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

        bool push(T item) {
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

}; //namespace khmer

#endif
