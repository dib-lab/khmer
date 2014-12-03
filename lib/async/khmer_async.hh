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
#include <time.h>
#include <boost/lockfree/queue.hpp>

#define VERBOSITY 0

using namespace boost::lockfree;

namespace khmer {

static inline timespec timediff(timespec start, timespec end)
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

#define TIMING 1
#if TIMING
#define TSTART() clock_gettime(CLOCK_MONOTONIC, &start_t);
#define TEND(dest) clock_gettime(CLOCK_MONOTONIC, &end_t); \
                   dest+= (timediff(start_t, end_t).tv_nsec / 1000000000.0);
#else
#define TSTART()
#define TEND(dest)
#endif

typedef khmer::read_parsers::Read Read;

typedef std::shared_ptr<Read> ReadPtr;

typedef capacity<50000> Cap;

typedef queue<khmer::HashIntoType, Cap> HashQueue;
typedef queue<const char *, Cap> CharQueue;
typedef queue<ReadPtr, Cap> ReadQueue;

class Hashtable;

// General deletion functor
template <class T> class Delete {
    void operator()(T ptr) { delete ptr; }    
};

// Destroy all objects on a queue
template <typename T> void flush_queue(queue<T, Cap>* q) {
    q->consume_all([](T ptr){delete ptr;});
}

inline char * copy_seq(ReadPtr read) {
    char * buffer = new char [read->sequence.size() + 1];
    memcpy(buffer, read->sequence.c_str(), read->sequence.size() + 1);
    return buffer;
}

// avoids needing to define two processors for paired
// and unpaired reads by pushing both types to the same
// input queue
class ReadBatch {
    
    protected:

        ReadPtr * batch;
        bool _size;

    public:

        ReadBatch (ReadPtr first, ReadPtr second) {
            batch = new ReadPtr[2];
            batch[0] = first;
            batch[1] = second;
            _size = true;
        }

        ReadBatch (ReadPtr read) {
            batch = new ReadPtr[1];
            batch[0] = read;
            _size = false;
        }

        ~ReadBatch () {
            delete [] batch;
        }

        ReadPtr first() {
            return batch[0];
        }

        ReadPtr second() {
            return batch[_size];
        }

        unsigned int size() {
            return (unsigned int) _size + 1;
        }
};

typedef ReadBatch * ReadBatchPtr;

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


};
#endif
