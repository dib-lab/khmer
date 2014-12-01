#ifndef KHMER_ASYNC_WRITERS_HH
#define KHMER_ASYNC_WRITERS_HH

#include "khmer_async.hh"
#include "async_models.hh"

namespace khmer {

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
class AsyncSequenceParser: public Async<ReadPtr> {
    
    protected:

        khmer::Hashtable * _ht;
        unsigned int _n_written;
        unsigned int _n_pushed;
        unsigned int _ksize;

    public:

        AsyncSequenceParser ():
                     khmer::Async<ReadPtr>()
        {
        }

        void start();
        virtual void consume(ReadQueue * q);
        bool push(ReadPtr &read);
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

}; //namespace khmer

#endif
