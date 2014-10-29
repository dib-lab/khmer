#include "hashtable.hh"
#include "async_hash.hh"
#include "cqueue.hh"

namespace khmer {

void AsyncHashWriter::start() {
    //std::cout << "Thread " << std::this_thread::get_id() << " spawning hash consumer" << std::endl;
    consume_hash_thread = new std::thread(&AsyncHashWriter::consume_hash, this, std::ref(hash_q));
}

void AsyncHashWriter::stop() {
    hash_q.push(_stop_sig);
    consume_hash_thread->join();
    //std::cout << "Thread " << std::this_thread::get_id() << " destroyed hash consumer" << std::endl;
}

void AsyncHashWriter::enqueue(HashIntoType khash) {
    hash_q.push(khash);
}

void AsyncHashWriter::consume_hash(CQueue<HashIntoType>& q) {
    HashIntoType khash;
    unsigned long long total_consumed = 0;
    //std::cout << "Thread " << std::this_thread::get_id() << " waiting on hashes" << std::endl;
    while(1) {
        khash = q.pop();
        if (khash == _stop_sig) {
            return;
        } else {
            _ht->count(khash);
            total_consumed++;
        }
        //if (total_consumed % 100000 == 0) {
        //    std::cout << "async consumed " << total_consumed << std::endl;
        //}
    }
}
}
