//
// Copyright (c) 2013 Juan Palacios juan.palacios.puyana@gmail.com
// Subject to the BSD 2-Clause License
// - see < http://opensource.org/licenses/BSD-2-Clause>
//

// Code borrowed from:
// https://github.com/juanchopanza/cppblog/tree/master/Concurrency/Queue

#ifndef CONCURRENT_QUEUE_
#define CONCURRENT_QUEUE_

#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>

namespace khmer {

template <typename T>
class CQueue
{
 public:

  T pop() 
  {
    std::unique_lock<std::mutex> mlock(mutex_);
    while (queue_.empty())
    {
      cond_.wait(mlock);
    }
    auto val = queue_.front();
    queue_.pop();
    return val;
  }

  void pop(T& item)
  {
    std::unique_lock<std::mutex> mlock(mutex_);
    while (queue_.empty())
    {
      cond_.wait(mlock);
    }
    item = queue_.front();
    queue_.pop();
  }

  void push(const T& item)
  {
    std::unique_lock<std::mutex> mlock(mutex_);
    queue_.push(item);
    mlock.unlock();
    cond_.notify_one();
  }
  CQueue()=default;
  CQueue(const CQueue&) = delete;            // disable copying
  CQueue& operator=(const CQueue&) = delete; // disable assignment
  
 private:
  std::queue<T> queue_;
  std::mutex mutex_;
  std::condition_variable cond_;
};
}
#endif
