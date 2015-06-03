#include <iostream>
#include <seqan/parallel.h>

using namespace seqan;

int main()
{
    // instantiate an empty Queue
    ConcurrentQueue<unsigned> queue;

    // start two threads
    SEQAN_OMP_PRAGMA(sections num_threads(2))
    {
        SEQAN_OMP_PRAGMA(section)
        {
            for (int i = 9999; i != 0; --i)
                appendValue(queue, i);
        }

        SEQAN_OMP_PRAGMA(section)
        {
            bool equal = true;
            for (int i = 9999; i != 0; --i)
                equal &= (i == popFront(queue));
            std::cout << (equal ? "SUCCESS" : "FAILURE") << std::endl;
        }
    }

    return 0;
}
