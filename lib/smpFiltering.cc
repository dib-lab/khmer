/* This file is a simple serial implementation of the folloing algorithm:

for read in file:
    if median(read) > C:
        discard();
    else:
        count(read);
        save_to_file(read);
*/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "threadedParsers.hh"
#include "counting.hh"
#include "primes.hh"
#include <math.h>

using namespace std;
using namespace khmer;
using namespace khmer:: threaded_parsers;

#define MAX_MEDIAN_COUNT 10 /* Only count if the median is < this */
#define TABLE_SIZE 2000000000

/* Forward Declarations */

int main(int argc, char **argv)
{
    /* Check commandline arguments */
    if (argc != 4)
    {
        cerr << "You must specify an input and output filename and the k-size" << endl;
        cerr << "\tUsage: " << argv[0] << " input filename output filename k-size" << endl;
        exit(-1);
    }

    unsigned int K;
    if (sscanf(argv[3], "%d", &K) != 1)
    {
        cerr << "Could not interpret K size of " << argv[3] << endl;
        exit(-1);
    }
    

    /* Initialize the counting hash */
    vector<HashIntoType> tableSizes;
    Primes pri(TABLE_SIZE);
    tableSizes.push_back(pri.get_next_prime());
    tableSizes.push_back(pri.get_next_prime());
    tableSizes.push_back(pri.get_next_prime());
    CountingHash h(K, tableSizes);

    /* Start parsing the fasta file */
    ThreadedIParserFactory* pf = ThreadedIParserFactory::get_parser(argv[1], 104857600);
    long long unsigned int totalCount = 0, keptCount = 0;

    #pragma omp parallel reduction(+:totalCount,keptCount)
    while (!pf->is_complete())
    {
        ThreadedIParser* p = pf->get_next_parser();

        /* Open the output file */
        ofstream outFile;
        string outputFileName(argv[2]);
        char numstr[21]; // enough to hold all numbers up to 64-bits
        sprintf(numstr, "%020lu", p->getEndPos());
        outputFileName += "_";
        outputFileName += numstr;
        outFile.open(outputFileName.c_str());
        if (! outFile.is_open())
        {
            cerr << "Failed to open file " << outputFileName;
            perror("");
            cerr << endl;
            exit(-1);
        }

        while (!p->is_complete())
        {
            Read r;
            BoundedCounterType medCount;
            float meanCount;
            float stdDev;
            r = p->get_next_read();
            totalCount++;
            h.get_median_count(r.seq, medCount, meanCount, stdDev);
            if (medCount < MAX_MEDIAN_COUNT)
            {
                keptCount++;
                /* Count it */
                KMerIterator kmers(r.seq.c_str(), K);
                while(!kmers.done())
                {
                    h.count(kmers.next());
                }

                outFile << ">" << r.name << endl;
                outFile << r.seq << endl;
            }
        }
        outFile.close();
    }


    printf("Total Count: %llu, Kept Count: %llu\n", totalCount, keptCount);
    printf("Hashtable Occupancy: %lf\n", (double)h.n_occupied() / TABLE_SIZE);
}
