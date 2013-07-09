/*#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdlib.h>             //used for c_str() function
#include <ctime>
#include <stdlib.h>
#include <time.h>   //both used to make sure we do have differnt random numbers
#include <set>
#include <vector>
*/
#include <time.h>
#include "outline_index.hh"

using namespace khmer;

int main(int argc,char *argv[])
{
    std::string  tagedKhmerFileName="", readsFileName="", queryFileName="";

    int optind=1,khmerSize=0,maxReadSize=0;

    while ((optind < argc) && (argv[optind][0]=='-')) {
        std::string sw = argv[optind];
        if (sw=="-tf") {
            optind++;
            tagedKhmerFileName=argv[optind];
        } else if (sw=="-rf") {
            optind++;
            readsFileName=argv[optind];
        } else if (sw=="-qf") {
            optind++;
            queryFileName = argv[optind];
        } else if(sw=="-ks") {
            optind++;
            khmerSize=atoi(argv[optind]);
        } else if(sw=="-rs") {
            optind++;
            maxReadSize=atoi(argv[optind]);
        }

        else if (sw=="-help") {
            std::cout<<"the following are the prgram arguments\n";
            std::cout<<"\t-tf  TagedKhmerFileName\n";
            std::cout<<"\t-rf  ReadsFileName\n";
            std::cout<<"\t-qf  QueryFileName\n";
            std::cout<<"\t-ks  Khmer Size\n";
            std::cout<<"\t-rs  Max Read Size\n";
            std::cout<<"\t-help disply this help\n";
            return 0;
        } else {
            std::cout << "Unknown switch: " << argv[optind] << std::endl;
        }

        optind++;
    }

    
    //creaate reads binary file
    std::cout<<"\ncreating reads binary file ....\n";
    std::string readsFileNameBin=readsFileName+".bin";
    convertFastaToBin(readsFileName,readsFileNameBin);
    std::cout<<"done...\n";
    
    //intilize the classes
    std::cout<<"\nload tagged khmers and sort them ...\n";
    std::vector<khmer::HashIntoType> sortedKhmerVector;
    unsigned int save_ksize;
    load_tagset(tagedKhmerFileName,sortedKhmerVector,save_ksize);
    std::cout<<"done..., the num of tags are:"<<sortedKhmerVector.size()<<std::endl;

    //build the index
    std::cout<<"\nbuild the index ...\n";
    //start the timer
    time_t start,end;
    double diff;
    time (&start);
    build_index(readsFileNameBin,sortedKhmerVector,save_ksize);
    time (&end);
    diff = difftime (end,start);
    printf ("It took you %.2lf seconds to build the tree.\n", diff );
    std::cout<<"done...\n";
    
    /*
    //query example
    std::cout<<"\nexact query search ...\n";
    //exactQuery(readsFileNameBin,queryFileName);

    //open the query file and save the tagged khmers in a vector
    std::cout<<"\tin Load_Queries...\n";
    std::fstream inQfile(queryFileName.c_str(),std::ios::in| std::ios::binary);
    assert(inQfile.is_open());
    khmer::HashIntoType khmer;
    std::vector<khmer::HashIntoType> qeuery_tagged_khmer;
    while (!(inQfile.eof())) {
        inQfile>>khmer;
        //std::cout<<khmer<<" ";
        qeuery_tagged_khmer.push_back(khmer);
    }

    //retrive the read ids that contains at lest one of the query tagged khmer
    std::vector<long> reads_ids;
    std::string indexfilename;
    indexfilename=readsFileNameBin+".index";
    retrieve_read_ids_by_tag(indexfilename,qeuery_tagged_khmer,reads_ids );
    //output the rertive id result
    std::cout<<"\toutput the rertive id result\n";
    for (int i=0; i< reads_ids.size(); i++) {
        std::cout<<reads_ids[i]<<" ";
    }
    std::cout<<std::endl;

    //retrieve reads by read ids
    std::vector<std::string> reads;
    retrieve_read_by_id(readsFileNameBin,reads_ids,reads);
    //output the retreieved reads
    std::cout<<"\toutput the retreieved reads\n";
    for (int i=0; i< reads.size(); i++) {
        std::cout<<reads[i]<<std::endl;
    }
    */
    return 0;
}
