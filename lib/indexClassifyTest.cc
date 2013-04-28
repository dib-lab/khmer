#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdlib.h>             //used for c_str() function
#include <ctime>
#include <stdlib.h>
//#include <time.h>   //both used to make sure we do have differnt random numbers
#include <set>
#include <vector>

#include "ktable.hh"
#include "khmer.hh"
#include "hashbits.hh"
#include "minilib.hh"

using namespace std;

int main(int argc,char *argv[])
{
	string  tagedKhmerFileName="", readsFileName="", queryFileName="";

	int optind=1,khmerSize=0,maxReadSize=0;
        
	while ((optind < argc) && (argv[optind][0]=='-')) {
	string sw = argv[optind];
	if (sw=="-tf") {
            	optind++;
            	tagedKhmerFileName=argv[optind];
	 }
        else if (sw=="-rf") {
  		optind++;
            	readsFileName=argv[optind];
        }
        else if (sw=="-qf") {
            	optind++;
	        queryFileName = argv[optind];
        }
	else if(sw=="-ks") {
                optind++;
                khmerSize=atoi(argv[optind]);
         }
	 else if(sw=="-rs") {
                optind++;
                maxReadSize=atoi(argv[optind]);
         }

        else if (sw=="-help") {
                cout<<"the following are the prgram arguments\n";
                cout<<"\t-tf  TagedKhmerFileName\n";
                cout<<"\t-rf  ReadsFileName\n";
                cout<<"\t-qf  QueryFileName\n";
		cout<<"\t-ks  Khmer Size\n";
		cout<<"\t-rs  Max Read Size\n";
                cout<<"\t-help disply this help\n";
		return 0;
        }
        else { cout << "Unknown switch: " << argv[optind] << endl;}
        
	optind++;	
	}
	
	//creaate reads binary file
	std::cout<<"\ncreating reads binary file ....\n";
	std::string readsFileNameBin=readsFileName+".bin";
	convertFastaToBin(readsFileName,readsFileNameBin);
	
	//intilize the classes
	std::cout<<"\nintilizaing classes ...\n";
	vector<khmer::HashIntoType> sortedKhmerVector;
	unsigned int save_ksize;
	//load_tagset(tagedKhmerFileName,sortedKhmerVector,save_ksize);
	loading_tagset(tagedKhmerFileName,sortedKhmerVector,save_ksize);
	
	//classfy the reads
	std::cout<<"\nclassify the reads into the classes ...\n";
	classifyReads(readsFileNameBin,sortedKhmerVector,save_ksize);
	
	//query example
	std::cout<<"\nexact query search ...\n";
	exactQuery(readsFileNameBin,queryFileName);
	
	return 1;
}

