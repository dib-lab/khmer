/*#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdlib.h>             //used for c_str() function
#include <ctime>
#include <stdlib.h>
//#include <time.h>   //both used to make sure we do have differnt random numbers
#include <set>
#include <vector>
*/
#include "indexClassify.hh"

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
                std::cout<<"the following are the prgram arguments\n";
                std::cout<<"\t-tf  TagedKhmerFileName\n";
                std::cout<<"\t-rf  ReadsFileName\n";
                std::cout<<"\t-qf  QueryFileName\n";
		std::cout<<"\t-ks  Khmer Size\n";
		std::cout<<"\t-rs  Max Read Size\n";
                std::cout<<"\t-help disply this help\n";
		return 0;
        }
        else { std::cout << "Unknown switch: " << argv[optind] << std::endl;}
        
	optind++;	
	}
	
/*
  std::map<char,int> mymap;

  mymap.insert ( std::pair<char,int>('a',100) );
  mymap.insert ( std::pair<char,int>('z',200) );

  std::pair<std::map<char,int>::iterator,bool> ret;
  ret = mymap.insert ( std::pair<char,int>('z',500) );
  if (ret.second==false) {
    std::cout << "element 'z' already existed";
    std::cout << " with a value of " << ret.first->second << '\n';
  }

  std::map<char,int>::iterator it = mymap.begin();
  mymap.insert (it, std::pair<char,int>('b',300));
  mymap.insert (it, std::pair<char,int>('c',400));
  std::map<char,int> anothermap;
  anothermap.insert(mymap.begin(),mymap.find('c'));
  std::cout << "mymap contains:\n";
  for (it=mymap.begin(); it!=mymap.end(); ++it)
  std::cout << it->first << " => " << it->second << '\n';

  std::cout << "anothermap contains:\n";
  for (it=anothermap.begin(); it!=anothermap.end(); ++it)
  std::cout << it->first << " => " << it->second << '\n';
  */

	
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
	build_index(readsFileNameBin,sortedKhmerVector,save_ksize);
	std::cout<<"done...\n";

	//query example
	std::cout<<"\nexact query search ...\n";
	//exactQuery(readsFileNameBin,queryFileName);
	
	//open the query file and save the tagged khmers in a vector
	std::cout<<"\tin Load_Queries...\n";
  	std::fstream inQfile(queryFileName.c_str(),std::ios::in| std::ios::binary);
  	assert(inQfile.is_open());
	khmer::HashIntoType khmer;
	std::vector<khmer::HashIntoType> qeuery_tagged_khmer;
	while (!(inQfile.eof())){
        inQfile>>khmer;
	//std::cout<<khmer<<" ";
	qeuery_tagged_khmer.push_back(khmer);
	}
	/*qeuery_tagged_khmer.push_back(5455492123);
	qeuery_tagged_khmer.push_back(5770913351);
	qeuery_tagged_khmer.push_back(1016590426926);
	qeuery_tagged_khmer.push_back(1007251177738);
	*/

	//retrive the read ids that contains at lest one of the query tagged khmer
	std::vector<long> reads_ids;
	std::string indexfilename;
        indexfilename=readsFileNameBin+".index";
	retrieve_read_ids_by_tag(indexfilename,qeuery_tagged_khmer,reads_ids );
	//output the rertive id result
	std::cout<<"\toutput the rertive id result\n";
	for (int i=0; i< reads_ids.size();i++)
		std::cout<<reads_ids[i]<<" ";
	std::cout<<std::endl;
	
	//retrieve reads by read ids
	std::vector<std::string> reads;
	retrieve_read_by_id(readsFileNameBin,reads_ids,reads);
	//output the retreieved reads
	std::cout<<"\toutput the retreieved reads\n";
	for (int i=0; i< reads.size(); i++)
		std::cout<<reads[i]<<std::endl;
	
	return 1;
}
