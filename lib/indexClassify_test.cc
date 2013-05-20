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
	/*std::cout<<"\nexact query search ...\n";
	exactQuery(readsFileNameBin,queryFileName);
	*/

	std::string infilename;
	infilename=readsFileNameBin+".index";
	std::vector<khmer::HashIntoType> qeuery_tagged_khmer;
	
	qeuery_tagged_khmer.push_back(5455492123);
	qeuery_tagged_khmer.push_back(5770913351);
	qeuery_tagged_khmer.push_back(1016590426926);
	qeuery_tagged_khmer.push_back(1007251177738);
	
	std::vector<long> reads_ids;
	
	retrieve_read_ids_by_tag(infilename,qeuery_tagged_khmer,reads_ids );
	for (int i=0; i< reads_ids.size();i++)
		std::cout<<reads_ids[i]<<" ";
	std::cout<<std::endl;

	std::vector<std::string> reads;

	retrieve_read_by_id(readsFileNameBin,reads_ids,reads);
	for (int i=0; i< reads.size(); i++)
		std::cout<<reads[i]<<std::endl;
	
	#if 0
  	int myints[] = {10,20,30,30,20,10,10,20};
  	std::vector<int> v(myints,myints+8);
	 std::sort (v.begin(), v.end());                
	std::vector<int>::iterator low,up;
 	low=std::lower_bound (v.begin(), v.end(), 20);
	up= std::upper_bound (v.begin(), v.end(), 20); 
	std::cout << "lower_bound at position " << (low- v.begin()) << '\n';
 	 std::cout << "upper_bound at position " << (up - v.begin()) << '\n';
	#endif
	return 1;
}
