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
    std::string  taggedkmerfilename="", readsfilename="", queryfilename="";

    unsigned int optind=1,kmer_size=0, density_size=0;

    while ((optind < argc) && (argv[optind][0]=='-')) {
        std::string sw = argv[optind];
        if (sw=="-tf") {
            optind++;
            taggedkmerfilename=argv[optind];
        } else if (sw=="-rf") {
            optind++;
            readsfilename=argv[optind];
        } else if (sw=="-qf") {
            optind++;
            queryfilename = argv[optind];
        } else if(sw=="-ks") {
            optind++;
            kmer_size=atoi(argv[optind]);
        } else if(sw=="-ds") {
            optind++;
            density_size=atoi(argv[optind]);
        }

        else if (sw=="-help") {
            std::cout<<"the following are the prgram arguments\n";
            std::cout<<"\t-tf  TagedKhmerFileName\n";
            std::cout<<"\t-rf  ReadsFileName\n";
            std::cout<<"\t-qf  QueryFileName\n";
            std::cout<<"\t-ks  k-mer Size\n";
            std::cout<<"\t-ds   density\n";
            std::cout<<"\t-help disply this help\n";
            return 0;
        } else {
            std::cout << "Unknown switch: " << argv[optind] << std::endl;
        }

        optind++;
    }

    //step 1: creaate reads binary file
    //std::cout<<"Creating reads binary file ....\n";
    std::string readsfilenamebin=readsfilename+".bin";
    //convertFastaToBin(readsfilename,readsfilenamebin);
    //std::cout<<"\tdone...\n";
    
    unsigned int ksize   = 0;
    unsigned int density = 0;

    //step 2: choose the seeds
    ksize = kmer_size ; density = density_size ;
    std::set<khmer::HashIntoType> all_seeds;
    //std::cout<<"Sampling k-mer every "<<density<<" positions\n";
    //all_seeds =  consume_fasta_and_tag( readsfilename, ksize, density);	
    //std::cout<<"\tdensity:"<<density<<"\t"<<"#num of seeds:"<<all_seeds.size()<<std::endl;
    /*for (std::set<khmer::HashIntoType>::iterator pi = all_seeds.begin(); pi != all_seeds.end();
         pi++) {
    	std::cout<<*pi<<" ";
    	}
    std::cout<<std::endl;
    */
    std::string density_str=""; std::ostringstream convert;     convert << density;     density_str = convert.str();
    std::string outfilename=readsfilename+".seedset.d"+density_str;
    //if the seeds is created then it is will be default tags file 
    taggedkmerfilename=outfilename;
    //std::cout<<"\tsaving the seeds set to file:"<<outfilename<<std::endl;
    //save_seedset(outfilename,all_seeds,ksize, density);
 
    //step 3: load the seeds
    std::cout<<"load seeds k-mers and sort them ...\n";
    std::vector<khmer::HashIntoType> sorted_kmer_vector;
    //unsigned int save_ksize=0;
    //unsigned int  _tag_density=0; 
    
    load_seedset(taggedkmerfilename,sorted_kmer_vector,ksize,density);

    /*std::cout<<"test if we correctly load the the seeds\n";
    for (int i=0; i<sorted_kmer_vector.size(); i++) {
        std::cout<<sorted_kmer_vector[i]<<" ";
    }
    std::cout<<std::endl;
    */

    //step 4: build the index
    //std::cout<<"Build the index ...\n";
    //start the timer
    time_t start,end;
    double diff;
    time (&start);
    //build_index(readsfilenamebin,density,sorted_kmer_vector,ksize);
    time (&end);
    diff = difftime (end,start);
    //printf ("\tIt took you %.2lf seconds to build the tree.\n", diff );
    //std::cout<<"\tdone...\n";
  
    //step 5: Query process
    std::cout<<"Start Load_Queries...\n";
    std::string indexfilename;	std::string statfilename;
    indexfilename=readsfilenamebin+".index.d"+density_str;
    statfilename=readsfilename+".stat.d"+density_str;
  #if 0
    //the query file is set of seq sampled from input file
    std::fstream inQfile(queryfilename.c_str(),std::ios::in);
    std::fstream out_stat_file(statfilename.c_str(),std::ios::out);
    assert(inQfile.is_open());	assert(out_stat_file.is_open());
    out_stat_file<<"q_seq_id  r_seq_id  score\n";
    
    std::string seq="";
    std::set<khmer::HashIntoType> qeuery_tagged_khmer;
    std::set<long> reads_ids;
    std::vector<std::string> reads;

    std::map<long,std::string> mymap;
    std::map<long,std::string>::iterator it;

    unsigned int q_id = 1, score =0 ;
    bool enable_scoring = true;
    

   unsigned long long total_q_seq = 0;
   using namespace khmer:: read_parsers;
   IParser * parser  = IParser::get_parser(queryfilename.c_str());
   Read read;
   while(!parser->is_complete())  {
	read = parser->get_next_read();
	std::transform(read.sequence.begin(), read.sequence.end(), read.sequence.begin(), ::toupper);
        total_q_seq++;
	seq=read.sequence;
	//std::cout<<"process q-seq-id"<<total_q_seq<<std::endl;
	extract_tags_from_seq(seq,ksize,indexfilename,qeuery_tagged_khmer);
	retrieve_read_ids_by_tag(indexfilename,qeuery_tagged_khmer,reads_ids );
  	/*std::cout<<"output the rertive id result\n";
        for (std::set<long>::iterator pi = reads_ids.begin(); pi != reads_ids.end();
         pi++) {
        	std::cout<<*pi<<" ";
        	}
	std::cout<<std::endl;
	*/
	if (enable_scoring){
        //retrieve reads by read ids
	retrieve_read_by_id(readsfilenamebin,reads_ids,mymap);
	//std::cout<<"output the retreieved reads\n";
	for (it=mymap.begin(); it!=mymap.end(); ++it) {
    		//std::cout << it->first << "\t" << it->second << '\n';
		score = 0 ;
                score = sim_measure ( seq, it->second, ksize ) ;
                //std::cout<<"score "<<score<<std::endl;
                out_stat_file << q_id << "\t" << it->first << "\t" << score << "\n" ;
			}
		}
	//clean the var for next query seq
	seq=""; 
	qeuery_tagged_khmer.clear();
        reads_ids.clear(); 
	q_id++;
	mymap.clear();
	} // iterat over q-seq
    std::cout<<"\tdone...\n";
    out_stat_file.close(); 
#endif
    bool enable_stat=true;
    //test the create_stat fucntion
    if (enable_stat)	{
    	std::cout<<"create and print statistics\n";
    	create_stat(statfilename,sorted_kmer_vector.size(),density);
 	}

    // the sampling process
    //samplefromfile();
    //samplefrombinary();
    
    return 0;
}
