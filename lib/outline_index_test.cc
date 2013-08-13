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

  
    //creaate reads binary file
    std::cout<<"\ncreating reads binary file ....\n";
    std::string readsfilenamebin=readsfilename+".bin";
    convertFastaToBin(readsfilename,readsfilenamebin);
    std::cout<<"done...\n";
  
  
    unsigned int ksize   = 0;
    unsigned int density = 0;

    std::cout<<"sampling k-mer postions\n";
    //std::cout<<"note: to sample every k-mer you must send d=1 if d=0 no k-mers will be chosen\n";
    ksize = kmer_size ; density = density_size ;
    std::set<khmer::HashIntoType> all_seeds;
    all_seeds =  consume_fasta_and_tag( readsfilename, ksize, density);	
    std::cout<<"d:"<<density<<"\t"<<"#num of seeds:"<<all_seeds.size()<<std::endl;
    for (std::set<khmer::HashIntoType>::iterator pi = all_seeds.begin(); pi != all_seeds.end();
         pi++) {
    	std::cout<<*pi<<" ";
    	}
    std::cout<<std::endl;
    std::string density_str=""; std::ostringstream convert;     convert << density;     density_str = convert.str();
    std::string outfilename=readsfilename+".seedset.d"+density_str;
    //if the seeds is created then it is will be default tags file 
    taggedkmerfilename=outfilename;
    std::cout<<"saving the seeds set to file:"<<outfilename<<std::endl;
    save_seedset(outfilename,all_seeds,ksize, density);
 
    //intilization
    std::cout<<"\nload seeds k-mers and sort them ...\n";
    std::vector<khmer::HashIntoType> sorted_kmer_vector;
    //unsigned int save_ksize=0;
    //unsigned int  _tag_density=0; 
    /*
    load_tagset(taggedkmerfilename,sorted_kmer_vector,ksize,density);
    std::cout<<"done..., the num of tags are:"<<sorted_kmer_vector.size()<<std::endl;
    std::cout<<"saved ksize is:"<<ksize<<"\n";
    std::cout<<"saved density is:"<<density<<"\n";
    test the printing procedure
    std::cout<<"print the tagset into file xxx.tagset.dx.str\n";
    print_tagset(taggedkmerfilename+"str", sorted_kmer_vector, ksize);
    sortedKhmerVector.clear();
    */
    //load seeds functions
    
    load_seedset(taggedkmerfilename,sorted_kmer_vector,ksize,density);
    std::cout<<"test if we correctly load the the seeds\n";
    for (int i=0; i<sorted_kmer_vector.size(); i++) {
        std::cout<<sorted_kmer_vector[i]<<" ";
    }
    std::cout<<std::endl;
    //std::cout<<"print the seedset into file xxx.seedset.dx.str.seed\n";
    //print_seedset(taggedkmerfilename+"str.",sorted_kmer_vector, ksize);


     //Iam not sure why this is not working any more!
    /*load_tagset(tagedKhmerFileName,sortedKhmerVector,save_ksize,_tag_density);
    std::cout<<"done..., the num of tags are:"<<sortedKhmerVector.size()<<std::endl;
    std::cout<<"	the saved k:"<<save_ksize<<std::cout<<std::endl;    
    std::cout<<"	the density:"<<_tag_density<<std::cout<<std::endl;
*/
    //build the index
    std::cout<<"\nbuild the index ...\n";
    //start the timer
    time_t start,end;
    double diff;
    time (&start);
    //build_index(readsFileNameBin,sortedKhmerVector,save_ksize);
    build_index(readsfilenamebin,density,sorted_kmer_vector,ksize);
    time (&end);
    diff = difftime (end,start);
    printf ("It took you %.2lf seconds to build the tree.\n", diff );
    std::cout<<"done...\n";
  
    
    
    //query example
    std::cout<<"\tin Load_Queries...\n";
    
    std::string indexfilename;	std::string statfilename;
    //std::string density_str="";	std::ostringstream convert;	convert << density;	density_str = convert.str(); 
    indexfilename=readsfilenamebin+".index.d"+density_str;
    statfilename=readsfilename+".stat.d"+density_str;// give info like desinty;
  
  //the query file is set of seq sampled from input file
    std::fstream inQfile(queryfilename.c_str(),std::ios::in);
    std::fstream out_stat_file(statfilename.c_str(),std::ios::out);
    assert(inQfile.is_open());	assert(out_stat_file.is_open());
    out_stat_file<<"q_seq_id r_seq_id score\n";
    
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
        total_q_seq++;
	seq=read.sequence;
	std::cout<<"process q-seq-id"<<total_q_seq<<std::endl;
	extract_tags_from_seq(seq,ksize,indexfilename,qeuery_tagged_khmer);
	retrieve_read_ids_by_tag(indexfilename,qeuery_tagged_khmer,reads_ids );
  	std::cout<<"output the rertive id result\n";
        for (std::set<long>::iterator pi = reads_ids.begin(); pi != reads_ids.end();
         pi++) {
        	std::cout<<*pi<<" ";
        	}
	
	std::cout<<std::endl;
	if (enable_scoring){
        //std::cout<<"\tcompute the similarities\n";
        //retrieve reads by read ids
	//retrieve_read_by_id(readsfilenamebin,reads_ids,reads);
	retrieve_read_by_id(readsfilenamebin,reads_ids,mymap);
	//output the retreieved reads
	std::cout<<"output the retreieved reads\n";
	for (it=mymap.begin(); it!=mymap.end(); ++it)
    		std::cout << it->first << "\t" << mymap[it->first]/*it->second*/ << '\n';
	
	/*
	for (int i=0; i< reads.size(); i++) {
		std::cout<<reads[i]<<std::endl;
		score = 0 ;
                score = sim_measure ( seq, reads[i], ksize ) ;
                std::cout<<"score "<<score<<std::endl;
		//note we are losing the actuall read id and just increment the vlaue
                out_stat_file << q_id << "\t" << i+1 << "\t" << score << "\n" ;
                        }
	*/
		}
	//clean the var for next query seq
	seq=""; 
	qeuery_tagged_khmer.clear();
        reads_ids.clear(); 
	q_id++;
	mymap.clear();
	} // whle

    out_stat_file.close(); 
 
    bool enable_stat=false;
    //test the create_stat fucntion
    if (enable_stat)	{
    	std::cout<<"create and print statistics\n";
    	create_stat(statfilename);
 	}
    // the sampling process
    //samplefrombinary();
    
    return 0;
}
