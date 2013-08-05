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
    //convertFastaToBin(readsFileName,readsFileNameBin);
    std::cout<<"done...\n";
    
    unsigned int ksize=0;
    unsigned int density=0;
    //intilization
    std::cout<<"\nload tagged khmers and sort them ...\n";
    std::vector<khmer::HashIntoType> sortedKhmerVector;
    //unsigned int save_ksize=0;
    //unsigned int  _tag_density=0; 

    load_tagset(tagedKhmerFileName,sortedKhmerVector,ksize,density);
    std::cout<<"done..., the num of tags are:"<<sortedKhmerVector.size()<<std::endl;
    std::cout<<"saved ksize is:"<<ksize<<"\n";
    std::cout<<"saved density is:"<<density<<"\n";
    //test the printing procedure
    std::cout<<"print the tagset into file xxx.tagset.dx.str\n";
    print_tagset(tagedKhmerFileName+"str", sortedKhmerVector, ksize);
    
    sortedKhmerVector.clear();
    //test the correctness of seeds functions
    load_seedset(tagedKhmerFileName,sortedKhmerVector,ksize,density);
    std::cout<<"print the seedset into file xxx.tagset.dx.str.seed\n";
    print_tagset(tagedKhmerFileName+"str.seed", sortedKhmerVector, ksize);


/*  //Iam not sure why this is not working any more!
    load_tagset(tagedKhmerFileName,sortedKhmerVector,save_ksize,_tag_density);
    std::cout<<"done..., the num of tags are:"<<sortedKhmerVector.size()<<std::endl;
    std::cout<<"	the saved k:"<<save_ksize<<std::cout<<std::endl;    std::cout<<"	the density:"<<_tag_density<<std::cout<<std::endl;*/
#if 0
    //build the index
    std::cout<<"\nbuild the index ...\n";
    //start the timer
    time_t start,end;
    double diff;
    time (&start);
    //build_index(readsFileNameBin,sortedKhmerVector,save_ksize);
    build_index(readsFileNameBin,density,sortedKhmerVector,ksize);
    time (&end);
    diff = difftime (end,start);
    printf ("It took you %.2lf seconds to build the tree.\n", diff );
    std::cout<<"done...\n";
  
    
    
    //query example
    std::cout<<"\nquery search ...\n";

    //open the query file and save the tagged khmers in a vector
    //this is an old version where we give only tagged k-mers
    std::cout<<"\tin Load_Queries...\n";
    /*std::fstream inQfile(queryFileName.c_str(),std::ios::in| std::ios::binary);
    assert(inQfile.is_open());
    khmer::HashIntoType khmer;
    std::vector<khmer::HashIntoType> qeuery_tagged_khmer;
    while (!(inQfile.eof())) {
        inQfile>>khmer;
        //std::cout<<khmer<<" ";
        qeuery_tagged_khmer.push_back(khmer);
    }
*/
    
    std::string indexfilename;
    std::string statfilename;
    std::string density_str="";	std::ostringstream convert;	convert << density;	density_str = convert.str(); 
    indexfilename=readsFileNameBin+".index.d"+density_str;
    statfilename=readsFileName+".d"+density_str+".stat";// give info like desinty;
    //the query file is set of seq sampled from input file
    std::fstream inQfile(queryFileName.c_str(),std::ios::in);
    std::fstream out_stat_file(statfilename.c_str(),std::ios::out);
    assert(inQfile.is_open());	assert(out_stat_file.is_open());
    out_stat_file<<"q_seq_id r_seq_id score\n";
    std::string seq="";
    std::vector<khmer::HashIntoType> qeuery_tagged_khmer;
    std::vector<long> reads_ids;
    std::vector<std::string> reads;
    unsigned int q_id=1;
    bool enable_scoring=true;
    unsigned int score,min, max, sum;
    unsigned int global_score=0,global_min=UINT_MAX, global_max=0, global_sum, global_cnt=0;
    while (!(inQfile.eof())) {
        if (q_id==11)	break; //run for specisfic num of queries  
	inQfile>>seq;
        //std::cout<<"\nseq id:"<<q_id<<" "<<seq<<std::cout<<std::endl;
        extract_tags_from_seq(seq,ksize,indexfilename,qeuery_tagged_khmer);
	//std::cout<<"the num of tagged k-mers in the q-seq is "<<qeuery_tagged_khmer.size()<<std::endl;
 	retrieve_read_ids_by_tag(indexfilename,qeuery_tagged_khmer,reads_ids );
	//output the rertive id result
	/*std::cout<<"\toutput the rertive id result\n";
	for (int i=0; i< reads_ids.size(); i++) {
        	std::cout<<reads_ids[i]<<" ";
    	}*/
	if (enable_scoring){
		//std::cout<<"\tcompute the similarities\n";
		//retrieve reads by read ids
		retrieve_read_by_id(readsFileNameBin,reads_ids,reads);
		//output the retreieved reads
		//std::cout<<"\toutput the retreieved reads\n";
		min=UINT_MAX; max=0; sum=0;
		for (int i=0; i< reads_ids.size(); i++) {
			score=0;
        		//std::cout<<reads[i]<<std::endl;
			score=sim_measure(seq,reads[i],ksize);
			//std::cout<<"score "<<score<<std::endl;
			out_stat_file<<q_id<<"\t"<<i+1<<"\t"<<score<<"\n";
			if (min>score)	min=score;
			if (max<score)	max=score;
			sum+=score;
    		}
		//std::cout<<"#reads:"<< reads_ids.size()<<" min:"<<min<<" max:"<<max<<" avg:"<<float(sum)/float(reads_ids.size())<<std::endl;
		global_cnt+=reads_ids.size(); global_score+=sum; if (global_min>min) global_min=min; if (global_max<max) global_max=max;
	}
	
	// clean the var for next query seq
	seq="";	qeuery_tagged_khmer.clear();
	reads_ids.clear(); q_id++;
    }           
    std::cout<<"global stat:\n";
    std::cout<<"global#reads:"<< (float)global_cnt/(float)q_id<<" global min:"<<global_min<<"global max:"<<global_max<<" global avg:"<<float(global_score)/float(global_cnt)<<std::endl;
    inQfile.close();
    out_stat_file.close();
     
 
/*
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
  #endif 
    // the sampling process
    //samplefrombinary();
    return 0;
}
