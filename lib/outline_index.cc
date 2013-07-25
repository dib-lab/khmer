//this is the codes for the mini lib function
#include "khmer.hh"

#include "outline_index.hh"
using namespace khmer;
//----------------------------------------
void khmer::convertFastaToBin(std::string readsFileName,std::string readsBinFileName)
{
    //std::cout<<"in convertFastaToBin...\n";
    std::fstream readBinFile;

    readBinFile.open(readsBinFileName.c_str(),std::ios::out|std::ios::binary);
    readBinFile.clear();

    long total_reads = 0;

    using namespace khmer:: read_parsers;

    IParser *                 parser  = IParser::get_parser(readsFileName.c_str());

    Read read;
    readNode readBin;
    std::string seq = "";
    //you may need to consider the other fields added in bleeding edge branch
    while(!parser->is_complete())  {
        read = parser->get_next_read();
        total_reads++;
        readBin.setPageNum(total_reads);
        //temp change the data the reads to upper cass
        //KHMER_EXTRA_SANITY_CHECKS flag in kmer iterators  must be set 
        std::transform(read.sequence.begin(), read.sequence.end(), read.sequence.begin(), ::toupper);
	 readBin.setName(read.name.size(),read.name);
        readBin.setSeq(read.sequence.size(),read.sequence);
        WriteToDiskRead(readBinFile,&readBin,total_reads);
    }
    std::cout<<"the num of reads consumed is :"<<total_reads<<std::endl;
    //write header of the bin file
    WriteToDiskHeader(readBinFile,&total_reads);
    readBinFile.close();
}
//------- Reads IO operations ----------
void khmer::WriteToDiskHeader(std::fstream& readsBinFile,long* numReads)
{
    //std::cout<<"WriteToDiskHeader\n";
    readsBinFile.seekp(0);
    readsBinFile.write((char*) numReads, sizeof(long));
}
void khmer::ReadFromDiskHeader(std::fstream& readsBinFile,long* numReads)
{
    //std::cout<<"ReadFromDiskHeader\n";
    readsBinFile.seekg(0);
    readsBinFile.read((char*)numReads,sizeof(long));
}

void khmer::WriteToDiskRead(std::fstream& readsBinFile,readNode* read,long pageNum)
{
    //std::cout<<"in WriteToDiskRead\n";
    readsBinFile.seekp(pageNum*blockSizeRead);
    readsBinFile.write((char*) read, sizeof(readNode));
}
void khmer::ReadFromDiskRead(std::fstream& readsBinFile,readNode* read,long pageNum)
{
    //std::cout<<"in ReadFromDiskRead\n";
    readsBinFile.seekg(pageNum*blockSizeRead);
    readsBinFile.read((char*)read,sizeof(readNode));
}

//------ tagSet IO operations ---------
void khmer::load_tagset(std::string infilename,std::vector<khmer::HashIntoType>& mykhmervector,unsigned int& save_ksize)
{
    //std::cout<<"in loading_tagset...\n";
    //std::cout<<"the file name sents is:"<<infilename<<std::endl;
    std::ifstream infile(infilename.c_str(), std::ios::binary);
    assert(infile.is_open());

    unsigned char version, ht_type;
    //unsigned int save_ksize = 0;
    unsigned int tagset_size = 0;
    unsigned int _tag_density = 0;
    infile.read((char *) &version, 1);
    infile.read((char *) &ht_type, 1);
    infile.read((char *) &save_ksize, sizeof(save_ksize));
    infile.read((char *) &tagset_size, sizeof(tagset_size));
    infile.read((char *) &_tag_density, sizeof(_tag_density));
    //std::cout<<"\nsave_ksize:"<<save_ksize<<" tagset_size:"<<tagset_size<<std::endl;

    khmer::HashIntoType * buf = new khmer::HashIntoType[tagset_size];
    infile.read((char *) buf, sizeof(khmer::HashIntoType) * tagset_size);

    //working with hashed hkmer
    for (int i=0; i<tagset_size; i++) {
        mykhmervector.push_back(buf[i]);
    }

    std::sort (mykhmervector.begin(),mykhmervector.end());
//intilize the classes
    /*for (std::vector<khmer::HashIntoType>::iterator it=mykhmervector.begin(); it!=mykhmervector.end(); ++it)
          std::cout<<*it<<" ";
    */
    delete buf;

    infile.close();
    return ;
}
//---- indexing  process ----
void khmer::build_index(std::string readsBinFileName, std::vector<khmer::HashIntoType>& sortedKhmerVector,unsigned int save_ksize)
{
    //std::cout<<"in building the index...\n";
    std::fstream readBinFile;
    readBinFile.open(readsBinFileName.c_str(),std::ios::in|std::ios::binary);
    long numReads=0;
    ReadFromDiskHeader(readBinFile,&numReads);
    std::cout<<"\tnum of reads in the file:"<<numReads<<std::endl;
    unsigned int numTkmer=sortedKhmerVector.size();
    std::cout<<"\tnum of tagged khmers:"<<numTkmer<<std::endl;
    
    int  classSize[numTkmer];		//the number of reads belong to this id
    long index;
    std::vector<long>* ptr[numTkmer];		//array of ptrs to set of read ids

    for (int i=0; i< numTkmer; i++) {
        classSize[i]=0;
        ptr[i]=new std::vector<long>;
    }

    readNode readBin;

    int seqLen;
    std::string  seq;
    std::string _seq="";
    for(int i=0; i<seqLen; i++) {
        _seq+=seq[i];
    }
    bool read_is_tagged;		//test if a read is successfully tagged
    //unsigned long error_cnt=0;	//keep track of the number of reads not tagged
    for (long i=0; i<numReads; i++) {
        //std::cout<<" read #:"<<i+1<<std::endl;

        //retreive a read number i
        ReadFromDiskRead(readBinFile,&readBin,i+1);
        seqLen=readBin.getSeqLength();
        seq=readBin.getSeq();
        for(int j=0; j<seqLen; j++) {
            _seq+=seq[j];
        }

        //extract the k-mers from this read and check if any is tagged one
        khmer::KMerIterator kmers(_seq.c_str(), save_ksize);
        khmer::HashIntoType kmer;
        std::vector<khmer::HashIntoType>::iterator low;
	//read_is_tagged=false;
        while(!kmers.done()) {
            kmer = kmers.next();
            low=std::lower_bound (sortedKhmerVector.begin(), sortedKhmerVector.end(), kmer);
            if (low != sortedKhmerVector.end() && *low == kmer) {
                index= low - sortedKhmerVector.begin();
                classSize[index]++;
		//std::cout<<"index:"<<index<<" i+1:"<<i+1<<std::endl;
                ptr[index]->push_back(i+1);
		read_is_tagged=true;// std::cout<<"read ID:"<<i+1<<"is tagged now\n";
            }
        } //while
	if (read_is_tagged==false){ 
		std::cout<<"ERROR read is not tagged\n";
		//error_cnt++;
	}

        // clean for the next read processing
        readBin.nullfy();
        _seq="";
    }
    std::cout<<"\tdone indexing Reads...\n";
    std::string outfilename2= readsBinFileName+".readDist";
    std::ofstream outfile2(outfilename2.c_str());
    std::cout<<"wrting the class size array into xxx.readDist file:\n";
    for (int i=0;i< numTkmer ;i++){ 
	outfile2<<i+1<<"\t"<<classSize[i]<<"\n";
	}
    outfile2.close();
    std::cout<<std::endl;
    
    std::cout<<"\tsaving the index information...\n";
    int indexSize=0;
    //save the index informaiton
    std::string outfilename= readsBinFileName+".index";
    std::ofstream outfile(outfilename.c_str(), std::ios::binary);

    //wrtie the size of taged k-mer sorted array
    outfile.write((const char *) &numTkmer, sizeof(numTkmer));
    //std::cout<<numTkmer<<std::endl;
    indexSize+=sizeof(numTkmer);

    //write the taged k-mer sorted array
    khmer::HashIntoType * buf = new khmer::HashIntoType[numTkmer];
    unsigned int i = 0;
    for (std::vector<khmer::HashIntoType>::iterator it = sortedKhmerVector.begin() ; it != sortedKhmerVector.end(); ++it,++i) {
        buf[i]= *it;/* std::cout<<*it<<" ";*/
    }
    std::cout<<std::endl;
    outfile.write((const char *) buf, sizeof(HashIntoType) * numTkmer);
    indexSize+=sizeof(HashIntoType) * numTkmer;
    delete buf;

    //create the array of acummulated sizes
    unsigned int sum=0;
    unsigned int * buff = new unsigned int [numTkmer];
    for (unsigned int i=0; i<numTkmer ; i++) {
        sum+=classSize[i];
        buff[i]=sum;
        //std::cout<<i<<":"<<sum<<"/";
    }
    outfile.write((const char *) buff, sizeof(unsigned int) * numTkmer);
    delete buff;
    indexSize+=sizeof(unsigned int) * numTkmer;
    //std::cout<<std::endl;

    //write set of associated arrays that map t-k-mer with a set of read ids
    long read_id;
    int cnt=0;
    for (long i=0; i< numTkmer ; i++) {
        //std::cout<<"\nT "<<i+1<<":";
        for (std::vector<long>::iterator it = ptr[i]->begin() ; it != ptr[i]->end(); ++it) {
            //std::cout << ' ' << *it;
            read_id=*it;
            outfile.write((const char *) &read_id, sizeof(long));
	    cnt++;
        }
        //	std::cout<<"/";
    }
    indexSize+=sizeof(long)*cnt;
    std::cout<<"the size of the index is:"<<indexSize<<std::endl;
    //std::cout<<std::endl;
    outfile.close();
    readBinFile.close();
}

//------ load index ------
//in this function we assume we can load all the index information into the memory

void khmer::load_index_header(std::string infilename,unsigned int& num_tagged_khmer,std::vector<khmer::HashIntoType>& sorted_khmer,std::vector<unsigned int>& accumulated_sizes)
{

    sorted_khmer.clear();
    accumulated_sizes.clear();

    std::ifstream infile(infilename.c_str(), std::ios::binary);
    assert(infile.is_open());

    unsigned int n_tagged_khmer = 0;
    infile.read((char*) &n_tagged_khmer,sizeof(unsigned int));
    num_tagged_khmer=n_tagged_khmer;

    HashIntoType * buf = new HashIntoType[num_tagged_khmer];
    infile.read((char *) buf, sizeof(HashIntoType) * num_tagged_khmer);
    for (unsigned int i = 0; i < num_tagged_khmer; i++) {
        sorted_khmer.push_back(buf[i]);
    }
    delete buf;

    unsigned int * buff = new unsigned int [num_tagged_khmer];
    infile.read((char *) buff, sizeof(unsigned int) * num_tagged_khmer);
    for (unsigned int i = 0; i < num_tagged_khmer; i++) {
        accumulated_sizes.push_back(buff[i]);
    }
    delete buff;
    infile.close();
}

//------ query -------
//given the index file and set of tagged kmers, retreieve the reads ids contanning this tagged khmers
void khmer::retrieve_read_ids_by_tag(std::string infilename,std::vector<khmer::HashIntoType>& qeuery_tagged_khmer,std::vector<long>& reads_ids )
{
    //std::cout<<"in retrieve_read_ids_by_tag...\n";
    // laod the index header
    unsigned int num_tagged_khmer=0;
    std::vector<khmer::HashIntoType> sorted_khmer;
    std::vector<unsigned int> accumulated_sizes;

    // open the index file to read the ID
    std::ifstream infile(infilename.c_str(), std::ios::binary);
    assert(infile.is_open());

    infile.read((char*) &num_tagged_khmer,sizeof(unsigned int));
    HashIntoType * buf1 = new HashIntoType[num_tagged_khmer];
    infile.read((char *) buf1, sizeof(HashIntoType) * num_tagged_khmer);
    for (unsigned int i = 0; i < num_tagged_khmer; i++) {
        sorted_khmer.push_back(buf1[i]);
    }
    delete buf1;
    unsigned int * buf2 = new unsigned int [num_tagged_khmer];
    infile.read((char *) buf2, sizeof(unsigned int) * num_tagged_khmer);
    for (unsigned int i = 0; i < num_tagged_khmer; i++) {
        accumulated_sizes.push_back(buf2[i]);
    }
    delete buf2;

    //save the starting location
    unsigned int start,index,offset,num_reads;
    start=sizeof(unsigned int)+num_tagged_khmer*sizeof(khmer::HashIntoType)+num_tagged_khmer*sizeof(unsigned int);
    std::vector<khmer::HashIntoType>::iterator low;
    khmer::HashIntoType kmer;
    // go through all the query tagged khmers
    //std::cout<<"strating query process...\n";
    for (int i=0; i<qeuery_tagged_khmer.size(); i++) {
        //std::cout<<"query # "<<i+1<<std::endl;
        kmer=qeuery_tagged_khmer[i];
        index=-1;
        //locate khmer of interest
        low=std::lower_bound (sorted_khmer.begin(), sorted_khmer.end(), kmer);
        if (low != sorted_khmer.end() && *low == kmer) {
            index= low - sorted_khmer.begin();
        }
        if (index==-1) {
            continue;    // query khmer is not in the tagged khmer set
        }
        if (index==0) {
            offset=0;
        } else {
            offset=accumulated_sizes[index-1];    //num of reads ids to be skipped
        }
        if (index==0) {
            num_reads=accumulated_sizes[index];
        } else {
            num_reads=accumulated_sizes[index]-accumulated_sizes[index-1];
        }
        //seeek the file
        infile.seekg(start+offset*sizeof(long));
        unsigned long * buf = new unsigned long [num_reads];
        infile.read((char *) buf, sizeof(long) * num_reads);
        for (unsigned int i = 0; i < num_reads; i++) {
            reads_ids.push_back(buf[i]);
        }
    }
    infile.close();
}

//given the read binary file and set of read ids, retrieve reads
void khmer::retrieve_read_by_id(std::string readsBinFileName, std::vector<long>& reads_ids, std::vector<std::string>& reads)
{
    //std::cout<<"in retrieve_read_by_id..\n";
    //open read binary file
    std::fstream readBinFile;
    readBinFile.open(readsBinFileName.c_str(),std::ios::in|std::ios::binary);

    reads.clear();
    readNode readBin;
    unsigned int seqLen;
    std::string  seq;
    std::string _seq="";
    for (unsigned int i=0; i<reads_ids.size(); i++) {
        ReadFromDiskRead(readBinFile,&readBin,reads_ids[i]);
        seqLen=readBin.getSeqLength();
        seq=readBin.getSeq();
        for(int j=0; j<seqLen; j++) {
            _seq+=seq[j];
        }
        reads.push_back(_seq);
        // clean for the next read processing
        readBin.nullfy();
        _seq="";
    }
    readBinFile.close();
}

//------ exact query ------
/*void khmer::exactQuery(std::string readsBinFileName,std::string queryFileName){
  std::cout<<"in Load_Queries...\n";
  std::fstream inQfile(queryFileName.c_str(),std::ios::in| std::ios::binary);
  assert(inQfile.is_open());
  std::fstream inRfile(readsBinFileName.c_str(),std::ios::in|std::ios::binary);
  assert(inRfile.is_open());

 khmer::HashIntoType khmer;
 std::stringstream ss;
 std::ifstream inputClassFile; std::string inputClassFileName;
 long location;	readNode readBin;char ch;
 //while not end of query file
 while (!(inQfile.eof())){
 	inQfile>>khmer;		ss.str("");	ss<<khmer;
	std::cout<<"qKhmer:"<<khmer<<std::endl;
 	inputClassFileName="classes/class";  inputClassFileName+=ss.str();
 	inputClassFile.open(inputClassFileName.c_str(),std::ios::in);
 	if (!inputClassFile.is_open()) std::cout<<"this tagged khmer is new!!!!\n";
 	else { location=-1;
		while (!(inputClassFile.eof())){
			//location=-1;
			inputClassFile>>location; std::cout<<"\tlocation:"<<location;
 			readBin.nullfy();
			ReadFromDiskRead(inRfile,&readBin,location);
 			std::cout<<"\tread:";	readBin.printSeq(); std::cout<<std::endl;
			inputClassFile>>ch;
		}
 	inputClassFile.close();
	}
 }
 inQfile.close();
 inRfile.close();

}
*/
//------- searching procedures -------------
void khmer::exhaustive_search(std::string readsBinFileName,std::string queryFileName){
	std::cout<<"exhaustive_search\n";
}
void khmer::approximate_search(std::string readsBinFileName,std::string queryFileName){
	std::cout<<"approximate_search\n";;
}

//------- Sampling Procedures ---------------
void khmer::samplefrombinary(){
    std::cout<<"in sampling from binary sequnce file\n";
    /*
    this function  create set of qureirs randomly from a given input file
 input: 1- binary file where each entry contins DNA sequnace
        2- size of the smaple, the number of the reads to be sampled
        3- fixed length of the sample, it is optinal in case we want to extract an equal leng
th substrings
        4- if length is not providede then random length will be chosen
        5- full flag is used to smaple the whole read not a substring from a read
    */
    std::string   readsBinFileName="", queryFileName="";
    int num_q=0;
    //mandatory inputs for all cases
    std::cout<<"enter the reads Binary file name:";	std::cin>>readsBinFileName;	std::cout<<std::endl;
    std::cout<<"enter the output query file name:";	std::cin>>queryFileName;	std::cout<<std::endl;
    std::cout<<"enter the number of qureies:",		std::cin>>num_q;		std::cout<<std::endl;
    
    std::fstream readBinFile;
    readBinFile.open(readsBinFileName.c_str(),std::ios::in|std::ios::binary);
    long numReads=0;
    ReadFromDiskHeader(readBinFile,&numReads);
    std::cout<<"\tnum of reads in the file:"<<numReads<<std::endl;
    std::fstream queryFile;
    queryFile.open(queryFileName.c_str(),std::ios::out);
    
    //some cases
    bool full_length=true;
    int n;
    std::cout<<"default sampling is the whole read,\n";
    std::cout<<" if you want sub-sting enter 0 otherwise enter 1:"; std::cin>>n;	std::cout<<std::endl;
    if (n==0) full_length=false;
   //we need to do the follwoing
   //randomly select k , the skips over the reads, to slect a read 
   //if we are in case one we are don
   //if we are in case two , then either fixed length of variable length
   //in fixed length case we have L, then select reandomly a pos in (s,e-L+1) then str=seq[pos,pos+l] 
   //if substring is vaialbe, then selct the pos in(s,e-1) and then select l in (1,e-pos+1) then str=seq[pos,pos+l)
   long location=0; readNode readBin; 
  if (full_length){
	std::cout<<"in smapling the entire reads mode\n";
	unsigned int k; k=rand()%100; std::cout<<"k:"<<k<<std::endl; //skipping paramater
	for (unsigned int i=0; i<num_q; i++){
		location+=(long)((i+k)%numReads)+1;
		std::cout<<location<<std::endl;
		readBin.nullfy();
                ReadFromDiskRead(readBinFile,&readBin,location);
                readBin.printSeq(queryFile); queryFile<<std::endl;
	
	}


	}
  if (!full_length){
	std::cout<<" in smapling sub strings mode\n";
	std::cout<<" the defualt is to sample variable length sub strings\n";
	std::cout<<" if you want to smaple with fixed length please enter the length othewise enter 0:";
	std::cin>>n; std::cout<<std::endl;
	if (n==0){
		std::cout<<"We are in smapling substring with variable length\n";
	}
	else	{
		std::cout<<" we are in sampling substring with fixed length equal to "<<n<<std::endl;
	}
	}

   readBinFile.close();
   queryFile.close();
}//end sampling procedure
