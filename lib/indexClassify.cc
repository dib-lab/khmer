//this is the codes for the mini lib function
#include "minilib.hh"

//----------------------------------------
void convertFastaToBin(std::string readsFileName,std::string readsBinFileName){
	std::cout<<"in convertFastaToBin...\n";
	std::fstream readBinFile;
	readBinFile.open(readsBinFileName.c_str(),std::ios::out|std::ios::binary); readBinFile.clear();
	
	long total_reads = 0;
	
	IParser* parser = IParser::get_parser(readsFileName.c_str());
  	Read read;
	readNode readBin;
	while(!parser->is_complete()) {
    		read = parser->get_next_read();
		//std::cout<<read.seq<<" "<<read.seq.size()<<std::endl;
		total_reads++;
		readBin.setPageNum(total_reads);
		readBin.setName(read.name.size(),read.name);
		readBin.setSeq(read.seq.size(),read.seq);
		WriteToDiskRead(readBinFile,&readBin,total_reads);
	} //end while 
	//write header of the bin file
	WriteToDiskHeader(readBinFile,&total_reads);
	readBinFile.close();
}
//------- Reads IO operations ----------
void WriteToDiskHeader(std::fstream& readsBinFile,long* numReads){
	//std::cout<<"WriteToDiskHeader\n";
	readsBinFile.seekp(0);
	readsBinFile.write((char*) numReads, sizeof(long));
}
void ReadFromDiskHeader(std::fstream& readsBinFile,long* numReads){
	//std::cout<<"ReadFromDiskHeader\n";
	readsBinFile.seekp(0);
	readsBinFile.read((char*)numReads,sizeof(long));
}

void WriteToDiskRead(std::fstream& readsBinFile,readNode* read,long pageNum){
	//std::cout<<"in WriteToDiskRead\n";
	readsBinFile.seekp(pageNum*blockSizeRead);
	readsBinFile.write((char*) read, sizeof(readNode));
}
void ReadFromDiskRead(std::fstream& readsBinFile,readNode* read,long pageNum) {
	//std::cout<<"in ReadFromDiskRead\n";
	readsBinFile.seekg(pageNum*blockSizeRead);
        readsBinFile.read((char*)read,sizeof(readNode));
}

//------ tagSet IO operations ---------
void load_tagset(std::string infilename,std::vector<khmer::HashIntoType>& mykhmervector, unsigned int& save_ksize)
{
  std::cout<<"in load_tagset...\n";
  std::ifstream infile(infilename.c_str(), std::ios::binary);
  assert(infile.is_open());

  //unsigned int save_ksize = 0;
  unsigned int tagset_size = 0;

  infile.read((char *) &save_ksize, sizeof(save_ksize));
  infile.read((char *) &tagset_size, sizeof(tagset_size));
  //std::cout<<"\nsave_ksize:"<<save_ksize<<" tagset_size:"<<tagset_size<<std::endl;
  
  khmer::HashIntoType * buf = new khmer::HashIntoType[tagset_size];
  infile.read((char *) buf, sizeof(khmer::HashIntoType) * tagset_size);
  
 //working with hashed hkmer
 for (int i=0; i<tagset_size; i++)
	mykhmervector.push_back(buf[i]);
 
 std::sort (mykhmervector.begin(),mykhmervector.end());
 //intilize the classes
  for (std::vector<khmer::HashIntoType>::iterator it=mykhmervector.begin(); it!=mykhmervector.end(); ++it)
   	openNewClass(*it);
 delete buf;
 
 infile.close();
 return ;
}

void loading_tagset(std::string infilename,std::vector<khmer::HashIntoType>& mykhmervector,unsigned int& save_ksize)
{
  std::cout<<"in loading_tagset...\n";
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
  std::cout<<"\nsave_ksize:"<<save_ksize<<" tagset_size:"<<tagset_size<<std::endl;

  khmer::HashIntoType * buf = new khmer::HashIntoType[tagset_size];
  infile.read((char *) buf, sizeof(khmer::HashIntoType) * tagset_size);

  //working with hashed hkmer
  for (int i=0; i<tagset_size; i++)
        mykhmervector.push_back(buf[i]);

 std::sort (mykhmervector.begin(),mykhmervector.end());
 //intilize the classes
  for (std::vector<khmer::HashIntoType>::iterator it=mykhmervector.begin(); it!=mykhmervector.end(); ++it)
        openNewClass(*it);
 delete buf;

 infile.close();
 return ;
}


//-----------------------
void openNewClass(khmer::HashIntoType newClass){
	//created and name the apropriate file
	//std::cout<<"in openNewClass...\n";
	std::stringstream ss;
        ss<<newClass;
        std::ofstream outputfile; std::string outputfilename;
        outputfilename="classes/class";	outputfilename+=ss.str();
        outputfile.open(outputfilename.c_str(),std::ios::out);
	outputfile.close();
}

//---- classifcation process ----
void classifyReads(std::string readsBinFileName, std::vector<khmer::HashIntoType>& sortedKhmerVector,unsigned int save_ksize){
	//std::cout<<"in classifyReads...\n";

	std::fstream readBinFile;
	readBinFile.open(readsBinFileName.c_str(),std::ios::in|std::ios::binary);
        long numReads=0;        
	ReadFromDiskHeader(readBinFile,&numReads);
        //std::cout<<"num of reads in the file:"<<numReads<<std::endl;

	readNode readBin;
	int seqLen; char* seq;
	for (long i=1; i<numReads+1; i++){
		ReadFromDiskRead(readBinFile,&readBin,i);
		seqLen=readBin.getSeqLength(); 
		seq=readBin.getSeq();
		classifyOneRead(i,seqLen,seq,sortedKhmerVector,save_ksize);
        	readBin.nullfy();
        }
	//std::cout<<"\t done classifyReads...\n";
        readBinFile.close();

}
void classifyOneRead(long location,int seqLen, char* seq,std::vector<khmer::HashIntoType>& sortedKhmerVector,unsigned int save_ksize){
  //std::cout<<"in classifyOneRead...\n";
  std::string _seq="";
  for(int i=0; i<seqLen; i++) _seq+=seq[i];
  const char * sp = _seq.c_str();
  khmer::KMerIterator kmers(sp, save_ksize);
  khmer::HashIntoType kmer;
  while(!kmers.done()) {
    kmer = kmers.next();
    //std::cout<<kmer<<std::endl;
    if (std::binary_search (sortedKhmerVector.begin(), sortedKhmerVector.end(), kmer)){
    //std::cout << "found!\n"; //else std::cout << "not found.\n";
    insertReadLocationToClassFile(kmer,location);}
 }
  //std::cout<<"\t done classifyOneRead...\n";
  return ;
}
//------ assinge read id into class -----
void insertReadLocationToClassFile(khmer::HashIntoType suffixName,long location){
	std::stringstream ss;
        ss<<suffixName;
        std::ofstream outputfile; std::string outputfilename;
        outputfilename="classes/class";  outputfilename+=ss.str();
        outputfile.open(outputfilename.c_str(),std::ios::app);
        if (!outputfile.is_open()) {std::cout<<"Error the data file cannot be open!!!!\n"; return;}
	outputfile<<location<<"\n";
        outputfile.close();
} 
//------ exact query ------
void exactQuery(std::string readsBinFileName,std::string queryFileName){
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

