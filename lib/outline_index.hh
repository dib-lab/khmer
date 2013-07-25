// this is a mini lib for indexing and query stuff
#ifndef IndexClassify_H
#define IndexClassify_H

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>		//used for binary search
#include <sstream>		//used for file name string stream
#include <string>               //used for files names
#include <stdlib.h>             //used for c_str() function
#include "hashtable.hh"
#define maxNameLength 80
#define maxSeqLength 200

namespace khmer
{

typedef long ReadID;				//we may need long long if we have many reads
typedef unsigned int khmerID;			//each khmer has a rank
typedef std::set<HashIntoType> TaggedKhmerSet; //if we decide to change taggedkhmer from vector to set data str
typedef std::set<ReadID> SetReadId;
typedef std::map<khmerID, SetReadId*> TagsToReadsMap;
//typedef std::map<khmerID,SetReadId*>;  //assing for each khmerID a set of reads Iet

class readNode
{
protected:
    long _pageNum;
    char _name[maxNameLength];
    char _seq[maxSeqLength];
    int _nameLength;
    int _seqLength;
    void init() {
        _pageNum=-1;
        _nameLength=_seqLength=0;
        for (int i=0; i<maxNameLength; i++) {
            _name[i]='*';
        }
        for (int i=0; i<maxSeqLength; i++) {
            _seq[i]='*';
        }
    }
public:
    //constructor
    readNode()	{
        init();
    }
    //destructor
    ~readNode() {
        ;
    }
    //accessors to get values
    long	getPageNum()	{
        return _pageNum;
    }
    int 	getNameLength()	{
        return _nameLength;
    }
    int 	getSeqLength()	{
        return _seqLength;
    }
    char *	getName()	{
        return _name;
    }
    char *	getSeq()	{
        return _seq;
    }

    //accessors to set values
    void setPageNum(long newPageNum)	{
        _pageNum=newPageNum;
    }
    void setNameLength(int newNameLength)	{
        _nameLength=newNameLength;
    }
    void setSeqLength(int newSeqLength)	{
        _seqLength=newSeqLength;
    }
    void setName(int newNameLength,std::string  newName)	{
        for (int i=0; i<newNameLength; i++) {
            _name[i]=newName[i];
        }
        this->setNameLength(newNameLength);
    }
    void setSeq(int newSeqLength,std::string  newSeq) {
        for (int i=0; i<newSeqLength; i++) {
            _seq[i]=newSeq[i];
        }
        this->setSeqLength(newSeqLength);
    }
    void nullfy()	{
        this->init();
    }
    // read IO
    void printPageNum()	{
        std::cout<<_pageNum;
    }
    void printName()	{
        for (int i=0; i<_nameLength; i++)	{
            std::cout<<_name[i];
        }
    }
    void printSeq()		{
        for (int i=0; i<_seqLength ; i++)	{
            std::cout<<_seq[i];
        }
    }
    void printPageNum(std::fstream& outFile) {
        outFile<<_pageNum;
    }
    void printName(std::fstream& outFile)	{
        for (int i=0; i<_nameLength; i++)	{
            outFile<<_name[i];
        }
    }
    void printSeq(std::fstream& outFile)	{
        for (int i=0; i<_seqLength ; i++)	{
            outFile<<_seq[i];
        }
    }

};

#define blockSizeRead sizeof(readNode)

//----------------------------------------------
void convertFastaToBin(std::string readsFileName,std::string readBinFileName);
//------------ Reads IO operations ---------------
void WriteToDiskHeader(std::fstream& readBinFile,long* numReads);
void ReadFromDiskHeader(std::fstream& readBinFile,long* numReads);

void WriteToDiskRead(std::fstream& readBinFile,readNode* read,long pageNum);
void ReadFromDiskRead(std::fstream& readBinFile,readNode* read,long pageNum);
//----------- tagSet IO operation ------------
void load_tagset(std::string infilename,std::vector<khmer::HashIntoType>& mykhmervector, unsigned int& k);
//---------- index ----------
void build_index(std::string readsBinFileName,std::vector<khmer::HashIntoType>& sortedKhmerVector,unsigned int  save_ksize );
void load_index_header(std::string infilename,unsigned int& num_tagged_khmer,std::vector<khmer::HashIntoType>& sorted_khmer,std::vector<unsigned int>& accumulated_sizes);
//--------- exact query -----------
void retrieve_read_ids_by_tag(std::string infilename,std::vector<khmer::HashIntoType>& qeuery_tagged_khmer,std::vector<long>& reads_ids );

void retrieve_read_by_id(std::string infilename, std::vector<long>& reads_ids, std::vector<std::string>& reads);

void exactQuery(std::string readsBinFileName,std::string queryFileName);
//------ searching Procedures ------
void exhaustive_search(std::string readsBinFileName,std::string queryFileName);
void approximate_search(std::string readsBinFileName,std::string queryFileName);
//------ sampling procedure --------
void samplefrombinary();


};

#endif //IndexClassify_H
