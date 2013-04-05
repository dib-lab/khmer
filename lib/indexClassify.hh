// this is a mini lib for indexing and query stuff
#ifndef MiniLib_H
#define MiniLib_H

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <string>               //used for files names
#include <stdlib.h>             //used for c_str() function

#include "parsers.hh"
#include "ktable.hh"
#include "hashbits.hh"
//using namespace std;

#define maxNameLength 80
#define maxSeqLength 200


class readNode{
protected:
	long _pageNum;
	char _name[maxNameLength];
	char _seq[maxSeqLength];
	int _nameLength;
	int _seqLength;
	void init(){
		_pageNum=-1;_nameLength=_seqLength=0;
		for (int i=0; i<maxNameLength; i++) _name[i]='*';
		for (int i=0; i<maxSeqLength; i++) _seq[i]='*';
	}
public:
	//constructor
	readNode()	{init();}
	//destructor
	~readNode(){ ;}
	//accessors to get values
	long	getPageNum()	{return _pageNum;}
	int 	getNameLength()	{return _nameLength;}
	int 	getSeqLength()	{return _seqLength;}
	char *	getName()	{return _name;}
	char *	getSeq()	{return _seq;}

	//accessors to set values
	void setPageNum(long newPageNum)	{_pageNum=newPageNum;}
	void setNameLength(int newNameLength)	{_nameLength=newNameLength;}
	void setSeqLength(int newSeqLength)	{_seqLength=newSeqLength;}
	void setName(int newNameLength,std::string  newName)	{
		for (int i=0; i<newNameLength; i++) _name[i]=newName[i];
		this->setNameLength(newNameLength);
	}
	void setSeq(int newSeqLength,std::string  newSeq){
		for (int i=0; i<newSeqLength; i++) _seq[i]=newSeq[i];
		this->setSeqLength(newSeqLength);
	}
	void nullfy()	{this->init();}
	// read IO
	void printPageNum()	{std::cout<<_pageNum;}
	void printName()	{for (int i=0; i<_nameLength; i++)	std::cout<<_name[i];}
	void printSeq()		{for (int i=0; i<_seqLength ; i++)	std::cout<<_seq[i];}
	void printPageNum(std::fstream& outFile){outFile<<_pageNum;}
	void printName(std::fstream& outFile)	{for (int i=0; i<_nameLength; i++)	outFile<<_name[i];}
	void printSeq(std::fstream& outFile)	{for (int i=0; i<_seqLength ; i++)	outFile<<_seq[i];}

};

#define blockSizeRead sizeof(readNode)

//----------------------------------------------
void convertFastaToBin(std::string readsFileName,std::string readBinFileName);
//------------ Reads IO operations ---------------
void WriteToDiskHeader(std::fstream& readBinFile,long* numReads);
void ReadFromDiskHeader(std::fstream& readBinFile,long* numReads);

void WriteToDiskRead(std::fstream& readBinFile,readNode* read,long pageNum);
void ReadFromDiskRead(std::fstream& readBinFile,readNode* read,long pageNum);
//----------- tagSet IO operation --------------
void load_tagset(std::string infilename, std::vector<khmer::HashIntoType>& tagedKhmerVector,unsigned int& k);
void loading_tagset(std::string infilename,std::vector<khmer::HashIntoType>& mykhmervector, unsigned int& k);
void openNewClass(khmer::HashIntoType newClass);
//--------- classification -----------
void classifyReads(std::string readsBinFileName,std::vector<khmer::HashIntoType>& sortedKhmerVector,unsigned intsave_ksize );
void classifyOneRead(long location,int seqLen, char* seq,std::vector<khmer::HashIntoType>& sortedKhmerVector,unsigned int save_ksize);
void insertReadLocationToClassFile(khmer::HashIntoType suffixName,long location);
//--------- exact query -----------
void exactQuery(std::string readsBinFileName,std::string queryFileName);
#endif //MiniLib_H
