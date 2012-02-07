#ifndef THREADED_PARSERS_H
#define THREADED_PARSERS_H

#include <iostream>
#include <string>
#include <string.h>
#include <fstream>
#include <assert.h>

#define THREADED_PARSER_CHUNK_SIZE	104857600

namespace khmer
{

    namespace threaded_parsers
    {

	struct Read
	{
	    std::string name;
	    std::string seq;
	    //std::string quality;
	};

	/* Base class for parsers, defines a set of methods that a parser must
	 * implement and a static function that returns the proper parser
	 * for different file types
	 */
	class ThreadedIParser
	{
	public:
	    virtual Read get_next_read() = 0;
	    virtual bool is_complete() = 0;
	    virtual ~ThreadedIParser() { }
	    virtual long int getEndPos() = 0;
	};

	class ThreadedFastaParser : public ThreadedIParser
	{
	private:
	    std::ifstream infile;
	    long int endPos;

	public:
	   ThreadedFastaParser(const std::string &inputfile, long int startPos, long int end);
	   ~ThreadedFastaParser() { infile.close();  }
	   Read get_next_read();
	   bool is_complete() { return (infile.tellg() >= endPos || infile.eof()); } 
	   long int getEndPos() { return endPos; }
	};

	class ThreadedFastqParser : public ThreadedIParser
	{
	private:
	   std::ifstream infile;
	   long int endPos;

	public:
	   ThreadedFastqParser(const std::string &inputfile, long int startPos, long int end);
	   ~ThreadedFastqParser() { infile.close(); }
	   Read get_next_read();
	   bool is_complete() { return infile.tellg() >= endPos || infile.eof(); }
	   long int getEndPos() { return endPos; }
	};

	/* Base class for parser factories. Each factory will return a parser
	 * object that will return all the reads within a chunk of the file
	 */
	class ThreadedIParserFactory
	{
	public:
	    virtual ThreadedIParser* get_next_parser() = 0;
	    virtual bool is_complete() = 0;
	    static ThreadedIParserFactory* get_parser(const std::string &inputfile, long int chunkSize);
	};

	class ThreadedFastaParserFactory : public ThreadedIParserFactory
	{
	private:
	    std::string filename;
	    long int curPos;
	    long int fileSize;
	    long int chunkSize;

	public:
	    ThreadedFastaParserFactory(const std::string &inputfile, 
		long int size);
	    ThreadedIParser* get_next_parser();
	    bool is_complete() { return curPos >= fileSize; }
	};

	class ThreadedFastqParserFactory : public ThreadedIParserFactory
	{
	private:
	    std::string filename;
	    long int curPos;
	    long int fileSize;
	    long int chunkSize;

	public:
	    ThreadedFastqParserFactory(const std::string &inputfile, 
		long int chunkSize);
	    ThreadedIParser* get_next_parser();
	    bool is_complete() { return curPos >= fileSize; }
	};
    
    };

};

#endif

// vim: set sts=4 sw=4:
