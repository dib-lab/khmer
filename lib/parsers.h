#ifndef PARSERS_H
#define PARSERS_H

#include <iostream>
#include <string>
#include <fstream>
#include <assert.h>

struct Read
{
   std::string name;
   std::string seq;
   //std::string quality;
};

class IParser
{
public:
   virtual Read get_next_read() = 0;
   virtual bool is_complete() = 0;
   virtual ~IParser() { }
   static IParser* get_parser(const std::string &inputfile);
};


class FastaParser : public IParser
{
private:
   std::ifstream infile;
   Read current_read;
public:
   FastaParser(const std::string &inputfile);
   ~FastaParser() { infile.close();  }
   Read get_next_read();
   bool is_complete() { return infile.eof(); } 
};

class FastqParser : public IParser
{
private:
   std::ifstream infile;
   Read current_read;
public:
   FastqParser(const std::string &inputfile);
   ~FastqParser() { infile.close(); }
   Read get_next_read();
   bool is_complete() { return infile.eof(); }
};

#endif
