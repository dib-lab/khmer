//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2014. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. 
//Contact: khmer-project@idyll.org
//
#ifndef PARSERS_H
#define PARSERS_H

#include <iostream>
#include <string>
#include <string.h>
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
   std::string next_name;
   bool one_read_left;
public:
   FastaParser(const std::string &inputfile);
   ~FastaParser() { infile.close();  }
   Read get_next_read();
   bool is_complete() { return !one_read_left && infile.eof(); } 
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
