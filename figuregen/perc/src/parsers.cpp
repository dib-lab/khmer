//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. 
// Contact: khmer-project@idyll.org
//
#include "parsers.h"

IParser* IParser::get_parser(const std::string &inputfile)
{
   std::string filename(inputfile);
   int found = filename.find_last_of(".");

   std::string type = filename.substr(found+1);

   if (type == "fq" || type == "fastq") {
      return new FastqParser(inputfile);
   } else if (type == "fa" || type == "fasta") {
      return new FastaParser(inputfile);
   } else {
      return new FastaParser(inputfile);
   }
}

FastaParser::FastaParser(const std::string &inputfile) : 
                         infile(inputfile.c_str())
{
   std::string line, seq;
   next_name = "";
   
   assert(infile.is_open());

   bool valid_read = 0;

   while (!valid_read)  {
      line = "";
      seq = "";
      if (next_name == "")  {
        getline(infile, current_read.name);
        assert(current_read.name[0] == '>');
        current_read.name = current_read.name.substr(1);
      }
      else  {
         current_read.name = next_name;
         next_name = "";
      }
      
      while(line[0] != '>' && !infile.eof()) {
         getline(infile, line);
         if (line[0] != '>') {
            seq += line;
         }
      }

      if ((int)seq.find('N') == -1)  {
         valid_read = 1;
      }

      if (line[0] == '>') {
         next_name = line.substr(1);
      } else {
         seq += line;
      }

   }

   one_read_left = false;

   if (infile.eof())
   {
      one_read_left = true;
   }

   current_read.seq = seq;
}

Read FastaParser::get_next_read()
{
   std::string line = "", seq = "";
   Read next_read = current_read;

   if (one_read_left)  {
      one_read_left = false;
      return next_read;
   }

   bool valid_read = 0;

   while (!valid_read)  {
      current_read.name = next_name;
      next_name = "";
      current_read.seq = "";

      getline(infile, line);

      while(line[0] != '>' && !infile.eof())
      {

         if (line[0] != '>')  {
            seq += line;
         }
         getline(infile, line);
      }

      if (line[0] == '>')  {
         next_name = line.substr(1);
      }

      if ((int)seq.find('N') == -1)  {
         valid_read = 1;
      }  else if (infile.eof())  {
         one_read_left = false;
         break;
      }

      current_read.seq = seq;
      seq = "";

      if (infile.eof()) {
         one_read_left = true;
      }

   }

   return next_read;
}



FastqParser::FastqParser(const std::string &inputfile) :
                         infile(inputfile.c_str())
{
   std::string line_three, quality_scores;

   assert(infile.is_open());

   bool valid_read = 0;

   while (!valid_read && !infile.eof())  {
      getline(infile, current_read.name);
      getline(infile, current_read.seq);
      getline(infile, line_three); 
      getline(infile, quality_scores);

      assert(current_read.name[0] == '@');
      assert(line_three[0] == '+' || line_three[0] == '#');
      assert(quality_scores.length() == current_read.seq.length());
   
      current_read.name = current_read.name.substr(1);

      if ((int)current_read.seq.find('N') == -1)  {
         valid_read = 1;
      }
   }
}

Read FastqParser::get_next_read()
{
   Read next_read = current_read;
   std::string line_three, quality_scores;

   bool valid_read = 0;

   while (!valid_read && !infile.eof())  {

      getline(infile, current_read.name);

      if (infile.eof())  {
         return next_read;
      }

      getline(infile, current_read.seq);
      getline(infile, line_three);
      getline(infile, quality_scores);
   
      assert(current_read.name[0] == '@');
      assert(line_three[0] == '+' || line_three[0] == '#');
      assert(quality_scores.length() == current_read.seq.length());

      current_read.name = current_read.name.substr(1);

      if ((int)current_read.seq.find('N') == -1)  {
         valid_read = 1;
      }
   }

   return next_read;
}
