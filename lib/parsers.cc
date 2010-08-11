#include "parsers.hh"

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
   std::string line, seq = "";

   assert(infile.is_open());

   getline(infile, current_read.name); 
   assert(current_read.name[0] == '>');
   current_read.name = current_read.name.substr(1);
   
   while(line[0] != '>' && !infile.eof()) {
      getline(infile, line);
      if (line[0] != '>') {
         seq += line;
      }
   }

   if (line[0] == '>') {
      next_name = line.substr(1);
   } else {
      seq += line;
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

   current_read.name = next_name;
   next_name = "";
   current_read.seq = "";

   while(line[0] != '>' && !infile.eof())
   {
      getline(infile, line);

      if (line[0] != '>') {
         seq += line;
      }
   }

   if (line[0] == '>') {
      next_name = line.substr(1);
   }

   current_read.seq = seq;

   if (infile.eof()) {
      one_read_left = true;
   }

   return next_read;
}

FastqParser::FastqParser(const std::string &inputfile) :
                         infile(inputfile.c_str())
{
   std::string line_three, quality_scores;

   assert(infile.is_open());

   getline(infile, current_read.name);
   getline(infile, current_read.seq);
   getline(infile, line_three); 
   getline(infile, quality_scores);

   assert(current_read.name[0] == '@');
   assert(line_three[0] == '+');
   assert(quality_scores.length() == current_read.seq.length());

   current_read.name = current_read.name.substr(1);
}

Read FastqParser::get_next_read()
{
   Read next_read = current_read;
   std::string line_three, quality_scores;

   getline(infile, current_read.name);

   if (infile.eof())  {
      return next_read;
   }

   getline(infile, current_read.seq);
   getline(infile, line_three);
   getline(infile, quality_scores);
   
   assert(current_read.name[0] == '@');
   assert(line_three[0] == '+');
   assert(quality_scores.length() == current_read.seq.length());

   current_read.name = current_read.name.substr(1);

   return next_read;
}

/*
int main()
{
   IParser* parser = IParser::get_parser("test.fa");
   
   while(!parser->is_complete())  {
      Read next_read = parser->get_next_read();
      std::cout << next_read.name << ": " << next_read.seq << std::endl;
   }

   delete parser;

   IParser* parser2 = IParser::get_parser("test.fq");
   while (!parser2->is_complete()) {
      Read next_read = parser->get_next_read();
      std::cout << next_read.name << ": " << next_read.seq << std::endl;
   }

   delete parser2;

   return 0;
}
*/
