#include "parsers.h"

FastaParser::FastaParser(const std::string &inputfile) : 
                         infile(inputfile.c_str())
{
   assert(infile.is_open());

   getline(infile, current_read.name);
   getline(infile, current_read.seq);

   assert(current_read.name[0] == '>');

   current_read.name = current_read.name.substr(1);
}

Read FastaParser::get_next_read()
{
   Read next_read = current_read;

   getline(infile, current_read.name);

   if (infile.eof())  {
      return next_read;
   }

   getline(infile, current_read.seq);

   assert(current_read.name[0] == '>' || infile.eof());

   current_read.name = current_read.name.substr(1);

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
void test(IParser* parser)
{
   while(!parser->is_complete())  {
      Read next_read = parser->get_next_read();
      std::cout << next_read.name << ": " << next_read.seq << std::endl;
   }
}

int main()
{
   FastaParser parser("test.fa");
   test(&parser);

   FastqParser parser2("test.fq");
   test(&parser2);

   return 0;
}
*/
