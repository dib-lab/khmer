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
   } else if (type == "gz") {
      if ((signed int)filename.find("fastq") != -1 || 
          (signed int)filename.find("fq") != -1) {
         return new FastqGzParser(inputfile);
      }
      else {
         return new FastaGzParser(inputfile);
      }
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
      infile.close();
      one_read_left = true;
   }

   current_read.seq = seq;
}

FastaGzParser::FastaGzParser(const std::string &inputfile)
{
   infile = gzopen(inputfile.c_str(), "rb");

   std::string line, seq;
   int MAXLINE = 1000;
   char tmp[1000];
   next_name = "";

   bool valid_read = 0;

   while (!valid_read)  {
      line = "";
      seq = "";
      if (next_name == "")  {
         gzgets(infile, tmp, MAXLINE);
         tmp[strlen(tmp)-1] = '\0';
         current_read.name = tmp;
         //getline(infile, current_read.name);
         assert(current_read.name[0] == '>');
         current_read.name = current_read.name.substr(1);
      }
      else  {
         current_read.name = next_name;
         next_name = "";
      }

      while(line[0] != '>' && !gzeof(infile)) {
         gzgets(infile, tmp, MAXLINE);
         tmp[strlen(tmp)-1] = '\0';
         line = tmp;
         //getline(infile, line);
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

   if (gzeof(infile))
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
         infile.close();
         one_read_left = false;
         break;
      }

      current_read.seq = seq;
      seq = "";

      if (infile.eof()) {
         infile.close();
         one_read_left = true;
      }

   }

   return next_read;
}

Read FastaGzParser::get_next_read()
{
   int MAXLINE = 1000;
   char tmp[MAXLINE];
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

      gzgets(infile, tmp, MAXLINE);
      tmp[strlen(tmp)-1] = '\0';
      line = tmp;

      while(line[0] != '>' && !gzeof(infile))
      {
         if (line[0] != '>')  {
            seq += line;
         }
         gzgets(infile, tmp, MAXLINE);
         tmp[strlen(tmp)-1] = '\0';
         line = tmp;
      }

      if (gzeof(infile)) {
         seq += line;
      }

      if (line[0] == '>')  {
         next_name = line.substr(1);
      }

      if ((int)seq.find('N') == -1)  {
         valid_read = 1;
      }  else if (gzeof(infile))  {
         one_read_left = false;
         break;
      }

      current_read.seq = seq;
      seq = "";

      if (gzeof(infile)) {
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

FastqGzParser::FastqGzParser(const std::string &inputfile)
{
   int MAXLINE = 1000;
   char tmp[MAXLINE];
   infile = gzopen(inputfile.c_str(), "rb");
   std::string line_three, quality_scores;

   one_read_left = 0;

   bool valid_read = 0;

   while (!valid_read && !gzeof(infile))  {
      //getline(infile, current_read.name);
      //getline(infile, current_read.seq);
      //getline(infile, line_three);
      //getline(infile, quality_scores);
      gzgets(infile, tmp, MAXLINE);
      tmp[strlen(tmp)-1] = '\0';
      current_read.name = tmp;
      gzgets(infile, tmp, MAXLINE);
      tmp[strlen(tmp)-1] = '\0';
      current_read.seq = tmp;
      gzgets(infile, tmp, MAXLINE);
      tmp[strlen(tmp)-1] = '\0';
      line_three = tmp;
      gzgets(infile, tmp, MAXLINE);
      tmp[strlen(tmp)-1] = '\0';
      quality_scores = tmp;

      //std::cout << current_read.name[0] << std::endl;

      assert(current_read.name[0] == '@');
      assert(line_three[0] == '+');
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
         infile.close();
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

Read FastqGzParser::get_next_read()
{
   int MAXLINE = 1000;
   char tmp[MAXLINE];
   Read next_read = current_read;
   std::string line_three, quality_scores;

   bool valid_read = 0;

   if (gzeof(infile)) {
      one_read_left = 0;
      return next_read;
   }

   while (!valid_read && !gzeof(infile))  {

      //getline(infile, current_read.name);
      gzgets(infile, tmp, MAXLINE);
      tmp[strlen(tmp)-1] = '\0';
      current_read.name = tmp;

      //if (gzeof(infile))  {
      //   one_read_left = 0;
      //   return next_read;
      //}

      //getline(infile, current_read.seq);
      //getline(infile, line_three);
      //getline(infile, quality_scores);
      gzgets(infile, tmp, MAXLINE);
      tmp[strlen(tmp)-1] = '\0';
      current_read.seq = tmp;
      gzgets(infile, tmp, MAXLINE);
      tmp[strlen(tmp)-1] = '\0';
      line_three = tmp;
      gzgets(infile, tmp, MAXLINE);
      tmp[strlen(tmp)-1] = '\0';
      quality_scores = tmp;

      assert(current_read.name[0] == '@');
      assert(line_three[0] == '+');
      assert(quality_scores.length() == current_read.seq.length());

      current_read.name = current_read.name.substr(1);

      if ((int)current_read.seq.find('N') == -1)  {
         valid_read = 1;
      }

      if (gzeof(infile) && valid_read) {
         one_read_left = 1;
      }
   }

   return next_read;
}


int main()
{
   IParser* parser = IParser::get_parser("test.fasta.gz");
   
   while(!parser->is_complete())  {
      Read next_read = parser->get_next_read();
      std::cout << next_read.name << ": " << next_read.seq << std::endl;
   }

   delete parser;

   IParser* parser2 = IParser::get_parser("test.fq.gz");
   while (!parser2->is_complete()) {
      Read next_read = parser2->get_next_read();
      std::cout << next_read.name << ": " << next_read.seq << std::endl;
   }

   delete parser2;

   return 0;
}
