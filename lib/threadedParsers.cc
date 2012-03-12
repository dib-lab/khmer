#include "threadedParsers.hh"

using namespace khmer;
using namespace khmer:: threaded_parsers;

/* This function returns the proper parser factory for a given
 * file type
 */
ThreadedIParserFactory* ThreadedIParserFactory::get_parser(const std::string &inputfile, long int chunkSize)
{
   std::string filename(inputfile);
   int found = filename.find_last_of(".");

   std::string type = filename.substr(found+1);

   if (type == "fq" || type == "fastq") {
      return new ThreadedFastqParserFactory(inputfile, chunkSize);
   } else if (type == "fa" || type == "fasta") {
      return new ThreadedFastaParserFactory(inputfile, chunkSize);
   } else {
      return new ThreadedFastaParserFactory(inputfile, chunkSize);
   }
}

/* Constructer for ThreadedFastaParserFactory */
ThreadedFastaParserFactory::ThreadedFastaParserFactory(
    const std::string &inputfile, long int size)
{
    /* Save off the filename and chunk size */
    filename = inputfile;
    chunkSize = size;
    curPos = 0;

    /* Get the file size */
    std::ifstream myFile;
    myFile.open(filename.c_str(), std::ifstream::in);
    myFile.seekg (0, std::ios::end);

    /* Subtract off one so that if there is only one
     * byte left (which is the case since we call unget)
     * we will still indicate that we have reached the end
     * of the file */
    fileSize = (long int)myFile.tellg() - 1;
    myFile.close();
}

ThreadedIParser* ThreadedFastaParserFactory::get_next_parser()
{
    /* Open the file */
    std::ifstream myFile;
    myFile.open(filename.c_str(), std::ifstream::in);
    long int myCurPos, endPos;

    do {
        /* Save off the current curPos so that we can use it in the compare
         * and swap later, to make this call threadsafe */
        myCurPos = curPos;

        /* Start the end position off at the current position + the chunk size */
        endPos = myCurPos + chunkSize;

        if (endPos > fileSize)
            endPos = fileSize;

        /* Read until we see
         * a '>', which indicates the start of the next read, or we reach EOF
         */
        myFile.seekg(endPos, std::ios::beg);

        while (!myFile.eof() && myFile.get() != '>');

        /* We get here as soon as we read a '>' */
        /* Put the '>' back */
        myFile.unget();
        endPos = myFile.tellg();

        if (endPos == -1)
            endPos = fileSize;

    } while(!__sync_bool_compare_and_swap(&curPos, myCurPos, endPos));
    

    myFile.close();
    return new ThreadedFastaParser(filename, myCurPos, endPos);
}


/* Constructer for ThreadedFastqParserFactory */
ThreadedFastqParserFactory::ThreadedFastqParserFactory(
    const std::string &inputfile, long int size)
{
    /* Save off the filename and chunk size */
    filename = inputfile;
    chunkSize = size;
    curPos = 0;

    /* Get the file size */
    std::ifstream myFile;
    myFile.open(filename.c_str(), std::ifstream::in);
    myFile.seekg (0, std::ios::end);
    /* Subtract off one so that if there is only one
     * byte left (which is the case since we call unget)
     * we will still indicate that we have reached the end
     * of the file */
    fileSize = (long int)myFile.tellg() - 1;
    myFile.close();
}

ThreadedIParser* ThreadedFastqParserFactory::get_next_parser()
{
    /* Open the file */
    std::ifstream myFile;
    myFile.open(filename.c_str(), std::ifstream::in);
    long int myCurPos, endPos;

    do {
        /* Save off the current curPos so that we can use it in the compare
         * and swap later, to make this call threadsafe */
        myCurPos = curPos;

        /* Start the end position off at the current position + the chunk size */
        endPos = myCurPos + chunkSize;

        if (endPos > fileSize)
            endPos = fileSize;

        /* Read until we see
         * a '@', which indicates the start of the next read, or we reach EOF
         * We have to be careful - a '@' can show up in the quality scores too
         * so watch for that
         */
        myFile.seekg(endPos, std::ios::beg);

        char prevChar = myFile.get();
        while (!myFile.eof())
        {
            char curChar = myFile.get();
            if (curChar == '@' && prevChar == '\n')
                break;
            prevChar = curChar;
        }

        /* We get here as soon as we read a '@' */
        /* Put the '@' back */
        if (!myFile.eof())
            myFile.unget();
        endPos = myFile.tellg();

        if (endPos == -1)
            endPos = fileSize;

    } while(!__sync_bool_compare_and_swap(&curPos, myCurPos, endPos));
    
    myFile.close();
    return new ThreadedFastqParser(filename, myCurPos, endPos);
}



ThreadedFastaParser::ThreadedFastaParser(const std::string &inputfile, 
    long int startPos, long int end) : 
                         infile(inputfile.c_str())
{

    assert(infile.is_open());

    /* Lay the start pointer down in the right location */
    infile.seekg(startPos, std::ios::beg);

    /* Set the end ptr */
    endPos = end;
}

Read ThreadedFastaParser::get_next_read()
{
   std::string line = "", seq = "";
   Read read;

   bool valid_read = 0;

   while (!valid_read)  {
      getline(infile, line);

      if (!infile.eof( ))
	assert(line[0] == '>');
      read.name = line.substr(1);

      /* Get all the lines that are sequence data */
      while(infile.get() != '>' && !infile.eof())
      {
         infile.unget();
         getline(infile, line);
         seq += line;
      }
      /* Put the character back if we are not at EOF */
      if (!infile.eof())
          infile.unget();

      if ((int)seq.find_first_of("Nn") == -1)  {
          valid_read = 1;
      }

      read.seq = seq;
      seq = "";

   }

   return read;
}

ThreadedFastqParser::ThreadedFastqParser(const std::string &inputfile, 
                long int startPos, long int end) :
                         infile(inputfile.c_str())
{
    assert(infile.is_open());

    /* Lay the start pointer down in the right location */
    infile.seekg(startPos, std::ios::beg);

    /* Set the end ptr */
    endPos = end;
}

Read ThreadedFastqParser::get_next_read()
{
   Read read;
   std::string line_three, quality_scores;

   bool valid_read = 0;

   while (!valid_read && !infile.eof())  {

      getline(infile, read.name);
      getline(infile, read.seq);
      getline(infile, line_three);
      getline(infile, quality_scores);
   
      assert(read.name[0] == '@');
      assert(line_three[0] == '+' || line_three[0] == '#');
      assert(quality_scores.length() == read.seq.length());

      read.name = read.name.substr(1);

      if ((int)read.seq.find_first_of("Nn") == -1)  {
         valid_read = 1;
      }
   }

   return read;
}

// vim: set sts=2 sw=2:
