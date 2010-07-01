#include "khmer.hh"
#include "hashtable.hh"
#include <iostream>

using namespace khmer;
using namespace std;

//
// filter_fasta_file: filter and trims a FASTA file into a new one
//

unsigned int Hashtable::filter_fasta_file(const std::string &inputfile,
					  const std::string &outputfile, 
					  int minLength, 
					  int threshold)
{
   string line;
   ifstream infile(inputfile.c_str());
   ofstream outfile;
   outfile.open(outputfile.c_str());
   int isRead = 0;
   string name;
   string seq;
   unsigned int n_kept = 0;

   if (infile.is_open())
   {
      while(!infile.eof())
      {
         getline(infile, line);
         if (line.length() == 0) {
            break;
	 }

         if (isRead)
         {
	    bool valid_read = true;
	    seq = line;

	    for (unsigned int i = 0; i < seq.length(); i++)  {
	      if (!is_valid_dna(seq[i])) {
		valid_read = false;
		break;
	      }
	    }

            if (valid_read && get_max_count(seq) >= minLength) {
               outfile << ">" << name << endl;
               outfile << seq << endl;

	       n_kept++;
            }
            name.clear();
            seq.clear();
         }
         else
         {
            name = line.substr(1, line.length()-1);
         }

         isRead = isRead? 0 : 1;
      }
   }
  
   infile.close();
   outfile.close();

   return n_kept;
}

//
// checkAndProcessRead: checks for non-ACGT characters
//

unsigned int Hashtable::checkAndProcessRead(const std::string &read,
                                            HashIntoType lower_bound,
                                            HashIntoType upper_bound)
{
   unsigned int i;
   bool is_valid = true;

   for (i = 0; i < read.length(); i++)  {
     if (!is_valid_dna(read[i])) {
         is_valid = false;
         break;
      }
   }

   if (is_valid) {
      return consume_string(read, lower_bound, upper_bound);
   }
   else  {
      return 0;
   }
}

//
// consume_fasta: consume a FASTA file of reads
//
unsigned int Hashtable::consume_fasta(const std::string &filename,
                              HashIntoType lower_bound,
                              HashIntoType upper_bound)
{
   string line;
   ifstream infile(filename.c_str());

   unsigned int n_consumed=0, n=0;

   string currName = "";
   string currSeq = "";

   if (infile.is_open())  {
      while(!infile.eof())  {
         getline(infile, line);

         if (line[0] == '>')  {
            if (currSeq != "")  {
               n++;
               if (n % 10000 == 0) { // @CTB remove me
                  cout << n << endl;
	       }

               n_consumed += checkAndProcessRead(currSeq, lower_bound,
						 upper_bound);
               currSeq = "";
            }
            currName = line.substr(1, line.length()-1);
         }
         else  {
            currSeq += line;
         }
      }
   }

   if (currSeq != "")  {
      n_consumed += checkAndProcessRead(currSeq, lower_bound, upper_bound);
   }

   infile.close();

   return n_consumed;
}

//
// consume_string: run through every k-mer in the given string, & hash it.
//

unsigned int Hashtable::consume_string(const std::string &s,
				       HashIntoType lower_bound,
				       HashIntoType upper_bound)
{
  const char * sp = s.c_str();
  unsigned int length = s.length();
  unsigned int n_consumed = 0;

  unsigned int mask = 0;
  for (unsigned int i = 0; i < _ksize; i++) {
    mask = mask << 2;
    mask |= 3;
  }

  HashIntoType h = 0, r = 0;
  bool bounded = true;

  if (lower_bound == upper_bound && upper_bound == 0) {
    bounded = false;
  }
  
  HashIntoType bin = _hash(sp, _ksize, &h, &r);

  if (!bounded || (bin >= lower_bound && bin < upper_bound)) {
    bin = bin % _tablesize;
    if (_counts[bin] != MAX_COUNT) {
      _counts[bin]++;
      n_consumed++;
    }
  }

  for (unsigned int i = _ksize; i < length; i++) {
    // left-shift the previous hash over
    h = h << 2;

    // 'or' in the current nt
    h |= twobit_repr(sp[i]);

    // mask off the 2 bits we shifted over.
    h &= mask;

    // now handle reverse complement
    r = r >> 2;
    r |= (twobit_comp(sp[i]) << (_ksize*2 - 2));

    bin = (h < r) ? h : r;

    if (!bounded || (bin >= lower_bound && bin < upper_bound)) {
      bin = bin % _tablesize;
      if (_counts[bin] != MAX_COUNT) {
	_counts[bin]++;
	n_consumed++;
      }
    }
  }
  return n_consumed;
}


BoundedCounterType Hashtable::get_min_count(const std::string &s,
					    HashIntoType lower_bound,
					    HashIntoType upper_bound)
{
  const unsigned int length = s.length();
  const char * sp = s.c_str();
  BoundedCounterType min_count = 255, count;

  unsigned int mask = 0;
  for (unsigned int i = 0; i < (unsigned int) _ksize; i++) {
    mask = mask << 2;
    mask |= 3;
  }

  HashIntoType h = 0, r = 0;
  bool bounded = true;

  if (lower_bound == upper_bound && upper_bound == 0) {
    bounded = false;
  }

  HashIntoType bin;
  
  bin = _hash(sp, _ksize, &h, &r);
  if (!bounded || (bin >= lower_bound && bin < upper_bound)) {
    min_count = this->get_count(bin);
  }

  for (unsigned int i = _ksize; i < length; i++) {
    // left-shift the previous hash over
    h = h << 2;

    // 'or' in the current nt
    h |= twobit_repr(sp[i]);

    // mask off the 2 bits we shifted over.
    h &= mask;

    // now handle reverse complement
    r = r >> 2;
    r |= (twobit_comp(sp[i]) << (_ksize*2 - 2));

    bin = (h < r) ? h : r;

    if (!bounded || (bin >= lower_bound && bin < upper_bound)) {
      count = this->get_count(bin);
    
      if (count < min_count) {
	min_count = count;
      }
    }
  }
  return min_count;
}

BoundedCounterType Hashtable::get_max_count(const std::string &s,
					    HashIntoType lower_bound,
					    HashIntoType upper_bound)
{
  const unsigned int length = s.length();
  const char * sp = s.c_str();
  BoundedCounterType max_count = 0, count;

  unsigned int mask = 0;
  for (unsigned int i = 0; i < (unsigned int) _ksize; i++) {
    mask = mask << 2;
    mask |= 3;
  }

  HashIntoType h = 0, r = 0;
  bool bounded = true;

  if (lower_bound == upper_bound && upper_bound == 0) {
    bounded = false;
  }

  HashIntoType bin = _hash(sp, _ksize, &h, &r);
  if (!bounded || (bin >= lower_bound && bin < upper_bound)) {
    max_count = this->get_count(bin);
  }

  for (unsigned int i = _ksize; i < length; i++) {
    // left-shift the previous hash over
    h = h << 2;

    // 'or' in the current nt
    h |= twobit_repr(sp[i]);

    // mask off the 2 bits we shifted over.
    h &= mask;

    // now handle reverse complement
    r = r >> 2;
    r |= (twobit_comp(sp[i]) << (_ksize*2-2));

    bin = (h < r) ? h : r;
    if (!bounded || (bin >= lower_bound && bin < upper_bound)) {
      count = this->get_count(bin);

      if (count > max_count) {
	max_count = count;
      }
    }
  }
  return max_count;
}

