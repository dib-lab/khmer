#include "khmer.hh"
#include "hashtable.hh"
#include <iostream>

using namespace khmer;
using namespace std;

//
// filter_fasta_file: filter and trims a FASTA file into a new one
//

void Hashtable::filter_fasta_file(const std::string &inputfile,
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

   if (infile.is_open())
   {
      while(!infile.eof())
      {
         getline(infile, line);
         if (line.length() == 0)
            break;

         if (isRead)
         {
            seq = line;

            if (get_max_count(seq) >= minLength) {
               outfile << ">" << name << endl;
               outfile << seq << endl;
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
}


//
// consume_fasta: consume a FASTA file of reads
//

void Hashtable::consume_fasta(const std::string &filename)
{
   string line;
   ifstream infile(filename.c_str());
   int isRead = 0, n =0;

   if (infile.is_open())
   {
     while (!infile.eof())
     {
       getline(infile, line);

       if (isRead) {
         n++;
         if (n % 10000 == 0)
           cout << n << endl;

         Hashtable::consume_string(line);
       }
       
       isRead = isRead? 0 : 1;
     }
  }
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


BoundedCounterType Hashtable::get_min_count(const std::string &s)
{
  const unsigned int length = s.length();
  const char * sp = s.c_str();
  BoundedCounterType min_count, count;

  unsigned int mask = 0;
  for (unsigned int i = 0; i < (unsigned int) _ksize; i++) {
    mask = mask << 2;
    mask |= 3;
  }

  HashIntoType h = 0, r = 0;
  
  _hash(sp, _ksize, &h, &r);

  if (h < r) {
    min_count = this->get_count(h);
  } else {
    min_count = this->get_count(r);  
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

    if (h < r) {
      count = this->get_count(h);
    } else {
      count = this->get_count(r);
    }
    
    if (count < min_count) {
      min_count = count;
    }
  }
  return min_count;
}

BoundedCounterType Hashtable::get_max_count(const std::string &s)
{
  const unsigned int length = s.length();
  const char * sp = s.c_str();
  BoundedCounterType max_count, count;

  unsigned int mask = 0;
  for (unsigned int i = 0; i < (unsigned int) _ksize; i++) {
    mask = mask << 2;
    mask |= 3;
  }

  HashIntoType h = 0, r = 0;

  _hash(sp, _ksize, &h, &r);

  if (h < r) {
    max_count = this->get_count(h);
  } else {
    max_count = this->get_count(r);
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

    if (h < r) {
      count = this->get_count(h);
    } else {
      count = this->get_count(r);    
    }

    if (count > max_count) {
      max_count = count;
    }
  }
  return max_count;
}

