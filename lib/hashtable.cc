#include <iostream>
#include <list>

#include "khmer.hh"
#include "hashtable.hh"

using namespace khmer;
using namespace std;

MinMaxTable * Hashtable::fasta_file_to_minmax(const std::string &inputfile,
					      unsigned int total_reads,
					      ReadMaskTable * readmask)
{
   string line;
   ifstream infile(inputfile.c_str());
   int isRead = 0;
   string name;
   string seq;
   unsigned int read_num = 0;

   MinMaxTable * mmt = new MinMaxTable(total_reads);

   if (infile.is_open()) {
     while(!infile.eof()) {
       getline(infile, line);
       if (line.length() == 0) {
	 break;
       }

       if (isRead) {
	 bool valid_read = true;
	 seq = line;
	 if (!readmask || readmask->get(read_num)) {

	   for (unsigned int i = 0; i < seq.length(); i++)  {
	     if (!is_valid_dna(seq[i])) {
	       valid_read = false;
	       break;
	     }
	   }

	   if (valid_read) {
	     BoundedCounterType minval = get_min_count(seq);
	     BoundedCounterType maxval = get_max_count(seq);

	     mmt->add_min(read_num, minval);
	     mmt->add_max(read_num, maxval);
	   }
	   name.clear();
	   seq.clear();
	 }

	 read_num += 1;
       }
       else {
	 name = line.substr(1, line.length()-1);
       }

       isRead = isRead ? 0 : 1;
     }
   }
  
   infile.close();

   return mmt;
}

//
// filter_fasta_file: filter and trims a FASTA file into a new one
//

ReadMaskTable * Hashtable::filter_fasta_file_max(const std::string &inputfile,
						 MinMaxTable &minmax,
						 BoundedCounterType threshold,
						 ReadMaskTable * old_readmask)
{
   string line;
   ifstream infile(inputfile.c_str());
   int isRead = 0;
   string name;
   string seq;
   unsigned int read_num = 0;
   ReadMaskTable * readmask = new ReadMaskTable(minmax.get_tablesize());

   if (old_readmask) {
     readmask->merge(*old_readmask);
   }

   if (infile.is_open()) {
     while(!infile.eof()) {
       getline(infile, line);
       if (line.length() == 0) {
	 break;
       }

       if (isRead) {
	 seq = line;
	 if (readmask->get(read_num)) {
	   BoundedCounterType maxval = minmax.get_max(read_num);

	   if (maxval < threshold) {
	     readmask->set(read_num, false);
	   }
	   name.clear();
	   seq.clear();
	 }

	 read_num += 1;
       }
       else {
	 name = line.substr(1, line.length()-1);
       }

       isRead = isRead ? 0 : 1;
     }
   }
  
   infile.close();

   return readmask;
}

unsigned int khmer::output_filtered_fasta_file(const std::string &inputfile,
					       const std::string &outputfile,
					       ReadMaskTable * readmask)
{
   string line;
   ifstream infile(inputfile.c_str());
   ofstream outfile;
   outfile.open(outputfile.c_str());
   int isRead = 0;
   string name;
   string seq;
   unsigned int n_kept = 0;
   unsigned int read_num = 0;

   if (infile.is_open()) {
     while(!infile.eof()) {
       getline(infile, line);
       if (line.length() == 0) {
	 break;
       }

       if (isRead) {
	 seq = line;

	 if (readmask->get(read_num)) {
	     outfile << ">" << name << endl;
	     outfile << seq << endl;

	     n_kept++;
	 }
	 name.clear();
	 seq.clear();

	 read_num++;
       } else {
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

void Hashtable::consume_fasta(const std::string &filename,
			      unsigned int &total_reads,
			      unsigned int &n_consumed,
			      HashIntoType lower_bound,
			      HashIntoType upper_bound,
			      ReadMaskTable ** orig_readmask,
			      bool update_readmask)
{
   total_reads = 0;
   n_consumed = 0;

   string line;
   ifstream infile(filename.c_str());

   if (!infile.is_open())  {
     return;
   }

   string currName = "";
   string currSeq = "";

   //
   // readmask stuff: were we given one? do we want to update it?
   // 

   ReadMaskTable * readmask = NULL;
   std::list<unsigned int> masklist;

   if (orig_readmask && *orig_readmask) {
     readmask = *orig_readmask;
   }

   //
   // iterate through the FASTA file & consume the reads.
   //

   while(1)  {
     getline(infile, line);

     if (line[0] == '>' || infile.eof())  {

       // do we have a sequence to process?
       if (currSeq != "")  {

	 // do we want to process it?
	 if (!readmask || readmask->get(total_reads)) {

	   // yep! process.

	   unsigned int this_n_consumed;
	   this_n_consumed = checkAndProcessRead(currSeq,
						 lower_bound,
						 upper_bound);

	   // was this an invalid sequence -> mark as bad?
	   if (this_n_consumed == 0 && update_readmask) {
	     if (readmask) {
	       readmask->set(total_reads, false);
	     } else {
	       masklist.push_back(total_reads);
	     }
	   } else {		// nope -- count it!
	     n_consumed += this_n_consumed;
	   }
	 }
	       
	 // reset the sequence info, increment read number, etc.
	 currSeq = "";
	 total_reads++;
	 if (total_reads % 10000 == 0) { // @CTB remove me!
	   cout << total_reads << endl;
	 }
       }

       // new sequence => new sequence name
       if (line[0] == '>') {
	 currName = line.substr(1, line.length()-1);
       }
     }
     else  {			// additional line for sequence
       currSeq += line;
     }
     
     // @ end of file? break out.
     if (infile.eof()) {
       break;
     }
   }

   infile.close();

   //
   // We've either updated the readmask in place, OR we need to create a
   // new one.
   //

   if (orig_readmask && update_readmask && readmask == NULL) {
     // allocate, fill in from masklist
     readmask = new ReadMaskTable(total_reads);

     list<unsigned int>::const_iterator it;
     for(it = masklist.begin(); it != masklist.end(); ++it) {
       readmask->set(*it, false);
     }
     *orig_readmask = readmask;
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
    }
    n_consumed++;
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
      }
      n_consumed++;
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

