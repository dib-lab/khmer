#ifndef HASHTABLE_HH
#define HASHTABLE_HH

#include "khmer.hh"
#include "seqfuncs.hh"
#include <fstream>
#include <string>

#define MAX_COUNT 255

namespace khmer {
  typedef unsigned char HashcountType;
  class Hashtable {
  protected:
    const unsigned int _ksize;
    const unsigned long long int _tablesize;

    HashcountType * _counts;

    void _allocate_counters() {
      _counts = new HashcountType[_tablesize];
      memset(_counts, 0, _tablesize * sizeof(HashcountType));
    }

  public:
    Hashtable(unsigned int ksize, unsigned long long int tablesize) :
      _ksize(ksize), _tablesize(tablesize) {
      _allocate_counters();
    }

    ~Hashtable() {
      delete _counts; _counts = NULL;
    }

    // accessor to get 'k'
    const unsigned int ksize() const { return _ksize; }

    // accessors to get table info
    const unsigned int n_entries() const { return _tablesize; }

    void count(const char * kmer) {
      unsigned long long int bin = _hash(kmer, _ksize) % _tablesize;
      if (_counts[bin] == MAX_COUNT) { return; }
      _counts[bin]++;
    }

    void count(unsigned int khash) {
      unsigned int bin = khash % _tablesize;
      if (_counts[bin] == MAX_COUNT) { return; }
      _counts[bin]++;
    }

    // get the count for the given k-mer.
    const unsigned long long int get_count(const char * kmer) const {
      unsigned long long int bin = _hash(kmer, _ksize) % _tablesize;
      return _counts[bin];
    }

    // get the count for the given k-mer hash.
    const unsigned int get_count(unsigned int khash) const {
      unsigned int bin = khash % _tablesize;
      return _counts[bin];
    }

    // count every k-mer in the string.
    void consume_string(const std::string &s);

    // count every k-mer in the FASTA file.
    void consume_fasta(const std::string &filename);

    // filter/trim through the given FASTA file.
    void filter_fasta_file(const std::string &inputfile,
                           const std::string &outputfile,
                           int minLength, 
                           int threshold);

    // @@CTB doc
    HashcountType get_min_count(const std::string &s);
    HashcountType get_max_count(const std::string &s);
  };

  class HashtableIntersect {
  protected:
    khmer::Hashtable * _kh1;
    khmer::Hashtable * _kh2;

  public:
    HashtableIntersect(unsigned int ksize, unsigned int tablesize1, unsigned int tablesize2)
    {
      _kh1 = new Hashtable(ksize, tablesize1);
      _kh2 = new Hashtable(ksize, tablesize2);
    }

    ~HashtableIntersect()
    {
      delete _kh1;
      delete _kh2;
    }

    // count every k-mer in the string.
    void consume_string(const std::string &s)
    {
      _kh1->consume_string(s);
      _kh2->consume_string(s);
    }

    HashcountType get_min_count(const std::string &s)
    {
      HashcountType kh1Min = _kh1->get_min_count(s);
      HashcountType kh2Min = _kh2->get_min_count(s);

      if (kh1Min < kh2Min)
        return kh1Min;
      else
        return kh2Min;
    }

    HashcountType get_max_count(const std::string &s)
    {
      HashcountType kh1Max = _kh1->get_max_count(s);
      HashcountType kh2Max = _kh2->get_max_count(s);

      if (kh1Max > kh2Max)
        return kh1Max;
      else
        return kh2Max;
    }
  };
};

#endif // HASHTABLE_HH
