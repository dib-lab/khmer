#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include "read_aligner.hh"
#include "khmer.hh"
#include "hashbits.hh"
#include "primes.hh"

using namespace khmer;

// testing Hashbits parameters 
const unsigned int ht_size = 1000000;
const unsigned int ht_count = 5;
const WordLength ksize = 30;

const unsigned int num_test_seqs = 1;
const std::string test_seqs[] { "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCTTTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTGTTGCAATCTTAACAACCTCTTTAC" };

//const std::string toalign = "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCTTTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTGTTGCAATCTTAACAACCTCTTTAC"; //perfect alignment
//const std::string toalign = "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAA"; //short perfect alignment
//const std::string toalign = "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTAGATGTTTGATTATCAA"; //short mismatch alignment
//const std::string toalign = "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTATTGATTATCAA"; //short insert alignment
//const std::string toalign = "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATTGTTTGATTATCAA"; //short delete alignment
//const std::string toalign = "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTATTGATTATCAA"; //short insert and delete alignment
//const std::string toalign = "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTATAGATTATCAA"; //short insert, delete, mismatch alignment
const std::string toalign = "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATAATTTTGCCGCTTTAACTGGGTCTAGTTTCTACTGCAAACTTTCCACCAACTAGTTTTTCTGCATCCTTTGTTGCAATCTTAACAACCTCTTTAC"; //hard alignment

int main(void) {
  Primes primetab( ht_size );
  ScoringMatrix sm = default_sm;
  std::vector<HashIntoType> ht_sizes;
  for ( unsigned int i = 0; i < ht_count; ++i )
  {
    ht_sizes.push_back( primetab.get_next_prime( ) );
  }
  CountingHash ht = CountingHash(ksize, ht_sizes);

  for(unsigned int index = 0;index < num_test_seqs;index++) {
    std::cout << "Loading test sequence " << index << ": " << test_seqs[index] << std::endl;
    ht.consume_string(test_seqs[index]);
  }

  std::cout << "Match: " << sm.match << ", mismatch: " << sm.mismatch << std::endl;
  for(unsigned int state = MM;state <= DD;state++) {
    std::cout << "Transition: " << state << ", prob: " << sm.tsc[state] << " " << pow(2, sm.tsc[state]) << std::endl;
  }

  std::cout << std::endl << test_seqs[0] << std::endl << toalign << std::endl << std::endl;

  ReadAligner aligner = ReadAligner(&ht);
  Alignment* alignment = aligner.align_test(toalign);
  //Alignment* alignment = aligner.align_test(test_seqs[0]);
  std::cout << "Bloom read:      " << test_seqs[0] << std::endl;
  std::cout << "Query read:      " << toalign << std::endl;
  std::cout << "Alignment Score: " << alignment->score << std::endl;
  std::cout << "Read alignment:  " << alignment->read_alignment << std::endl;
  std::cout << "Graph alignment: " << alignment->graph_alignment << std::endl;
}
