/*
  This line intentionally left blank to bother mr-c

  TODO:
  - Free allocated memory
    - subalign will have to return an alignment object
    - or take an alignment object to fill out (*)
  - Document python return tuple (named tuple?)
  - Decide whether probabilities should be modifible without recompiling
  - Implement trusted kmer detection (cutoff vs probability)
  - Model error vs trusted kmers
      - Another state? (probably not...)
      - Another emission probability?
      - Have to make sure they aren't less probable than a mismatch, or insert/mismatch/insert
 */

#ifndef READ_ALIGNER_HH
#define READ_ALIGNER_HH

#include "khmer.hh"
#include "counting.hh"

#include <limits>
#include <algorithm>
#include <set>
#include <vector>
#include <queue>
#include <memory>

#define READ_ALIGNER_DEBUG 0

namespace khmer {

  enum State { MATCH, INSERT_READ, INSERT_GRAPH };
  enum Nucl {A, C, G, T};
  static const char nucl_lookup[4] = {'A', 'C', 'G', 'T'};
  
  struct AlignmentNode {
    AlignmentNode* prev;
    Nucl base;
    unsigned int seq_idx;
    State state;
    HashIntoType fwd_hash;
    HashIntoType rc_hash;

    double score;
    double f_score;
    double h_score;
    bool trusted;

    unsigned int length;
    
    AlignmentNode(AlignmentNode* _prev, Nucl _emission, int _seq_idx, State _state, HashIntoType _fwd_hash, HashIntoType _rc_hash, unsigned int _length)
      :prev(_prev), base(_emission), seq_idx(_seq_idx), state(_state), fwd_hash(_fwd_hash), rc_hash(_rc_hash), length(_length) {}
    
    bool operator== (const AlignmentNode& rhs) const {
      return (seq_idx == rhs.seq_idx) && (state == rhs.state) &&
	uniqify_rc(fwd_hash, rc_hash) == uniqify_rc(rhs.fwd_hash, rhs.rc_hash);
    }

    bool operator< (const AlignmentNode& rhs) const {
      return f_score < rhs.f_score;
    }
  };

  class AlignmentNodeCompare {
  public:
    bool operator()(AlignmentNode* o1, AlignmentNode* o2) {
      if (o1->f_score < o2->f_score) {
	return true;
      } else {
	return false;
      }
    }
  };
  
  typedef std::priority_queue<AlignmentNode*, std::vector<AlignmentNode*>, AlignmentNodeCompare> NodeHeap;


  struct ScoringMatrix {
    const double trusted_match;
    const double trusted_mismatch;
    const double untrusted_match;
    const double untrusted_mismatch;

    const double* tsc;

    ScoringMatrix(double trusted_match, double trusted_mismatch, double untrusted_match, double untrusted_mismatch, double* trans): trusted_match(trusted_match), trusted_mismatch(trusted_mismatch), untrusted_match(untrusted_match), untrusted_mismatch(untrusted_mismatch), tsc(trans) {}
  };

    
  struct Alignment {
    std::string graph_alignment;
    std::string read_alignment;
    std::string trusted;
    double score;
    bool truncated;
  };

  // Constants for state transitions
  enum Transition { MM, MI, MD, IM, II, DM, DD, disallowed };
  // log probabilities for state transitions
  static double trans_default[] = { log2(.9), log2(.05), log2(.05),
				    log2(.95), log2(.05),
				    log2(.95), log2(.05)};
  
  
  class ReadAligner {
  private:

    Alignment* ExtractAlignment(AlignmentNode*, bool forward, const std::string&);
    void Enumerate(NodeHeap&, std::vector<AlignmentNode*>& all_nodes, AlignmentNode*, bool, const std::string&);
    Alignment* Subalign(AlignmentNode*, unsigned int, bool, const std::string&);
    
    // These variables are required to use the _revhash and hash macros
    // might as well just compute them once
    const HashIntoType bitmask;
    const unsigned int rc_left_shift;

    khmer::CountingHash* m_ch;
    ScoringMatrix m_sm;

    unsigned int m_trusted_cutoff;
    double m_bits_theta;
    
    HashIntoType comp_bitmask(WordLength k) {
      HashIntoType ret = 0;
      for (unsigned int i = 0; i < k; i++) {
	ret = (ret << 2) | 3;
      }
      return ret;
    }
    
  public:
    Alignment* Align(const std::string&);

    ReadAligner(khmer::CountingHash* ch, unsigned int trusted_cutoff, double bits_theta)
      : bitmask(comp_bitmask(ch->ksize())), rc_left_shift(ch->ksize() * 2 - 2), m_ch(ch), m_sm(log2(.985), log2(.01), log2(.004), log2(.001), trans_default), m_trusted_cutoff(trusted_cutoff), m_bits_theta(bits_theta) {}
  };
}

#endif // READ_ALIGNER_HH
