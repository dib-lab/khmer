#ifndef READ_ALIGNER_HH
#define READ_ALIGNER_HH

#include "khmer.hh"
#include "counting.hh"

#include <algorithm>
#include <set>
#include <vector>

namespace khmer {

  enum States { match, insert, deletion };
  enum Nucl {A, C, G, T};
  static const char nucl_lookup[4] = {'A', 'C', 'G', 'T'};
  
  struct AlignmentNode {
    AlignmentNode* prev;
    double score;
    double f_score;
    double h_score;
    unsigned int seq_idx;

    HashIntoType fwd_hash;
    HashIntoType rc_hash;
    States state;
    Nucl base;

    AlignmentNode(AlignmentNode* _prev, Nucl _emission, int _seq_idx, States _state, HashIntoType _fwd_hash, HashIntoType _rc_hash)
      :prev(_prev), base(_emission), seq_idx(_seq_idx), state(_state), fwd_hash(_fwd_hash), rc_hash(_rc_hash) {}
    
    bool operator== (const AlignmentNode& rhs) const {
      return (seq_idx == rhs.seq_idx) && (state == rhs.state) && (base == rhs.base) &&
	uniqify_rc(fwd_hash, rc_hash) == uniqify_rc(rhs.fwd_hash, rhs.rc_hash);
    }

    bool operator< (const AlignmentNode& rhs) const {
      return true;
    }
  };
  typedef std::priority_queue<AlignmentNode*> NodeHeap;


  struct ScoringMatrix {
    const double match;
    const double mismatch;
    const double gap;

    ScoringMatrix(double match, double mismatch, double gap): match(match), mismatch(mismatch), gap(gap) {}
  };

    
  struct Alignment {
    std::string alignment;
    double score;
  };

  //pam 1 = 1.99 2 = 1.97
  static ScoringMatrix sm();
  
  class ReadAligner {
  private:
    static const float errorOffset = 20.0f;
    void enumerate(const NodeHeap&, khmer::AlignmentNode*, unsigned char, const std::string&);
    
    khmer::CountingHash * ch;
    ScoringMatrix* sm;
    int maxErrorRegion;

    const HashIntoType bitmask;
    const unsigned int rc_left_shift;

    HashIntoType comp_bitmask(WordLength k) {
      HashIntoType ret = 0;
      for (unsigned int i = 0; i < k; i++) {
	ret = (ret << 2) | 3;
      }
      return ret;
    }
    
  public:
    AlignmentNode* subalign(AlignmentNode*, unsigned int, const std::string&);
    Alignment align(const std::string&, const std::string&, int);
    
    ReadAligner(khmer::CountingHash* _ch,  int maxErrorReg=-1)
      : bitmask(_ch->ksize()), rc_left_shift(_ch->ksize() * 2 - 2), ch(_ch){
      maxErrorRegion = maxErrorReg;
    }
  };
}

#endif // READ_ALIGNER_HH
