#ifndef READ_ALIGNER_HH
#define READ_ALIGNER_HH

#include "khmer.hh"
#include "counting.hh"

#include <algorithm>
#include <set>
#include <vector>
#include <queue>

namespace khmer {

  enum State { match, insertion, deletion };
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

    AlignmentNode(AlignmentNode* _prev, Nucl _emission, int _seq_idx, State _state, HashIntoType _fwd_hash, HashIntoType _rc_hash)
      :prev(_prev), base(_emission), seq_idx(_seq_idx), state(_state), fwd_hash(_fwd_hash), rc_hash(_rc_hash) {}
    
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
    const double match;
    const double mismatch;
    const double* tsc;

    ScoringMatrix(double match, double mismatch, double* trans): match(match), mismatch(mismatch), tsc(trans) {}
  };

    
  struct Alignment {
    std::string graph_alignment;
    std::string read_alignment;
    double score;
  };

  static const char* trans_labels[] = {"mm", "mi", "md", "im", "ii", "dm", "dd", "disallowed"};
  static const char* state_labels[] = {"match", "insert", "delete"};
  //pam 1 = 1.99 2 = 1.97
  enum Transition { MM, MI, MD, IM, II, DM, DD, disallowed };
  static double trans_default[] = { log2(.96), log2(.01), log2(.03), log2(.95), log2(.05), log2(.95), log2(.05)};
  static const ScoringMatrix default_sm(log2(.99), log2(.01), trans_default);
  
  class ReadAligner {
  private:
    static const float errorOffset = 20.0f;

    Alignment* extract_alignment(khmer::AlignmentNode*, const std::string&, const std::string&);
    void enumerate(NodeHeap&, khmer::AlignmentNode*, bool, const std::string&);
    
    khmer::CountingHash * ch;
    const ScoringMatrix* sm;
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
    AlignmentNode* subalign(AlignmentNode*, unsigned int, bool, const std::string&);
    Alignment* align(const std::string&, const std::string&, int);
    Alignment* align_test(const std::string&);

    ReadAligner(khmer::CountingHash* _ch,  int maxErrorReg=-1)
      : ch(_ch), bitmask(comp_bitmask(_ch->ksize())), rc_left_shift(_ch->ksize() * 2 - 2), sm(&default_sm) {
      maxErrorRegion = maxErrorReg;
    }
  };
}

#endif // READ_ALIGNER_HH
