//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE. Contact: ctb@msu.edu
//

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

namespace khmer
{

enum State { MATCH, INSERT_READ, INSERT_GRAPH,
             MATCH_UNTRUSTED, INSERT_READ_UNTRUSTED, INSERT_GRAPH_UNTRUSTED
           };

// Constants for state transitions
enum Transition { MM, MIr, MIg, MMu, MIru, MIgu,
                  IrM, IrIr, IrMu, IrIru,
                  IgM, IgIg, IgMu, IgIgu,
                  MuM, MuIr, MuIg, MuMu, MuIru, MuIgu,
                  IruM, IruIr, IruMu, IruIru,
                  IguM, IguIg, IguMu, IguIgu,
                  disallowed
                };

/*
Ig_t-Ig_t       0.2294619
Ig_t-Ig_u       0.0021453
Ig_t-M_t        0.7611255
Ig_t-M_u        0.0072673
Ig_u-Ig_t       0.0431328
Ig_u-Ig_u       0.1821200
Ig_u-M_t        0.1384551
Ig_u-M_u        0.6362921
Ir_t-Ir_t       0.4647955
Ir_t-Ir_u       0.0096792
Ir_t-M_t        0.5196194
Ir_t-M_u        0.0059060
Ir_u-Ir_t       0.0036995
Ir_u-Ir_u       0.5885548
Ir_u-M_t        0.1434529
Ir_u-M_u        0.2642928
M_t-Ig_t        0.0000334
M_t-Ig_u        0.0000003
M_t-Ir_t        0.0000735
M_t-Ir_u        0.0000017
M_t-M_t 0.9848843
M_t-M_u 0.0150068
M_u-Ig_t        0.0001836
M_u-Ig_u        0.0004173
M_u-Ir_t        0.0000262
M_u-Ir_u        0.0033370
M_u-M_t 0.0799009
M_u-M_u 0.9161349

*/
// log probabilities for state transitions
static double trans_default[] = { log2(0.9848843), log2(0.0000735), log2(0.0000334), log2(0.0150068), log2(0.0000017), log2(0.0000003),  // M_t
                                  log2(0.5196194), log2(0.4647955), log2(0.0059060), log2(0.0096792),                        // Ir_t
                                  log2(0.7611255), log2(0.2294619), log2(0.0072673), log2(0.0021453),                        // Ig_t
                                  log2(0.0799009), log2(0.0000262), log2(0.0001836), log2(0.9161349), log2(0.0033370), log2(0.0004173),  // M_u
                                  log2(0.1434529), log2(0.0036995), log2(0.2642928), log2(0.5885548),                        // Ir_u
                                  log2(0.1384551), log2(0.0431328), log2(0.6362921), log2(0.1821200),                        // Ig_u
                                };
/*{ log2(.80), log2(.045), log2(.045), log2(.06), log2(.025), log2(.025),
                                  log2(.875), log2(.045), log2(.055), log2(.025),
                                  log2(.875), log2(.045), log2(.055), log2(.025),
			    log2(.80), log2(.045), log2(.045), log2(.06), log2(.025), log2(.025),
                                  log2(.875), log2(.045), log2(.055), log2(.025),
                                  log2(.875), log2(.045), log2(.055), log2(.025),
};*/

enum Nucl {A, C, G, T};
static const char nucl_lookup[4] = {'A', 'C', 'G', 'T'};
static const double background_prob = 0;//log2(.99);

struct AlignmentNode {
    AlignmentNode* prev;
    Nucl base;
    size_t seq_idx;
    State state;
    Transition trans;
    HashIntoType fwd_hash;
    HashIntoType rc_hash;

    double score;
    double f_score;
    double h_score;
    bool trusted;

    size_t num_indels;

    size_t length;

    AlignmentNode(AlignmentNode* _prev, Nucl _emission, size_t _seq_idx,
                  State _state, Transition _trans, HashIntoType _fwd_hash,
                  HashIntoType _rc_hash, size_t _length)
        :prev(_prev), base(_emission), seq_idx(_seq_idx),
         state(_state), trans(_trans), fwd_hash(_fwd_hash),
         rc_hash(_rc_hash), score(0), f_score(0), h_score(0), trusted(false),
         num_indels(0), length(_length) {}

    bool operator== (const AlignmentNode& rhs) const
    {
        return (seq_idx == rhs.seq_idx) && (state == rhs.state) &&
               uniqify_rc(fwd_hash, rc_hash) == uniqify_rc(rhs.fwd_hash, rhs.rc_hash)
               && trans == rhs.trans;
    }

    bool operator< (const AlignmentNode& rhs) const
    {
        return f_score < rhs.f_score;
    }
};

class AlignmentNodeCompare
{
public:
    bool operator()(AlignmentNode* o1, AlignmentNode* o2)
    {
        if (o1->f_score < o2->f_score) {
            return true;
        } else {
            return false;
        }
    }
};

typedef std::priority_queue<AlignmentNode*,
        std::vector<AlignmentNode*>,
        AlignmentNodeCompare> NodeHeap;

struct ScoringMatrix {
    const double trusted_match;
    const double trusted_mismatch;
    const double untrusted_match;
    const double untrusted_mismatch;

    const double* tsc;

    ScoringMatrix(double trusted_match, double trusted_mismatch,
                  double untrusted_match, double untrusted_mismatch,
                  double* trans)
        : trusted_match(trusted_match), trusted_mismatch(trusted_mismatch),
          untrusted_match(untrusted_match),
          untrusted_mismatch(untrusted_mismatch), tsc(trans) {}
};


struct Alignment {
    std::string graph_alignment;
    std::string read_alignment;
    std::string trusted;
    double score;
    bool truncated;
};


class ReadAligner
{
private:

    Alignment* ExtractAlignment(AlignmentNode*,
                                bool forward, const std::string&);

    void Enumerate(NodeHeap&, std::vector<AlignmentNode*>& all_nodes,
                   AlignmentNode*, bool, const std::string&);
    Alignment* Subalign(AlignmentNode*, size_t, bool, const std::string&);

#if READ_ALIGNER_DEBUG
    void WriteNode(AlignmentNode* curr);
#endif

    // These variables are required to use the _revhash and hash macros
    // might as well just compute them once
    const HashIntoType bitmask;
    const size_t rc_left_shift;

    khmer::CountingHash* m_ch;
    ScoringMatrix m_sm;

    size_t m_trusted_cutoff;
    double m_bits_theta;

    HashIntoType comp_bitmask(WordLength k)
    {
        HashIntoType ret = 0;
        for (size_t i = 0; i < k; i++) {
            ret = (ret << 2) | 3;
        }
        return ret;
    }

public:
    Alignment* Align(const std::string&);

    ReadAligner(khmer::CountingHash* ch,
                BoundedCounterType trusted_cutoff, double bits_theta)
        : bitmask(comp_bitmask(ch->ksize())),
          rc_left_shift(ch->ksize() * 2 - 2),
          m_ch(ch), m_sm(
              log2(.955), log2(.04), log2(.004),
              log2(.001), trans_default),
          m_trusted_cutoff(trusted_cutoff),
          m_bits_theta(bits_theta)
    {
#if READ_ALIGNER_DEBUG
        std::cerr << "Trusted cutoff: " << m_trusted_cutoff
                  << " bits theta: " << bits_theta
                  << " trusted match: " << m_sm.trusted_match
                  << " untrusted match: " << m_sm.untrusted_match
                  << " trusted mismatch: " << m_sm.trusted_mismatch
                  << " untrusted mismatch: " << m_sm.untrusted_mismatch
                  << std::endl;
#endif
    }
};
}

#endif // READ_ALIGNER_HH
