/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2013-2015, Michigan State University.
Copyright (C) 2015-2016, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)

Contact: khmer-project@idyll.org
*/
#include <ctype.h>
#include <algorithm>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <utility>

#include "oxli/hashtable.hh"
#include "oxli/oxli_exception.hh"
#include "oxli/read_aligner.hh"

namespace oxli
{

Alignment * _empty_alignment()
{
    Alignment* ret = new Alignment;
    ret->score = -std::numeric_limits<double>::infinity();
    ret->read_alignment = "";
    ret->graph_alignment = "";
    ret->truncated = true;
    return ret;
}

static Nucl _ch_to_nucl(char base)
{
    base = toupper(base);

    Nucl e = A;
    switch(base) {
    case 'A':
        e = A;
        break;
    case 'C':
        e = C;
        break;
    case 'G':
        e = G;
        break;
    case 'T':
    case 'U':
        e = T;
        break;
    }
    return e;
}

struct del_alignment_node_t {
    del_alignment_node_t& operator()(AlignmentNode* p)
    {
        delete p;
        return *this;
    }
};

del_alignment_node_t del_alignment_node()
{
    return del_alignment_node_t();
}

/*
  Compute the null model (random sequence) log odds probability
  for a given length
 */
double GetNull(size_t length)
{
    return log2(.25) * length + log2(1.0 / (length + 1));
}


/*
  Turn two states in to a transition, or disallowed if the
  transition isn't modelled
 */
Transition get_trans(State s1, State s2)
{
    if (s1 == MATCH) {
        if (s2 == MATCH) {
            return MM;
        } else if (s2 == INSERT_GRAPH) {
            return MIg;
        } else if (s2 == INSERT_READ) {
            return MIr;
        } else if (s2 == MATCH_UNTRUSTED) {
            return MMu;
        } else if (s2 == INSERT_GRAPH_UNTRUSTED) {
            return MIgu;
        } else if (s2 == INSERT_READ_UNTRUSTED) {
            return MIru;
        }
    } else if (s1 == INSERT_GRAPH) {
        if (s2 == MATCH) {
            return IgM;
        } else if (s2 == INSERT_GRAPH) {
            return IgIg;
        } else if (s2 == MATCH_UNTRUSTED) {
            return IgMu;
        } else if (s2 == INSERT_GRAPH_UNTRUSTED) {
            return IgIgu;
        }
    } else if (s1 == INSERT_READ) {
        if(s2 == MATCH) {
            return IrM;
        } else if (s2 == INSERT_READ) {
            return IrIr;
        } else if(s2 == MATCH_UNTRUSTED) {
            return IrMu;
        } else if (s2 == INSERT_READ_UNTRUSTED) {
            return IrIru;
        }
    } else if (s1 == MATCH_UNTRUSTED) {
        if (s2 == MATCH) {
            return MuM;
        } else if (s2 == INSERT_GRAPH) {
            return MuIg;
        } else if (s2 == INSERT_READ) {
            return MuIr;
        } else if (s2 == MATCH_UNTRUSTED) {
            return MuMu;
        } else if (s2 == INSERT_GRAPH_UNTRUSTED) {
            return MuIgu;
        } else if (s2 == INSERT_READ_UNTRUSTED) {
            return MuIru;
        }
    } else if (s1 == INSERT_GRAPH_UNTRUSTED) {
        if (s2 == MATCH) {
            return IguM;
        } else if (s2 == INSERT_GRAPH) {
            return IguIg;
        } else if (s2 == MATCH_UNTRUSTED) {
            return IguMu;
        } else if (s2 == INSERT_GRAPH_UNTRUSTED) {
            return IguIgu;
        }
    } else if (s1 == INSERT_READ_UNTRUSTED) {
        if(s2 == MATCH) {
            return IruM;
        } else if (s2 == INSERT_READ) {
            return IruIr;
        } else if(s2 == MATCH_UNTRUSTED) {
            return IruMu;
        } else if (s2 == INSERT_READ_UNTRUSTED) {
            return IruIru;
        }
    }

    return disallowed;
}

void ReadAligner::Enumerate(
    NodeHeap& open,
    std::vector<AlignmentNode*>& all_nodes,
    AlignmentNode* curr,
    bool forward,
    const std::string& seq
)
{
    size_t next_seq_idx;
    size_t remaining;
    double hcost;
    double sc;
    Transition trans;
    HashIntoType fwd = curr->fwd_hash;
    HashIntoType rc = curr->rc_hash;
    HashIntoType next_fwd, next_rc;

    AlignmentNode* next;

    if (forward) {
        next_seq_idx = curr->seq_idx + 1;
        remaining = seq.size() - next_seq_idx;
        if (next_seq_idx >= seq.size()) {
            return;
        }
    } else {
        next_seq_idx = curr->seq_idx - 1;
        remaining = next_seq_idx;
    }

    double match_sc;
    double mismatch_sc;
    int start_state;
    int end_state;

    // loop for MATCHes and INSERT_READs
    for (int i = A; i <= T; i++) {
        unsigned char next_nucl = nucl_lookup[i];

        if(forward) {
            next_fwd = next_f(fwd, next_nucl);
            next_rc = next_r(rc, next_nucl);
        } else {
            next_fwd = prev_f(fwd, next_nucl);
            next_rc = prev_r(rc, next_nucl);
        }

        HashIntoType hash = uniqify_rc(next_fwd, next_rc);
        BoundedCounterType kmerCov = m_ch->get_count(hash);

        if (kmerCov == 0) {
            continue;
        } else if (kmerCov < m_trusted_cutoff) {
            start_state = MATCH_UNTRUSTED;
            end_state = INSERT_GRAPH_UNTRUSTED;
            //match_sc = m_sm.untrusted_match;
            //mismatch_sc = m_sm.untrusted_mismatch;
            match_sc = m_sm.trusted_match;
            mismatch_sc = m_sm.trusted_mismatch;
        } else {
            start_state = MATCH;
            end_state = INSERT_GRAPH;
            match_sc = m_sm.trusted_match;
            mismatch_sc = m_sm.trusted_mismatch;
        }

        for(int next_state_iter = start_state;
                next_state_iter <= end_state;
                next_state_iter++) {

            State next_state = static_cast<State>(next_state_iter);
            trans = get_trans(curr->state, next_state);
            hcost = m_sm.tsc[get_trans(next_state, MATCH)]
                    + (m_sm.tsc[MM] + m_sm.trusted_match)
                    * ((remaining == 0) ?
                       0 : (remaining - 1));

            if(trans == disallowed) {
                continue;
            }

            if(next_state == MATCH || next_state == MATCH_UNTRUSTED) {
                if(next_nucl == seq[next_seq_idx]) {
                    sc = match_sc;
                } else {
                    sc = mismatch_sc;
                }
            } else {
                sc = background_prob;
            }

            if(next_state == MATCH || next_state == MATCH_UNTRUSTED) {
                next = new AlignmentNode(curr, (Nucl)i,
                                         next_seq_idx, (State)next_state, trans,
                                         next_fwd, next_rc, curr->length + 1);
                next->num_indels = curr->num_indels;
            } else if(next_state == INSERT_READ || next_state == INSERT_READ_UNTRUSTED) {
                next = new AlignmentNode(curr, (Nucl)i,
                                         next_seq_idx, (State)next_state, trans,
                                         curr->fwd_hash, curr->rc_hash,
                                         curr->length + 1);
                next->num_indels = curr->num_indels + 1;
            } else if(next_state == INSERT_GRAPH || next_state == INSERT_GRAPH_UNTRUSTED) {
                next = new AlignmentNode(curr, (Nucl)i,
                                         curr->seq_idx, (State)next_state, trans,
                                         next_fwd, next_rc, curr->length);
                next->num_indels = curr->num_indels + 1;
            }

            next->score = curr->score + sc + m_sm.tsc[trans];
            next->trusted = (kmerCov >= m_trusted_cutoff);
            next->cov = kmerCov;
            next->h_score = hcost;
            next->f_score = next->score + next->h_score;

            // TODO(fishjord) make max indels tunable)
            if (next->num_indels < 3
                    && next->score - GetNull(next->length) > next->length * m_bits_theta) {
                open.push(next);
                all_nodes.push_back(next);
            } else {
                delete next;
            }
        }
    }
}

#if READ_ALIGNER_DEBUG
void ReadAligner::WriteNode(AlignmentNode* curr)
{
    std::cerr << "curr: " << curr << " "
              << curr->prev << " " << " state=" << curr->state << " "
              << m_ch->unhash_dna(curr->fwd_hash) << " "
              << m_ch->unhash_dna(curr->rc_hash)
              << " cov="
              << m_ch->get_count(uniqify_rc(curr->fwd_hash, curr->rc_hash))
              << " emission=" << curr->base
              << " seqidx=" << curr->seq_idx
              << " score=" << curr->score
              << " fscore=" << curr->f_score
              << " bits_saved=" << curr->score - GetNull(curr->length)
              << std::endl;
}
#endif

Alignment* ReadAligner::Subalign(AlignmentNode* start_vert,
                                 size_t seqLen,
                                 bool forward,
                                 const std::string& seq)
{
    std::vector<AlignmentNode*> all_nodes;
    NodeHeap open;
    std::map<AlignmentNode, unsigned int> closed;
    open.push(start_vert);

    AlignmentNode* curr = NULL;
    AlignmentNode* best = NULL;
    std::map<AlignmentNode, unsigned int>::iterator tmp;

    unsigned int times_closed = 0;

    while (!open.empty()) {
        curr = open.top();

#if READ_ALIGNER_DEBUG
        WriteNode(curr);
        if (curr->prev == NULL) {
            std::cerr << "\tprev = null" << std::endl;
        } else {
            std::cerr << "\tprev = ";
            WriteNode(curr->prev);
        }
#endif
        open.pop();

        if(best == NULL ||
                (best->score - GetNull(best->length) <
                 curr->score - GetNull(curr->length))) {
            best = curr;
        }

        if (curr->seq_idx == seqLen-1 ||
                curr->seq_idx == 0) {
            best = curr;
            break;
        }

        tmp = closed.find(*curr);
        if(tmp == closed.end()) {  //Hasn't been closed yet
            //do nothing
            times_closed = 0;
        } else if (tmp->first.score > curr->score) { //Better than what we've closed
            times_closed = tmp->second;
            closed.erase(tmp);
        } else if (tmp->first.score == curr->score) { //Same as what we've closed
            times_closed = tmp->second;
            closed.erase(tmp);
        } else {
            continue;
        }

        if (times_closed > 200) {
            continue;
        }

        closed[*curr] = times_closed + 1;

        Enumerate(open, all_nodes, curr, forward, seq);
    }

    Alignment* ret = ExtractAlignment(best, forward, seq);
    std::for_each(all_nodes.begin(), all_nodes.end(), del_alignment_node());

    return ret;
}

Alignment* ReadAligner::ExtractAlignment(AlignmentNode* node,
        bool forward,
        const std::string& read)
{
    Alignment* ret = new Alignment;
    if(node == NULL) {
        ret->score = 0;
        ret->read_alignment = "";
        ret->graph_alignment = "";
        ret->trusted = "";
        ret->truncated = true;
        return ret;
    }

    if (!(node->seq_idx < read.length())) {
        delete ret;
        throw oxli_exception();
    }
    std::string read_alignment = "";
    std::string graph_alignment = "";
    std::string trusted = "";
    std::vector<BoundedCounterType> covs;
    size_t farthest_seq_idx = node->seq_idx;
    ret->score = node->score;
    ret->truncated = (node->seq_idx != 0)
                     && (node->seq_idx != read.length() - 1);
#if READ_ALIGNER_DEBUG
    std::cerr << "Alignment end: " << node->prev << " "
              << node->base << " " << node->seq_idx << " "
              << node->state << " " << node->score << std::endl;
#endif

    char read_base;
    char graph_base;

#if READ_ALIGNER_DEBUG
    std::cerr << "graph_base" << "\t" << "read_base" << "\t"
              << "score\th_score\tf_score\tlength\tstate"
              << "\ttrusted?\tseq_idx\tfwd_hash\trc_hash" << std::endl;
#endif

    while(node != NULL && node->prev != NULL) {
        if(node->state == MATCH || node->state == MATCH_UNTRUSTED) {
            graph_base = toupper(nucl_lookup[node->base]);
            read_base = read[node->seq_idx];
        } else if(node->state == INSERT_READ || node->state == INSERT_READ_UNTRUSTED) {
            graph_base = '-';
            read_base = tolower(read[node->seq_idx]);
        } else if(node->state == INSERT_GRAPH
                  || node->state == INSERT_GRAPH_UNTRUSTED) {
            graph_base = tolower(nucl_lookup[node->base]);
            read_base = '-';
        } else {
            graph_base = '?';
            read_base = '?';
        }
#if READ_ALIGNER_DEBUG
        std::cerr << graph_base << "\t" << read_base << "\t"
                  << node->score << "\t" << node->h_score << "\t"
                  << node->f_score << "\t" << node->length << "\t"
                  << node->state << "\t" << node->trusted << "\t"
                  << node->seq_idx << "\t"
                  << m_ch->unhash_dna(node->fwd_hash) << "\t"
                  << m_ch->unhash_dna(node->rc_hash) << std::endl;
#endif

        if(forward) {
            graph_alignment = graph_base + graph_alignment;
            read_alignment = read_base + read_alignment;
            trusted = ((node->trusted)? "T" : "F") + trusted;
            covs.insert(covs.begin(), node->cov);
        } else {
            graph_alignment = graph_alignment + graph_base;
            read_alignment = read_alignment + read_base;
            trusted = trusted + ((node->trusted)? "T" : "F");
        }

        node = node->prev;
    }
    ret->graph_alignment = graph_alignment;
    ret->read_alignment = read_alignment;
    ret->trusted = trusted;
    ret->covs = covs;

    if(ret->truncated) {
        std::string new_graph_alignment;
        if (forward) {
            new_graph_alignment = graph_alignment +
                                  read.substr(farthest_seq_idx + 1, std::string::npos);
        } else {
            new_graph_alignment = read.substr(0, node->seq_idx)
                                  + graph_alignment;
        }
        ret->graph_alignment = new_graph_alignment;
    }

    return ret;

}

struct SearchStart {
    size_t kmer_idx;
    size_t k_cov;
    std::string kmer;
};

Alignment* ReadAligner::Align(const std::string& read)
{
    WordLength k = m_ch->ksize();
    size_t num_kmers = read.length() - k + 1;

    SearchStart start;
    start.k_cov = 0;
    start.kmer_idx = 0;

    for (size_t i = 0; i < num_kmers; i++) {
        std::string kmer = read.substr(i, k);

        size_t kCov = m_ch->get_count(kmer.c_str());
        if(kCov > start.k_cov) {
            start.kmer_idx = i;
            start.k_cov = kCov;
            start.kmer = kmer;
        }
    }

    if(start.k_cov == 0) {
        return _empty_alignment();
    }

    HashIntoType fhash = 0, rhash = 0;
    fhash = m_ch->hash_dna_top_strand(start.kmer.c_str());
    rhash = m_ch->hash_dna_bottom_strand(start.kmer.c_str());

#if READ_ALIGNER_DEBUG
    std::cerr << "Starting kmer: " << start.kmer << " "
              << m_ch->unhash_dna(fhash) << " "
              << m_ch->unhash_dna(rhash)
              << " cov: " << start.k_cov << " idx: " << start.kmer_idx << ", "
              << start.kmer_idx + k - 1
              << " emission: " << start.kmer[k - 1] << std::endl;
#endif

    Nucl e = _ch_to_nucl(start.kmer[k - 1]);
    AlignmentNode startingNode = AlignmentNode(NULL,
                                 e, start.kmer_idx + k - 1,
                                 MATCH, MM, fhash, rhash, k);
    startingNode.f_score = 0;
    startingNode.h_score = 0;
    Alignment* forward = NULL;
    Alignment* reverse = NULL;
    size_t final_length = 0;

    if(start.k_cov >= m_trusted_cutoff) {
        startingNode.score = k * m_sm.trusted_match + k * m_sm.tsc[MM];
    } else {
        startingNode.score = k * m_sm.untrusted_match + k * m_sm.tsc[MM];
    }

    forward = Subalign(&startingNode, read.length(), true, read);
    final_length = forward->read_alignment.length() + k;

    startingNode.seq_idx = start.kmer_idx;
    reverse = Subalign(&startingNode, read.length(), false, read);
    final_length += reverse->read_alignment.length();

    Alignment* ret = new Alignment;

    // We've actually counted the starting node score
    // twice, so we need to adjust for that
    ret->score = reverse->score + forward->score - startingNode.score;
    ret->read_alignment = reverse->read_alignment +
                          start.kmer + forward->read_alignment;
    ret->graph_alignment = reverse->graph_alignment +
                           start.kmer + forward->graph_alignment;
    ret->score = ret->score - GetNull(final_length);
    ret->truncated = forward->truncated || reverse->truncated;

#if READ_ALIGNER_DEBUG
    fprintf(stderr,
            "FORWARD\n\tread_aln:%s\n\tgraph_aln:%s\n\tscore:%f\n\ttrunc:%d\n",
            forward->read_alignment.c_str(), forward->graph_alignment.c_str(),
            forward->score, forward->truncated);
    fprintf(stderr,
            "REVERSE\n\tread_aln:%s\n\tgraph_aln:%s\n\tscore:%f\n\ttrunc:%d\n",
            reverse->read_alignment.c_str(), reverse->graph_alignment.c_str(),
            reverse->score, reverse->truncated);
#endif

    delete forward;
    delete reverse;

    return ret;
}

Alignment* ReadAligner::AlignForward(const std::string& read)
{
    WordLength k = m_ch->ksize();

    // start with seed at position 0
    SearchStart start;
    start.kmer = read.substr(0, k);
    start.kmer_idx = 0;
    start.k_cov = m_ch->get_count(start.kmer.c_str());

    if(start.k_cov == 0) {
        return _empty_alignment();
    }

    HashIntoType fhash = 0, rhash = 0;
    fhash = m_ch->hash_dna_top_strand(start.kmer.c_str());
    rhash = m_ch->hash_dna_bottom_strand(start.kmer.c_str());

#if READ_ALIGNER_DEBUG
    std::cerr << "Starting kmer: " << start.kmer << " "
              << m_ch->unhash_dna(fhash) << " "
              << m_ch->unhash_dna(rhash)
              << " cov: " << start.k_cov << " idx: " << start.kmer_idx << ", "
              << start.kmer_idx + k - 1
              << " emission: " << start.kmer[k - 1] << std::endl;
#endif

    Nucl e = _ch_to_nucl(start.kmer[k - 1]);
    AlignmentNode startingNode = AlignmentNode(NULL,
                                 e, start.kmer_idx + k - 1,
                                 MATCH, MM, fhash, rhash, k);
    startingNode.f_score = 0;
    startingNode.h_score = 0;
    Alignment* forward = NULL;
    size_t final_length = 0;

    if(start.k_cov >= m_trusted_cutoff) {
        startingNode.score = k * m_sm.trusted_match + k * m_sm.tsc[MM];
    } else {
        startingNode.score = k * m_sm.untrusted_match + k * m_sm.tsc[MM];
    }

    forward = Subalign(&startingNode, read.length(), true, read);
    final_length = forward->read_alignment.length() + k;

    Alignment* ret = new Alignment;

    ret->score = forward->score;
    ret->read_alignment = start.kmer + forward->read_alignment;
    ret->graph_alignment = start.kmer + forward->graph_alignment;
    ret->score = ret->score - GetNull(final_length);
    ret->truncated = forward->truncated;

    ret->covs = forward->covs;
    ret->covs.insert(ret->covs.begin(), start.k_cov);
    for (WordLength i = 0; i < k - 1; i++) {
        ret->covs.push_back(0);
    }

#if READ_ALIGNER_DEBUG
    fprintf(stderr,
            "FORWARD\n\tread_aln:%s\n\tgraph_aln:%s\n\tscore:%f\n\ttrunc:%d\n",
            forward->read_alignment.c_str(), forward->graph_alignment.c_str(),
            forward->score, forward->truncated);
#endif

    delete forward;
    return ret;
}

ScoringMatrix ReadAligner::getScoringMatrix()
{
    return m_sm;
}

}
