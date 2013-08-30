/*
  file: read_aligner.cc
  author: Jordan Fish
  purpose: align...reads?
*/
#include "read_aligner.hh"

namespace khmer {

  /*
    Compute the null model (random sequence) log odds probability for a given length
   */
  double GetNull(int length)
  {
    return log2(.25) * length + log2(1.0 / (length + 1));
  }

  /*
    Turn two states in to a transition, or disallowed if the transition isn't modeled
   */
  Transition get_trans(State s1, State s2) {
    if(s1 == MATCH) {
      if(s2 == MATCH) {
	return MM;
      } else if(s2 == INSERT_GRAPH) {
	return MD;
      } else if(s2 == INSERT_READ) {
	return MI;
      }
    } else if(s1 == INSERT_GRAPH) {
      if(s2 == MATCH) {
	return DM;
      } else if(s2 == INSERT_GRAPH) {
	return DD;
      }
    } else if(s1 == INSERT_READ) {
      if(s2 == MATCH) {
	return IM;
      } else if(s2 == INSERT_READ) {
	return II;
      }
    }

    return disallowed;
  }

  void ReadAligner::Enumerate(
		       NodeHeap& open,
		       AlignmentNode* curr,
		       bool forward,
		       const std::string& seq
		       ) {
    int next_seq_idx;
    int remaining;
    int kmerCov;
    double hcost;
    double sc;
    Transition trans;
    HashIntoType fwd = curr->fwd_hash;
    HashIntoType rc = curr->rc_hash;
    HashIntoType next_fwd, next_rc, hash;
    bool isCorrect;
    unsigned char next_nucl;
    
    AlignmentNode* next;
    
    if (forward) {
      next_seq_idx = curr->seq_idx + 1;
      remaining = seq.size() - next_seq_idx;
    } else {
      next_seq_idx = curr->seq_idx - 1;
      remaining = next_seq_idx;
    }

    //std::cerr << "Enumerating: " << _revhash(fwd, ch->ksize()) << " " << _revhash(rc, ch->ksize()) << std::endl;
    
    // loop for MATCHes and INSERT_READs
    for (int i = A; i <= T; i++) {
      next_nucl = nucl_lookup[i];
    
      if(forward) {
	next_fwd = next_f(fwd, next_nucl);
	next_rc = next_r(rc, next_nucl);
      } else {
	next_fwd = prev_f(fwd, next_nucl);
	next_rc = prev_r(rc, next_nucl);
      }
      
      hash = uniqify_rc(next_fwd, next_rc);
      kmerCov = m_ch->get_count(hash);

      if (kmerCov == 0) {
	continue;
      }

      for(int next_state_iter = MATCH;next_state_iter <= INSERT_GRAPH;next_state_iter++) {
	State next_state = static_cast<State>(next_state_iter);
	trans = get_trans(curr->state, next_state);
	hcost = m_sm.tsc[get_trans(next_state, MATCH)] + m_sm.tsc[MM] * (remaining - 1);
	
	if(trans == disallowed) {
	  continue;
	}

	if(next_state == MATCH) {
	  if(next_nucl == seq[next_seq_idx]) {
	    sc = m_sm.match;
	  } else {
	    sc = m_sm.mismatch;
	  }
	} else {
	  sc = 0;
	}

	if(next_state == MATCH) {
	  next = new AlignmentNode(curr, (Nucl)i, next_seq_idx, (State)next_state,
				   next_fwd, next_rc, curr->length + 1);
	} else if(next_state == INSERT_GRAPH) {
	  next = new AlignmentNode(curr, (Nucl)i, curr->seq_idx, (State)next_state,
				   next_fwd, next_rc, curr->length);
	} else if(next_state == INSERT_READ) {
	  next = new AlignmentNode(curr, (Nucl)i, next_seq_idx, (State)next_state,
				   curr->fwd_hash, curr->rc_hash, curr->length + 1);
	}

	next->score = curr->score + sc + m_sm.tsc[trans];
	isCorrect = true; //kmerCov >= 1; //isCorrectKmer(nextKmerCov, lambdaOne, lambdaTwo);
	next->h_score = hcost;
	next->f_score = next->score + next->h_score;
	//std::cerr << "\t" << i << " " << next_nucl << " " << next->score << " " << next->h_score << " " << next->f_score << " " << state_labels[next->state] << " " << trans_labels[trans] << " " << hash << " " << _revhash(hash, m_ch->ksize()) << " " << _revhash(next_fwd, m_ch->ksize()) << " " << _revhash(next_rc, ch->ksize()) << " " << kmerCov << std::endl;
	
	if (isCorrect && next->score - GetNull(next->length) > next->length) {
	  open.push(next);
	} else {
	  delete next;
	}
      }
    }
  }

  AlignmentNode* ReadAligner::Subalign(AlignmentNode* startVert,
				   unsigned int seqLen,
				   bool forward,
				   const std::string& seq) {
    
    NodeHeap open;
    std::set<AlignmentNode> closed;
    open.push(startVert);
    AlignmentNode* curr;
    AlignmentNode* best = NULL;

    while (!open.empty()) {
      curr = open.top();
      //std::cerr << "curr: " << curr->prev << " " << _revhash(curr->fwd_hash, m_ch->ksize()) << " " << _revhash(curr->rc_hash, m_ch->ksize()) << " " << curr->base << " " << curr->seq_idx << " " << state_labels[curr->state] << " " << curr->score << " " << curr->f_score << std::endl;
      open.pop();

      if(best == NULL || (best->score - GetNull(best->length) > curr->score - GetNull(curr->length))) {
	best = curr;
      }
      
      if (curr->seq_idx == seqLen-1 ||
	  curr->seq_idx == 0) {
	return curr;
      }

      if(set_contains(closed, *curr)) {
	continue;
      }

      closed.insert(*curr);
      
      Enumerate(open, curr, forward, seq);
    }
    
    return curr;
  }

  Alignment* ReadAligner::ExtractAlignment(AlignmentNode* node, const std::string& read, const std::string& kmer) {
    Alignment* ret = new Alignment;
    if(node == NULL) {
      ret->score = 0;
      ret->read_alignment = "";
      ret->graph_alignment = "";
      ret->truncated = true;
      return ret;
    }

    assert(node->seq_idx < read.length());
    assert(node->seq_idx > 0);
    std::string read_alignment = "";
    std::string graph_alignment = "";
    ret->score = node->score - GetNull(read.length());
    ret->truncated = (node->seq_idx != 0) && (node->seq_idx != read.length() - 1);
    //std::cerr << "Alignment end: " << node->prev << " " << node->base << " " << node->seq_idx << " " << node->state << " " << node->score << std::endl;

    char read_base;
    char graph_base;

    while(node != NULL && node->prev != NULL) {
      if(node->state == MATCH) {
	graph_base = toupper(nucl_lookup[node->base]);
	read_base = read[node->seq_idx];
      } else if(node->state == INSERT_READ) {
	graph_base = '-';
	read_base = tolower(read[node->seq_idx]);
      } else if(node->state == INSERT_GRAPH) {
	graph_base = tolower(nucl_lookup[node->base]);
	read_base = '-';
      } else {
	graph_base = '?';
	read_base = '?';
      }
      graph_alignment = graph_base + graph_alignment;
      read_alignment = read_base + read_alignment;
      node = node->prev;
    }
    ret->graph_alignment = kmer + graph_alignment;
    ret->read_alignment = kmer + read_alignment;
    
    return ret;
  }
  
  Alignment* ReadAligner::Align(const std::string& read) {
    int k = m_ch->ksize();
    for (unsigned int i = 0; i < read.length() - k + 1; i++) {
      std::string kmer = read.substr(i, k);
      
      int kCov = m_ch->get_count(kmer.c_str());
      if(kCov > 0) {
	HashIntoType fhash, rhash;
	_hash(kmer.c_str(), k, fhash, rhash);
	//std::cerr << "Starting kmer: " << kmer << " " << _revhash(fhash, m_ch->ksize()) << " " << _revhash(rhash, m_ch->ksize()) << " idx: " << i << ", " << i + k - 1 << " emission: " << kmer[k - 1] << std::endl;
	char base = toupper(kmer[k - 1]);
	Nucl e;
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
	case 'T': case 'U':
	  e = T;
	  break;
	}
	AlignmentNode start(NULL, e, i + k - 1, MATCH, fhash, rhash, k);
	start.score = k * m_sm.match + k * m_sm.tsc[MM];
	AlignmentNode* end = Subalign(&start, read.length(), true, read);
	return ExtractAlignment(end, read, kmer);
      }
    }
    return NULL;
  }

}
