/*
  file: read_aligner.cc
  author: Jordan Fish
  purpose: align...reads?
*/
#include "read_aligner.hh"

namespace khmer {

  struct del_alignment_node_t
  {
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
		       std::vector<AlignmentNode*>& all_nodes,
		       AlignmentNode* curr,
		       bool forward,
		       const std::string& seq
		       ) {
    int next_seq_idx;
    int remaining;
    unsigned int kmerCov;
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
	    sc = (kmerCov >= m_trusted_cutoff)? m_sm.trusted_match : m_sm.untrusted_match;
	  } else {
	    sc = (kmerCov >= m_trusted_cutoff)? m_sm.trusted_mismatch : m_sm.untrusted_mismatch;
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
	next->trusted = (kmerCov >= m_trusted_cutoff);
	isCorrect = true; //kmerCov >= 1; //isCorrectKmer(nextKmerCov, lambdaOne, lambdaTwo);
	next->h_score = hcost;
	next->f_score = next->score + next->h_score;
	
	if (isCorrect && next->score - GetNull(next->length) > next->length * m_bits_theta) {
	  open.push(next);
	  all_nodes.push_back(next);
	} else {
	  delete next;
	}
      }
    }
  }

  Alignment* ReadAligner::Subalign(AlignmentNode* start_vert,
				   unsigned int seqLen,
				   bool forward,
				   const std::string& seq) {
    std::vector<AlignmentNode*> all_nodes;
    NodeHeap open;
    std::set<AlignmentNode> closed;
    open.push(start_vert);
    //all_nodes.insert(start_vert);
    
    AlignmentNode* curr = NULL;
    AlignmentNode* best = NULL;

    while (!open.empty()) {
      curr = open.top();
      #if READ_ALIGNER_DEBUG
      std::cerr << "curr: " << curr->prev << " "
		<< _revhash(curr->fwd_hash, m_ch->ksize()) << " " << _revhash(curr->rc_hash, m_ch->ksize())
		<< " cov=" << m_ch->get_count(uniqify_rc(curr->fwd_hash, curr->rc_hash))
		<< " emission=" << curr->base
		<< " seqidx=" << curr->seq_idx
		<< " score=" << curr->score
	        << " fscore=" << curr->f_score
	        << " bits_saved=" << curr->score - GetNull(curr->length) << " " << ((best != NULL)? best->score - GetNull(best->length) : -1) << std::endl;
      #endif
      open.pop();

      if(best == NULL ||
	 (best->score - GetNull(best->length) < curr->score - GetNull(curr->length))) {
	best = curr;
      }
      
      if (curr->seq_idx == seqLen-1 ||
	  curr->seq_idx == 0) {
	best = curr;
	break;
      }

      if(set_contains(closed, *curr)) {
	continue;
      }

      closed.insert(*curr);
      
      Enumerate(open, all_nodes, curr, forward, seq);
    }

    Alignment* ret = ExtractAlignment(best, forward, seq);
    std::for_each(all_nodes.begin(), all_nodes.end(), del_alignment_node());
    
    return ret;
  }

  Alignment* ReadAligner::ExtractAlignment(AlignmentNode* node, bool forward, const std::string& read) {
    Alignment* ret = new Alignment;
    if(node == NULL) {
      ret->score = 0;
      ret->read_alignment = "";
      ret->graph_alignment = "";
      ret->trusted = "";
      ret->truncated = true;
      return ret;
    }

    assert(node->seq_idx < read.length());
    assert(node->seq_idx >= 0);
    std::string read_alignment = "";
    std::string graph_alignment = "";
    std::string trusted = "";
    std::string cov = "";
    ret->score = node->score;
    ret->truncated = (node->seq_idx != 0) && (node->seq_idx != read.length() - 1);
    #if READ_ALIGNER_DEBUG
    std::cerr << "Alignment end: " << node->prev << " " << node->base << " " << node->seq_idx << " " << node->state << " " << node->score << std::endl;
    #endif

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
      if(forward) {
	graph_alignment = graph_base + graph_alignment;
	read_alignment = read_base + read_alignment;
	trusted = ((node->trusted)? "T" : "F") + trusted;
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
    
    return ret;
  }

  struct SearchStart {
    int kmer_idx;
    unsigned int k_cov;
    std::string kmer;
  };
  
  Alignment* ReadAligner::Align(const std::string& read) {
    int k = m_ch->ksize();
    int num_kmers = read.length() - k + 1;
    int mid_kmer = num_kmers / 2;

    SearchStart start;
    start.k_cov = 0;
    
    for (unsigned int i = 0; i < num_kmers; i++) {
      std::string kmer = read.substr(i, k);
      
      unsigned int kCov = m_ch->get_count(kmer.c_str());
      if(kCov > 0) {
	if(start.k_cov < m_trusted_cutoff || abs(mid_kmer - i) < abs(mid_kmer - start.kmer_idx)) {
	  start.kmer_idx = i;
	  start.k_cov = kCov;
	  start.kmer = kmer;
	}
      }
    }

    if(start.k_cov > 0) {
      HashIntoType fhash = 0, rhash = 0;
      _hash(start.kmer.c_str(), k, fhash, rhash);
      #if READ_ALIGNER_DEBUG
      std::cerr << "Starting kmer: " << start.kmer << " " << _revhash(fhash, m_ch->ksize()) << " " << _revhash(rhash, m_ch->ksize()) << " cov: " << start.k_cov << " idx: " << start.kmer_idx << ", " << start.kmer_idx + k - 1 << " emission: " << start.kmer[k - 1] << std::endl;
      #endif
      char base = toupper(start.kmer[k - 1]);
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
      case 'T': case 'U':
	e = T;
	break;
      }
      
      AlignmentNode startingNode = AlignmentNode(NULL, e, start.kmer_idx + k - 1, MATCH, fhash, rhash, k);
      startingNode.f_score = 0;
      startingNode.h_score = 0;
      Alignment* forward = NULL;
      Alignment* reverse = NULL;
      int final_length = 0;
      
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
      ret->score = reverse->score + forward->score - startingNode.score; //We've actually counted the starting node score twice, so we need to adjust for that
      ret->read_alignment = reverse->read_alignment + start.kmer + forward->read_alignment;
      ret->graph_alignment = reverse->graph_alignment + start.kmer + forward->graph_alignment;
      ret->score = ret->score -  GetNull(final_length);
      ret->truncated = forward->truncated || reverse->truncated;

      #if READ_ALIGNER_DEBUG
      fprintf(stderr, "FORWARD\n\tread_aln:%s\n\tgraph_aln:%s\n\tscore:%f\n\ttrunc:%d\n", forward->read_alignment.c_str(), forward->graph_alignment.c_str(), forward->score, forward->truncated);
      fprintf(stderr, "REVERSE\n\tread_aln:%s\n\tgraph_aln:%s\n\tscore:%f\n\ttrunc:%d\n", reverse->read_alignment.c_str(), reverse->graph_alignment.c_str(), reverse->score, reverse->truncated);
      #endif
      
      delete forward;
      delete reverse;
      return ret;
    } else {

      Alignment* ret = new Alignment;
      ret->score = -std::numeric_limits<double>::infinity();
      ret->read_alignment = "";
      ret->graph_alignment = "";
      ret->truncated = true;
      return ret;
    }
  }
}
