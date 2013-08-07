#include "read_aligner.hh"

namespace khmer {
// http://www.wesleysteiner.com/professional/del_fun.html

  void ReadAligner::enumerate(
		       const NodeHeap& open,
		       AlignmentNode* curr,
		       unsigned char forward,
		       const std::string& seq
		       ) {
    int next_seq_idx;
    int remaining;
    int kmerCov;
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
    
    // loop for matches and insertions
    for (int i = A; i <= T; i++) {
      next_nucl = nucl_lookup[i];
    
      if(forward) {
	next_fwd = next_f(fwd, next_nucl);
	next_rc = next_r(rc, next_nucl);
      } else {
	next_fwd = next_f(fwd, next_nucl);
	next_rc = next_r(rc, next_nucl);
      }
      
      hash = uniqify_rc(next_fwd, next_rc);
      kmerCov = ch->get_count(hash);
    
      if (kmerCov == 0) {
	continue;
      }
    
      // match
      next = new AlignmentNode(curr, (Nucl)i, next_seq_idx, match,
			       next_fwd, next_rc);
      next->score = curr->score + sm->score(next_nucl, seq[next_seq_idx]);
      isCorrect = kmerCov >= 2; //isCorrectKmer(nextKmerCov, lambdaOne, lambdaTwo);
      if (!isCorrect) {
	next->score += errorOffset;
      }
      next->h_score = sm->trans[MM] * remaining;
      next->f_score = next->score + next->h_score;
      
      if (isCorrect) {
	open.push(next);
      } else {
	delete next;
      }
    }
  }

  AlignmentNode* ReadAligner::subalign(AlignmentNode* startVert,
				   unsigned int seqLen,
				   bool forward,
				   const std::string& seq) {
    
    NodeHeap open();
    open.push(startVert);
    AlignmentNode* curr;

    while (!open.empty()) {
      curr = open.top();
      open.pop();
      
      if (curr->stateNo == seqLen-1 ||
	  curr->stateNo == 0) {
	return curr;
      }
      
      enumerate(open, forward, seq);
    }
    
    return NULL;
  }

  Alignment ReadAligner::align(const std::string& seq,
			   const std::string& kmer,
			   int index) {
    AlignmentNode* leftStart = new AlignemntNode(NULL,
						 kmer[0],
						 index,
						 'm',
						 Kmer(kmer));
    AlignmentNode* rightStart = new Node(NULL,
					 kmer[kmer.length()-1],
					 index + kmer.length()-1,
					 'm',
					 Kmer(kmer));
    AlignmentNode* leftGoal = subalign(leftStart,
				       seq.length(),
				       0,
				       leftClosed,
				       leftOpen,
				       seq);
    AlignmentNode* rightGoal = subalign(rightStart,
					seq.length(),
					1,
					rightClosed,
					rightOpen,
					seq);
    
    std::string align = extractString(leftGoal, 0, &readDels) +
      kmer +
      extractString(rightGoal, 1, &readDels);
    
    return Alignment(readDels, align);
  }

  Alignment ReadAligner::alignRead(const std::string& read) {
    std::vector<unsigned int> markers;
    bool toggleError = 1;
    
    unsigned int longestErrorRegion = 0;
    unsigned int currentErrorRegion = 0;
    
    unsigned int k = ch->ksize();

    std::set<CandidateAlignment> alignments;
    CandidateAlignment best = CandidateAlignment();

    std::string graphAlign = "";

    for (unsigned int i = 0; i < read.length() - k + 1; i++) {
      std::string kmer = read.substr(i, k);

      assert(kmer.length() == k);
      
      int kCov = ch->get_count(kmer.c_str());
      
      bool isCorrect = isCorrectKmer(kCov, lambdaOne, lambdaTwo);
      
      if (isCorrect && currentErrorRegion) {
	currentErrorRegion = 0;
      }
      
      if (!isCorrect) {
	currentErrorRegion++;
	
	if (currentErrorRegion > longestErrorRegion) {
	  longestErrorRegion = currentErrorRegion;
	}
      }
      
      if (toggleError && isCorrect) {
	markers.push_back(i);
	toggleError = 0;
      } else if (!toggleError && !isCorrect) {
	markers.push_back(i-1);
	toggleError = 1;
      }
    }
    
    // couldn't find a seed k-mer
    if (markers.size() == 0) {
      //std::cout << "Couldn't find a seed k-mer." << std::endl;
      return best;
    }

   // exceeded max error region parameter
    if (longestErrorRegion > maxErrorRegion && maxErrorRegion >= 0) {
      return best;
    }

   // read appears to be error free
    if (markers.size() == 1 && markers[0] == 0) {
      std::map<int,int> readDels;
      CandidateAlignment retAln = CandidateAlignment(readDels, read);
      return retAln;
    }
   
    unsigned int startIndex = 0;

    if (markers[0] != 0) {
      unsigned int index = markers[0];
      Alignment aln = align(ch,
			    read.substr(0, index+k),
			    read.substr(index, k),
			    index);
      
      graphAlign += aln.alignment.substr(0,aln.alignment.length()-k); 
      startIndex++;

      if (markers.size() > 1) {
	graphAlign += read.substr(index, markers[1]-index);
      } else {
	graphAlign += read.substr(index);
      }
    } else {
      graphAlign += read.substr(0, markers[1]-markers[0]);
      startIndex++;
    }

    for (unsigned int i = startIndex; i < markers.size(); i+=2) {
      unsigned int index = markers[i];

      if (i == markers.size()-1) {
	Alignment aln = align(ch,
			      read.substr(index),
			      read.substr(index, k),
			      0);
	graphAlign += aln.alignment.substr(0,aln.alignment.length());
	break;
      } else {
	Alignment aln = align(ch, 
			      read.substr(index, markers[i+1]-index+k),
			      read.substr(index, k),
			      0);
	size_t kmerInd = aln.alignment.rfind(read.substr(markers[i+1], k));
	if (kmerInd == std::string::npos) {
	  return best;
	} else {
	  graphAlign += aln.alignment.substr(0, kmerInd);
	}
      }

      // add next correct region to alignment
      if (i+1 != markers.size()-1) {
	graphAlign += read.substr(markers[i+1], markers[i+2]-markers[i+1]);
      } else {
	graphAlign += read.substr(markers[i+1]);
      }
    }
   
    return NULL;
  }
}
