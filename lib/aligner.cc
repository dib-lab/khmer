#include "aligner.hh"

// http://www.wesleysteiner.com/professional/del_fun.html
template<class T>
struct del_fun_t {
   del_fun_t& operator()(T* p) {
      delete p;
      return *this;
   }
};

template<class T>
del_fun_t<T> del_fun() {
   return del_fun_t<T>();
}

std::set<Node*>::iterator node_set_find(std::set<Node*>& a,
                    Node* val) {
   std::set<Node*>::iterator it;

   for (it = a.begin(); it != a.end(); it++) {
      if (**it == *val) {
         return it;
      }
   }

   return a.end();
}

std::vector<Node*>::iterator node_vector_find(std::vector<Node*>& a,
                       Node* val) {
   std::vector<Node*>::iterator it;

   for (it = a.begin(); it != a.end(); it++) {
      if (**it == *val) {
         return it;
      }
   }

   return a.end();
}

Node * Aligner::subalign(Node * startVert,
                        unsigned int seqLen,
                        unsigned char forward,
                        std::set<Node*>& closed,
                        std::vector<Node*>& open,
                        const std::string& seq) {

   make_heap(open.begin(), open.end());
   open.push_back(startVert);
   std::push_heap(open.begin(), open.end(), NodeCompare());

   while (!open.empty()) {
      Node * curr = open.front();
      std::pop_heap(open.begin(), open.end(), NodeCompare());
      open.pop_back();

      closed.insert(curr);

      if (curr->stateNo == seqLen-1 ||
          curr->stateNo == 0) {
         return curr;
      }

      std::queue<Node*> nodes = curr->enumerate(ch, sm, forward, seq,
                                                lambdaOne,
                                                lambdaTwo);

      while (!nodes.empty()) {
         Node * next = nodes.front();
         nodes.pop();

         std::vector<Node*>::iterator where = node_vector_find(open, next);
         std::set<Node*>::iterator in_closed = node_set_find(closed, next);

         if (in_closed == closed.end() &&
             (where == open.end() ||
              next->gval < (*where)->gval)) {
            open.push_back(next);
            std::push_heap(open.begin(), open.end(), NodeCompare());
           

         } else {
            delete next;
         }
      } 
   }

   return NULL;
}

std::string Aligner::extractString(Node* goal,
                                   unsigned char forward,
                                   std::map<int,int>* readDeletions) {
   std::string ret;

   while (goal->parent != NULL) {
      char b = goal->emission;

      ret += b;

      if (goal->state == 'i' && readDeletions != NULL) {
         if (readDeletions->count(goal->stateNo) == 0) {
            (*readDeletions)[goal->stateNo] = 1;
         } else {
            (*readDeletions)[goal->stateNo]++;
         }
      }

      goal = goal->parent;
   }


   // reverse the string if we are going in the forward direction
   if (forward) {
      std::string tmp;

      for (int i = ret.length()-1; i >= 0; i--) {
         tmp += ret[i];
      }

      ret = tmp;
   }
   return ret;
}

CandidateAlignment Aligner::align(khmer::CountingHash * ch,
                                    const std::string& seq,
                                    const std::string& kmer,
                                    int index) {
   std::set<Node*> leftClosed;
   std::set<Node*> rightClosed;
   std::vector<Node*> leftOpen;
   std::vector<Node*> rightOpen;

   Node * leftStart = new Node(NULL,
                         kmer[0],
                         index,
                         'm',
                         Kmer(kmer));
   Node * rightStart = new Node(NULL,
                          kmer[kmer.length()-1],
                          index + kmer.length()-1,
                          'm',
                          Kmer(kmer));
   Node * leftGoal = subalign(leftStart,
                              seq.length(),
                              0,
                              leftClosed,
                              leftOpen,
                              seq);
   Node * rightGoal = subalign(rightStart,
                               seq.length(),
                               1,
                               rightClosed,
                               rightOpen,
                               seq);

   if (leftGoal == NULL || rightGoal == NULL) {
      for_each(leftOpen.begin(), leftOpen.end(), del_fun<Node>());
      for_each(rightOpen.begin(), rightOpen.end(), del_fun<Node>());
      for_each(leftClosed.begin(), leftClosed.end(), del_fun<Node>());
      for_each(rightClosed.begin(), rightClosed.end(), del_fun<Node>());
      return CandidateAlignment();
   }

   std::map<int,int> readDels;

   std::string align = extractString(leftGoal, 0, &readDels) +
                       kmer +
                       extractString(rightGoal, 1, &readDels);

   /*
   // score up the alignment
   double score = 0;
   int readIndex = 0;
   int tmpDels = 0;
   for (int i = 0; i < (int)align.length(); i++) {
      if (tmpDels > 0) {
         score += sm->score(align[i], '-');
         tmpDels--;
      } else if (readDels.count(i) == 0) {
         score += sm->score(align[i], seq[readIndex]);
         readIndex++;
      } else {
         score += sm->score(align[i], seq[readIndex]);
         tmpDels = readDels[i];
         readIndex++;
      }
   }
   */
  
   // memory cleanup!
   for_each(leftOpen.begin(), leftOpen.end(), del_fun<Node>());
   for_each(rightOpen.begin(), rightOpen.end(), del_fun<Node>());
   for_each(leftClosed.begin(), leftClosed.end(), del_fun<Node>());
   for_each(rightClosed.begin(), rightClosed.end(), del_fun<Node>()); 

   return CandidateAlignment(readDels, align);
}

void Aligner::printErrorFootprint(const std::string& read) {
   unsigned int k = ch->ksize();

   for (unsigned int i = 0; i < read.length() - k + 1; i++) {
      std::string kmer = read.substr(i, k);

      assert(kmer.length() == k);

      int kCov = ch->get_count(kmer.c_str());

      bool isCorrect = isCorrectKmer(kCov, lambdaOne, lambdaTwo);

      std::cout << isCorrect;
   }

   std::cout << std::endl;
}

CandidateAlignment Aligner::alignRead(const std::string& read) {
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
#if (0) // Jason's original code
   if (longestErrorRegion > maxErrorRegion && maxErrorRegion >= 0) {
#else
   if (longestErrorRegion > maxErrorRegion && maxErrorRegion < UINT_MAX) {
#endif
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
      CandidateAlignment aln = align(ch,
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
         CandidateAlignment aln = align(ch,
                                        read.substr(index),
                                        read.substr(index, k),
                                        0);
         graphAlign += aln.alignment.substr(0,aln.alignment.length());
         break;
      } else {
         CandidateAlignment aln = align(ch, 
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
   
   std::map<int,int> readDels;
   CandidateAlignment retAln = CandidateAlignment(readDels, graphAlign);
   return retAln;
}
