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

      std::queue<Node*> nodes = curr->enumerate(ch, sm, forward, seq);

      while (!nodes.empty()) {
         Node * next = nodes.front();
         nodes.pop();

         std::vector<Node*>::iterator where = node_vector_find(open, next);
         std::set<Node*>::iterator in_closed = node_set_find(closed, next);

         if (in_closed == closed.end() &&
             (where == open.end() ||
              next->gval > (*where)->gval)) {
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

   // score up the alignment
   int score = 0;
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
  
   // memory cleanup!
   for_each(leftOpen.begin(), leftOpen.end(), del_fun<Node>());
   for_each(rightOpen.begin(), rightOpen.end(), del_fun<Node>());
   for_each(leftClosed.begin(), leftClosed.end(), del_fun<Node>());
   for_each(rightClosed.begin(), rightClosed.end(), del_fun<Node>()); 

   return CandidateAlignment(readDels, align, score);
}

CandidateAlignment Aligner::alignRead(const std::string& read) {
   unsigned int k = ch->ksize();

   std::set<CandidateAlignment> alignments;
   CandidateAlignment best = CandidateAlignment();

   for (unsigned int i = 0; i < read.length() - k + 1; i++) {
      std::string kmer = read.substr(i, k);

      assert(kmer.length() == k);

      if (ch->get_count(kmer.c_str())) {
         if (best.alignment.find(kmer) == std::string::npos) {
            CandidateAlignment aln = align(ch, read, kmer, i);
            if (aln.alignment != "")  {
               alignments.insert(aln);
               if (aln.score > best.score || best.score == 1) {
                  best = aln;
               }
            }
         }
      }
   }

/*
   // find the best alignment...
   CandidateAlignment * best = NULL;
   std::set<CandidateAlignment>::iterator it;
   for (it=alignments.begin(); it != alignments.end(); it++) {
      CandidateAlignment curr = (*it);
      if (best == NULL) {
         best = &curr;
      } else if (best->score < curr.score) {
         best = &curr;
      }
   }
*/
   if (best.score < 1) {
      return best;
   } else {
      return CandidateAlignment();
   }
}
