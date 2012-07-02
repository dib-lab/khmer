#include "readaligner.hh"

AStarSearchNode * assn_set_find(std::set<AStarSearchNode*>* a,
                               AStarSearchNode* val) {

   std::set<AStarSearchNode*>::iterator it;
   
   for (it = a->begin(); it != a->end(); it++) {
      if (**it == *val) {
         return *it;
      }
   }

   return NULL;
}

AStarSearchNode * assn_vector_find(std::vector<AStarSearchNode*>* a,
                                  AStarSearchNode* val) {

   std::vector<AStarSearchNode*>::iterator it;

   for (it = a->begin(); it != a->end(); it++) {
      if (**it == *val) {
         return (*it);
      }
   }

   return NULL;
}

AStarSearchNode* ReadAligner::subalign(AStarSearchNode* startVert,
                                       NodeEnumerator* enumerator,
                                       int seqLen,
                                       std::set<AStarSearchNode*>* closed,
                                       std::string seq) {
   std::vector<AStarSearchNode*> open;

   open.push_back(startVert);
   std::push_heap(open.begin(), open.end(), ASSNCompare());

   AStarSearchNode* curr;

   while (!open.empty()) {
      curr = open.front();
      std::pop_heap(open.begin(), open.end(), ASSNCompare());
      open.pop_back();
     
      closed->insert(curr);

      int max = 0;

      if (curr->stateNo > max) {
         max = curr->stateNo;
         //std::cout << closed->size() << " " << open.size() << " " << curr->deletes << " " << curr->kmer << " " << curr->stateNo << std::endl;
      }

      if (curr->stateNo == seqLen -1 ||
          curr->stateNo == 0) {
         std::cout << "returning curr " << closed->size() << " " << open.size() << std::endl;
         return curr;
      }

      std::queue<AStarSearchNode*> nodes = enumerator->enumerateNodes(curr, ch);
      
      while (!nodes.empty()) {
         AStarSearchNode* next = nodes.front();
         nodes.pop();

         AStarSearchNode* where = assn_vector_find(&open, next);
         AStarSearchNode* in_closed = assn_set_find(closed, next);
         if (in_closed == NULL &&
               (where == NULL ||
               next->score > (where)->score)) {
         
            int inClosed = (in_closed == NULL)? 0 : 1;
            int inOpen = (where == NULL)? 0 : 1;

            //std::cout << next->kmer << " " << next->state << " " << next->stateNo << " " << inClosed << " " << inOpen << std::endl;
            open.push_back(next);
            std::push_heap(open.begin(), open.end(), ASSNCompare());
         }
      }
   }   
   std::cout << "returning NULL " << closed->size() << std::endl;
   return NULL;
}

std::string ReadAligner::extractString(AStarSearchNode* goal,
                                       char forward,
                                       std::map<int,int>* readDeletions) {
   std::string ret;
   
   while (goal->discoveredFrom != NULL) {
      char b = goal->emission;
      
      ret += b;

      if (goal->state == 'i' && readDeletions != NULL) {
         if (readDeletions->count(goal->stateNo) == 0) {
            (*readDeletions)[goal->stateNo] = 1;
         } else {
            (*readDeletions)[goal->stateNo]++;
         }
      }
  
      goal = goal->discoveredFrom;
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

CandidateAlignment* ReadAligner::align(khmer::CountingHash * ch,
                                      std::string seq,
                                      std::string kmer,
                                      int index) {
   std::set<AStarSearchNode*> leftClosed;
   std::set<AStarSearchNode*> rightClosed;
   AStarSearchNode * leftStart = new AStarSearchNode(NULL,
                                                     kmer[0],
                                                     index,
                                                     'm',
                                                     kmer);
   AStarSearchNode * rightStart = new AStarSearchNode(NULL,
                                                      kmer[kmer.length()-1],
                                                      index + kmer.length()-1,
                                                      'm',
                                                      kmer);
   AStarSearchNode * leftGoal = subalign(leftStart, 
                                         new NodeEnumerator(0, seq, sm),
                                         seq.length(), 
                                         &leftClosed,
                                         seq);
   //rightStart->snps = leftGoal->snps;
   AStarSearchNode * rightGoal = subalign(rightStart,
                                          new NodeEnumerator(1, seq, sm),
                                          seq.length(),
                                          &rightClosed,
                                          seq);

   if (leftGoal == NULL || rightGoal == NULL) {
      return NULL;
   }

   std::map<int,int> readDels;

   std::string align = extractString(leftGoal, 0, &readDels) + 
                       kmer + 
                       extractString(rightGoal, 1, &readDels);

   int score = 0;
   int readIndex = 0;
   int tmpDels = 0;
   for (int i = 0; i < (int)align.length(); i++) {
      if (tmpDels > 0) {
         score += sm->score(align[i], '-');
         tmpDels--;
      }  else if (readDels.count(i) == 0) {
         score += sm->score(align[i], seq[readIndex]);
         readIndex++;
      }  else  {
         score += sm->score(align[i], seq[readIndex]);
         tmpDels = readDels[i];
         readIndex++;
      }
   }

   clearSet(&leftClosed);
   clearSet(&rightClosed);

   CandidateAlignment* aln = new CandidateAlignment(readDels, align, score);

   return aln;
}

CandidateAlignment ReadAligner::alignRead(std::string read) {
   int k = ch->ksize();   

   std::set<CandidateAlignment*> alignments;

   for (int i = 0; i < (int)read.length() - k + 1; i++) {
      std::string kmer = read.substr(i, k);  

      std::cout << kmer << std::endl;

      assert(kmer.length() == k);
   
      if (ch->get_count(kmer.c_str())) {
         CandidateAlignment* aln = align(ch, read, kmer, i);
         if (aln != NULL) {
            alignments.insert(aln);
         }
      }
   }

   // find the best alignment...
   CandidateAlignment* best = NULL;
   std::set<CandidateAlignment*>::iterator it;
   for (it=alignments.begin(); it != alignments.end(); it++) {
      CandidateAlignment* curr = *it;
      if (best == NULL) {
         best = curr;  
      } else if (best->score < curr->score) {
         best = curr;
      }
   }

   if (best != NULL) {
      CandidateAlignment ret = CandidateAlignment(best->readDeletions,
                                                  best->alignment,
                                                  best->score);
      return ret;
   }  else  {
      std::map<int,int> emptyMap;
      CandidateAlignment ret = CandidateAlignment(emptyMap, "", 0);
      return ret;
   }
}

void ReadAligner::clearSet(std::set<AStarSearchNode*>* s) {
   std::set<AStarSearchNode*>::iterator it;

   for (it=s->begin(); it != s->end(); it++) {
      delete (*it);
   }
}
