// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Stephan Aiche <aiche@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_FIND_PEX_H
#define SEQAN_HEADER_FIND_PEX_H

// uncomment this for verbose debug output
//#define SEQAN_DEBUG_PEX

namespace SEQAN_NAMESPACE_MAIN
{

struct Hierarchical;
struct NonHierarchical;


template <typename TVerification, typename TMultiFinder = WuManber>
struct Pex;

typedef Pex<Hierarchical,AhoCorasick>      PexHierarchical;
typedef Pex<NonHierarchical,AhoCorasick>   PexNonHierarchical;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TVerification, typename TMultiFinder>
struct FindBeginPatternSpec< Pattern<TNeedle, Pex<TVerification , TMultiFinder > > >:
    DefaultFindBeginPatternSpec<>
{
};

//////////////////////////////////////////////////////////////////////////////

/*!
 * @mfn Pattern#PexMultiFinder
 * @headerfile <seqan/find.h>
 * @brief Determines the multiple exact string matching algorithm used by the Pex algorithm.
 *
 * @signature PexMultiFinder<Pattern<TNeedle,Pex<TVerification, TMultiFinder> > >::Type;
 *
 * @tparam TMultiFinder The specification for the multiple exact string matching algorithm that should be used
 *                      with the Pex algorithm.
 *
 * @return Type Pattern type of the multiple exact string matching algorithm for the specified Pattern.
 *
 * For a description of Pattern usage, see @link Pattern @endlink.
 *
 * Overload this metafunction if you want to use something else for verification then
 * <tt>Pattern&lt;String&lt;Segment&l;tNeedle&gt; &gt;, TMultiFinder&gt;</tt>.
 */

template<typename T>
struct PexMultiFinder;

template<typename TNeedle, typename TVerification, typename TMultiFinder>
struct PexMultiFinder< Pattern<TNeedle, Pex<TVerification , TMultiFinder > > >
{
  typedef Pattern<String<Segment<TNeedle> > , TMultiFinder> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TPosition,typename TScore,typename TVerifier,typename TNeedle>
struct PexRange_{
  TPosition start,end;
  TScore error;
  TVerifier verifier;
};

//////////////////////////////////////////////////////////////////////////////

/*!
 * @class PexPattern
 * @extends Pattern
 * @headerfile <seqan/find.h>
 * @brief Provides a fast approximate string matching filter that splits the needle into several pieces that are
 *        searched with a multiple exact string matching algorithm and later verified.
 *
 * @signature template <typename TNeedle, typename TVerification, typename TMultiFinder>
 *            class Pattern<TNeedle, Pex<TVerification, TMultiFinder> >;
 *
 * @tparam TNeedle       The needle type. Type: @link ContainerConcept @endlink
 * @tparam TVerification Determines if the hierarchical verification proposed by Navarro and Beaza-Yates is used or
 *                       not.
 * @tparam TMultiFinder  Specifies the algorithm for the multiple exact string matching algorithm.  Type: AhoCorasick.
 *
 * There are two defaults available: <tt>PexHierarchical</tt> and <tt>PexNonHiearchical</tt> (e.g.
 * <tt>Pattern&lt;CharString&gt;, PexHierarchical&gt;</tt> that both use the @link AhoCorasickPattern Aho-Corasick
 * @endlink algorithm for the multiple exact string matching.
 */

/*!
 * @class HierarchicalPexPattern
 * @extends PexPattern
 * @headerfile <seqan/find.h>
 * @brief By using this specialization, the hierarchircal verification is enabled.
 *
 * @signature template <typename TNeedle, typename TMultiFinder>
 *            class Pattern<TNeedle, Pex<Hierarchical, TMultiFinder> >;
 *
 * @tparam TNeedle      The needle type. Type: @link ContainerConcept @endlink
 * @tparam TMultiFinder Specifies the algorithm for th emultiple exact string matching algorithm.
 */

/*!
 * @class NonHierarchicalPexPattern
 * @extends PexPattern
 * @headerfile <seqan/find.h>
 * @brief By using this specialization, the hierarchircal verification is disabled.
 *
 * @signature template <typename TNeedle, typename TMultiFinder>
 *            class Pattern<TNeedle, Pex<NonHierarchical, TMultiFinder> >;
 *
 * @tparam TNeedle      The needle type. Type: @link ContainerConcept @endlink
 * @tparam TMultiFinder Specifies the algorithm for th emultiple exact string matching algorithm.
 */

template <typename TNeedle, typename TVerification, typename TMultiFinder>
class Pattern<TNeedle, Pex<TVerification, TMultiFinder > >:
    public FindBegin_<Pattern<TNeedle, Pex<TVerification, TMultiFinder > > >
{
 public:
   typedef typename Position<TNeedle>::Type TPosition;
   typedef unsigned TScore;
   typedef Pattern<TNeedle, MyersUkkonen > TVerifier;
   typedef typename PexMultiFinder<
                       Pattern<TNeedle, Pex<TVerification,TMultiFinder > >
                                  >::Type TMFinder;

   // the maximal accepted error
   TScore limit;
   // reference to the needle
   Holder<TNeedle> data_host;
   // pattern object for the multi pattern search
   TMFinder multiPattern;
   // needles for the multi pattern search
   String<Segment<TNeedle> >  splitted_needles;

   // data store for the verification tree respectively the splitted needle
   std::map<unsigned, PexRange_<TPosition,TScore,TVerifier,TNeedle> > range_table;
   // map leafs of the tree to parts of the needle
   std::map<unsigned, unsigned> leaf_map;

   // store the infixes for the verifiers
   String<Segment<TNeedle> > segment_store;

   // track position where the last occurrence was found
   unsigned lastFPos;
   unsigned lastFNdl;

   // indicator to track if we already found an occurrence
   bool findNext,patternNeedsInit;

   unsigned needleLength;

   Pattern() :
       limit(1), lastFPos(0), lastFNdl(0), findNext(false), patternNeedsInit(true)
   {}

   template <typename TNeedle2>
   Pattern(TNeedle2 const & ndl) :
       limit(1), lastFPos(0), lastFNdl(0), findNext(false), patternNeedsInit(true)
   {
     setHost(*this, ndl);
   }

   template <typename TNeedle2>
   Pattern(TNeedle2 const & ndl, int _limit = -1) :
       limit(-_limit), lastFPos(0), lastFNdl(0), findNext(false), patternNeedsInit(true)
   {
     setHost(*this, ndl);
   }
};

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TNeedle2, typename TVerification, typename TMultiFinder>
void setHost (Pattern<TNeedle, Pex<TVerification,TMultiFinder > > & me, TNeedle2 const & needle)
{
  // initialisation of the find-tree etc. will be done when patternInit
  // is called to assure that we already know the scoreLimit
  me.data_host = needle;
  me.needleLength = length(needle);
  me.findNext = false;
  me.patternNeedsInit = true;
}

template <typename TNeedle, typename TNeedle2, typename TVerification, typename TMultiFinder>
void setHost (Pattern<TNeedle, Pex<TVerification,TMultiFinder > > & me, TNeedle2 & needle)
{
  setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TVerification, typename TMultiFinder>
inline typename Host<Pattern<TNeedle, Pex<TVerification,TMultiFinder > > >::Type &
host(Pattern<TNeedle, Pex<TVerification,TMultiFinder > > & me)
{
SEQAN_CHECKPOINT
  return value(me.data_host);
}

template <typename TNeedle, typename TVerification, typename TMultiFinder>
inline typename Host<Pattern<TNeedle, Pex<TVerification,TMultiFinder > > const>::Type &
host(Pattern<TNeedle, Pex<TVerification,TMultiFinder > > const & me)
{
SEQAN_CHECKPOINT
  return value(me.data_host);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TMultiFinder>
int _getRoot(Pattern<TNeedle, Pex<NonHierarchical, TMultiFinder > > & me)
{
SEQAN_CHECKPOINT
  return length(me.splitted_needles);
}

template <typename TNeedle, typename TMultiFinder>
int _getRoot(Pattern<TNeedle, Pex<Hierarchical, TMultiFinder > > &)
{
SEQAN_CHECKPOINT
  return 1;
}

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn PexPattern#getScore
 * @headerfile <seqan/find.h>
 * @brief Score of the last found match in approximate searching.
 *
 * @signature TScoreValue getScore(pattern);
 *
 * @param[in] pattern A pex pattern that can be used for approximate searching.
 *
 * @return TScoreValue The score of the last match found using <tt>pattern</tt>.  If no match was found then the value
 *                     is undefined.
 */
template <typename TNeedle, typename TVerification, typename TMultiFinder>
int getScore(Pattern<TNeedle, Pex<TVerification,TMultiFinder > > & me)
{
SEQAN_CHECKPOINT
  return getScore(me.range_table[_getRoot(me)].verifier);
}

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn PexPattern#scoreLimit
 * @headerfile <seqan/find.h>
 * @brief The minimal score a match must reach in approximate searching.
 *
 * @signature TScoreValue scoreLimit(pattern);
 *
 * @param[in] pattern The pattern to query.
 *
 * @return TScoreValue The score limit value.
 */
template <typename TNeedle, typename TVerification, typename TMultiFinder>
inline int
scoreLimit(Pattern<TNeedle, Pex<TVerification,TMultiFinder > > const & me)
{
SEQAN_CHECKPOINT
  return - (int) me.limit;
}


//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn PexPattern#setSoreLimit
 * @headerfile <seqan/find.h>
 * @brief Set the minimal score a match must reach in approximate serach.
 *
 * @signature void setScoreLimit(pattern, limit);
 *
 * @param[in,out] pattern The pattern to set the limit for.
 * @param[in]     limit   The limit score value to set.
 *
 * @return TScoreValue The score limit value.
 */


template <typename TNeedle, typename TScoreValue,typename TVerification, typename TMultiFinder>
inline void
setScoreLimit(Pattern<TNeedle, Pex<TVerification,TMultiFinder > > & me,
              TScoreValue _limit)
{
SEQAN_CHECKPOINT
  me.patternNeedsInit = true;
  me.limit = (- _limit);
}

//////////////////////////////////////////////////////////////////////////////
//   PexNonHierarchical -- functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TFinder, typename TMultiFinder>
void _patternInit(Pattern<TNeedle, Pex<NonHierarchical, TMultiFinder > > &me, TFinder &)
{
SEQAN_CHECKPOINT
  typedef typename Position<TNeedle>::Type TPosition;
  typedef unsigned TScore;
  typedef Pattern<TNeedle,MyersUkkonen> TVerifier;

/*
  // split pattern
  unsigned k = me.limit + 1;
  unsigned seg_len = me.needleLength / k; //std::floor(me.needleLength/k);

  clear(me.splitted_needles);
  clear(me.range_table);
  clear(me.segment_store);
  unsigned s = 0;
  unsigned c = 0;
  unsigned i = 0;
  while(s < me.needleLength)
  {
    PexRange_<TPosition,TScore,TVerifier,TNeedle> pr;
    pr.start = s;
    pr.end = (c == me.limit ? me.needleLength : s + seg_len);
    pr.error = 0;

    insert(me.range_table,i,pr);
    appendValue(me.splitted_needles,infix(value(me.data_host),pr.start,pr.end));
    s += (c == me.limit ? me.needleLength : seg_len);
    ++c;
    ++i;
  }
*/
  //split pattern (improved)
  unsigned k = me.limit + 1;

  clear(me.splitted_needles);
  clear(me.range_table);
  clear(me.segment_store);
  unsigned int pos = 0;
  for (unsigned int i = 0; i < k; ++i)
  {
    PexRange_<TPosition,TScore,TVerifier,TNeedle> pr;
    pr.start = pos;
    pr.end = (pos = me.needleLength * (i + 1) / k);
    pr.error = 0;

    insert(me.range_table, i, pr);
    appendValue(me.splitted_needles, infix(value(me.data_host), pr.start, pr.end));
  }
//*/

  me.lastFPos = 0;
  me.lastFNdl = 0;

  // insert complete needle in range table to use the verifier
  appendValue(me.segment_store,infix(value(me.data_host),0,me.needleLength));
  PexRange_<TPosition,TScore,TVerifier,TNeedle> pr;
  pr.start = 0;
  pr.end = me.needleLength;
  pr.error = me.limit;
  setHost(pr.verifier,me.segment_store[0]);
  setScoreLimit(pr.verifier, - static_cast<int>(me.limit));
  insert(me.range_table,length(me.splitted_needles),pr);

  // init multipattern finder
  setHost(me.multiPattern,me.splitted_needles);

  me.patternNeedsInit = false;
  me.findNext = false;
  _findBeginInit(me, needle(me));

#ifdef SEQAN_DEBUG_PEX
  std::cout << " -------------------------------------------------  " << std::endl;
  std::cout << "                   PATTERN INIT                     " << std::endl;
  std::cout << "Needle:   " << value(me.data_host) << std::endl;
  std::cout << "|Needle|: " << me.needleLength << std::endl;
  std::cout << "seg_len:  " << me.needleLength / k << std::endl;
  std::cout << "limit:    " << me.limit << std::endl;
  std::cout << "k:        " << k << std::endl;
  std::cout << "computed following needles for multipattern search: " << std::endl;
  for(unsigned i = 0;i < length(me.splitted_needles);++i)  std::cout << me.splitted_needles[i] << std::endl;
  std::cout << " -------------------------------------------------  " << std::endl;
#endif

}

//////////////////////////////////////////////////////////////////////////////

template <typename TFinder, typename TNeedle, typename TMultiFinder>
inline bool find (TFinder & finder, Pattern<TNeedle, Pex<NonHierarchical, TMultiFinder > > & me)
{
SEQAN_CHECKPOINT

  typedef typename Host<TFinder>::Type    THost;
  typedef Segment<THost>                  THostSegment;
  typedef Finder<THostSegment>            THSFinder;
  TFinder mf(finder);
  unsigned startPos;

  if (empty(finder))
  {
    _finderSetNonEmpty(finder);
  }
  if(me.patternNeedsInit)
  {
     _patternInit(me, finder);
  }

  if(me.findNext){
    startPos = position(finder);
    int start = me.lastFPos - me.range_table[me.lastFNdl].start - me.limit;
    int end   = me.lastFPos + (me.needleLength - me.range_table[me.lastFNdl].start) + me.limit;

    // adjust start and end if they point over the edges of host(finder)
    start = (start < 0 ? 0 : start);
    end = (end > static_cast<int>(length(host(finder))) ? static_cast<int>(length(host(finder))) : end);

    THostSegment s(infix(host(finder),start,end));
    THSFinder f(s);

    while(find(f,me.range_table[_getRoot(me)].verifier))
    {
      unsigned nP = start + position(f);
      if(nP > startPos){
    // compute new position
    unsigned offset = nP - position(finder);
    finder += offset;
    me.findNext = true;
    _setFinderEnd(finder);
    return true;
      }
    }
    // reset mf finder to old position
    unsigned mf_offset = position(finder) - me.lastFPos;
    mf -= mf_offset;
  }
  me.findNext = false;
  startPos = position(finder);

  while(find(mf,me.multiPattern))
  {
    int s = position(mf) - me.range_table[position(me.multiPattern)].start - me.limit;
    int e   = position(mf) + (me.needleLength - me.range_table[position(me.multiPattern)].start) + me.limit;

    // adjust start and end if they point over the edges of host(finder)
    s = (s < 0 ? 0 : s);
    e = (e > static_cast<int>(length(host(finder))) ? static_cast<int>(length(host(finder))) : e);

    THostSegment i(infix(host(mf),s,e));
    THSFinder f(i);
    while(find(f,me.range_table[_getRoot(me)].verifier))
    {
      unsigned nP = s + position(f);
      if(nP > startPos){
    // compute new position
    unsigned offset = nP - position(finder);
    finder += offset;
    me.lastFPos = position(mf);
    me.lastFNdl = position(me.multiPattern);
    me.findNext = true;
    _setFinderEnd(finder);
    return true;
      }
    }
  }
  // set finder to end position
  unsigned t = length(host(finder))- position(finder);
  finder += t;

  return false;
}

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//   PexHierarchical -- functions
//////////////////////////////////////////////////////////////////////////////
/*
template <typename TNeedle, typename TMultiFinder>
void _createTree(Pattern<TNeedle, Pex<Hierarchical, TMultiFinder > > &me, unsigned start, unsigned end,
         unsigned k, unsigned parent, unsigned direction ,unsigned idx, unsigned plen)
{
  //create tree like proposed in Navarro & Raffinot
  // direction == 0 .. choose left child in the tree
  // direction == 1 .. choose right child in the tree

#ifdef SEQAN_DEBUG_PEX
  std::cout << "called _createTree:" << std::endl;
  std::cout << "  start: " << start << std::endl;
  std::cout << "  end  : " << end << std::endl;
  std::cout << "  seq  : " << infix(value(me.data_host),start,end + 1) << std::endl;
  std::cout << "  k    : " << k << std::endl;
  std::cout << "  paren: " << parent << std::endl;
  std::cout << "  direc: " << direction << std::endl;
  std::cout << "  idx  : " << idx << std::endl;
  std::cout << "  plen : " << plen << std::endl;
  std::cout << " ----------------------------- " << std::endl;
#endif
  typedef typename Position<TNeedle>::Type TPosition;
  typedef unsigned TScore;
  typedef Pattern<TNeedle,MyersUkkonen> TVerifier;

  PexRange_<TPosition,TScore,TVerifier,TNeedle> pr;
  pr.start = start;
  pr.end = end;
  pr.error = k;

  appendValue(me.segment_store,infix(value(me.data_host),pr.start,pr.end + 1));
  setScoreLimit(pr.verifier, - static_cast<int>(pr.error));
  setHost(pr.verifier, me.segment_store[length(me.segment_store) - 1]);

  unsigned left = k/2 + 1; //std::ceil(static_cast<double>(k + 1)/2);
  unsigned cur_idx = (parent << 1) + direction;

  // insert pr into the tree
  insert(me.range_table,cur_idx,pr);

  if(k == 0){
    appendValue(me.splitted_needles,infix(value(me.data_host),pr.start,pr.end + 1));
#ifdef SEQAN_DEBUG_PEX
    std::cout << "inserted : " << me.splitted_needles[length(me.splitted_needles) - 1] << " into splitted needles" << std::endl;
    std::cout << "assign to leaf_map " << length(me.splitted_needles) - 1 << " value " << cur_idx << std::endl;
    std::cout << " ----------------------------- " << std::endl;
#endif
    me.leaf_map[length(me.splitted_needles) - 1] = cur_idx;
  }else{
    // recusivly create the rest of the tree
//    _createTree(me, start, start + left * plen - 1, std::floor(static_cast<double>(left * k)/ static_cast<double>(k + 1)),cur_idx,0,idx,plen);
    _createTree(me, start, start + left * plen - 1, left * k / (k+1),cur_idx,0,idx,plen);
//    _createTree(me,  start + left * plen, end, std::floor(static_cast<double>((k + 1 - left)*k)/ static_cast<double>(k + 1)),cur_idx,1,idx + left,plen);
    _createTree(me,  start + left * plen, end, (k + 1 - left)*k / (k+1),cur_idx,1,idx + left,plen);
  }
}
*/
template <typename TNeedle, typename TMultiFinder>
void _createTree(Pattern<TNeedle, Pex<Hierarchical, TMultiFinder > > &me,
                 unsigned start, unsigned end,
                 unsigned k,
                 unsigned parent,
                 unsigned direction,
                 unsigned idx,
                 unsigned plen)
{
  typedef typename Position<TNeedle>::Type TPosition;
  typedef unsigned TScore;
  typedef Pattern<TNeedle,MyersUkkonen> TVerifier;

  PexRange_<TPosition,TScore,TVerifier,TNeedle> pr;
  pr.start = start;
  pr.end = end;
  pr.error = k;

  appendValue(me.segment_store,infix(value(me.data_host),pr.start,pr.end + 1));
  setScoreLimit(pr.verifier, - static_cast<int>(pr.error));
  setHost(pr.verifier, me.segment_store[length(me.segment_store) - 1]);

  unsigned cur_idx = (parent << 1) + direction;

  // insert pr into the tree
  insert(me.range_table,cur_idx,pr);

  if(k == 0)
  {
    appendValue(me.splitted_needles,infix(value(me.data_host),pr.start,pr.end + 1));
    me.leaf_map[length(me.splitted_needles) - 1] = cur_idx;
  }
  else
  {
      unsigned int lower_2power = 1 << log2(k+1);
      unsigned int len = end - start+1;
      unsigned int right_k = lower_2power/2-1;
      unsigned int left_k = k - right_k-1;
      unsigned int left_len = len * (left_k+1) / (k+1);
      _createTree(me, start, start + left_len-1, left_k, cur_idx, 0, idx,plen);
      _createTree(me, start + left_len, end, right_k, cur_idx, 1, idx + (left_k+1),plen);
  }
}

template <typename TNeedle, typename TFinder, typename TMultiFinder>
void _patternInit(Pattern<TNeedle, Pex<Hierarchical, TMultiFinder > > &me, TFinder &)
{
  unsigned k = me.limit + 1;
  unsigned plen = me.needleLength / k; //std::floor(static_cast<double>(me.needleLength)/static_cast<double>(k));

  // reset
  clear(me.splitted_needles);
  clear(me.range_table);
  clear(me.leaf_map);
  clear(me.segment_store);

  // build the verification tree
  _createTree(me, 0, me.needleLength - 1,me.limit, 0, 1 , 0, plen);

  me.lastFPos = 0;
  me.lastFNdl = 0;
  setHost(me.multiPattern, me.splitted_needles);
  me.patternNeedsInit = false;
  me.findNext = false;

  _findBeginInit(me, needle(me));

#ifdef SEQAN_DEBUG_PEX
  std::cout << " -------------------------------------------------  " << std::endl;
  std::cout << "                   PATTERN INIT                     " << std::endl;
  std::cout << "Needle:   " << value(me.data_host) << std::endl;
  std::cout << "|Needle|: " << me.needleLength << std::endl;
  std::cout << "limit:    " << me.limit << std::endl;
  std::cout << "k:        " << k << std::endl;
  std::cout << "computed following needles for multipattern search: " << std::endl;
  for(unsigned i = 0;i < length(me.splitted_needles);++i)  std::cout << me.splitted_needles[i] << std::endl;
  std::cout << " -------------------------------------------------  " << std::endl;
#endif
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFinder, typename TNeedle, typename TMultiFinder>
inline bool find (TFinder & finder, Pattern<TNeedle, Pex<Hierarchical, TMultiFinder > > & me)
{
SEQAN_CHECKPOINT

  typedef typename Host<TFinder>::Type    THost;
  typedef Segment<THost>                  THostSegment;
  typedef Finder<THostSegment>            THSFinder;
  TFinder mf(finder);
  unsigned startPos;

  if (empty(finder))
  {
    _finderSetNonEmpty(finder);
  }
  if(me.patternNeedsInit)
  {
     _patternInit(me, finder);
  }
  if(me.findNext){
    // we found an occurrence
    startPos = position(finder);
    unsigned pnode = _getRoot(me); // use root
    unsigned in = me.range_table[me.leaf_map[me.lastFNdl]].start;

    int p1 = me.lastFPos - (in - me.range_table[pnode].start) - me.range_table[pnode].error;
    int p2 = me.lastFPos + (me.range_table[pnode].end - in + 1) + me.range_table[pnode].error;

    // adjust start and end if they point over the edges of host(finder)
    p1 = (p1 < 0 ? 0 : p1);
    p2 = (p2 > static_cast<int>(length(host(finder))) ? static_cast<int>(length(host(finder))) : p2);
    THostSegment i(infix(host(mf),p1,p2));
    THSFinder f(i);

    while(find(f,me.range_table[pnode].verifier))
    {
      unsigned nP = p1 + position(f);
      if(nP > startPos)
    {
      // compute new position
      unsigned offset = nP - position(finder);
      finder += offset;
      me.findNext = true;
      _setFinderEnd(finder);
      return true;
    }
    }
    // reset mf finder to old position
    unsigned mf_offset = position(finder) - me.lastFPos;
    mf -= mf_offset;
  }
  me.findNext = false;
  startPos = position(finder);

  while(find(mf,me.multiPattern))
  {
    // get found leaf
    unsigned node = me.leaf_map[position(me.multiPattern)];
    unsigned in = me.range_table[node].start;
    node = node >> 1;
    bool cand = true;

    while( cand && node != 1) // stop when reaching root
    {
      int p1 = position(mf) - (in - me.range_table[node].start) - me.range_table[node].error;
      int p2 = position(mf) + (me.range_table[node].end - in + 1) + me.range_table[node].error;

      // adjust start and end if they point over the edges of host(finder)
      p1 = (p1 < 0 ? 0 : p1);
      p2 = (p2 > static_cast<int>(length(host(finder))) ? static_cast<int>(length(host(finder))) : p2);
      THostSegment i(infix(host(mf),p1,p2));
      THSFinder f(i);
      cand = find(f,me.range_table[node].verifier);
      node = node >> 1;
    }
    // if we verfied till here .. verify the complete pattern
    if(cand){
      // we found an occurrence
      node = _getRoot(me); // use root
      int p1 = position(mf) - (in - me.range_table[node].start) - me.range_table[node].error;
      int p2 = position(mf) + (me.range_table[node].end - in + 1) + me.range_table[node].error;

      // adjust start and end if they point over the edges of host(finder)
      p1 = (p1 < 0 ? 0 : p1);
      p2 = (p2 > static_cast<int>(length(host(finder))) ? static_cast<int>(length(host(finder))) : p2);
      THostSegment i(infix(host(mf),p1,p2));
      THSFinder f(i);
      while(find(f,me.range_table[node].verifier))
      {
    unsigned nP = p1 + position(f);
    if(nP > startPos)
    {
      // compute new position
      unsigned offset = nP - position(finder);
      finder += offset;
      me.lastFPos = position(mf);
      me.lastFNdl = position(me.multiPattern);
      me.findNext = true;
      _setFinderEnd(finder);
      return true;
    }
      }
    }
  }
  // nothing more to find -> set finder to end position
  unsigned t = length(host(finder))- position(finder);
  finder += t;

  return false;
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_..

