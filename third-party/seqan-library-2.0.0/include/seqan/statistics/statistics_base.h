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

#ifndef SEQAN_STATISTICS_STATISTICS_BASE_H_
#define SEQAN_STATISTICS_STATISTICS_BASE_H_

namespace seqan
{

template <typename TAlgorithm, typename TFloat, typename TAlphabet>
void _numOccurrences(TFloat &nW, String<TAlphabet>& haystack, StringSet<String<TAlphabet> >& needle, TAlgorithm const &);

/*
 * @fn _zscore
 * @headerfile <seqan/statistics.h>
 * @brief Auxiliary function to compure the z-score index for a set of patterns w.r.t. a set of text strings and a
 *        MarkovMovel.
 *
 * @signature TFloat _zscore(W, X, M, tag);
 *
 * @param[in] W   The StringSet to use as the words.
 * @param[in] X   The StringSet to use as the text.
 * @param[in] M   The MarkovModel to use.
 * @param[in] tag The tag to select the algorithm.
 *
 * @return TFloat The z-score, TFloat is the TFloat from the MarkovModel M.
 */

template <typename TAlgorithm, typename TFloat,  typename TStringSet, typename TAlphabet, typename TSpec>
TFloat _zscore(TStringSet W,  TStringSet& X, MarkovModel<TAlphabet, TFloat, TSpec> & M, TAlgorithm const &)
{


    TFloat z_score=0;
    TFloat nW=0;
    //compute occurrances
    for(unsigned int i=0; i< length(X); i++)
    {
        String<TAlphabet> temp = getValueById(X, i);
        _numOccurrences(nW, temp, W, TAlgorithm());
    }

    //compute expectation
    TFloat E = expectation(W, X, M);
//std::cout<<"\nE:"<<E;
    //compute variance
    TFloat V = _computeVariance(W, X, M, E);
//std::cout<<"\nV:"<<V;
    //compute z-score
    z_score=(nW-E)/sqrt(V);

    return z_score;

}

/*
 * @fn _numOccurences
 * @headerfile <seqan/statistics.h>
 * @brief Auxiliary function to compute the number of occurences of a set of patterns in a set of text strings.
 *
 * @signature void _numOccurences(W, haystack, needle, algoTag);
 *
 * @param[in,out] W        A counter, incremented by the number of occurences.
 * @param[in]     haystack The text strings.
 * @param[in]     needle   The set of patterns.
 * @param[in]     algoTag  The tag to select the online text search algorithm with.
 */

//Fixed to  AhoCorasick in original code, reason???
template <typename TAlgorithm, typename TFloat, typename TAlphabet>
void _numOccurrences(TFloat &nW, String<TAlphabet> &haystack, StringSet<String<TAlphabet> > &needle, TAlgorithm const &)
{
    SEQAN_CHECKPOINT;
    Finder<String<TAlphabet> > finder(haystack);
    Pattern<StringSet<String<TAlphabet> >, TAlgorithm> pattern(needle);
    while (find(finder, pattern))
    {
        nW++;
    }
}

/*
 * @fn _computeExpecttation
 * @brief Auxiliary function to compute the expectation for a set of patterns w.r.t. a text string  and a MarkovModel.
 *
 * @signature TFloat _computeExpectation(mm, W, n);
 *
 * @param[in] mm The MarkovModel.
 * @param[in] W  The set of patterns, StringSet.
 * @param[in] n  The length of the string, unsigned.
 *
 * @return TFloat the expectation valuefor W w.r.t. a string and M.
 */

template <typename TAlphabet, typename TFloat, typename TSpec>
TFloat _computeExpectation(MarkovModel<TAlphabet, TFloat, TSpec> &mm,
                     StringSet<String<TAlphabet> > &W, unsigned int n)
{
    TFloat E=0;
    for (unsigned int i=0; i<length(W); i++){
        String<TAlphabet> temp = getValueById(W, i);
        E += (n - length(temp) + 1)*mm.emittedProbability(temp);
    }
    return E;
}


/*
 * @fn _computeVariance
 * @headerfile <seqan/statistics.h>
 * @brief Auxiliary function to compute the variance for a set of patterns w.r.t. a set of text strings and a MarkovModel.
 *
 * @signature TFloat _computeVariance(W, X, M);
 *
 * @param[in] W The set of words.
 * @param[in] X The text set.
 * @param[in] M The MarkovModel to use.
 * @param[in] E
 *
 * @return TFloat The variance for W w.r.t. X and M.
 *
 * @section Remarks
 *
 * If the alphabet is Dna, then the suitable correction factors are computed.
 */

// TODO(holtgrew): W could become const-ref.

template <typename TFloat, typename TAlphabet, typename TSpec>
TFloat _computeVariance( StringSet<String<TAlphabet> > W,  StringSet<String<TAlphabet> > &X, MarkovModel<TAlphabet, TFloat, TSpec> &M, TFloat &E)
{
    //V=B+2C-E^2
    TFloat V = E;

    //C=D+A

    //compute A and D

    TFloat A = 0;
    TFloat D = 0;
    TFloat tmpA, eQPPPe, eQPPQPPe;
    unsigned int sizeW=length(W);
    unsigned int n;

    String <TFloat> pStar;
    resize(pStar, sizeW, 0);

    Shape<TAlphabet, SimpleShape> orderShape;
    resize(orderShape, M.order);

    for(unsigned int j=0; j<sizeW; j++){
        String<TAlphabet> string =getValueById(W, j);

        int row = hash(orderShape,begin(string));
        TFloat p = 1;
        for(unsigned int i=1; i<length(string)-M.order+1; i++)
        {
            int column=hash(orderShape,begin(string)+i);
            p*=value(M.transition,row,column);
            row = column;
        }
        value(pStar, j) = p;
    }



    for(unsigned int z=0; z<length(X); z++){

      for(unsigned int i=0; i<length(X); i++){

         n = length(getValueById(X, i));

         for(unsigned int j=0; j<sizeW; j++){

            String<TAlphabet> Wj =getValueById(W, j);

            TFloat q = (TFloat) (n-(2*length(Wj))+2);

            for(unsigned int k=0; k<sizeW; k++){

                tmpA=value(pStar,j)*value(pStar,k);

                unsigned int jfirst, jlast, kfirst;

                jfirst = hash(orderShape,begin(Wj));

                jlast = hash(orderShape,end(Wj)-M.order);

                kfirst = hash(orderShape,begin(getValueById(W, k)));

                eQPPPe = value(M._qppp, jlast,kfirst);

                eQPPQPPe = value(M._qppqpp, jlast,kfirst);

                tmpA  *= value(M.stationaryDistribution, jfirst) * ((q*(q+1)/2)* value(M.stationaryDistribution, kfirst) - (q-1)*eQPPPe - eQPPQPPe);

                A += tmpA;
            }
         }
      }

      // Compute D
      D+= _overlapExpectation(W,M,length(getValueById(X, z)));
    }




    //Compute Variance
    V += (2*A) + (2*D) -  std::pow((double) E, (int) 2);

    //return V;
    return V;
}

/*
 * @fn _overlapExpectation
 * @headerfile <seqan/statistics.h>
 * @brief Auxiliary function necessary when correction factors have to be computed.
 *
 * @signature TFloat _overlapExpectation(W, M, n);
 *
 * @param[in] W The set of words.
 * @param[in] M The MarkovModel to use for the computation.
 * @param[in] n The length of the text.
 *
 * @return TFloat The expectation value for overlapping, TFloat is the TFloat from the type of M.
 */

template <typename TFloat, typename TAlphabet, typename TSpec>
TFloat _overlapExpectation(StringSet<String<TAlphabet> > W, MarkovModel<TAlphabet, TFloat, TSpec> &M, unsigned int n)
{
    TFloat E_overlap = 0;
    unsigned int sizeW = length(W);
    for(unsigned int i=0; i<sizeW; i++)
    {
        String<TAlphabet> patt1 = getValueById(W, i);
        unsigned int size1 = length(patt1);
        for(unsigned int j=0; j<sizeW; j++)
        {
            String<TAlphabet> patt2 = getValueById(W, j);
            unsigned int k=1;
            unsigned int size2 = length(patt2);
            if(size1>size2)
            {
                k = size1 - size2 + 1;
            }
            for(; k<size1; k++)
            {
                if(isEqual(infix(patt1,begin(patt1)+k,end(patt1)),infix(patt2,begin(patt2),begin(patt2)+k-1)))
                {
                    String<TAlphabet> temp = infix(patt1, begin(patt1),begin(patt1)+k-1);
                    append(temp,infix(patt2,begin(patt2),end(patt2)));
                    E_overlap += (n - size1 + 1)*M.emittedProbability(temp);
                }
            }
        }
    }
    return E_overlap;
}

/*
 * @fn _addReverseComplemenents
 * @headerfile <seqan/statistics.h>
 * @brief Computes the reverse complemenets of a set of strings in the input.
 *
 * @signature void _addReverseComplements(ss);
 *
 * @param[in,out] ss A String set to expand.
 */

template <typename TAlphabet>
void _addReveseComplements(StringSet<String<TAlphabet> > &stringSet)
{
    unsigned int num= length(stringSet);

    for(unsigned int i=0; i< num; i++){
           DnaStringReverseComplement mycom(getValueById(stringSet, i));
         appendValue(stringSet, mycom);
    }
}


///////////////////////////////////////////////////////////////////////
// Extern functions to be provided by SeqAn
///////////////////////////////////////////////////////////////////////

typedef Dna TDnaAlphabet;
typedef String<TDnaAlphabet> TDnaSequence;

/*!
 * @fn zscore
 * @headerfile <seqan/statistics.h>
 * @brief Computes the z-score index for a set of patterns w.r.t. a set of text strings and a MarkovModel.
 *
 * @signature TFloat zscore(W, X, M, algoTag);
 *
 * @param[in] W       The StringSet of pattern strings.
 * @param[in] X       The StringSet of text strings.
 * @param[in] M       The MarkovModel object.
 * @param[in] algoTag The algorithm to exploit to compute the number of occurrences of patterns in the text strings
 *                    (see @link AhoCorasickPattern @endlink etc.).
 *
 * @return TFloat The z-score for W w.r.t. X and M, TFloat is the TFloat from the type of M.
 */

template <typename TAlgorithm, typename TFloat, typename TSpec, typename TStringSet, typename TAlphabet>
TFloat zscore(TStringSet W,  TStringSet &X, MarkovModel<TAlphabet, TFloat, TSpec> &M, TAlgorithm const & algorithmTag)
{
    ensureAuxMatrices(M);
       return _zscore(W,X,M, algorithmTag);
}

template <typename TAlgorithm, typename TFloat, typename TSpec, typename TDnaSequence>
TFloat zscore(StringSet<TDnaSequence> W,  StringSet<TDnaSequence> &X, MarkovModel<Dna, TFloat, TSpec> &M, TAlgorithm const &)
{
   //add-reverse complements
   _addReveseComplements(W);

    ensureAuxMatrices(M);

   TFloat z_score=0;
   TFloat nW=0;
   //compute occurrences
   for(unsigned int i=0; i < length(X); i++)
   {
         String<Dna> temp = getValueById(X, i);
        _numOccurrences(nW, temp, W, TAlgorithm());
    }

    //compute expectation
    TFloat E = expectation(W, X, M);
    //std::cout<<"\nE: "<<E<<"\n";
    //compute variance
    TFloat V = _computeVariance(W, X, M, E);
    //std::cout<<"\nV: "<<V<<"\n";
    //compute correction factor
    TFloat correction = 0;

    unsigned int n;
    unsigned int sizeW= length(W);

    for(unsigned int j=0; j<length(X); j++){

         n = length(getValueById(X, j));

        for(unsigned int i=0; i<sizeW; i++)
        {
            String<Dna> patt = getValueById(W, i);
            DnaStringReverseComplement revpatt(patt);
            String<Dna> revc= revpatt;
            if (isEqual(patt,revc))
            {
                correction += (n-length(patt)+1)*M.emittedProbability(revc);
            }
        }
    }

    V+= correction;

    //compute z-score
    z_score=(nW-E)/sqrt(V);
    //std::cout<<"\nnW: "<<nW<<"\n";
    //std::cout<<"\nZ: "<<z_score<<"\n";
    return z_score;
}

/*!
 * @fn variance
 * @headerfile <seqan/statistics.h>
 * @brief Computes the variance for a set of patterns w.r.t a set of text strings and a MarkovModel.
 *
 * @signature TFloat variance(W, X, M);
 *
 * @param[in] W The StringSet of word strings.
 * @param[in] X The StringSet of text strings.
 * @param[in] M The MarkovModel to use.
 *
 * @return TFloat The variance for W w.r.t. X and M, TFloat is the TFloat from the type of M.
 */

template <typename TFloat, typename TAlphabet, typename TSpec>
TFloat variance(StringSet<String<TAlphabet> > &W, StringSet<String<TAlphabet> >& X, MarkovModel<TAlphabet, TFloat, TSpec> & M)
{
   TFloat E = expectation(W, X, M);

   return _computeVariance(W,X,M,E);
}

//Special case for DNA sequences, reverse complement sequences are added
template <typename TFloat, typename TSpec>
TFloat variance(StringSet<String<Dna> > W, StringSet<String<Dna> > &X, MarkovModel<Dna, TFloat, TSpec> & M)
{

   //add-reverse complements
    _addReveseComplements(W);

    TFloat E = expectation(W, X, M);

    TFloat var =  _computeVariance(W,X,M,E);

    //compute correction factor
    TFloat correction = 0;

    unsigned int n;
    unsigned int sizeW= length(W);


   for(unsigned int j=0; j<length(X); j++){

         n = length(getValueById(X, j));

        for(unsigned int i=0; i<sizeW; i++)
        {
            String<Dna> patt = getValueById(W, i);
            DnaStringReverseComplement revpatt(patt);
            String<Dna> revc= revpatt;
            if (isEqual(patt,revc))
            {
                correction += (n-length(patt)+1)*M.emittedProbability(revc);
            }
        }
    }
    var+=correction;

  return var;
}

/*!
 * @fn expectation
 * @headerfile <seqan/statistics.h>
 * @brief Computes the expectation for a set of patterns w.r.t a set of text strings and a MarkovModel.
 *
 * @signature TFloat expectation(W, X, M);
 *
 * @param[in] W The StringSet of word strings.
 * @param[in] X The StringSet of text strings.
 * @param[in] M The MarkovModel to use.
 *
 * @return TFloat The expectation for W w.r.t. X and M, TFloat is the TFloat from the type of M.
 */

// TODO(holtgrew): Add const?

template <typename TAlphabet, typename TFloat, typename TSpec>
TFloat expectation(StringSet<String<TAlphabet> > & W, StringSet<String<TAlphabet> > &X, MarkovModel<TAlphabet, TFloat, TSpec> &M)
{
    unsigned int n;
    TFloat E = 0;

    for(unsigned int i=0; i<length(X); i++){
         n = length(getValueById(X, i));
        E += _computeExpectation(M, W, n);
    }

    return E;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STATISTICS_STATISTICS_BASE_H_
