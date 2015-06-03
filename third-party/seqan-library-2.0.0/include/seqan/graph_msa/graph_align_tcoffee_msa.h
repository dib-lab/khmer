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
// Author: Tobias Rausch <rausch@embl.de>
// ==========================================================================

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_MSA_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_MSA_H

namespace seqan
{

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class MsaOptions
// --------------------------------------------------------------------------

/*!
 * @class MsaOptions
 * @headerfile <seqan/graph_msa.h>
 * @brief Configuration for @link globalMsaAlignment @endlink
 *
 * @signature template <typename TAlphabet, typename TScore>
 *            struct MsaOptions;
 *
 * @tparam TAlphabet The alphabet type to use.
 * @tparam TScore    The @link Score @endlink type to use.
 */

template <typename TAlphabet, typename TScore>
struct MsaOptions
{
public:
    /*!
     * @var bool MsaOptions::rescore
     * @brief Whether or not to rescore segment matches after refinment (default: <tt>true</tt>).
     */
    bool rescore;

    /*!
     * @var unsigned MsaOptions::outputFormat;
     * @brief Output Format (<tt>0</tt> - FASTA, <tt>1</tt> MSF, default: <tt>0</tt>).
     */
    unsigned outputFormat;

    /*!
     * @var TScore MsaOptions::sc;
     * @brief The @link Score @endlink object to use.
     */
    TScore sc;

    /*!
     * @var unsigned MsaOptions::build;
     * @brief Methods for computing guide tre.
     *
     * 0 Neighbor-joining, 1 UPGMA single linkage, 2 UPGMA complete linkage,
     * 3 UPGMA average linkage, 4 UPGMA weighted average linkage.
     */
    unsigned build;

    /*!
     * @var TUnsignedString MsaOptions::method;
     * @brief All methods to compute segment matches, of Type <tt>String&lt;unsigned&gt;</tt>.
     *
     * 0 global alignment, 1 local alignments, 2 overlap alignments, 3 longest common subsequence
     */
    String<unsigned> method;

    // TODO(holtgrew): Document me!
    // Various input and output file names
    String<std::string> alnfiles;       // External alignment files
    String<std::string> libfiles;       // T-Coffee library files
    String<std::string> blastfiles;     // Blast match files
    String<std::string> mummerfiles;    // MUMmer files
    std::string outfile;                // Output file name
    std::string seqfile;                // Sequence file name
    std::string infile;                 // Alignment file for alignment evaluation
    std::string treefile;               // Guide tree file

    /*!
     * @fn MsaOptions::MsaOptions
     * @brief Default constructor.
     * @signature MsaOptions::MsaOptions();
     */
    MsaOptions() : rescore(true), outputFormat(0), build(0)
    {}
};

// --------------------------------------------------------------------------
// Function evaluateAlignment()
// --------------------------------------------------------------------------

// TODO(holtgrew): Document me!

template <typename TAlphabet, typename TScore>
void evaluateAlignment(MsaOptions<TAlphabet, TScore> const & msaOpt)
{
    typedef typename Value<TScore>::Type TScoreValue;
    typedef String<TAlphabet> TSequence;
    typedef typename Size<TSequence>::Type TSize;
    StringSet<TSequence, Owner<> > origStrSet;
    StringSet<String<char> > names;

    // Read the sequences
    std::fstream strm;
    strm.open(msaOpt.infile.c_str(), std::ios_base::in | std::ios_base::binary);
    read(strm, origStrSet, names, FastaAlign());
    strm.close();

    // Make a dependent StringSet
    typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
    TDepSequenceSet strSet(origStrSet);

    // Read the alignment
    typedef String<Fragment<> > TFragmentString;
    String<TScoreValue> scores;
    TFragmentString matches;
    std::fstream strm_lib;
    strm_lib.open(msaOpt.infile.c_str(), std::ios_base::in | std::ios_base::binary);
    read(strm_lib, matches, scores, names, FastaAlign());
    strm_lib.close();

    // Build the alignment graph
    typedef Graph<Alignment<TDepSequenceSet, TSize> > TGraph;
    TGraph g(strSet);
    buildAlignmentGraph(matches, g, FrequencyCounting());

    // Print the scoring information
    TScoreValue gop = msaOpt.sc.data_gap_open;
    TScoreValue gex = msaOpt.sc.data_gap_extend;
    std::cout << "Scoring parameters:" << std::endl;
    std::cout << "*Gap opening: " << gop << std::endl;
    std::cout << "*Gap extension: " << gex << std::endl;
    std::cout << "*Scoring matrix: " << std::endl;
    TSize alphSize = ValueSize<TAlphabet>::VALUE;
    std::cout << "   ";
    for (TSize col = 0; col < alphSize; ++col)
        std::cout << TAlphabet(col) << ',';
    std::cout << std::endl;
    for (TSize row = 0; row < alphSize; ++row)
    {
        for (TSize col = 0; col < alphSize; ++col)
        {
            if (col == 0)
                std::cout << TAlphabet(row) << ": ";
            std::cout << score(msaOpt.sc, TAlphabet(row), TAlphabet(col));
            if (col < alphSize - 1)
                std::cout << ',';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // Print the alignment information
    TSize numGapEx = 0;
    TSize numGap = 0;
    TSize numPairs = 0;
    TSize alignLen = 0;
    String<TSize> pairCount;
    String<char> mat;
    if (convertAlignment(g, mat))
    {
        TScoreValue alignScore = alignmentEvaluation(g, msaOpt.sc, numGapEx, numGap, numPairs, pairCount, alignLen);
        std::cout << "Alignment Score: " << alignScore << std::endl;
        std::cout << "Alignment Length: " << alignLen << std::endl;
        std::cout << "#Match-Mismatch pairs: " << numPairs << std::endl;
        std::cout << "Score contribution by match-mismatch pairs: " << (alignScore - (((TScoreValue) numGap * gop) + ((TScoreValue) numGapEx * gex))) << std::endl;
        std::cout << "#Gap extensions: " << numGapEx << std::endl;
        std::cout << "Score contribution by gap extensions: " << ((TScoreValue) numGapEx * gex) << std::endl;
        std::cout << "#Gap openings: " << numGap << std::endl;
        std::cout << "Score contribution by gap openings: " << ((TScoreValue) numGap * gop) << std::endl;
        std::cout << std::endl;
        std::cout << "#Pairs: " << std::endl;
        std::cout << "   ";
        for (TSize col = 0; col < alphSize; ++col)
            std::cout << TAlphabet(col) << ',';
        std::cout << std::endl;
        for (TSize row = 0; row < alphSize; ++row)
        {
            for (TSize col = 0; col < alphSize; ++col)
            {
                if (col == 0)
                    std::cout << TAlphabet(row) << ": ";
                std::cout << value(pairCount, row * alphSize + col);
                if (col < alphSize - 1)
                    std::cout << ',';
            }
            std::cout << std::endl;
        }
    }
    else
    {
        std::cout << "No valid alignment!" << std::endl;
    }
}

// --------------------------------------------------------------------------
// Function globalMsaAlignment()
// --------------------------------------------------------------------------

template <typename TStrSpec, typename TSpec, typename TList, typename TScore, typename TSegmentMatches, typename TScores>
void _appendSegmentMatches(StringSet<String<AminoAcid, TStrSpec>, Dependent<TSpec> > const & str,
                             TList const & pList,
                             TScore const &,
                             TSegmentMatches & matches,
                             TScores & scores)
{
    Blosum62 local_score(-1, -8);
    appendSegmentMatches(str, pList, local_score, matches, scores, LocalPairwiseLibrary());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TStrSpec, typename TSpec, typename TList, typename TScore, typename TSegmentMatches, typename TScores>
void _appendSegmentMatches(StringSet<String<TValue, TStrSpec>, Dependent<TSpec> > const & str,
                             TList const & pList,
                             TScore const & score_type,
                             TSegmentMatches & matches,
                             TScores & scores)
{
    appendSegmentMatches(str, pList, score_type, matches, scores, LocalPairwiseLibrary());
}

/*!
 * @fn globalMsaAlignment
 * @headerfile <seqan/graph_msa.h>
 * @brief Compute a global multiple sequence alignment.
 *
 * @signature void globalMsaAlignment(align, score);
 * @signature void globalMsaAlignment(gAlign, score[, options]);
 *
 * @param[in,out] gAlign  An @link AlignmentGraph @endlink containing two or more sequences.
 * @param[in,out] align   A @link Align @endlink object with two or more sequences to align.
 * @param[in]     score   The @link Score @endlink to use for computing the alignment.
 * @param[in]     options The @link MsaOptions @endlink to use for the configuration.
 *
 * The resulting alignment is stored in <tt>align</tt>/<tt>gAlign</tt>.
 */

template <typename TStringSet, typename TCargo, typename TSpec, typename TStringSet1, typename TNames, typename TAlphabet, typename TScore>
void globalMsaAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> > & gAlign,
                        TStringSet1 & sequenceSet,
                        TNames & sequenceNames,
                        MsaOptions<TAlphabet, TScore> const & msaOpt)
{
    typedef typename Value<TScore>::Type TScoreValue;
    typedef typename Size<TStringSet>::Type TSize;
    typedef Graph<Alignment<TStringSet, TSize> > TGraph;
    //typedef typename Id<TGraph>::Type TId;
    typedef double TDistanceValue;

    // Initialize alignment object
    clear(gAlign);
    assignStringSet(gAlign, sequenceSet);

    // Some alignment constants
    TStringSet & seqSet = stringSet(gAlign);
    TSize nSeq = length(seqSet);
    TSize threshold = 30;

    // Select all possible pairs for global and local alignments
    String<TSize> pList;
    selectPairs(seqSet, pList);

    // Set-up a distance matrix
    typedef String<TDistanceValue> TDistanceMatrix;
    TDistanceMatrix distanceMatrix;

    // Containers for segment matches and corresponding scores
    typedef String<Fragment<> > TFragmentString;
    TFragmentString matches;
    typedef String<TScoreValue> TScoreValues;
    TScoreValues scores;

    // Include segment matches from subalignments
    if (!empty(msaOpt.alnfiles))
    {
        typedef typename Iterator<String<std::string> const, Standard>::Type TIter;
        TIter begIt = begin(msaOpt.alnfiles, Standard());
        TIter begItEnd = end(msaOpt.alnfiles, Standard());
        for (; begIt != begItEnd; goNext(begIt))
        {
            std::ifstream strm_lib;
            strm_lib.open((*begIt).c_str(), std::ios_base::in | std::ios_base::binary);
            read(strm_lib, matches, scores, sequenceNames, FastaAlign());
            strm_lib.close();
        }
    }

    // Include computed segment matches
    if (!empty(msaOpt.method))
    {
        typedef typename Iterator<String<unsigned int> const, Standard>::Type TIter;
        TIter begIt = begin(msaOpt.method, Standard());
        TIter begItEnd = end(msaOpt.method, Standard());
        for (; begIt != begItEnd; goNext(begIt))
        {
            if (*begIt == 0)
                appendSegmentMatches(seqSet, pList, msaOpt.sc, matches, scores, distanceMatrix, GlobalPairwiseLibrary());
            else if (*begIt == 1)
                _appendSegmentMatches(seqSet, pList, msaOpt.sc, matches, scores);
            else if (*begIt == 2)
            {
                Nothing noth;
                appendSegmentMatches(seqSet, pList, msaOpt.sc, matches, scores, noth, AlignConfig<true, true, true, true>(), GlobalPairwiseLibrary());
            }
            else if (*begIt == 3)
                appendSegmentMatches(seqSet, pList, matches, scores, LcsLibrary());
        }
    }

    // Include a T-Coffee library
    if (!empty(msaOpt.libfiles))
    {
        typedef typename Iterator<String<std::string> const, Standard>::Type TIter;
        TIter begIt = begin(msaOpt.libfiles, Standard());
        TIter begItEnd = end(msaOpt.libfiles, Standard());
        for (; begIt != begItEnd; goNext(begIt))
        {
            std::ifstream strm_lib;
            strm_lib.open((*begIt).c_str(), std::ios_base::in | std::ios_base::binary);
            read(strm_lib, matches, scores, sequenceNames, TCoffeeLib());
            strm_lib.close();
        }
    }

    // Include MUMmer segment matches
    if (!empty(msaOpt.mummerfiles))
    {
        typedef typename Iterator<String<std::string> const, Standard>::Type TIter;
        TIter begIt = begin(msaOpt.mummerfiles, Standard());
        TIter begItEnd = end(msaOpt.mummerfiles, Standard());
        for (; begIt != begItEnd; goNext(begIt))
        {
            std::ifstream strm_lib;
            strm_lib.open((*begIt).c_str(), std::ios_base::in | std::ios_base::binary);
            read(strm_lib, matches, scores, seqSet, sequenceNames, MummerLib());
            strm_lib.close();
        }
    }

    // Include BLAST segment matches
    if (!empty(msaOpt.blastfiles))
    {
        typedef typename Iterator<String<std::string> const, Standard>::Type TIter;
        TIter begIt = begin(msaOpt.blastfiles, Standard());
        TIter begItEnd = end(msaOpt.blastfiles, Standard());
        for (; begIt != begItEnd; goNext(begIt))
        {
            std::ifstream strm_lib;
            strm_lib.open((*begIt).c_str(), std::ios_base::in | std::ios_base::binary);
            read(strm_lib, matches, scores, sequenceNames, BlastLib());
            strm_lib.close();
        }
    }

    // Use these segment matches for the initial alignment graph
    TGraph g(seqSet);
    if (!msaOpt.rescore)
        buildAlignmentGraph(matches, scores, g, FractionalScore());
    else
        buildAlignmentGraph(matches, scores, g, msaOpt.sc, ReScore());
    clear(matches);
    clear(scores);

    // Guide tree
    Graph<Tree<TDistanceValue> > guideTree;
    if (!empty(msaOpt.treefile))
    {
        std::fstream strm_tree;
        strm_tree.open(msaOpt.treefile.c_str(), std::ios_base::in | std::ios_base::binary);
        read(strm_tree, guideTree, sequenceNames, NewickFormat());  // Read newick tree
        strm_tree.close();
    }
    else
    {
        // Check if we have a valid distance matrix
        if (empty(distanceMatrix))
            getDistanceMatrix(g, distanceMatrix, KmerDistance());
        // Get distance matrix values for a precision of 10 decimal digits.
        for (unsigned i = 0; i < length(distanceMatrix); ++i)
            distanceMatrix[i] = static_cast<__int64>(distanceMatrix[i] * 1e10) / 1e10;
        if (msaOpt.build == 0)
            njTree(distanceMatrix, guideTree);
        else if (msaOpt.build == 1)
            upgmaTree(distanceMatrix, guideTree, UpgmaMin());
        else if (msaOpt.build == 2)
            upgmaTree(distanceMatrix, guideTree, UpgmaMax());
        else if (msaOpt.build == 3)
            upgmaTree(distanceMatrix, guideTree, UpgmaAvg());
        else if (msaOpt.build == 4)
            upgmaTree(distanceMatrix, guideTree, UpgmaWeightAvg());
    }
    clear(distanceMatrix);

    // Triplet extension
    if (nSeq < threshold)
        tripletLibraryExtension(g);
    else
        tripletLibraryExtension(g, guideTree, threshold / 2);

    // Progressive Alignment
    progressiveAlignment(g, guideTree, gAlign);

    clear(guideTree);
    clear(g);

    //TStringSet& str = stringSet(gAlign);
    //for(TSize i = 0; i<length(str);++i) {
    //    for(TSize j=0;j<length(str[i]);++j) {
    //        std::cout << ordValue(str[i][j]) << ',';
    //    }
    //    std::cout << std::endl;
    //}
}

template <typename TStringSet, typename TCargo, typename TSpec, typename TScore>
void globalMsaAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> > & gAlign,
                        TScore const & scoreObject)
{
    //typedef typename Value<TStringSet>::Type TString;
    //typedef typename Value<TString>::Type TAlphabet;
    TStringSet sequenceSet = stringSet(gAlign);
    String<String<char> > sequenceNames;
    resize(sequenceNames, length(sequenceSet), String<char>("tmpName"));
    MsaOptions<AminoAcid, TScore> msaOpt;
    msaOpt.sc = scoreObject;
    appendValue(msaOpt.method, 0);  // Global pairwise
    appendValue(msaOpt.method, 1);  // Local pairwise
    globalMsaAlignment(gAlign, sequenceSet, sequenceNames, msaOpt);
}

template <typename TSource, typename TSpec, typename TScore>
void globalMsaAlignment(Align<TSource, TSpec> & align,
                        TScore const & scoreObject)
{
    typedef StringSet<TSource, Dependent<> > TStringSet;
    TStringSet sequenceSet = stringSet(align);
    Graph<Alignment<TStringSet, void, WithoutEdgeId> > gAlign(sequenceSet);
    globalMsaAlignment(gAlign, scoreObject);
    convertAlignment(gAlign, align);
}

// TODO(holtgrew): Remove the following?

//////////////////////////////////////////////////////////////////////////////
// Just two testing functions
//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet, typename TMatches>
void
_debugMatches(TStringSet & str,
              TMatches & matches)
{
    typedef typename Id<TStringSet>::Type TId;
    typedef typename Size<TStringSet>::Type TSize;

    // Print all the matches
    std::cout << "The sequences:" << std::endl;
    for (TSize i = 0; i < length(str); ++i)
    {
        std::cout << positionToId(str, i) << ':' << str[i] << std::endl;
    }
    std::cout << "The matches:" << std::endl;
    for (TSize i = 0; i < length(matches); ++i)
    {
        TId tmp_id1 = sequenceId(matches[i], 0);
        std::cout << tmp_id1 << ',' << fragmentBegin(matches[i], tmp_id1) << ',';
        for (TSize j = fragmentBegin(matches[i], tmp_id1); j < fragmentBegin(matches[i], tmp_id1) + fragmentLength(matches[i], tmp_id1); ++j)
        {
            std::cout << str[idToPosition(str, tmp_id1)][j];
        }
        TId tmp_id2 = sequenceId(matches[i], 1);
        std::cout << ',' << tmp_id2 << ',' << fragmentBegin(matches[i], tmp_id2) << ',';
        for (TSize j = fragmentBegin(matches[i], tmp_id2); j < fragmentBegin(matches[i], tmp_id2) + fragmentLength(matches[i], tmp_id2); ++j)
        {
            std::cout << str[idToPosition(str, tmp_id2)][j];
        }
        std::cout << std::endl;

        SEQAN_ASSERT(sequenceId(matches[i], 0) != sequenceId(matches[i], 1));
        SEQAN_ASSERT(fragmentBegin(matches[i], tmp_id1) < length(str[idToPosition(str, tmp_id1)]));
        SEQAN_ASSERT(fragmentBegin(matches[i], tmp_id1) + fragmentLength(matches[i], tmp_id1) <= length(str[idToPosition(str, tmp_id1)]));
        SEQAN_ASSERT(fragmentBegin(matches[i], tmp_id2) < length(str[idToPosition(str, tmp_id2)]));
        SEQAN_ASSERT(fragmentBegin(matches[i], tmp_id2) + fragmentLength(matches[i], tmp_id2) <= length(str[idToPosition(str, tmp_id2)]));
        SEQAN_ASSERT(fragmentLength(matches[i], tmp_id2) == fragmentLength(matches[i], tmp_id1));
    }
}

template <typename TGraph>
void _debugRefinedMatches(TGraph & g)
{
    typedef typename Id<TGraph>::Type TId;
    //typedef typename Size<TGraph>::Type TSize;
    typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;

    std::cout << "Refined matches" << std::endl;
    TEdgeIterator it_tmp(g);
    for (; !atEnd(it_tmp); ++it_tmp)
    {
        TId id1 = sequenceId(g, sourceVertex(it_tmp));
        TId id2 = sequenceId(g, targetVertex(it_tmp));
        std::cout << id1 << ',' << fragmentBegin(g, sourceVertex(it_tmp)) << ',';
        std::cout << label(g, sourceVertex(it_tmp));
        std::cout << ',' << id2 << ',' << fragmentBegin(g, targetVertex(it_tmp)) << ',';
        std::cout << label(g, targetVertex(it_tmp));
        std::cout << " (" << cargo(*it_tmp) << ")";
        std::cout << std::endl;

        SEQAN_ASSERT(sequenceId(g, sourceVertex(it_tmp)) != sequenceId(g, targetVertex(it_tmp)));
        SEQAN_ASSERT(fragmentBegin(g, sourceVertex(it_tmp)) < length((stringSet(g))[idToPosition((stringSet(g)), id1)]));
        SEQAN_ASSERT(fragmentBegin(g, sourceVertex(it_tmp)) + fragmentLength(g, sourceVertex(it_tmp)) <= length((stringSet(g))[idToPosition((stringSet(g)), id1)]));
        SEQAN_ASSERT(fragmentBegin(g, targetVertex(it_tmp)) < length((stringSet(g))[idToPosition((stringSet(g)), id2)]));
        SEQAN_ASSERT(fragmentBegin(g, targetVertex(it_tmp)) + fragmentLength(g, targetVertex(it_tmp)) <= length((stringSet(g))[idToPosition((stringSet(g)), id2)]));
        SEQAN_ASSERT(fragmentLength(g, sourceVertex(it_tmp)) == fragmentLength(g, targetVertex(it_tmp)));

    }
}

}  // namespace seqan

#endif  // #ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_MSA_H
