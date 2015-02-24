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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_STORE_ALL_H
#define SEQAN_HEADER_STORE_ALL_H

//#include <stdio.h>

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Contig Store Configuration
//////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew): Document.

template <typename TSpec = void>
struct FragmentStoreConfig
{
    typedef String<Dna5Q>    TReadSeq;
    typedef String<Dna5Q>    TContigSeq;

    typedef double            TMean;
    typedef double            TStd;
    typedef signed char        TMappingQuality;

    typedef void                    TReadStoreElementSpec;
    typedef Owner<ConcatDirect<> >    TReadSeqStoreSpec;
    typedef void                    TMatePairStoreElementSpec;
    typedef void                    TLibraryStoreElementSpec;
    typedef void                    TContigStoreElementSpec;
    typedef void                    TContigFileSpec;
    typedef void                    TAlignedReadStoreElementSpec;
    typedef Owner<ConcatDirect<> >    TAlignedReadTagStoreSpec;
    typedef void                    TAnnotationStoreElementSpec;

    typedef Alloc<>                    TReadNameSpec;
    typedef Owner<ConcatDirect<> >    TReadNameStoreSpec;
};

//////////////////////////////////////////////////////////////////////////////
// Fragment Store
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class FragmentStore
 * @headerfile <seqan/store.h>
 * @brief Multi-container to store contigs, reads, multiple read alignments and genome annotations.
 *
 * @signature template <[typename TSpec[, typename TConfig]]>
 *            class FragmentStore;
 *
 * @tparam TSpec   The specialializing type.  Default: <tt>void</tt>.
 * @tparam TConfig The configuration struct.  Default: <tt>FragmentStoreConfig&lt;TSpec&gt;</tt>.
 *
 * The FragmentStore is a data structure specifically designed for read mapping, genome assembly or gene
 * annotation. These tasks typically require lots of data structures that are related to each other like: reads,
 * mate-pairs, reference genome; pairwise alignments; genome annotation.
 *
 * The FragmentStore subsumes all these data structures in an easy to use interface.  It represents a multiple alignment
 * of millions of reads or mate-pairs against a reference genome consisting of multiple contigs.  Additionally, regions
 * of the reference genome can be annotated with features like 'gene', 'mRNA', 'exon', 'intro' or custom features.  The
 * FragmentStore supports I/O functions to read/write a read alignment in Sam or Amos format and to read/write
 * annotations in Gff/Gtf format.
 *
 * The FragmentStore can be compared with a database where each table (called "store") is implemented as a String member
 * of the FragmentStore class.  The rows of each table (implemented as structs) are referred by their ids which are
 * their positions in the string and not stored explicitly.  The only exception is the alignedReadStore whose elements
 * of type AlignedReadStoreElement contain an id-member as they may be rearranged in arbitrary order, e.g. by increasing
 * genomic positions or by readId.  Many stores have an associated name store to store element names.  Each name store
 * is a StringSet that stores the element name at the position of its id.  All stores are present in the FragmentStore
 * and empty if unused.  The concrete types, e.g. the position types or read/contig alphabet, can be easily changed by
 * defining a custom config struct which is a template parameter of the FragmentStore class.
 *
 * @section Examples
 *
 * Load read alignments and a reference genome and display the multiple alignment in a genomic range:
 *
 * @include demos/store/store_example.cpp
 *
 * @code{.console}
 * ATTTAAGAAATTACAAAATATAGTTGAAAGCTCTAACAATAGACTAAACCAAGCAGAAGAAAGAGGTTCAGAACTTGAAGACAAGTCTCTTATGAATTAA
 * ATTTAA  AATTACAAAATATAGTTGAAAGCTCTAACAATAGA   AACCAAGCAGAAGAAAGAGGTTCAGAACTTGAAGA  AGTCTCTTATGAATTAA
 * ATTTA GAAATTACAAAATATAGTTGAAAGCTCTAACAATA ACTAAACCAAGCAGAAGAAAGAGGTTCAGAACTTG AGACAAGTCTCTTATGAATTAA
 * attta GAAATTACAAAATATAGTTGAAAGCTCTAACAATAG    AACCAAGCAGAAGAAAGAGGCTCAGAACTTGAAGA  AGTCTCTTATGAATTAA
 * ATTTAA   ATTACAAAATATAGTTGAAAGATCTAACAATAGAC    CCAAGCAGAAGAAAGAGGTTCAGAACTTGAAGACAA     TTATGAATTAA
 * ATTTAAGAA TTACAAAATATAGTTGAAAGCTCTAACAATAGACT     AAGCAGAAGAAAGAGGTTCAGAACTTGAAGACAAG     TATGAATTAA
 * ATTTAAGAAA  ACAAAATATAGTTGAAAGCTCTAACAATAGACTAA     GCAGAAGAAAGAGGTTCAGAACTTGAAGACAAGTC    ATGAATTAA
 * ATTTAAGAAA  ACAAAATATAGTTGAAAGCTCTAACAATAGACTAA      CAGAAGAAAGAGGTTCAGAACTTGAAGACAAGTCT    TGAATTAA
 * ATTTAAGAAA  ACAAAATATAGTTGAAAGCTCTAACAATAGACTAA      CAGAAGAAAGAGGTTCANANNNTGANGACAAGTCT    TGAATTAA
 * ATTTAAGAAATT CAAAATATAGTTGAAAGCTCTAACAATAGACTAAA       GAAGAAAGAGGTTCAGAACTTGAAGACAAGTCTCT   GAATTAA
 * ATTTAAGAAAT   AAAATATAGTTGAAAGCTCTAACAATAGACTAAAC       AAGAAAGAGGTTCAGAACTTGAAGACAAGTCTCGT  GAATTAA
 * ATTTAAGAAAT   AAAATATAGTTGAAAGCTCTAACAATAGACTAAAC       AAGAAAGAGGTTCAGAACTTGAAGACAAGTCTCTT   AATTAA
 * ATTTAAGAAAT    AAATATAGTTGAAAGCTCTAACAATAGACTAAACC        GAAAGAGGTTCAGAACTTGAAGACAAGTCTCTTATG
 * ATTTAAGAAATT   AAATATAGTTGAAAGCTCTAACAATAGACTAAACC          AAGAGGTTCAGAACTTGAAGACAAGTCTCTTATGA
 * ATTTAAGAAATT    AATATAGTTGAAAGCTCTAACAATAGACTAAACCAA        AAGAGGTTCAGAACTTGAAGACAAGTCTCTTATGA
 * ATTTAAGAAATTACA  ATATAGTTGAAAGCTCTAACAATAGACTAAACCAA          GAGGTTCAGAACTTGAAGACAAGTCTCTTATGAAT
 * ATTTAAGAAATTACAA   ATAGTTGAAAGCTCTAACAATAGACTAAACCAAGC        GAGGTTCAGAACTTGAAGACAAGTCTCTTATGAAT
 * ATTTAAGAAATTACAAAATA AGTTGAAAGCTCTAACAATAGACTAAACCAAGCAG       AGGTTCAGAACTTGAAGACAAGTCTCTTATGAATT
 * ATTTAAGAAATTACAAAATAT  TTGAAAGCTCTAACAATAGACTAAACCAAGCAGAA      GGTTCAGAACTTGAAGACAAGTCTCTTATGAATTA
 * ATTTAAGAAATTACAAAATATA   GAAAGCTCTAACAATAGACTAAACCAAGCAGAAGAAAGAG TTCAGAACTTGAAGACAAGTCTCTTATGAATTAA
 * ATTTAAGAAATTACAAAATATAGTTGAA    CTAACAATAGACTAAACCAAGCAGAAGAAAGAGTT      CTTGAAGACAAGTCTCTTATGAATTAA
 * ATTTAAGAAATTACAAAATATAGTTGAAA   CTAACAATAGACTAAACCAAGCAGAAGAAAGAGGTT      TTGAAGACAAGTCTCTTATGAATTAA
 * ATTTAAGAAATTACAAAATATAGTTGAAAG   TAACAATAGACTAAACCAAGCAGAAGAAAGAGGTT       TGAAGACAAGTCTCTTATGAATTAA
 * ATTTAAGAAATTACAAAATATAGTTGAAAGCTCT ACAATAGACTAAACCAAGCAGAAGAAAGAGGTTCA     TGAAGACAAGTCTCTTATGAATTAA
 *   TTAAGAAATTACAAAATATAGTTGAAAGCTCTAAC    GACTAAACCAAGCAGAAGAAAGAGGTTCAGAACTT AAGACAAGTCTCTTATGAATTAA
 *    TAAGAAATTACAAAATATAGTTGAAAGCTCTAACAATAGA                     GGTTCAGAACTTGAAGACAAGTCTCTTATGAATTA
 *           TTACAAAATATAGTTGAAAGCTCTAACAATAGACT                   GGTTCAGAACTTGAAGACAAGTCTCTTATGAATTA
 *                    ATAGTTGAAAGCTCTAACAATAGACTAAACCAAGC           GTTCAGAACTTGAAGACAAGTCTCTTATGAATTAA
 *                           AAAGCTCTAACAATAGACTAAACCAAGCAGAAGAA      TCAGAACTTGAAGACAAGTCTCTTATGAATTAA
 *                           AAAGCTCTAACAATAGACTAAACCAAGCAGAAGAA               NAAGACAAGTCTCTTATGAATTAA
 *                            AAGCTCTAACAATAGACTAAACCAAGCAGAAGAAA              GAAGACAAGTCTCTTATGAATTAA
 *                                  TAACAATAGACTAAACCAAGCAGAAGAAAGAGGTT               AGTCTCTTATGAATTAA
 *                                  TAACAATAGACTAAACCAAGCAGAAGAAAGAGGTT                GTCTCTTATGAATTAA
 *                                   AACAATAGACTAAACCAAGCAGAAGAAAGAGGTTC
 *                                   AACAATAGACTAAACCAAGCAGAAGAAAGAGGTTC
 *                                      AATAGACTAAACCAAGCAGAAGAAAGAGGTTCAGA
 *                                      AATAGACTAAACCAAGCAGAAGAAAGAGGTTCAGA
 * @endcode
 *
 * @section Remarks
 *
 * The following figures visualize the relations between the different stores:
 *
 * <img width="597" height="363" src="FragmentStore.png" title="Stores that are involved in the representation of a multiple read alignment." />
 *
 * <img width="597" height="363" src="AnnotationStore.png" title="Stores that are involved in the representation of a genome alignment." />
 */

/*!
 * @typedef FragmentStore::TReadStore
 * @brief Type of the @link FragmentStore::readStore @endlink member.
 *
 * @signature typedef (..) TFragmentStore::TReadStore;
 *
 *
 * @typedef FragmentStore::TReadSeqStore
 * @brief Type of the @link FragmentStore::readSeqStore @endlink member.
 *
 * @signature typedef (..) TFragmentStore::TReadSeqStore;
 *
 *
 * @typedef FragmentStore::TMatePairStore
 * @brief Type of the @link FragmentStore::matePairStore @endlink member.
 *
 * @signature typedef (..) TFragmentStore::TMatePairStore;
 *
 *
 * @typedef FragmentStore::TContigFileStore
 * @brief Type of the @link FragmentStore::contigFileStore @endlink member.
 *
 * @signature typedef (..) TFragmentStore::TContigFileStore;
 *
 *
 * @typedef FragmentStore::TContigStore
 * @brief Type of the @link FragmentStore::contigStore @endlink member.
 *
 * @signature typedef (..) TFragmentStore::TContigStore;
 *
 *
 * @typedef FragmentStore::TAlignedReadStore
 * @brief Type of the @link FragmentStore::alignedReadStore @endlink member.
 *
 * @signature typedef (..) TFragmentStore::TAlignedReadStore;
 *
 *
 * @typedef FragmentStore::TAnnotationStore
 * @brief Type of the @link FragmentStore::annotationStore @endlink member.
 *
 * @signature typedef (..) TFragmentStore::TAnnotationStore;
 *
 *
 * @typedef FragmentStore::TAlignQualityStore
 * @brief Type of the @link FragmentStore::alignQualityStore @endlink member.
 *
 * @signature typedef (..) TFragmentStore::TAlignQualityStore;
 *
 *
 * @typedef FragmentStore::TAlignedReadTagStore
 * @brief Type of the @link FragmentStore::alignedReadTagStore @endlink member.
 *
 * @signature typedef (..) TFragmentStore::TAlignedReadTagStore;
 *
 *
 * @typedef FragmentStore::TReadNameStore
 * @brief Type of the @link FragmentStore::readNameStore @endlink member.
 *
 * @signature typedef (..) TFragmentStore::TReadNameStore;
 *
 *
 * @typedef FragmentStore::TMatePairNameStore
 * @brief Type of the @link FragmentStore::matePairNameStore @endlink member.
 *
 * @signature typedef (..) TFragmentStore::TMatePairNameStore;
 *
 *
 * @typedef FragmentStore::TLibraryNameStore
 * @brief Type of the @link FragmentStore::libraryNameStore @endlink member.
 *
 * @signature typedef (..) TFragmentStore::TLibraryNameStore;
 *
 *
 * @typedef FragmentStore::TContigNameStore
 * @brief Type of the @link FragmentStore::contigNameStore @endlink member.
 *
 * @signature typedef (..) TFragmentStore::TContigNameStore;
 *
 *
 * @typedef FragmentStore::TAnnotationNameStore
 * @brief Type of the @link FragmentStore::annotationNameStore @endlink member.
 *
 * @signature typedef (..) TFragmentStore::TAnnotationNameStore;
 *
 *
 * @typedef FragmentStore::TAnnotationTypeStore
 * @brief Type of the @link FragmentStore::annotationTypeStore @endlink member.
 *
 * @signature typedef (..) TFragmentStore::TAnnotationTypeStore;
 *
 *
 * @typedef FragmentStore::TAnnotationKeyStore
 * @brief Type of the @link FragmentStore::annotationKeyStore @endlink member.
 *
 * @signature typedef (..) TFragmentStore::TAnnotationKeyStore;
 */

/*!
 * @var FragmentStore::TReadStore FragmentStore::readStore;
 * @brief A @link String @endlink that maps from readId to matePairId.
 *
 * The value type is @link ReadStoreElement @endlink.
 *
 *
 * @var FragmentStore::TReadSeqStore FragmentStore::readSeqStore;
 * @brief @link StringSet @endlink that maps from readId to readSeq.
 *
 *
 * @var FragmentStore::TMatePairStore FragmentStore::matePairStore;
 * @brief @link String @endlink that maps from matePairId to (readId[2], libId).
 *
 * The value type is @link MatePairStoreElement @endlink.
 *
 *
 * @var FragmentStore::TLibraryStore FragmentStore::libraryStore;
 * @brief @link String @endlink that maps from libId to (mean, std).
 *
 * Value type is @link LibraryStoreElement @endlink.
 *
 *
 * @var FragmentStore::TContigFileStore FragmentStore::contigFileStore;
 * @brief @link String @endlink that maps from contigId to (contigSeq, contigGaps, contigFileId).
 *
 * Value type is @link ContigFile @endlink.
 *
 *
 * @var FragmentStore::TContigStore FragmentStore::contigStore;
 * @brief @link String @endlink that maps from contigId to (contigSeq, contigGaps, contigFileId).
 *
 * Value type is @link ContigStoreElement @endlink.
 *
 *
 * @var FragmentStore::TAlignedReadStore FragmentStore::alignedReadStore;
 * @brief @link String @endlink that stores (alignId, readId, contigId, pairMatchId, beginPos, endPos, gapAnchors).
 *
 * The value type is @link AlignedReadStoreElement @endlink.
 *
 * @section Remarks
 *
 * You can sort <tt>alignedReadStore</tt> using @link sortAlignedReads @endlink.  After sorting, you can use the
 * functions @link lowerBoundAlignedReads @endlink and @link upperBoundAlignedReads @endlink to perform a binary search,
 * e.g. for accessing only a subrange.
 *
 *
 * @see lowerBoundAlignedReads
 * @see upperBoundAlignedReads
 * @see sortAlignedReads
 *
 *
 * @var FragmentStore::TAnnotationStore FragmentStore::annotationStore
 * @brief String that maps from annoId to (contigId, typeId, beginPos, endPos, parentId, lastChildId, nextSiblingId, values).
 *
 * The value type is @link AnnotationStoreElement @endlink.
 *
 * Instead of accessing this store directly, consider to use the high-level interface provided by @link
 * AnnotationTreeIterator @endlink.
 *
 * @var FragmentStore::TAlignQualityStore FragmentStore::alignQualityStore;
 * @brief @link String @endlink that maps from alignId to (pairScore, score, errors).
 *
 * The value type is @link AlignQualityStoreElement @endlink.
 *
 *
 * @var FragmentStore::TAlignedReadTagStore FragmentStore::alignedReadTagStore;
 * @brief @link StringSet @endlink that maps from alignId to alignTag.
 *
 *
 * @var FragmentStore::TReadNameStore FragmentStore::readNameStore;
 * @brief @link StringSet @endlink that maps from readId to readName.
 *
 *
 * @var FragmentStore::TContigNameStore FragmentStore::matePairNameStore;
 * @brief @link StringSet @endlink that maps from contigId to matePairName.
 *
 *
 * @var FragmentStore::TLibraryNameStore FragmentStore::libraryNameStore;
 * @brief A @link StringSet @endlink that maps from libId to libName.
 *
 *
 * @var FragmentStore::TContigNameStore FragmentStore::contigNameStore;
 * @brief @link StringSet @endlink that maps from contigId to contigName.
 *
 *
 * @var FragmentStore::TAnnotationNameStore FragmentStore::annotationNameStore;
 * @brief @link StringSet @endlink that maps from annoId to annoName;
 *
 *
 * @var FragmentStore::TAnnotationTypeStore FragmentStore::annotationTypeStore;
 * @brief @link StringSet @endlink that maps from typeId to the type name of an annotation, e.g. "gene" or "exon".
 *        typeId is a member of the @link AnnotationStoreElement @endlink.
 *
 * @section Remarks
 *
 * There are @link FragmentStore::PredefinedAnnotationTypes predefined type ids @endlink for commonly used types
 * e.g. <tt>ANNO_GENE</tt> or <tt>ANNO_EXON</tt> which can be used to set the @link AnnotationStoreElement::typeId
 * @endlink directly as a fast alternative to @link AnnotationTreeIterator#getType @endlink and @link
 * AnnotationTreeIterator#setType @endlink.
 *
 *
 * @var FragmentStore::TAnnotationKeyStore FragmentStore::annotationKeyStore;
 * @brief @link StringSet @endlink that maps from keyId to the name of a key.  The keyId is used to address @link
 *        AnnotationStoreElement::values @endlink of an annotation.
 */

/*!
 * @enum FragmentStore::PredefinedAnnotationTypes
 * @headerfile <seqan/store.h>
 * @brief The @link FragmentStore @endlink predefines some commonly used @link AnnotationStoreElement::typeId @endlink
 *        values.
 *
 * These types can be used to compare or set the @link AnnotationStoreElement::typeId @endlink directly as a
 * fast alternative to @link AnnotationTreeIterator#getType @endlink and @link AnnotationTreeIterator#setType
 * @endlink.
 *
 * @signature enum FragmentStore::PredefinedAnnotationTypes;
 *
 * @val FragmentStore::PredefinedAnnotationTypes FragmentStore::ANNO_ROOT;
 * @brief The root node ("&lt;node&gt;").
 *
 * @val FragmentStore::PredefinedAnnotationTypes FragmentStore::ANNO_GENE;
 * @brief A gene ("gene").
 *
 * @val FragmentStore::PredefinedAnnotationTypes FragmentStore::ANNO_MRNA;
 * @brief An mRNA sequence, aka transcript ("mRNA").
 *
 * @val FragmentStore::PredefinedAnnotationTypes FragmentStore::ANNO_CDS;
 * @brief A coding region ("CDS");
 *
 * @val FragmentStore::PredefinedAnnotationTypes FragmentStore::ANNO_EXON;
 * @brief An exon ("exon").
 *
 * @val FragmentStore::PredefinedAnnotationTypes FragmentStore::ANNO_FIVE_PRIME_UTR;
 * @brief A 5' untranslated region ("five_prime_UTR").
 *
 * @val FragmentStore::PredefinedAnnotationTypes FragmentStore::ANNO_INTRON;
 * @brief An intron ("intron").
 *
 * @val FragmentStore::PredefinedAnnotationTypes FragmentStore::ANNO_THREE_PRIME_UTR;
 * @brief A 3' untranslated region ("three_prime_UTR").
 */


template <typename TSpec = void, typename TConfig = FragmentStoreConfig<TSpec> >
class FragmentStore
{
private:
    typedef typename TConfig::TReadStoreElementSpec            TReadStoreElementSpec;
    typedef typename TConfig::TReadSeqStoreSpec                TReadSeqStoreSpec;
    typedef typename TConfig::TReadNameSpec                    TReadNameSpec;
    typedef typename TConfig::TReadNameStoreSpec            TReadNameStoreSpec;
    typedef typename TConfig::TMatePairStoreElementSpec        TMatePairStoreElementSpec;
    typedef typename TConfig::TLibraryStoreElementSpec        TLibraryStoreElementSpec;
    typedef typename TConfig::TContigStoreElementSpec        TContigStoreElementSpec;
    typedef typename TConfig::TContigFileSpec                TContigFileSpec;
    typedef typename TConfig::TAlignedReadStoreElementSpec    TAlignedReadStoreElementSpec;
    typedef typename TConfig::TAlignedReadTagStoreSpec        TAlignedReadTagStoreSpec;
    typedef typename TConfig::TAnnotationStoreElementSpec    TAnnotationStoreElementSpec;

public:
    typedef typename TConfig::TMean                    TMean;
    typedef typename TConfig::TStd                    TStd;
    typedef typename TConfig::TMappingQuality        TMappingQuality;

    typedef typename TConfig::TReadSeq                TReadSeq;
    typedef typename TConfig::TContigSeq            TContigSeq;

    typedef typename Position<TReadSeq>::Type        TRSeqPos_;
    typedef typename Position<TContigSeq>::Type     TCSeqPos_;
    typedef typename MakeSigned_<TRSeqPos_>::Type    TReadPos;
    typedef typename MakeSigned_<TCSeqPos_>::Type    TContigPos;

    typedef GapAnchor<TReadPos>                        TReadGapAnchor;
    typedef GapAnchor<TContigPos>                    TContigGapAnchor;

    typedef StringSet<CharString>                    TNameStore;

    typedef AnnotationStoreElement< TContigPos, TAnnotationStoreElementSpec >    TAnnotationStoreElement;
    typedef typename TAnnotationStoreElement::TId                                TAnnotationStoreElementId;

    typedef String< ReadStoreElement< TReadStoreElementSpec > >                                                TReadStore;
    typedef String< MatePairStoreElement< TMatePairStoreElementSpec > >                                        TMatePairStore;
    typedef String< LibraryStoreElement< TMean, TStd, TLibraryStoreElementSpec > >                            TLibraryStore;
    typedef String< ContigStoreElement< TContigSeq, TContigGapAnchor, TContigStoreElementSpec > >            TContigStore;
    typedef String< ContigFile< TContigFileSpec > >                                                            TContigFileStore;
    typedef String< AlignedReadStoreElement< TContigPos, TReadGapAnchor, TAlignedReadStoreElementSpec > >    TAlignedReadStore;
    typedef String< AlignQualityStoreElement< TMappingQuality >    >                                            TAlignQualityStore;
    typedef StringSet<CharString, TAlignedReadTagStoreSpec>                                                    TAlignedReadTagStore;
    typedef String< TAnnotationStoreElement >                                                                TAnnotationStore;
    typedef String< IntervalTree< TContigPos, TAnnotationStoreElementId > >                                    TIntervalTreeStore;
    typedef StringSet<TReadSeq, TReadSeqStoreSpec>                                                            TReadSeqStore;

    typedef StringSet< String<char, TReadNameSpec>, TReadNameStoreSpec>                                        TReadNameStore;
    typedef TNameStore                                                                                        TMatePairNameStore;
    typedef TNameStore                                                                                        TLibraryNameStore;
    typedef TNameStore                                                                                        TContigNameStore;
    typedef TNameStore                                                                                        TAnnotationNameStore;
    typedef TNameStore                                                                                        TAnnotationTypeStore;
    typedef TNameStore                                                                                        TAnnotationKeyStore;

    // main containers
    TReadStore            readStore;              // readId       -> matePairId
    TMatePairStore        matePairStore;          // matePairId   -> readId0, readId1, libraryId
    TLibraryStore        libraryStore;           // libraryId    -> libSizeMean, libSizeStd
    TContigStore        contigStore;            // contigId     -> contigSeq, contigGaps, contigFileId
    TContigFileStore    contigFileStore;        // contigFileId -> fileName, firstContigId
    TAlignedReadStore    alignedReadStore;       //              -> id, readId, contigId, pairMatchId (not matePairId!), beginPos, endPos, gaps
    TAnnotationStore    annotationStore;        // annoId       -> parentId, contigId, beginPos, endPos
    TIntervalTreeStore    intervalTreeStore_F;    // treeId (same as contigId)    -> intervalTree (F: forward strand)
    TIntervalTreeStore    intervalTreeStore_R;    //                         (R: reverse complement strand)

                                            // REMARKS:
                                            // 1)
                                            //    beginPos <= endPos     forward strand
                                            //    beginPos >  endPos     backward strand (reverse complement)
                                            // 2)
                                            //    The alignedReadStore can arbitrarily be resorted. The unique identifier id should
                                            //    be used to address additional information for each alignedRead in additional tables.

    // we store the read sequences in a seperate stringset to reduce the memory overhead
    TReadSeqStore        readSeqStore;

    // extra Sam fields
    TAlignQualityStore        alignQualityStore;
    TAlignedReadTagStore    alignedReadTagStore;

    // retrieve the names of reads, mate-pairs, libraries, contigs, annotations by their ids
    TReadNameStore            readNameStore;
    TMatePairNameStore        matePairNameStore;
    TLibraryNameStore        libraryNameStore;
    TContigNameStore        contigNameStore;
    TAnnotationNameStore    annotationNameStore;
    TAnnotationTypeStore    annotationTypeStore;
    TAnnotationKeyStore        annotationKeyStore;

    NameStoreCache<TReadNameStore, CharString>            readNameStoreCache;
    NameStoreCache<TContigNameStore, CharString>        contigNameStoreCache;
    NameStoreCache<TAnnotationNameStore, CharString>    annotationNameStoreCache;
    NameStoreCache<TAnnotationTypeStore, CharString>    annotationTypeStoreCache;
    NameStoreCache<TAnnotationKeyStore, CharString>        annotationKeyStoreCache;

    enum {
        ANNO_ROOT,
        ANNO_DELETED,
        ANNO_GENE,
        ANNO_MRNA,
        ANNO_CDS,
        ANNO_EXON,
        ANNO_FIVE_PRIME_UTR,
        ANNO_INTRON,
        ANNO_THREE_PRIME_UTR,
        ANNO_PREDEFINED
    };

    FragmentStore():
        readNameStoreCache(readNameStore),
        contigNameStoreCache(contigNameStore),
        annotationNameStoreCache(annotationNameStore),
        annotationTypeStoreCache(annotationTypeStore),
        annotationKeyStoreCache(annotationKeyStore)
    {
        // ATTENTION: The order of these keywords must correspond to the order of the enums above.
        appendName(annotationTypeStore, "<root>", annotationTypeStoreCache);
        appendName(annotationTypeStore, "<deleted>", annotationTypeStoreCache);
        appendName(annotationTypeStore, "gene", annotationTypeStoreCache);
        appendName(annotationTypeStore, "mRNA", annotationTypeStoreCache);
        appendName(annotationTypeStore, "CDS", annotationTypeStoreCache);
        appendName(annotationTypeStore, "exon", annotationTypeStoreCache);
        appendName(annotationTypeStore, "five_prime_UTR", annotationTypeStoreCache);
        appendName(annotationTypeStore, "intron", annotationTypeStoreCache);
        appendName(annotationTypeStore, "three_prime_UTR", annotationTypeStoreCache);
        _storeClearAnnotations(*this);
    }

    void swap(FragmentStore & other)
    {
        seqan::swap(alignQualityStore, other.alignQualityStore);
        seqan::swap(alignedReadStore, other.alignedReadStore);
        seqan::swap(alignedReadTagStore, other.alignedReadTagStore);
        seqan::swap(annotationKeyStore, other.annotationKeyStore);
        seqan::swap(annotationNameStore, other.annotationNameStore);
        seqan::swap(annotationStore, other.annotationStore);
        seqan::swap(annotationTypeStore, other.annotationTypeStore);
        seqan::swap(contigFileStore, other.contigFileStore);
        seqan::swap(contigNameStore, other.contigNameStore);
        seqan::swap(contigStore, other.contigStore);
        seqan::swap(intervalTreeStore_F, other.intervalTreeStore_F);
        seqan::swap(intervalTreeStore_R, other.intervalTreeStore_R);
        seqan::swap(libraryNameStore, other.libraryNameStore);
        seqan::swap(libraryStore, other.libraryStore);
        seqan::swap(matePairNameStore, other.matePairNameStore);
        seqan::swap(matePairStore, other.matePairStore);
        seqan::swap(readNameStore, other.readNameStore);
        seqan::swap(readSeqStore, other.readSeqStore);
        seqan::swap(readStore, other.readStore);

        refresh(annotationKeyStoreCache);
        refresh(annotationNameStoreCache);
        refresh(annotationTypeStoreCache);
        refresh(contigNameStoreCache);
        refresh(readNameStoreCache);

        refresh(other.annotationKeyStoreCache);
        refresh(other.annotationNameStoreCache);
        refresh(other.annotationTypeStoreCache);
        refresh(other.contigNameStoreCache);
        refresh(other.readNameStoreCache);
    }
};

template <typename TSpec, typename TConfig>
inline void swap(FragmentStore<TSpec, TConfig> & lhs, FragmentStore<TSpec, TConfig> & rhs)
{
    lhs.swap(rhs);
}

//////////////////////////////////////////////////////////////////////////////

template < typename TSpec, typename TConfig >
struct VertexDescriptor< FragmentStore<TSpec, TConfig> >
{
    typedef FragmentStore<TSpec, TConfig>        TFragmentStore_;
    typedef typename Id<TFragmentStore_>::Type    Type;
};

//////////////////////////////////////////////////////////////////////////////
template < typename TConfig, typename TSpec, typename TIterSpec >
inline typename Iterator<FragmentStore<TSpec, TConfig>, TIterSpec >::Type
begin(FragmentStore<TSpec, TConfig> &store, TIterSpec const)
{
    return Iter<FragmentStore<TSpec, TConfig>, TIterSpec>(store);
}

template < typename TConfig, typename TSpec, typename TIterSpec >
inline typename Iterator<FragmentStore<TSpec, TConfig> const, TIterSpec >::Type
begin(FragmentStore<TSpec, TConfig> const &store, TIterSpec const)
{
    return Iter<FragmentStore<TSpec, TConfig> const, TIterSpec>(store);
}

//////////////////////////////////////////////////////////////////////////////
template < typename TConfig, typename TSpec, typename TIterSpec >
inline typename Iterator<FragmentStore<TSpec, TConfig>, TIterSpec >::Type
end(FragmentStore<TSpec, TConfig> &store, TIterSpec const)
{
    return Iter<FragmentStore<TSpec, TConfig>, TIterSpec>(store, MinimalCtor());
}

template < typename TConfig, typename TSpec, typename TIterSpec >
inline typename Iterator<FragmentStore<TSpec, TConfig> const, TIterSpec >::Type
end(FragmentStore<TSpec, TConfig> const &store, TIterSpec const)
{
    return Iter<FragmentStore<TSpec, TConfig> const, TIterSpec>(store, MinimalCtor());
}

//////////////////////////////////////////////////////////////////////////////

template < typename TConfig, typename TSpec >
inline void
_storeClearAnnotations(FragmentStore<TSpec, TConfig> & me)
{
    typedef FragmentStore<TSpec, TConfig>                TFragmentStore;
    typedef typename TFragmentStore::TAnnotationStore    TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type        TAnnotation;

    resize(me.annotationStore, 1);
    resize(me.annotationNameStore, 1);
    resize(me.annotationTypeStore, (unsigned)TFragmentStore::ANNO_PREDEFINED);

    TAnnotation root;
    root.typeId = 0;
    me.annotationStore[0] = root;
    me.annotationNameStore[0] = "<root>";
}

template < typename TConfig, typename TSpec >
inline void
_storeRemoveTempAnnoNames(FragmentStore<TSpec, TConfig> & me)
{
    typedef FragmentStore<TSpec, TConfig>                        TFragmentStore;
    typedef typename TFragmentStore::TAnnotationStore            TAnnotationStore;
    typedef typename TFragmentStore::TAnnotationNameStore        TAnnotationNameStore;
    typedef typename TFragmentStore::TAnnotationTypeStore        TAnnotationTypeStore;

    typedef typename GetValue<TAnnotationNameStore>::Type        TName;
    typedef typename GetValue<TAnnotationTypeStore>::Type        TType;
    typedef typename Iterator<TAnnotationStore, Standard>::Type    TAnnoIter;
    typedef typename Position<TAnnotationStore>::Type            TPosition;

    TAnnoIter it = end(me.annotationStore, Standard()) - 1;
    TAnnoIter itBegin = begin(me.annotationStore, Standard());
    TPosition pos = it - itBegin;
    for (; itBegin <= it; --it, --pos)
    {
        if ((*it).typeId >= length(me.annotationTypeStore))
            continue;
        TName name = me.annotationNameStore[pos];
        TType type = me.annotationTypeStore[(*it).typeId];
        if (length(name) > length(type) + 3)
            if (prefix(name, 2) == "__" && infix(name, 2, 2 + length(type)) == type && name[2 + length(type)] == '_')
                assignValue(me.annotationNameStore, pos, "");
    }
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TConfig, typename TId>
inline typename GetValue<typename FragmentStore<TSpec, TConfig>::TAnnotationNameStore>::Type
getAnnoName(FragmentStore<TSpec, TConfig> const & store, TId id)
{
    return store.annotationNameStore[id];
}

template <typename TSpec, typename TConfig, typename TId>
inline typename GetValue<typename FragmentStore<TSpec, TConfig>::TAnnotationTypeStore>::Type
getAnnoType(FragmentStore<TSpec, TConfig> const & store, TId id)
{
    return store.annotationTypeStore[id];
}

template <typename TSpec, typename TConfig, typename TId>
inline CharString
getAnnoUniqueName(FragmentStore<TSpec, TConfig> const & store, TId id)
{
    if (id < length(store.annotationNameStore) && !empty(getAnnoName(store, id)))
        return getAnnoName(store, id);

    std::stringstream tmp;
    tmp << "__" << getAnnoType(store, store.annotationStore[id].typeId) << '_' << id;
    return tmp.str();
}


//////////////////////////////////////////////////////////////////////////////
// append functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// _storeAppendRead
//
// adds a new entry to the read store if neccessary. Otherwise it writes the
// correct Id in the variable using to qname to identify it
// If needed a mate pair entry is created.
// Returns true, if the read hasn't been in the read store before and was appended.

template <typename TSpec, typename TConfig, typename TId, typename TName, typename TString, typename TFlag, typename TContext>
inline bool
_storeAppendRead (
    FragmentStore<TSpec, TConfig> & fragStore,
    TId & readId,
    TName const & qname,
    TString const & readSeq,
    TFlag const flag,
    TContext &)
{
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;
    typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;

    // search for readId by name (could be me or my mate)
    bool found = getIdByName(fragStore.readNameStore, qname, readId, fragStore.readNameStoreCache);

    // if naming scheme is xx/1, xx/2 or xx/L, xx/R try to look up my mate
    if (!found && (flag & 1) == 1 && length(qname) >= 2 && qname[length(qname) - 2] == '/')
    {
        CharString mate;

        char tag = back(qname);
        if (tag == '1' || tag == '2')
        {
            mate = qname;
            back(mate) = (tag == '1')? '2': '1';
        }
        else if (tag == 'L' || tag == 'R')
        {
            mate = qname;
            back(mate) = (tag == 'L')? 'R': 'L';
        }
        found = getIdByName(fragStore.readNameStore, mate, readId, fragStore.readNameStoreCache);
    }

    if (found)
    {
        if ((flag & 1) == 1)
        {
            // if the read is in the store and paired
            // check the mate pair store if it is the same mate of the pair
            // assuming that only one flag 0x040 or 0x0080 is 1
            int inPair = 1 - ((flag & 0x40) >> 6);    // bit 7 is set => inPair = 0
                                                    // else inPair = 1 (even if bits 6 and 7 are not set)

            TId matePairId = fragStore.readStore[readId].matePairId;
            if (matePairId != TMatePairElement::INVALID_ID)
            {
                readId = fragStore.matePairStore[matePairId].readId[inPair];
                if (readId == TMatePairElement::INVALID_ID)
                {
                    // create new entry in read and read name store
                    // set sequence and mate pair ID in new read store element
                    readId = appendRead(fragStore, readSeq, matePairId);
                    // add the identifier to the read name store
                    appendName(fragStore.readNameStore, qname, fragStore.readNameStoreCache);
                    // set the ID in the mate pair store
                    fragStore.matePairStore[matePairId].readId[inPair] = readId;
                    return true;
                }
                // else == I am already in the store
            }
            // else == my mate said he has no mate (do nothing)
        }
        return false;
    }

    // if the read name is not in the store
    // create new entry in read and read name store
    readId = length(fragStore.readStore);

    // if the read is paired
    if ((flag & 1) == 1)
    {
        TMatePairElement mateElem;
        // set the first or second read ID in the mate pair element
        TId matePairId = length(fragStore.matePairStore);
        mateElem.readId[(flag & 0x80) >> 7] = readId;
        // get a new mate pair ID and add the new mate pair element
        appendValue(fragStore.matePairStore, mateElem);
        // set the new mate pair ID in the read element
        appendRead(fragStore, readSeq, matePairId);
    }
    // if read is not paired
    else
        appendRead(fragStore, readSeq);

    appendName(fragStore.readNameStore, qname, fragStore.readNameStoreCache);
    return true;
}

//////////////////////////////////////////////////////////////////////////////
// _storeAppendContig
//
// adds a new entry to the contig store if neccessary. Otherwise it writes the
// correct Id in the variable using to qname to identify it

template <typename TSpec, typename TConfig, typename TId, typename TName>
inline void
_storeAppendContig (
    FragmentStore<TSpec, TConfig> & fragStore,
    TId & contigId,
    TName & rName)
{
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;
    typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;

    if (!getIdByName(fragStore.contigNameStore, rName, contigId, fragStore.contigNameStoreCache))
    {
        // if the contig is not in the store yet
        // set the ID on the last entry after appending
        contigId = length(fragStore.contigStore);
        // append contig store
        appendName(fragStore.contigNameStore, rName, fragStore.contigNameStoreCache);
        appendValue(fragStore.contigStore, TContigElement());
//        std::cout << "added contig:" << rName << std::endl;
    }
}

//////////////////////////////////////////////////////////////////////////////
// _annotationAppendAnnotation
//
// adds a new entry to the read store if neccessary. Otherwise it writes the
// correct Id in the variable using to qname to identify it

template <typename TSpec, typename TConfig, typename TId, typename TName, typename TTypeId>
inline void
_storeAppendAnnotationName (
    FragmentStore<TSpec, TConfig> & fragStore,
    TId & annotationId,
    TName & annotationName,
    TTypeId typeId)
{
    SEQAN_ASSERT_EQ(length(fragStore.annotationStore), length(fragStore.annotationNameStore));
    if (!empty(annotationName) && getIdByName(fragStore.annotationNameStore, annotationName, annotationId, fragStore.annotationNameStoreCache))
    {
        do
        {
            // allow different annotations to have the same name (but different typeId)
            if (typeId == maxValue<TTypeId>() || fragStore.annotationStore[annotationId].typeId == typeId)
                return;
            ++annotationId;
        } while (annotationId < length(fragStore.annotationNameStore) && fragStore.annotationNameStore[annotationId] == annotationName);
    }
    // if the annotation is not in the store yet
    // set the ID on the last entry after appending
    annotationId = length(fragStore.annotationNameStore);
    // append to annotationName store
    appendName(fragStore.annotationNameStore, annotationName, fragStore.annotationNameStoreCache);
    // we also need to append an annotation to store the typeId in case of duplicate annotation names
    resize(fragStore.annotationStore, length(fragStore.annotationStore) + 1);
    back(fragStore.annotationStore).typeId = typeId;
}

template <typename TSpec, typename TConfig, typename TId, typename TName>
inline void
_storeAppendAnnotationName (
    FragmentStore<TSpec, TConfig> & fragStore,
    TId & annotationId,
    TName & annotationName)
{
    _storeAppendAnnotationName(fragStore, annotationId, annotationName, maxValue<TId>());
}

template <typename TSpec, typename TConfig, typename TId, typename TName>
inline void
_storeAppendType (
    FragmentStore<TSpec, TConfig> & fragStore,
    TId & typeId,
    TName & annotationType)
{
    if (!getIdByName(fragStore.annotationTypeStore, annotationType, typeId, fragStore.annotationTypeStoreCache))
    {
        // if the annotation type name is not in the store yet
        // set the ID on the last entry after appending
        typeId = length(fragStore.annotationTypeStore);
        // append to annotationType store
        if (!empty(annotationType))
            appendName(fragStore.annotationTypeStore, annotationType, fragStore.annotationTypeStoreCache);
//        std::cout << "added type:" << annotationType << std::endl;
    }
}

template <typename TSpec, typename TConfig, typename TId, typename TName>
inline void
_storeAppendKey (
    FragmentStore<TSpec, TConfig> & fragStore,
    TId & keyId,
    TName & annotationKey)
{
    if (!getIdByName(fragStore.annotationKeyStore, annotationKey, keyId, fragStore.annotationKeyStoreCache))
    {
        // if the key name is not in the store yet
        // set the ID on the last entry after appending
        keyId = length(fragStore.annotationKeyStore);
        // append to annotationKey store
        if (!empty(annotationKey))
            appendName(fragStore.annotationKeyStore, annotationKey, fragStore.annotationKeyStoreCache);
//        std::cout << "added key:" << annotationKey << std::endl;
    }
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TConfig, typename TAnnotation, typename TKey, typename TValue>
inline void
annotationAssignValueByKey (
    FragmentStore<TSpec, TConfig> & fragStore,
    TAnnotation & annotation,
    TKey const & key,
    TValue const & value)
{
    typedef typename TAnnotation::TValues    TValues;
    typedef typename Size<TValues>::Type    TKeyId;

    TKeyId keyId = 0;
    _storeAppendKey(fragStore, keyId, key);
    if (length(annotation.values) <= keyId)
        resize(annotation.values, keyId + 1);
    assignValue(annotation.values, keyId, value);
}

template <typename TSpec, typename TConfig, typename TAnnotation, typename TKey, typename TValue>
inline bool
annotationGetValueByKey (
    TValue & value,
    FragmentStore<TSpec, TConfig> & fragStore,
    TAnnotation const & annotation,
    TKey const & key)
{
    typedef typename TAnnotation::TValues    TValues;
    typedef typename Size<TValues>::Type    TKeyId;

    TKeyId keyId = 0;
    if (!getIdByName(fragStore.annotationKeyStore, key, keyId, fragStore.annotationKeyStoreCache))
        return false;

    if (keyId >= length(annotation.values))
        return false;

    if (empty(annotation.values[keyId]))
        return false;

    assign(value, annotation.values[keyId]);
    return true;
}

template <typename TSpec, typename TConfig, typename TAnnotation, typename TKey>
inline typename Id<TAnnotation>::Type
annotationGetValueIdByKey (
    FragmentStore<TSpec, TConfig> & fragStore,
    TAnnotation const & annotation,
    TKey const & key)
{
    typedef typename TAnnotation::TValues           TValues;
    typedef typename Size<TValues>::Type            TKeyId;

    static const typename Value<typename TAnnotation::TValues const>::Type emptyString;

    TKeyId keyId = 0;
    if (getIdByName(fragStore.annotationKeyStore, key, keyId, fragStore.annotationKeyStoreCache))
        if (keyId < length(annotation.values) && !empty(annotation.values[keyId]))
            return keyId;
    return TAnnotation::INVALID_ID;
}

/*!
 * @fn FragmentStore#clearContigs
 * @brief REvmoes all contigs from a FragmentStore.
 *
 * @signature void clearContigs(store);
 *
 * @param[in,out] store The FragmentStore to remove all contigs from.
 *
 * This function clears the @link FragmentStore::contigStore @endlink and @link FragmentStore::contigNameStore @endlink.
 */

template <typename TSpec, typename TConfig>
inline void
clearContigs(FragmentStore<TSpec, TConfig> &me)
{
    clear(me.contigStore);
    clear(me.contigNameStore);
    refresh(me.contigNameStoreCache);
}


//////////////////////////////////////////////////////////////////////////////
// Read Store Accessors
//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn FragmentStore::clearReads
 * @brief Removes all reds from a FragmentStore.
 *
 * @signature void clearReads(store);
 *
 * @param[in,out] store The FragmentStore to remove all reads from.
 *
 * Clears the @link FragmentStore::readStore @endlink, @link FragmentStore::readSeqStore @endlink, and @link
 * FragmentStore::readNameStore @endlink.
 */

template <typename TSpec, typename TConfig>
inline void
clearReads(FragmentStore<TSpec, TConfig> &me)
{
    clear(me.readStore);
    clear(me.readSeqStore);
    clear(me.readNameStore);
    refresh(me.readNameStoreCache);
}

/*!
 * @fn FragmentStore#appendRead
 * @brief Append a read to a FragmentStore.
 *
 * @signature TSize appendRead(store, read[, matePairId]);
 * @signature TSize appendRead(store, read, name[, matePairId]);
 *
 * @param[in,out] store      The FragmentStore to append the read to.
 * @param[in]     read       The read sequence.  Type: @link ContainerConcept @endlink.
 * @param[in]     name       The name of the read.  Type: @link CharString @endlink.
 * @param[in]     matePairId ID of the mate-pair that this read is part of.  Default:
 *                           <tt>FragmentStore::INVALID_ID</tt> which corresponds to an unmated read.
 *
 * @return TSize The readId of the newly appended read.  TSize is the size type of the @link FragmentStore::readStore
 *               @endlink.
 *
 * This funciton appends a single read to the @link FragmentStore::readStore @endlink and @link
 * FragmentStore::readSeqStore @endlink.
 *
 * @see FragmentStore#getRead
 */


template <typename TSpec, typename TConfig, typename TRead, typename TId>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TReadStore>::Type
appendRead(
    FragmentStore<TSpec, TConfig> &me,
    TRead const &read,
    TId matePairId)
{
    SEQAN_ASSERT_EQ(length(me.readStore), length(me.readSeqStore));

    typedef typename FragmentStore<TSpec, TConfig>::TReadStore TReadStore;
    typename Value<TReadStore>::Type r;
    r.matePairId = matePairId;

    appendValue(me.readStore, r, Generous());
    appendValue(me.readSeqStore, read, Generous());
    return length(me.readStore) - 1;
}

template <typename TSpec, typename TConfig, typename TRead, typename TId>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TReadStore>::Type
appendRead(
    FragmentStore<TSpec, TConfig> &me,
    TRead const &read,
    CharString const &name,
    TId matePairId)
{
    SEQAN_ASSERT_EQ(length(me.readStore), length(me.readSeqStore));

    typedef typename FragmentStore<TSpec, TConfig>::TReadStore TReadStore;
    typename Value<TReadStore>::Type r;
    r.matePairId = matePairId;

    appendValue(me.readStore, r, Generous());
    appendValue(me.readSeqStore, read, Generous());
    appendValue(me.readNameStore, name, Generous());
    return length(me.readStore) - 1;
}

template <typename TSpec, typename TConfig, typename TRead>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TReadStore>::Type
appendRead(
    FragmentStore<TSpec, TConfig> &me,
    TRead const &read)
{
    typedef typename FragmentStore<TSpec, TConfig>::TReadStore TReadStore;
    typedef typename Value<TReadStore>::Type TReadStoreElement;

    return appendRead(me, read, TReadStoreElement::INVALID_ID);
}

template <typename TSpec, typename TConfig, typename TRead>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TReadStore>::Type
appendRead(
    FragmentStore<TSpec, TConfig> &me,
    TRead const &read,
    CharString const &name)
{
    typedef typename FragmentStore<TSpec, TConfig>::TReadStore TReadStore;
    typedef typename Value<TReadStore>::Type TReadStoreElement;

    return appendRead(me, read, name, TReadStoreElement::INVALID_ID);
}

/*!
 * @fn FragmentStore#getRead
 * @brief Returns the read with the given readId.
 *
 * @signature TRead getRead(store, id);
 *
 * @param[in] store The FragmentStore to query for the read.
 * @param[in] id    The id of the read.
 *
 * @return TRead The entry from the @link FragmentStore::readStore @endlink.  TRead is the value type of the @link
 *               FragmentStore::readStore @endlink.
 */

template <typename TSpec, typename TConfig, typename TId>
inline typename Value<typename FragmentStore<TSpec, TConfig>::TReadSeqStore>::Type
getRead(
    FragmentStore<TSpec, TConfig> &me,
    TId id)
{
    return value(me.readSeqStore, id);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn FragmentStore#appendAlignedRead
 * @brief Appends an aligned read entyr to a fragment store.
 *
 * @signature TSize appendAlignedRead(store, readId, contigId, beginPos, endPos[, pairMatchId]);
 *
 * @param[in,out] store       The FragmentStore to append to.
 * @param[in]     readId      The id of the read to append an alignment for.
 * @param[in]     contigId    The id of the contig of the alignment.
 * @param[in]     beginPos    The begin position of the alignment.
 * @param[in]     endPos      The end position of the alignment.
 * @param[in]     pairMatchId The id of the alignedRead pair.  Default:
 *                            <tt>FragmentStore::INVALID_ID</tt> which corresponds to an unmated read.
 *
 * @return TSize The alignedReadId of the alignment.
 *
 * @section Remarks
 *
 * This function appends a single read alignment to the @link FragmentStore::alignedReadStore @endlink.  Note that this
 * really only adds a match.  To generate a global alignment out of all these matches, use @link
 * FragmentStore#convertMatchesToGlobalAlignment @endlink.
 *
 * @see FragmentStore#appendRead
 */

template <typename TSpec, typename TConfig, typename TReadId, typename TContigId, typename TPos, typename TPairMatchId>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TAlignedReadStore>::Type
appendAlignedRead(
        FragmentStore<TSpec, TConfig> & store,
        TReadId const & readId,
        TContigId const & contigId,
        TPos const & beginPos,
        TPos const & endPos,
        TPairMatchId const & pairMatchId)
{
    SEQAN_CHECKPOINT;
    typedef typename FragmentStore<TSpec, TConfig>::TAlignedReadStore TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type TAlignedReadStoreElement;

    TAlignedReadStoreElement element(length(store.alignedReadStore), readId, contigId, beginPos, endPos);
    element.pairMatchId = pairMatchId;
    appendValue(store.alignedReadStore, element);

    return back(store.alignedReadStore).id;
}

template <typename TSpec, typename TConfig, typename TReadId, typename TContigId, typename TPos>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TAlignedReadStore>::Type
appendAlignedRead(
        FragmentStore<TSpec, TConfig> & store,
        TReadId const & readId,
        TContigId const & contigId,
        TPos const & beginPos,
        TPos const & endPos)
{
    SEQAN_CHECKPOINT;
    typedef typename FragmentStore<TSpec, TConfig>::TAlignedReadStore TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type TAlignedReadStoreElement;
    return appendAlignedRead(store, readId, contigId, beginPos, endPos, TAlignedReadStoreElement::INVALID_ID);
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn FragmentStore#appendMatePair
 * @brief Appends the two reads of a mate pair to a FragmentStore.
 *
 * @signature TSize appendMatePair(store, readSeq1, readSeq2[, name1, name2]);
 *
 * @param[in,out] store    The FragmentStore to append the mate pair to.
 * @param[in]     readSeq1 The read sequence of the first read.
 * @param[in]     readSeq2 The read sequence of the second read.
 * @param[in]     name1    The read name of the first read.
 * @param[in]     name2    The read name of the first read.
 *
 * @return TSize The matePairId of the newly appended mate pair.  TSize is the size type of the @link
 *               FragmentStore::matePairStore @endlink.
 *
 * @section Remarks
 *
 * This function appends two reads to the @link FragmentStore::readStore @endlink and @link FragmentStore::readSeqStore
 * @endlink and a mate pair entry for both them to the @link FragmentStore::matePairStore @endlink.  If names are given,
 * they are appended to @link FragmentStore::readNameStore @endlink.
 */

template <typename TSpec, typename TConfig, typename TRead>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TMatePairStore>::Type
appendMatePair(
    FragmentStore<TSpec, TConfig> &me,
    TRead const &read1,
    TRead const &read2)
{
    typedef FragmentStore<TSpec, TConfig>            TFragmentStore;
    typedef typename TFragmentStore::TReadStore        TReadStore;
    typedef typename TFragmentStore::TMatePairStore    TMatePairStore;

    typename Value<TReadStore>::Type r;
    typename Value<TMatePairStore>::Type mp;
    r.matePairId = length(me.matePairStore);
    mp.readId[0] = length(me.readStore);
    mp.readId[1] = length(me.readStore) + 1;

    appendValue(me.readStore, r, Generous());
    appendValue(me.readStore, r, Generous());
    appendValue(me.matePairStore, mp, Generous());
    appendValue(me.readSeqStore, read1, Generous());
    appendValue(me.readSeqStore, read2, Generous());
    return length(me.matePairStore) - 1;
}

template <typename TSpec, typename TConfig, typename TRead>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TMatePairStore>::Type
appendMatePair(
    FragmentStore<TSpec, TConfig> &me,
    TRead const &read1,
    TRead const &read2,
    CharString const &name1,
    CharString const &name2)
{
    SEQAN_ASSERT_EQ(length(me.readStore), length(me.readSeqStore));

    typedef FragmentStore<TSpec, TConfig>            TFragmentStore;
    typedef typename TFragmentStore::TReadStore        TReadStore;
    typedef typename TFragmentStore::TMatePairStore    TMatePairStore;

    typename Value<TReadStore>::Type r;
    typename Value<TMatePairStore>::Type mp;
    r.matePairId = length(me.matePairStore);
    mp.readId[0] = length(me.readStore);
    mp.readId[1] = length(me.readStore) + 1;

    appendValue(me.readStore, r, Generous());
    appendValue(me.readStore, r, Generous());
    appendValue(me.matePairStore, mp, Generous());
    appendValue(me.readSeqStore, read1, Generous());
    appendValue(me.readSeqStore, read2, Generous());
    appendValue(me.readNameStore, name1, Generous());
    appendValue(me.readNameStore, name2, Generous());
    return length(me.matePairStore) - 1;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn FragmentStore#compactAlignedReads
 * @brief Remove invalid aligned reads and rename the alignId's sequentially beginning with 0.
 *
 * @signature TSize compactAlignedReads(store);
 *
 * @param[in,out] store The FragmentStore to compact the aligned reads of.
 *
 * @returns TSize The new size of the @link FragmentStore::alignedReadStore @endlink.  TSize is the size type of the
 *                @link FragmentStore::alignedReadStore @endlink.
 *
 * @section Remarks
 *
 * This function removes all entries from the @link FragmentStore::alignedReadStore @endlink whose alignId is equal to
 * <tt>INVALID_ID</tt> as well as orphan entries in @link FragmentStore::alignQualityStore @endlink.  Afterwards, the
 * alignIds are renamed sequentially beginning with 0.  This function can be used to remove alignments which are flagged
 * by previously setting their id to INVALID_ID.
 */

// 1. remove aligned reads with invalid ids
// 2. rename ids beginning with 0
template <typename TSpec, typename TConfig>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TAlignedReadStore>::Type
compactAlignedReads(FragmentStore<TSpec, TConfig> &me)
{
    typedef FragmentStore<TSpec, TConfig>                            TFragmentStore;
    typedef typename TFragmentStore::TAlignedReadStore                TAlignedReadStore;
    typedef typename TFragmentStore::TAlignQualityStore                TAlignQualityStore;

    typedef typename Value<TAlignedReadStore>::Type                    TAlignedRead;
    typedef typename Id<TAlignedRead>::Type                            TId;
    typedef typename Size<TAlignQualityStore>::Type                    TAQSize;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type    TAlignedReadIter;
    typedef typename Iterator<TAlignQualityStore, Standard>::Type    TAlignQualityIter;

    sortAlignedReads(me.alignedReadStore, SortId());

    TAlignedReadIter itAR = begin(me.alignedReadStore, Standard());
    TAlignedReadIter itARend = end(me.alignedReadStore, Standard());
    TAlignQualityIter itAQ = begin(me.alignQualityStore, Standard());
    TAlignQualityIter itAQbegin = itAQ;
    TAQSize aqSize = length(me.alignQualityStore);
    TId newId = 0;

    for (; itAR != itARend; ++itAR, ++newId)
    {
        TId id = (*itAR).id;
        if (id == TAlignedRead::INVALID_ID) break;    // we assume that invalid ids are at the end of the AlignedReadStore
        if (id < aqSize)
        {
            *itAQ = *(itAQbegin + id);
            ++itAQ;
        }
        (*itAR).id = newId;
    }

    resize(me.alignedReadStore, newId, Exact());
    resize(me.alignQualityStore, itAQ - itAQbegin, Exact());
    return newId;
}

/*!
 * @fn FragmentStore#compactPairMatchIds
 * @brief Renames pairMatchId sequentially beginning with 0.
 *
 * @signature TSize compactPairMatchIds(store);
 *
 * @param[in,out] store The FragmentStore to compact pair match ids of.
 *
 * @return TSize The number of pair matches.  TSize is the size type of @link FragmentStore::alignedReadStore @endlink.
 *
 * @section Remarks
 *
 * This function renames the pairMatchId in the @link FragmentStore::alignedReadStore @endlink sequentially beginning
 * with 0.  Two read alignments can be identified to be pair match if they have the same pairMatchId.  Please note that
 * paired reads not necessarily have to map as a pair match, e.g. if they are on different ocntigs or have the same
 * orientation or a wrong insert size.
 */

// rename pair match ids beginning with 0, returns the number of pair matches
template <typename TSpec, typename TConfig>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TAlignedReadStore>::Type
compactPairMatchIds(FragmentStore<TSpec, TConfig> &me)
{
    typedef FragmentStore<TSpec, TConfig>                            TFragmentStore;
    typedef typename TFragmentStore::TAlignedReadStore                TAlignedReadStore;

    typedef typename Value<TAlignedReadStore>::Type                    TAlignedRead;
    typedef typename Id<TAlignedRead>::Type                            TId;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type    TAlignedReadIter;

    sortAlignedReads(me.alignedReadStore, SortPairMatchId());

    TAlignedReadIter itAR = begin(me.alignedReadStore, Standard());
    TAlignedReadIter itARend = end(me.alignedReadStore, Standard());
    if (itAR == itARend) return 0;

    TId lastId = (*itAR).pairMatchId;
    TId newId = 0;
    for (; itAR != itARend; ++itAR)
    {
        TId id = (*itAR).pairMatchId;
        if (id == TAlignedRead::INVALID_ID) break;    // we assume that invalid ids are at the end of the AlignedReadStore
        if (lastId < id)
        {
            lastId = id;
            ++newId;
        }
        (*itAR).pairMatchId = newId;
    }
    return newId + 1;
}

/*!
 * @fn FragmentStore#calculateInsertSizes
 * @brief Calcualtes a string wtih insert sizes for each pair match.
 *
 * @signature void calculateInsertSizes(insertSizes, store);
 *
 * @param[out] insertSizes A @link String @endlink of insert sizes.  This string is accordingly resized and can be
 *                         addressed by the pairMatchId.
 * @param[in]  store       The @link FragmentStore @endlink to compute the insert sizes for.
 *
 * @section Remarks
 *
 * This function calls @link FragmentStore#compactPairMatchIds @endlink first and calcualte the insert size for every
 * pair match.  The insert size of a pair match is the outer distance between the two matches.
 */

template <typename TLibSizeString, typename TSpec, typename TConfig>
inline void
calculateInsertSizes(TLibSizeString &insertSizes, FragmentStore<TSpec, TConfig> &me)
{
    typedef FragmentStore<TSpec, TConfig>                            TFragmentStore;
    typedef typename TFragmentStore::TAlignedReadStore                TAlignedReadStore;

    typedef typename Value<TAlignedReadStore>::Type                    TAlignedRead;
    typedef typename Id<TAlignedRead>::Type                            TId;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type    TAlignedReadIter;
    typedef typename TFragmentStore::TContigPos                        TGPos;

    TAlignedReadIter it = begin(me.alignedReadStore, Standard());
    TAlignedReadIter itEnd = end(me.alignedReadStore, Standard());

    resize(insertSizes, compactPairMatchIds(me), Exact());
    TId lastId = TAlignedRead::INVALID_ID;
    TGPos leftMatePos = 0;
    for (; it != itEnd; ++it)
    {
        TId id = (*it).pairMatchId;
        if (id == TAlignedRead::INVALID_ID) break;    // we assume that invalid ids are at the end of the AlignedReadStore
        if (id != lastId) {
            leftMatePos = (*it).beginPos;
            lastId = id;
        } else {
            if ((*it).beginPos < leftMatePos)
                insertSizes[id] = leftMatePos - (*it).beginPos;
            else
                insertSizes[id] = (*it).beginPos - leftMatePos;
        }
    }
}

/*!
 * @fn FragmentStore#getMateNo
 * @brief Returns the mate number for a read given a readId.
 *
 * @signature int getMateNo(store, readId);
 *
 * @param[in] store  The FragmentStore with the read.
 * @param[in] readId The readId.
 *
 * @return int The mate number (0 for the first mate, 1 for the second mate) of the read in its mate pair or -1 if the
 *             read is not paired.
 */

template <typename TSpec, typename TConfig, typename TId>
inline int
getMateNo(FragmentStore<TSpec, TConfig> const & me, TId readId)
{
    typedef FragmentStore<TSpec, TConfig>            TFragmentStore;
    typedef typename TFragmentStore::TReadStore        TReadStore;
    typedef typename TFragmentStore::TMatePairStore    TMatePairStore;

    typedef typename Value<TReadStore>::Type        TRead;
    typedef typename Value<TMatePairStore>::Type    TMatePair;

    if (readId != TRead::INVALID_ID)
    {
        TRead const &r = me.readStore[readId];
        if (r.matePairId != TRead::INVALID_ID)
        {
            TMatePair const &mp = me.matePairStore[r.matePairId];
            if (mp.readId[0] == readId) return 0;
            if (mp.readId[1] == readId) return 1;
        }
    }
    return -1;
}

/*!
 * @fn FragmentStore#calculateMateIndices
 * @brief Calculates a string that maps the readId of a read to the index of its mate in the alignedReadStore.
 *
 * @signature void calculateMateIndices(mateIndices, store);
 *
 * @param[out] mateIndices A @link String @endlink with the resulting mate indices.  This string is accordingly resized
 *                         and can be addressed by the readId.
 * @param[in]  store       The @link FragmentStore @endlink.
 *
 * Entries of reads without a mate contain <tt>INVALID_ID</tt>.
 */

// calculate index of the other mate for each pair match
template <typename TMateIndexString, typename TSpec, typename TConfig>
inline void
calculateMateIndices(TMateIndexString & mateIndices, FragmentStore<TSpec, TConfig> const & me)
{
    typedef FragmentStore<TSpec, TConfig> const                     TFragmentStore;
    typedef typename TFragmentStore::TAlignedReadStore const        TAlignedReadStore;

    typedef typename Value<TAlignedReadStore>::Type                 TAlignedRead;
    typedef typename Id<TAlignedRead>::Type                         TId;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type    TAlignedReadIter;

    TAlignedReadIter it = begin(me.alignedReadStore, Standard());
    TAlignedReadIter itEnd = end(me.alignedReadStore, Standard());

    for (TId idx = 0; it != itEnd; ++it, ++idx)
    {
        TId id = (*it).pairMatchId;
        if (id == TAlignedRead::INVALID_ID) continue;
        if (length(mateIndices) < 2*id + 2)
            resize(mateIndices, 2*id + 2, TAlignedRead::INVALID_ID, Generous());
        SEQAN_ASSERT_NEQ(getMateNo(me, (*it).readId), -1);
        mateIndices[2*id + 1 - getMateNo(me, (*it).readId)] = idx;
    }
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @class AlignedReadLayout
 * @headerfile <seqan/store.h>
 * @brief Stores a two dimensional visible layout of a multi-read alignment.
 *
 * @signature struct AlignedReadLayout;
 *
 * @var TContigRows AlignedReadLayout::contigRows;
 * @brief 2D multi-read layout.
 *
 * Stores for a contig and row the ids of aligned reads from left to right.  <tt>contigRows[contigId][row]</tt> stores
 * the alignId of all aligned reads from left to right assigned to the same row.  <tt>row</tt> is the row of the
 * alignment in the multiple sequence alignment and contigId the id of the reference contig.
 */


struct AlignedReadLayout
{
    typedef String<unsigned>    TRow;
    typedef String<TRow>        TRows;
    typedef String<TRows>        TContigRows;

    TContigRows contigRows;            // rows string, each row is a string of ids of alignedReads from left to right
    String<Pair<int> > mateCoords;    // coords of mate pair
};

/*!
 * @fn AlignedReadLayout#layoutAlignment
 * @brief Calculates a visible layout of aligned reads.
 *
 * @signature void layoutAlignment(layout, store);
 *
 * @param[out] layout The resulting layout structure.
 * @param[in]  store  The FragmentStore.
 *
 * For each contig, thisf unction layouts all reads in rows from up to down reusing empty row spaces.
 *
 * @see AlignedReadLayout#printAlignment
 */

template <typename TSpec, typename TConfig>
void layoutAlignment(AlignedReadLayout &layout, FragmentStore<TSpec, TConfig> &store)
{
    typedef FragmentStore<TSpec, TConfig>                            TFragmentStore;

    typedef typename TFragmentStore::TAlignedReadStore                TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type                    TAlignedRead;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type    TAlignedReadIter;

    typedef typename Id<TAlignedRead>::Type                            TId;
    typedef typename TFragmentStore::TContigPos                        TContigPos;

    typedef typename AlignedReadLayout::TRows                        TRows;
    typedef typename AlignedReadLayout::TRow                        TRow;
    typedef typename Iterator<TRows>::Type                            TRowsIter;

    // sort matches by increasing begin positions
    sortAlignedReads(store.alignedReadStore, SortBeginPos());
    sortAlignedReads(store.alignedReadStore, SortContigId());

    clear(layout.contigRows);
    TAlignedReadIter it = begin(store.alignedReadStore, Standard());
    TAlignedReadIter itEnd = end(store.alignedReadStore, Standard());

    for (TId id = 0; it != itEnd; ++it, ++id)
    {
        if ((*it).contigId == TAlignedRead::INVALID_ID) continue;
        if (length(layout.contigRows) <= (*it).contigId)
            resize(layout.contigRows, (*it).contigId + 1);

        TRowsIter lit = begin(layout.contigRows[(*it).contigId], Standard());
        TRowsIter litEnd = end(layout.contigRows[(*it).contigId], Standard());

        TContigPos beginPos = _min((*it).beginPos, (*it).endPos);

        for (; lit != litEnd; ++lit)
        {
            if (empty(*lit)) break;
            TAlignedRead &align = store.alignedReadStore[back(*lit)];
            if (_max(align.beginPos, align.endPos) < beginPos)            // maybe <= would be better
                break;                                                    // but harder to differ between reads borders
        }

        if (lit == litEnd)
        {
            TRow s;
            appendValue(s, id);
            appendValue(layout.contigRows[(*it).contigId], s);
        } else
            appendValue(*lit, id);
    }
}

template <typename TStream, typename TContigGaps, typename TReadGaps, typename TAlignedRead, typename TLine>
inline void _printRead(
    TStream &stream,
    AlignedReadLayout &,
    TContigGaps &contigGaps,
    TReadGaps &readGaps,
    TAlignedRead &alignedRead,
    TLine)
{
//    stream << readGaps;

    typedef typename Iterator<TContigGaps, Standard>::Type TContigIterator;
    typedef typename Iterator<TReadGaps, Standard>::Type TIterator;

    TContigIterator cit = begin(contigGaps, Standard()) + (_min(alignedRead.beginPos, alignedRead.endPos) + beginPosition(readGaps));
    TIterator it = begin(readGaps, Standard());
    TIterator itEnd = end(readGaps, Standard());

    char identityChar = (alignedRead.beginPos < alignedRead.endPos)? '.' : ',';

    for (; it != itEnd; ++it, ++cit)
    {
        if (isGap(it))
            stream << '*';
        else if (convert<Dna5>(*cit) == convert<Dna5>(*it))
            stream << identityChar;
        else
            stream << convert<Dna5>(*it);
    }
}

template <typename TStream, typename TContigGaps, typename TContigName>
inline void _printContig(
    TStream &stream,
    AlignedReadLayout &,
    TContigGaps &contigGaps,
    TContigName const &)
{
    stream << contigGaps;
}

/*!
 * @fn AlignedReadLayout#printAlignment
 * @brief Prints a window of the visible layout of reads into a std::outstream.
 *
 * @signature void printAlignment(stream, layout, store, contigID, posBegin, posEnd, lineBegin, lineEnd);
 *
 * @param[in,out] stream    The std::ostream to print to.
 * @param[in]     layout    The @link AlignedReadLayout @endlink computed earlier in @link
 *                          AlignedReadLayout#layoutAlignment @endlink.
 * @param[in]     store     The FragmentStore that this layout belongs to.
 * @param[in]     contigID  The id of the contig for the alignments to print.
 * @param[in]     posBegin  The begin position of the window.
 * @param[in]     posEnd    The end position of the window.
 * @param[in]     lineBegin The first line of the window.
 * @param[in]     lineEnd   The end line of the window.
 *
 * The window coordinates (beginPos, ...) may be chosen bigger than the lyaout is.  The empty space is then filled with
 * spaces.
 */

template <typename TStream, typename TSpec, typename TConfig, typename TContigId, typename TPos, typename TNum>
void printAlignment(
    TStream &stream,
    AlignedReadLayout &layout,
    FragmentStore<TSpec, TConfig> &store,
    TContigId contigId,
    TPos posBegin, TPos posEnd,
    TNum lineBegin, TNum lineEnd)
{
    typedef FragmentStore<TSpec, TConfig>                            TFragmentStore;

    typedef typename TFragmentStore::TAlignedReadStore                TAlignedReadStore;
    typedef typename TFragmentStore::TContigStore                    TContigStore;

    typedef typename Value<TContigStore>::Type                        TContig;
    typedef typename Value<TAlignedReadStore>::Type                    TAlignedRead;
    typedef typename TFragmentStore::TReadSeq                        TReadSeq;
    typedef typename Id<TAlignedRead>::Type                            TId;
    typedef typename TFragmentStore::TContigPos                        TContigPos;
    typedef typename TContig::TContigSeq                            TContigSeq;

    typedef AlignedReadLayout::TRows                                TRows;
    typedef typename Value<TRows>::Type                                TRow;
    typedef typename Size<TRows>::Type                                TRowsSize;
    typedef typename Iterator<TRows>::Type                            TRowsIter;
    typedef typename Iterator<TRow>::Type                            TRowIter;

    typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >    TContigGaps;
    typedef Gaps<CharString, AnchorGaps<typename TAlignedRead::TGapAnchors> >    TReadGaps;

    typedef typename DirectionIterator<TStream, Output>::Type       TOutIter;

    TContigGaps    contigGaps;
    if ((TId)contigId < length(store.contigStore))
    {
        set(contigGaps.data_source, store.contigStore[contigId].seq);
        set(contigGaps.data_gaps, store.contigStore[contigId].gaps);
//        TContigGaps    contigGaps(store.contigStore[contigId].seq, store.contigStore[contigId].gaps);
        setClippedBeginPosition(contigGaps, posBegin);
        setClippedEndPosition(contigGaps, posEnd);
        _printContig(stream, layout, contigGaps, store.contigNameStore[contigId]);
    };
    TOutIter iter = directionIterator(stream, Output());
    writeValue(iter, '\n');

    if ((TId)contigId >= length(layout.contigRows))
        return;

    if ((TRowsSize)lineEnd > length(layout.contigRows[contigId])) lineEnd = length(layout.contigRows[contigId]);
    if ((TRowsSize)lineBegin >= (TRowsSize)lineEnd) return;

    TRowsIter lit = begin(layout.contigRows[contigId], Standard()) + lineBegin;
    TRowsIter litEnd = begin(layout.contigRows[contigId], Standard()) + lineEnd;
    TReadSeq readSeq;
    CharString readSeqString;
    clearClipping(contigGaps);

    for (TNum line = 1; lit < litEnd; ++lit, ++line)
    {
        TRowIter itEnd = end(*lit, Standard());
        TRowIter left = begin(*lit, Standard());
        TRowIter right = itEnd;
        TRowIter mid = itEnd;

        while (left < right)
        {
            mid = left + (right - left) / 2;
            TAlignedRead &align = store.alignedReadStore[*mid];

            if (align.contigId < (TId)contigId || (align.contigId == (TId)contigId && (TPos)_max(align.beginPos, align.endPos) <= posBegin))
                left = mid + 1;    // what we search is in the right part
            else
                right = mid;    //            ...           left part
        }

        TPos cursor = posBegin;
        for (; mid < itEnd; ++mid)
        {
//            if (*mid >= lastRead) continue;
            TAlignedRead &align = store.alignedReadStore[*mid];
            if (align.contigId != (TId)contigId) break;

            TReadGaps   readGaps(readSeqString, align.gaps);
            TContigPos  left = align.beginPos;
            TContigPos    right = align.endPos;
            TContigPos    cBegin = _min(left, right);
            TContigPos    cEnd = _max(left, right);

            if ((TPos)cEnd <= posBegin) continue; // shouldn't occur
            if (posEnd <= (TPos)cBegin) break;

            readSeq = store.readSeqStore[align.readId];
            if (left > right)
            {
                reverseComplement(readSeq);
                readSeqString = readSeq;
                toLower(readSeqString);
            } else
                readSeqString = readSeq;

            if ((TPos)cBegin < posBegin)
                setClippedBeginPosition(readGaps, posBegin - (TPos)cBegin);
            else
            {
                TOutIter iter = directionIterator(stream, Output());
                for (; cursor < (TPos)cBegin; ++cursor)
                    writeValue(iter, ' ');
            }

            if (posEnd < (TPos)cEnd)
                setClippedEndPosition(readGaps, posEnd - (TPos)cBegin);

            _printRead(stream, layout, contigGaps, readGaps, align, line);
            cursor = cEnd;
        }
        TOutIter iter = directionIterator(stream, Output());
        writeValue(iter, '\n');
    }
}

/*!
 * @fn FragmentStore#convertMatchesToGlobalAlignment
 * @brief Converts all matches to a multiple global alignment in gap-space.
 *
 * @signature void convertMatchesToGlobalAlignment(store, score, shrinkMatches);
 *
 * @param[in,out] store The fragment store. Types: FragmentStore
 * @param[in]     score A score object used by @link globalAlignment @endlink in this function.
 * @param[in]     shrinkMatches States whether the matches should be shrinked. Types: True, False
 *
 * @section Remarks
 *
 * Before calling this function all <tt>gaps</tt> structures in @link FragmentStore::alignedReadStore @endlink and @link
 * FragmentStore::contigStore @endlink must be empty, i.e. there are no gaps in the alignments.  This function iterates
 * over entries in the @link FragmentStore::alignedReadStore @endlink and semi-global aligns each read to its contig
 * segments given by begin and end position.  Gaps introduced by these pair-wise alignments are then inserted to the
 * affected contig and reads correspondingly.
 *
 * The invariant that positions in the @link FragmentStore::alignedReadStore @endlink are in gap-space holds before
 * (there were no gaps in alignments) and after calling this functions.
 *
 * If the @link FragmentStore::alignQualityStore @endlink of the @link FragmentStore @endlink is empty when
 * <tt>convertMatchesToGlobalAlignment()</tt> is called then the @link FragmentStore::alignQualityStore @endlink is
 * filled with the edit distance scores.
 */

template <typename TSpec, typename TConfig, typename TScore, typename TShrinkMatches>
void convertMatchesToGlobalAlignment(FragmentStore<TSpec, TConfig> &store, TScore const & score, TShrinkMatches const &)
{
    typedef FragmentStore<TSpec, TConfig>                            TFragmentStore;

    //typedef typename TFragmentStore::TReadStore                        TReadStore;
    //typedef typename TFragmentStore::TReadSeqStore                    TReadSeqStore;
    typedef typename TFragmentStore::TAlignedReadStore                TAlignedReadStore;
    typedef typename TFragmentStore::TContigStore                    TContigStore;

    //typedef typename Value<TReadStore>::Type                        TRead;
    typedef typename Value<TContigStore>::Type                        TContig;
    typedef typename Value<TAlignedReadStore>::Type                    TAlignedRead;

    typedef typename TFragmentStore::TReadSeq                        TReadSeq;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type    TAlignedReadIter;
    typedef typename Id<TAlignedRead>::Type                            TId;
    typedef typename TFragmentStore::TContigPos                        TContigPos;

    typedef typename TContig::TContigSeq                            TContigSeq;
    typedef Align<TReadSeq, ArrayGaps>                                TAlign;
    typedef Gaps<TReadSeq, ArrayGaps>                                TGaps;

    typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >    TContigGaps;
    typedef Gaps<TReadSeq, AnchorGaps<typename TAlignedRead::TGapAnchors> >    TReadGaps;
    typedef typename Iterator<TContigGaps>::Type                            TContigIter;
    typedef typename Iterator<TReadGaps>::Type                                TReadIter;

    // Number of alignments with a quality score (edit distance) before computing global alignment.  We will fill the
    // store for all alignments if there are no scores yet.
    unsigned alignQualityStoreLengthPre = length(store.alignQualityStore);
    if (alignQualityStoreLengthPre == 0u)
        resize(store.alignQualityStore, length(store.alignedReadStore));

    // sort matches by increasing begin positions
    sortAlignedReads(store.alignedReadStore, SortBeginPos());
    sortAlignedReads(store.alignedReadStore, SortContigId());

    TReadSeq readSeq;
    TId lastContigId = TAlignedRead::INVALID_ID;
    TAlignedReadIter it = begin(store.alignedReadStore, Standard());
    TAlignedReadIter itEnd = end(store.alignedReadStore, Standard());
    TAlignedReadIter firstOverlap = begin(store.alignedReadStore, Standard());

//    TAlignedReadIter theIt = it;
//    for (; theIt != itEnd; ++theIt)
//    {
//        if (theIt->id == 29971)
//            break;
//    }
//    TReadSeq theReadSeq = store.readSeqStore[theIt->readId];
//    if (theIt->beginPos > theIt->endPos)
//        reverseComplement(theReadSeq);

//    TContigPos    cBeginPrev = 0;
    for (; it != itEnd;)
    {
        TContigPos    left = (*it).beginPos;
        TContigPos    right = (*it).endPos;
        TContigPos    cBegin = _min(left, right);
//    if (cBegin<cBeginPrev)
//    {
//        std::cout<<"SORTING ERROR"<<std::endl;
//    }
//    cBeginPrev=cBegin;
        TContigPos    cEnd = _max(left, right);
        TContigGaps    contigGaps(store.contigStore[(*it).contigId].seq, store.contigStore[(*it).contigId].gaps);
        TReadGaps    readGaps(readSeq, (*it).gaps);

        readSeq = store.readSeqStore[(*it).readId];
        if (left > right)
            reverseComplement(readSeq);

        // 1. Calculate pairwise alignment
        TAlign align;
        resize(rows(align), 2);
        assignSource(row(align, 0), infix(store.contigStore[it->contigId].seq, cBegin, cEnd));
        assignSource(row(align, 1), readSeq);
//        int ud = store.alignQualityStore[it->id].errors;
//        int ld = -ud;
//        if (IsSameType<TShrinkMatches, True>::VALUE)
//            globalAlignment(align, score, AlignConfig<true, false, false, true>(), ld, ud, Gotoh());
//        else
//            globalAlignment(align, score, ld, ud);

        int qualValue = 0;
        if (IsSameType<TShrinkMatches, True>::VALUE)
            qualValue = globalAlignment(align, score, AlignConfig<true, false, false, true>(), Gotoh());
        else
            qualValue = globalAlignment(align, score);
        // Update quality score entry if there were no scores before convertMatchesToGlobal() call.
        if (alignQualityStoreLengthPre > 0u && it->id >= alignQualityStoreLengthPre)
        {
            store.alignQualityStore[it->id].errors = qualValue / scoreMismatch(score);
            store.alignQualityStore[it->id].score = -qualValue / scoreMismatch(score);
        }
//    if (cBegin > 5349500 && cBegin < 5349700)
//    {
//        if (theIt->id != 29971)
//        {
//            theIt = begin(store.alignedReadStore, Standard());
//            for (; theIt != itEnd; ++theIt)
//            {
//                if (theIt->id == 29971)
//                break;
//            }
//        }
//
//        TReadGaps readGaps(theReadSeq, theIt->gaps);
//        TContigGaps    contigGaps(store.contigStore[theIt->contigId].seq, store.contigStore[theIt->contigId].gaps);
//        if ((*theIt).beginPos <= (*theIt).endPos)
//        {
//            setBeginPosition(contigGaps, (*theIt).beginPos);
//            setEndPosition(contigGaps, (*theIt).endPos);
//        } else
//        {
//            setBeginPosition(contigGaps, (*theIt).endPos);
//            setEndPosition(contigGaps, (*theIt).beginPos);
//        }
////              __int64 pos = positionGapToSeq(contigGaps, _min(theIt->beginPos, theIt->endPos)) + 1;
////        __int64 mpos = 0;
////        std::cout << "it->id == " << it->id << std::endl;
////        std::cout << "cBegin == " << cBegin << std::endl;
////        std::cout << contigGaps << std::endl;
////        std::cout << readGaps << std::endl;
////        std::cout << std::endl;
////        std::cout << std::endl;
////        std::cout << std::endl;
////        std::cout << std::endl;
//    }

//        if (infix(store.readNameStore[it->readId], 0, length("SRR049254.14375884")) == "SRR049254.14375884")
//        {
//            std::cerr << "cBegin == " << cBegin << "  id == " << (*it).id << std::endl;
//            std::cerr << align << std::endl;
//        }

        // If there is a shorter alignment then there are gaps in the beginning of the read alignment row.  In this
        // case, we move the current alignment further forward in the alignedReadStore and continue to work with the
        // next one.
        if (isGap(row(align, 1), 0))
        {
//            std::cerr << "cBegin == " << cBegin << " beginPosition(row(align, 1)) == " << beginPosition(row(align, 1)) << std::endl;
//            std::cerr << "id == " << it->id << ", read id == " << it->readId << std::endl;
//            std::cerr << align << std::endl;
            // Update aligned read element.
            cBegin += toViewPosition(row(align, 1), 0);
            bool reverse = it->beginPos > it->endPos;
            it->beginPos = cBegin;
            it->endPos = cEnd;
            if (reverse)
                std::swap(it->beginPos, it->endPos);
            // Move read element.
            TAlignedReadIter itEnd = end(store.alignedReadStore, Standard());
            TAlignedReadIter itContigEnd = upperBoundAlignedReads(it + 1, itEnd, it->contigId, SortContigId());
            TAlignedReadIter itTarget = lowerBoundAlignedReads(it + 1, itContigEnd, cBegin, SortBeginPos());
            if (itTarget != it + 1)  // Do not swap with self.
            {
                typename Value<TAlignedRead>::Type tmp = *it;
                arrayMove(it + 1, itTarget, it);
                *(itTarget - 1) = tmp;
//                if (_min((itTarget-1)->beginPos, (itTarget-1)->endPos) > _min((*itTarget).beginPos, (*itTarget).endPos))
//                {
//                std::cout << "MOVE ERROR" <<std::endl;
//                }
//                if (_min((itTarget-2)->beginPos, (itTarget-2)->endPos) > _min((itTarget-1)->beginPos, (itTarget-1)->endPos))
//                {
//                std::cout << "MOVE ERROR2" <<std::endl;
//                }
            }
            continue;
        }
        SEQAN_ASSERT_EQ(toViewPosition(row(align, 1), 0), 0u);

        // 2. Skip non-overlapping matches
        cBegin = positionSeqToGap(contigGaps, cBegin);
        if (lastContigId != (*it).contigId)
        {
            firstOverlap = it;
            lastContigId = (*it).contigId;
        } else
            while (firstOverlap != it && _max((*firstOverlap).beginPos, (*firstOverlap).endPos) <= cBegin)
                ++firstOverlap;

        // 3. Iterate over alignment
        setClippedBeginPosition(contigGaps, cBegin);

        TContigIter cIt = begin(contigGaps);
        TReadIter rIt = begin(readGaps);
        typename Iterator<TGaps>::Type it1 = begin(row(align, 0));
        typename Iterator<TGaps>::Type it2 = begin(row(align, 1));

        unsigned beginLocalContigGaps = toViewPosition(row(align, 0), 0);
        //std::cerr << "CONTIG\t" << contigGaps << "\n";
        // Heuristic (hack) for gaps in the beginning, so the following
        // does not happen:
        //
        // contig: XXXX------AAA
        //             CCC---AAA
        //                CCCAAA
        //
        // But instead, the C's are all aligned.
        if (beginLocalContigGaps > 0u) {
            unsigned i = beginLocalContigGaps;
            do {
                goPrevious(cIt);
                i -= 1;
                cBegin -= 1;
            } while (isGap(cIt) && i > 0u);
            if (!isGap(cIt)) {
                goNext(cIt);
                cBegin += 1;
            }
        }

        // We will clip off trailing gaps in the read row.
        unsigned charsEndPos = toViewPosition(row(align, 1), length(source(row(align, 1))));
        setClippedEndPosition(row(align, 0), charsEndPos);

//        std::cerr << "firstOverlap - theIt == " << firstOverlap - theIt << std::endl;
//        if (firstOverlap - theIt > 0)
//        {
//            if (_min(firstOverlap->beginPos, firstOverlap->endPos) < _max(theIt->beginPos, theIt->endPos))
//                std::cerr << "INVARIANT DOES NOT HOLD" << std::endl;
//        }

        for (; !atEnd(cIt) && !atEnd(it1); goNext(cIt), goNext(rIt))
        {
            bool isGapContig = isGap(cIt);
            if (isGapContig != isGap(it1))
            {
                if (isGapContig)
                {
                    // *** gap in contig of the global alignment ***
                    // copy exisiting contig gap
                    insertGaps(rIt, 1);
                    continue;
                } else
                {
                    // *** gap in contig of the pairwise alignment ***
                    // insert padding gaps in contig and reads
                    TContigPos insPos = cIt.current.gapPos;
                    insertGaps(cIt, 1);
                    for (TAlignedReadIter j = firstOverlap; j != it; ++j)
                    {
                        TContigPos rBegin = _min(j->beginPos, j->endPos);
                        TContigPos rEnd = _max(j->beginPos, j->endPos);
                        if (rBegin < insPos && insPos < rEnd)
                        {
//                    if (j->id == 29971)
//                    {
//                        std::cerr << "OK computer" << std::endl;
//                    }
                            if (rBegin < insPos)
                            {
                                TReadGaps gaps(store.readSeqStore[j->readId], j->gaps);
                                insertGap(gaps, insPos - rBegin);
                            } else
                            {
                                // shift beginPos if insertion was at the front of the read
                                if (j->beginPos < j->endPos)
                                    ++j->beginPos;
                                else
                                    ++j->endPos;
                            }
                            // shift endPos as the alignment was elongated or shifted
                            if (j->beginPos < j->endPos)
                                ++j->endPos;
                            else
                                ++j->beginPos;
                        } else if (insPos <= rBegin) {
                            ++j->endPos;
                            ++j->beginPos;
                        }
                    }
                }
            }
            if (isGap(it2))
            {
                // *** gap in read of pairwise alignment ***
                // copy gaps from alignment
                insertGaps(rIt, 1);
            }
            goNext(it1);
            goNext(it2);
        }

        // store new gap-space alignment borders
        cEnd = cBegin + length(readGaps);
        if (left < right)
        {
            (*it).beginPos = cBegin;
            (*it).endPos = cEnd;
        } else
        {
            (*it).beginPos = cEnd;
            (*it).endPos = cBegin;
        }

        ++it;  // Go to next alignment.
/*
//        if (interesting)
        {
            String<String<unsigned> > layout;
            layoutAlignment(layout, store, (*it).contigId);
            std::cout << store.readNameStore[(*it).readId] << std::endl;
            std::cout << readGaps << '\t' << cBegin << '\t' << cEnd << std::endl << std::endl;
            printAlignment(std::cout, layout, store, (*it).contigId, (int)cBegin-20, (int)cEnd+20, 0, 40, 1 + (it - begin(store.alignedReadStore, Standard())));
//            getc(stdin);
        }
*/
//        if (store.readNameStore[(*it).readId] == "read3305")
//            return;
    }

    // AlignedReadLayout layout;
    // layoutAlignment(layout, store);
    // std::cerr << "(int)length(store.contigStore[0].gaps) == " << (int)length(store.contigStore[0].gaps) << '\n';
    // printAlignment(std::cout, Raw(), layout, store, 0, -10, (int)(length(store.contigStore[0].seq) * 1.1), 0, 40);
}

/*!
 * @fn FragmentStore#convertPairWiseToGlobalAlignment
 * @brief Converts pairwise alignments to a multiple global alignment.
 *
 * @signature void convertPairWiseToGlobalAlignment(store, pairwiseContigGaps);
 *
 * @param[in,out] store The fragment store. Types: FragmentStore
 * @param[in,out] pairwiseContigGaps
 *                      A @link String @endlink of anchored contig gaps for every pairwise alignment.
 *
 * @section Remarks
 *
 * Before calling this function the <tt>gaps</tt> structures in the @link FragmentStore::contigStore @endlink must be
 * empty, i.e. there are no gaps in the contig. The pairwise alignment gaps of the reads are stored in the <tt>gaps</tt>
 * structure in the @link FragmentStore::alignedReadStore @endlink, whereas the pairwise alignment gaps of the contig
 * are stored in the <tt>pairwiseContigGaps</tt> string.
 *
 * After calling this functions all positions in the @link FragmentStore::alignedReadStore @endlink are in gap-space.
 */


template <typename TFragmentStore>
struct LessConvertPairWiseToGlobalAlignment:
    public std::binary_function<
        typename Value<typename TFragmentStore::TAlignedReadStore>::Type,
        typename Value<typename TFragmentStore::TAlignedReadStore>::Type,
        bool>
{
    typedef typename TFragmentStore::TAlignedReadStore      TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type         TAlignedRead;
    typedef typename Size<TAlignedReadStore>::Type          TSize;

    typedef typename TFragmentStore::TContigStore           TContigStore;
    typedef typename Value<TContigStore>::Type              TContig;
    typedef typename TFragmentStore::TContigPos             TContigPos;
    typedef typename TContig::TContigSeq                    TContigSeq;
    typedef typename TContig::TGapAnchors                   TContigAnchors;
    typedef Gaps<TContigSeq, AnchorGaps<TContigAnchors> >   TContigGaps;

    TFragmentStore &store;
    TSize offset;

    LessConvertPairWiseToGlobalAlignment(TFragmentStore &store_, TSize offset_):
        store(store_),
        offset(offset_) {}

    inline bool
    operator() (TAlignedRead const& a1, TAlignedRead const& a2) const
    {
//        if (a1.contigId < a2.contigId) return true;
//        if (a1.contigId > a2.contigId) return false;
//
//        TContigPos pos1 = _min(a1.beginPos, a1.endPos);
//        TContigPos pos2 = _min(a2.beginPos, a2.endPos);
//        TContig & contig = store.contigStore[a1.contigId];
//
//        if (a1.id < offset)
//        {
//            TContigGaps    contigGaps(contig.seq, contig.gaps);
//            pos1 = positionGapToSeq(contigGaps, pos1);
//        }
//
//        if (a2.id < offset)
//        {
//            TContigGaps    contigGaps(contig.seq, contig.gaps);
//            pos2 = positionGapToSeq(contigGaps, pos2);
//        }
//
        if (a1.id <  offset && a2.id >= offset) return true;
        if (a1.id >= offset && a2.id <  offset) return false;

        if (a1.contigId < a2.contigId) return true;
        if (a1.contigId > a2.contigId) return false;

        TContigPos pos1 = _min(a1.beginPos, a1.endPos);
        TContigPos pos2 = _min(a2.beginPos, a2.endPos);
        return pos1 < pos2;
    }
};


template <typename TAlignedReadIter, typename TPos, typename TSize, typename TReadSeqStore>
void _alignedReadsInsertGaps(
    TAlignedReadIter alignedReadsBegin,
    TAlignedReadIter alignedReadsEnd,
    TPos insPos,
    TSize blkLen,
    TReadSeqStore &readSeqStore)
{
    typedef typename Value<TReadSeqStore>::Type                             TReadSeq;
    typedef typename Value<TAlignedReadIter>::Type                          TAlignedRead;
    typedef Gaps<TReadSeq, AnchorGaps<typename TAlignedRead::TGapAnchors> > TReadGaps;
    typedef __int64                                                         TContigPos;

    for (TAlignedReadIter j = alignedReadsBegin; j != alignedReadsEnd; ++j)
    {

        TContigPos rBegin = _min(j->beginPos, j->endPos);
        TContigPos rEnd = _max(j->beginPos, j->endPos);
        if (rBegin < insPos && insPos < rEnd)
        {
            if (rBegin < insPos)
            {
                TReadGaps gaps(readSeqStore[j->readId], j->gaps);
                insertGaps(gaps, insPos - rBegin, blkLen);
            }
            else
            {
                // shift beginPos if insertion was at the front of the read
                if (j->beginPos < j->endPos)
                    j->beginPos += blkLen;
                else
                    j->endPos += blkLen;
            }
            // shift endPos as the alignment was elongated or shifted
            if (j->beginPos < j->endPos)
                j->endPos += blkLen;
            else
                j->beginPos += blkLen;
        }
    }
}

// The following function merges a pairwise alignment (contigGapsPW, readGaps) into the multiple alignment of
// (contigGapsGlobal, firstOverlappingAlignedRead, firstOverlappingAlignedRead+1, ... , currentAlignedRead - 1)
// It distinguishes two cases:
//
// 1. Contig gaps that are in the global but not in the pairwise alignment:
//    Insert these gaps into the read of the pairwise alignment
//
// 2. Contig gaps that are in the pairwise but not in the multiple alignment:
//    Insert these gaps into the global contig and all not yet considered overlapping reads
//    (firstOverlappingAlignedRead ... currentAlignedRead - 1)
//
template <typename TContigGapsGlobal, typename TContigGapsPW, typename TReadGaps, typename TAlignedReadIter, typename TReadSeqStore>
void _twoWayMergeAlignments(
    TContigGapsGlobal &contigGapsGlobal,
    TContigGapsPW &contigGapsPW,
    TReadGaps &readGaps,
    TAlignedReadIter overlappingAlignedReadsBegin1,
    TAlignedReadIter overlappingAlignedReadsEnd1,
    TAlignedReadIter overlappingAlignedReadsBegin2,
    TAlignedReadIter overlappingAlignedReadsEnd2,
    TReadSeqStore &readSeqStore)
{
    // gap iterators
    typedef typename Iterator<TContigGapsGlobal>::Type  TContigGlobalIter;
    typedef typename Iterator<TContigGapsPW>::Type      TContigPWIter;
    typedef typename Iterator<TReadGaps>::Type          TReadIter;
    typedef __int64                                     TContigPos;

    TContigGlobalIter cIt = begin(contigGapsGlobal);
    TContigPWIter pIt = begin(contigGapsPW);
    TReadIter rIt = begin(readGaps);

    /*
    std::cout << "contigGlobal\t" << contigGapsGlobal << std::endl;
    std::cout << "contigPW    \t" << contigGapsPW << std::endl;
    std::cout << "readPW      \t" << readGaps << std::endl;
    */

    __int64 blkLen = 0;
    for (; !atEnd(rIt); goFurther(rIt, blkLen), goFurther(cIt, blkLen))
    {
        bool isGapContig = isGap(cIt);
        bool isGapLocalContig = isGap(pIt);
        blkLen = _min(blockLength(cIt), blockLength(pIt));
        SEQAN_ASSERT_GT(blkLen, 0u);
//          SEQAN_ASSERT_LT(blkLen, length(contigGapsGlobal));

        if (isGapContig != isGapLocalContig)
        {
            if (isGapContig)
            {
                // CASE 1.
                //
                // *** gap in contig of the global alignment ***
                // copy exisiting contig gap
                insertGaps(rIt, blkLen);
                continue;
            }
            else
            {
                // CASE 2.
                //
                // *** gap in contig of the pairwise alignment ***
                // insert padding gaps in contig and reads
                TContigPos insPos = cIt.current.gapPos;
                insertGaps(cIt, blkLen);
                _alignedReadsInsertGaps(overlappingAlignedReadsBegin1, overlappingAlignedReadsEnd1, insPos, blkLen, readSeqStore);
                _alignedReadsInsertGaps(overlappingAlignedReadsBegin2, overlappingAlignedReadsEnd2, insPos, blkLen, readSeqStore);
            }
        }

        // fast forward the whole block
        goFurther(pIt, blkLen);
    }
}

template <typename TSpec, typename TConfig, typename TContigGapsString>
void convertPairWiseToGlobalAlignment(FragmentStore<TSpec, TConfig> &store, TContigGapsString &gaps)
{
    typedef FragmentStore<TSpec, TConfig>                            TFragmentStore;

    // stores
    //typedef typename TFragmentStore::TReadStore                        TReadStore;
    typedef typename TFragmentStore::TReadSeqStore                    TReadSeqStore;
    typedef typename TFragmentStore::TAlignedReadStore                TAlignedReadStore;
    typedef typename TFragmentStore::TContigStore                    TContigStore;
    //typedef typename TFragmentStore::TContigSeq                      TContigSeq;

    // store elements
    //typedef typename Value<TReadStore>::Type                        TRead;
    typedef typename Value<TReadSeqStore>::Type                        TReadSeq;
    typedef typename Value<TContigStore>::Type                        TContig;
    typedef typename Value<TAlignedReadStore>::Type                    TAlignedRead;

    typedef typename Iterator<TAlignedReadStore, Standard>::Type    TAlignedReadIter;
    typedef typename Id<TAlignedRead>::Type                            TId;
    typedef typename TFragmentStore::TContigPos                        TContigPos;

    // gap structures
    //typedef Gaps<TContigSeq/*Nothing*/, AnchorGaps<typename TContig::TGapAnchors> >            TContigGapsGlobal;
    //typedef Gaps<TContigSeq/*Nothing*/, AnchorGaps<typename Value<TContigGapsString>::Type> >    TContigGapsPW;
    typedef Gaps</*TContigSeq*/Nothing, AnchorGaps<typename TContig::TGapAnchors> >             TContigGapsGlobal;
    typedef Gaps</*TContigSeq*/Nothing, AnchorGaps<typename Value<TContigGapsString>::Type> >    TContigGapsPW;
    typedef Gaps<TReadSeq, AnchorGaps<typename TAlignedRead::TGapAnchors> >                     TReadGaps;

    typename Size<TAlignedReadStore>::Type offsetRecord = length(store.alignedReadStore) - length(gaps);

    // sort matches by increasing begin positions
    LessConvertPairWiseToGlobalAlignment<TFragmentStore> lessFunctor(store, offsetRecord);
    sortAlignedReads(store.alignedReadStore, lessFunctor);

    TReadSeq readSeq;
    TId lastContigId = TAlignedRead::INVALID_ID;
    TAlignedReadIter it = begin(store.alignedReadStore, Standard());
    TAlignedReadIter itEnd = end(store.alignedReadStore, Standard());

    TAlignedReadIter firstOverlapOld = begin(store.alignedReadStore, Standard());
    TAlignedReadIter lastOverlapOld = firstOverlapOld;
    TAlignedReadIter firstNew = begin(store.alignedReadStore, Standard()) + offsetRecord;
    TAlignedReadIter firstOverlapNew = firstNew;

    for (; it != itEnd; ++it)
    {
        if ((*it).id < offsetRecord)
            continue;

        TContigPos    left = it->beginPos;
        TContigPos    right = it->endPos;
        TContigPos    cBegin = _min(left, right);
        TContigPos    cEnd = _max(left, right);

        // 1. Initialize gap structures
        TContigGapsGlobal   contigGapsGlobal(/*store.contigStore[it->contigId].seq, */store.contigStore[it->contigId].gaps);
        /*
        TContigSeq contigInfix = infix(store.contigStore[it->contigId].seq, cBegin, cEnd);
        if (left > right)
            reverseComplement(contigInfix);
        */
        TContigGapsPW       contigGapsPW(/*contigInfix, */gaps[it->id - offsetRecord]);
        TReadGaps           readGaps(store.readSeqStore[it->readId], it->gaps);

        SEQAN_ASSERT(dependent(contigGapsGlobal.data_gaps));
        SEQAN_ASSERT(dependent(readGaps.data_gaps));


        // 1111111111111
        // 22222222
        //   33333
        //     44444444444
        //       555555
        //  |cBegin---cEnd)

        TId contigId = it->contigId;
        cBegin = positionSeqToGap(contigGapsGlobal, cBegin);

        // 2a. Skip non-overlapping old alignments
        while (firstOverlapOld != firstNew &&
               (firstOverlapOld->contigId < contigId || _max(firstOverlapOld->beginPos, firstOverlapOld->endPos) <= cBegin))
            ++firstOverlapOld;

        if (lastOverlapOld < firstOverlapOld)
            lastOverlapOld = firstOverlapOld;

        while (lastOverlapOld != firstNew &&
               (lastOverlapOld->contigId < contigId || _min(lastOverlapOld->beginPos, lastOverlapOld->endPos) < cEnd))
            ++lastOverlapOld;

        // 2b. Skip non-overlapping new matches
        if (lastContigId != contigId)
        {
            firstOverlapNew = it;
            lastContigId = contigId;
        }
        else
        {
            while (firstOverlapNew != it && _max(firstOverlapNew->beginPos, firstOverlapNew->endPos) <= cBegin)
                ++firstOverlapNew;
        }

        // 3. Iterate over alignment
        setClippedBeginPosition(contigGapsGlobal, cBegin);

        _twoWayMergeAlignments(
                contigGapsGlobal, contigGapsPW, readGaps,
                firstOverlapOld, lastOverlapOld, firstOverlapNew, it,
                store.readSeqStore);

        // store new gap-space alignment borders
        cEnd = cBegin + length(readGaps);
        if (left < right)
        {
            (*it).beginPos = cBegin;
            (*it).endPos = cEnd;
        }
        else
        {
            (*it).beginPos = cEnd;
            (*it).endPos = cBegin;
        }
    }
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
