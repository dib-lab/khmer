// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

#include <iostream>

#ifndef SEQAN_HEADER_STORE_IO_SAM_H
#define SEQAN_HEADER_STORE_IO_SAM_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// Parsing functions
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// _parseReadCigar
    
    template <typename TFile, typename TCigarString, typename TChar>
    inline void
    _parseReadCigar(TFile & file, TCigarString & cigar, TChar & c)
    {
//IOREV _nodoc_ _hasCRef_ could be simplified by using other _is or _parse calls
        typedef typename Value<TCigarString>::Type  TCigarElement;
        typedef typename TCigarElement::TCount      TCount;

        clear(cigar);
        
        // if the CIGAR is not set and '*'
        if (c == '*')
        {
            c = _streamGet(file);
            return;
        }
        
        while (!_streamEOF(file)) 
        {
            TCount count = _parseReadNumber(file, c);
            if (c >= 'a' && c <= 'z')
                c = c + 'A' - 'a';
            appendValue(cigar, TCigarElement(c, count));
            
            c = _streamGet(file);
            if (c == ' ' || c == '\t' || c == '\n') break;
        }
    }
    
//////////////////////////////////////////////////////////////////////////////
// _parseReadSamIdentifier
    
    template<typename TFile, typename TString, typename TChar>
    inline void
    _parseReadSamIdentifier(TFile & file, TString & str, TChar& c)
    {
//IOREV _nodoc_ _hasCRef_ _duplicate_ same as generic _parseUntilWhitespace?
        if (c == ' ' || c == '\t' || c == '\n') return;
        appendValue(str, c);
        while (!_streamEOF(file)) 
        {
            c = _streamGet(file);
            if (c == ' ' || c == '\t' || c == '\n') return;
            appendValue(str, c);
        }
    }
    
//////////////////////////////////////////////////////////////////////////////
// _parseIsDna
    
    template<typename TChar>
    inline bool
    _parseIsDna(TChar const & c)
    {
//IOREV _bug_ according to Dna5 doc c is already uppercase, so it should be (c == x) || (c - 'A' + 'a' == x) OR JUST: return tolower(c) == tolower(x);
        char x = TChar(Dna5(c));
        return (c == x) || (c + 'A' - 'a' == x);
    }
    
//////////////////////////////////////////////////////////////////////////////
//_parseReadDnaSeq
    
    template<typename TFile, typename TString, typename TChar>
    inline void
    _parseReadDnaSeq(TFile & file, TString & str, TChar & c)
    {
//IOREV _nodoc_ _hasCRef_ _duplicate_ same as generic _parseUntilWhitespace?
        TChar first = c;
        if (!_streamEOF(file)) 
            c = _streamGet(file);

        if (!_parseIsDna(first))
            return;
        appendValue(str, first, Generous());
        
        for (; !_streamEOF(file) && _parseIsDna(c); c = _streamGet(file))
            appendValue(str, c, Generous());
    }
        
//////////////////////////////////////////////////////////////////////////////
// _parseIsPhredQual
    
    template <typename TChar>
    inline bool
    _parseIsPhredQual(TChar c)
    {
//IOREV _nodoc_ what does the title mean? same as return isprint(c) && c != ' '
        return c >= '!' && c <= '~';
    }
    
//////////////////////////////////////////////////////////////////////////////
// _parseReadSeqQual
//
    
    template<typename TFile, typename TQualString, typename TChar>
    inline void
    _parseReadSeqQual(TFile & file, TQualString & str, TChar & c)
    {
//IOREV _nodoc_ unclearn what this actually does
        typedef typename Size<TQualString>::Type                TSize;
        typedef typename Iterator<TQualString, Standard>::Type  TIter;
        
        if (!_parseIsPhredQual(c)) return;

        TIter itBegin = begin(str, Standard());
        TIter it = itBegin; 
        TSize rest = length(str);
        
        do {
            int q = c - '!';
            if (!_streamEOF(file)) 
                c = _streamGet(file);
            else
                if (rest > 1)
                    rest = 1;
            
            if (q == '*' - '!' && !_parseIsPhredQual(c) && it == itBegin)
                return;
            
            if (rest != 0)
            {
                assignQualityValue(*it, q);
                ++it;
                --rest;
            }
        } while (_parseIsPhredQual(c) && !_streamEOF(file));
    }
    
//////////////////////////////////////////////////////////////////////////////
// _parseReadCharsUntilEndOfLine
//
// Reads all symbols till the next '\n' and writes them in the CharString str
// the c is the first character after the '\n'.
    
    template<typename TFile, typename TChar>
    inline void
    _parseReadCharsUntilEndOfLine(TFile & file, String<char> & str, TChar& c)
    {
//IOREV replace with generic _parseLine() function that ignores '\r' (this calls adds a lonely '\r' to the end of str if file is windows-formatted)
        // read all chars till '\n'
        while (c != '\n')
        {
            appendValue(str, c, Generous());
            if (_streamEOF(file)) return;
            c = _streamGet(file);
        }
        
        // read the first char after the '\n'
        if (!_streamEOF(file))
            c = _streamGet(file);
    }

//////////////////////////////////////////////////////////////////////////////
// appendAlignment
    
    template<typename TSpec, typename TConfig, typename TId, typename TPos, typename TGaps>
    inline typename Size<typename FragmentStore<TSpec, TConfig>::TAlignedReadStore>::Type
    appendAlignment(
        FragmentStore<TSpec, TConfig> & fragStore, 
        TId readId, 
        TId contigId, 
        TPos beginPos, 
        TPos endPos, 
        TGaps const & gaps)
    {
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
        typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;
        
        TId id = length(fragStore.alignedReadStore);
        TAlignedElement alignedElem = TAlignedElement(id, readId, contigId, beginPos, endPos, gaps);
        appendValue(fragStore.alignedReadStore, alignedElem);
        
        return id;
    }
    


//////////////////////////////////////////////////////////////////////////////
// read functions for Sam
//////////////////////////////////////////////////////////////////////////////

    
//////////////////////////////////////////////////////////////////////////////
// _generatePairMatchIds
//
    template <typename TPos, typename TId>
    struct MatchMateInfo_
    {
        TId     readId;
        TId     contigId;
        TId     pairMatchId;
        TId     matePairId;//:(sizeof(TId)*8-1);
        bool    reversed;
        TPos    beginPos;
    };
    
    template <typename TFragStore>
    struct AlignedMateLess_
    {
        TFragStore &fragStore;
        
        AlignedMateLess_(TFragStore &fragStore_) :
            fragStore(fragStore_) {}

        template <typename TAlignedRead>
        inline bool 
        operator() (TAlignedRead const& a, TAlignedRead const& b) const 
        {
            if (a.contigId < b.contigId) return true;
            if (a.contigId > b.contigId) return false;

            typename TFragStore::TContigPos posA = _min(a.beginPos, a.endPos);
            typename TFragStore::TContigPos posB = _min(b.beginPos, b.endPos);            
            if (posA < posB) return true;
            if (posA > posB) return false;
            
            bool reversedA = (a.beginPos > a.endPos);
            bool reversedB = (b.beginPos > b.endPos);
            if (reversedA != reversedB) return reversedB;

            typedef typename TFragStore::TMatePairStore     TMatePairStore;
            typedef typename Value<TMatePairStore>::Type    TMatePair;
            typename TMatePair::TId matePairIdB = TMatePair::INVALID_ID;

            if (a.readId >= length(fragStore.readStore))
                return false;
                
            if (b.readId < length(fragStore.readStore))
                matePairIdB = fragStore.readStore[b.readId].matePairId;

            return (fragStore.readStore[a.readId].matePairId < matePairIdB);
        }
    };

    struct MatchMateInfoLess_
    {
        template <typename TMInfo>
        inline bool 
        operator() (TMInfo const &a, TMInfo const &b) const 
        {
            if (a.contigId < b.contigId) return true;
            if (a.contigId > b.contigId) return false;
            if (a.beginPos < b.beginPos) return true;
            if (a.beginPos > b.beginPos) return false;
            if (a.reversed != b.reversed) return b.reversed;
            return (a.matePairId < b.matePairId);
        }
    };
    
    template <typename TAlignedRead, typename TMInfo, typename TFragStore>
    inline int 
    _compareAlignedReadAndMateInfo(TAlignedRead const &a, TMInfo const &b, TFragStore const &fragStore)
    {
        if (a.contigId < b.contigId) return -1;
        if (a.contigId > b.contigId) return 1;

        typename TFragStore::TContigPos posA = _min(a.beginPos, a.endPos);
        if (posA < b.beginPos) return -1;
        if (posA > b.beginPos) return 1;
        
        bool reversedA = (a.beginPos > a.endPos);
        if (!reversedA && b.reversed) return -1;
        if (reversedA && !b.reversed) return 1;

        typedef typename TFragStore::TMatePairStore     TMatePairStore;
        typedef typename Value<TMatePairStore>::Type    TMatePair;
        typename TMatePair::TId matePairIdA = TMatePair::INVALID_ID;

        if (a.readId < length(fragStore.readStore))
            matePairIdA = fragStore.readStore[a.readId].matePairId;
            
        if (matePairIdA < b.matePairId) return -1;
        if (matePairIdA > b.matePairId) return 1;
        return 0;
    }

    template<typename TSpec, typename TConfig, typename TMatchMateInfos>
    inline void 
    _generatePairMatchIds (
        FragmentStore<TSpec, TConfig> & fragStore,
        TMatchMateInfos & matchMateInfos)
    {
        typedef FragmentStore<TSpec, TConfig>                           TFragmentStore;        
        typedef typename TFragmentStore::TAlignedReadStore              TAlignedReadStore;
        typedef typename Value<TAlignedReadStore>::Type                 TAlignedRead;
        typedef typename Iterator<TAlignedReadStore, Standard>::Type    TIter;    
        typedef typename Iterator<TMatchMateInfos, Standard>::Type      TMIter;    
                
        TIter it = begin(fragStore.alignedReadStore, Standard());
        TIter itEnd = end(fragStore.alignedReadStore, Standard());
        TMIter mit = begin(matchMateInfos, Standard());
        TMIter mitEnd = end(matchMateInfos, Standard());

        if (it == itEnd || mit == mitEnd) return;

        // sort the aligned read store by: begin position, contig name
        std::sort(it,  itEnd,  AlignedMateLess_<TFragmentStore>(fragStore));
        std::sort(mit, mitEnd, MatchMateInfoLess_());

        while (true)
        {
            // skip already aligned reads
            while (it->pairMatchId != TAlignedRead::INVALID_ID)
                if (++it == itEnd) return;

            int cmp = _compareAlignedReadAndMateInfo(*it, *mit, fragStore);

            if (cmp == 0)   // both are equal -> link them
                (*it).pairMatchId = (*mit).pairMatchId;

            if (cmp >= 0)   // MateInfo is less or equal
            {
                if (++mit == mitEnd) return;
//                if (cmp>0)
//                    std::cout << "mateInfo:   contigId="<<mit->contigId<<"  beginPos="<<mit->beginPos<<"  reversed="<<mit->reversed<<"  matePairId="<<mit->matePairId<<std::endl;
            }

            if (cmp <= 0)   // AlignedRead is less or equal
            {
                if (++it == itEnd) return;
//                if (cmp<0)
//                    std::cout << "alignedR:   contigId="<<it->contigId<<"  beginPos="<<_min(it->beginPos,it->endPos)<<"  reversed="<<(it->beginPos > it->endPos)<<"  matePairId="<<fragStore.readStore[it->readId].matePairId<<std::endl;
            }
        }
    }    

//////////////////////////////////////////////////////////////////////////////
// read

///.Function.read.param.tag.type:Tag.File Format.tag.Sam
    
    struct FragStoreImportFlags
    {
        bool importRead:1;
        bool importReadSeq:1;
        bool importReadName:1;
        bool importReadAlignment:1;
        bool importReadAlignmentQuality:1;
        bool importReadAlignmentTags:1;

        FragStoreImportFlags():
            importRead(true),
            importReadSeq(true),
            importReadName(true),
            importReadAlignment(true),
            importReadAlignmentQuality(true),
            importReadAlignmentTags(true)
        {}
    };

    inline void
    clear(FragStoreImportFlags & flags)
    {
        flags.importRead = false;
        flags.importReadSeq = false;
        flags.importReadName = false;
        flags.importReadAlignment = false;
        flags.importReadAlignmentQuality = false;
        flags.importReadAlignmentTags = false;
    }

    template<typename TFile, typename TSpec, typename TConfig>
    inline void 
    read(
        TFile & file,
        FragmentStore<TSpec, TConfig> & fragStore,
        Sam,
        FragStoreImportFlags const & importFlags)
    {
//IOREV not sure if recordreading or batchreading
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
        typedef typename TFragmentStore::TContigPos TContigPos;
        typedef typename Id<TFragmentStore>::Type TId;
        
        // data structure to temporarily store the gaps that need to be inserted in the contig sequences
        typedef MatchMateInfo_<TContigPos, TId> TMatchMateInfo;
        typedef String<TMatchMateInfo> TMatchMateInfos;
        typedef StringSet<String<typename TFragmentStore::TContigGapAnchor>, Owner<ConcatDirect<> > > TContigAnchorGaps;

        // data structure to temporarily store information about match mates
        TMatchMateInfos matchMateInfos;
        TContigAnchorGaps contigAnchorGaps;
        
        if (_streamEOF(file)) return;

        // get first character from the stream
        char c = _streamGet(file);
        
        // Read in header section
        _readHeader(file, fragStore, c, Sam());
        
        // Read in alignments section
        _readAlignments(file, fragStore, contigAnchorGaps, matchMateInfos, c, Sam(), importFlags);
        
        if (importFlags.importReadAlignment)
        {
            // set the match mate IDs using the information stored in matchMateInfos
            _generatePairMatchIds(fragStore, matchMateInfos);

            convertPairWiseToGlobalAlignment(fragStore, contigAnchorGaps);
        }
    }

    template<typename TFile, typename TSpec, typename TConfig>
    inline void 
    read(
        TFile & file,
        FragmentStore<TSpec, TConfig> & fragStore,
        Sam)
    {
        read(file, fragStore, Sam(), FragStoreImportFlags());
    }

//////////////////////////////////////////////////////////////////////////////
// _readHeader

    template<typename TFile, typename TSpec, typename TConfig, typename TChar>
    inline void 
    _readHeader (
        TFile & file,
        FragmentStore<TSpec, TConfig> &,
        TChar & c,
        Sam)
    {
//IOREV _stub_ this isn't implemented yet
        // skip header for now
        while (c == '@')
            _parseSkipLine(file, c);
    }

    
//////////////////////////////////////////////////////////////////////////////
// _readAlignments
//
// reads in alignement sections from a Sam file

    template <typename TFragmentStore>
    struct FragStoreSAMContext
    {
        typedef typename Id<TFragmentStore>::Type                                   TId;
        typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type    TAlignedElement;
        typedef typename TAlignedElement::TGapAnchors                               TReadGapAnchors;
        typedef String<typename TFragmentStore::TContigGapAnchor>                   TContigAnchorGaps;

        TId                 readId;
        TId                 contigId;
        TReadGapAnchors     readGapAnchors;
        TContigAnchorGaps   contigGapAnchors;
    };


    template<typename TFile, typename TSpec, typename TConfig, typename TContigAnchorGaps, typename TMatchMateInfos, typename TChar>
    inline void 
    _readAlignments (
        TFile & file,
        FragmentStore<TSpec, TConfig> & fragStore,
        TContigAnchorGaps & contigAnchorGaps,   
        TMatchMateInfos & matchMateInfos,
        TChar & c,
        Sam,
        FragStoreImportFlags const & importFlags)
    {
//IOREV _nodoc_ docusmentation in code, but unclear
        // create dummy entries in Sam specific aligned read quality store and aligned read tag store
        // is needed so the ID in the aligned store can be use to access the other stores
        // even if there exists previous entries without
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
        typedef typename TFragmentStore::TAlignQualityStore TAlignQualityStore;
        typedef typename TFragmentStore::TReadSeqStore TReadSeqStore;
        typedef typename Size<TReadSeqStore>::Type TReadSeqStoreSize;
        typedef typename Value<TAlignQualityStore>::Type TAlignQuality;
        
        TAlignQuality q;
        q.score = maxValue(q.score);
        int diff = length(fragStore.alignedReadStore) - length(fragStore.alignQualityStore);
        for(int i = 0; i < diff; ++i)
            appendValue(fragStore.alignQualityStore, q, Generous());
        
        diff = length(fragStore.alignedReadStore) - length(fragStore.alignedReadTagStore);
        for(int i = 0; i < diff; ++i)
            appendValue(fragStore.alignedReadTagStore, "", Generous());
        
        // read in alignments
        FragStoreSAMContext<TFragmentStore> contextSAM;
        refresh(fragStore.contigNameStoreCache);
        refresh(fragStore.readNameStoreCache);

        while (!_streamEOF(file))
            _readOneAlignment(file, fragStore, contigAnchorGaps, matchMateInfos, c, contextSAM, importFlags);

        if (importFlags.importReadSeq)
        {
            TReadSeqStoreSize emptyReads = 0;
            for(TReadSeqStoreSize i = 0; i < length(fragStore.alignedReadStore); ++i)
                if (empty(fragStore.readSeqStore[fragStore.alignedReadStore[i].readId]))
                {
                    ++emptyReads;
    //                std::cerr << "Read sequence empty for " << fragStore.readNameStore[fragStore.alignedReadStore[i].readId] << std::endl;
                }
            if (emptyReads != 0)
                std::cerr << "Warning: " << emptyReads << " read sequences are empty." << std::endl;
        }
    }
    
        
    template <typename TReadSeq, typename TCigar, typename TPos, typename TId, typename TFragmentStore>
    inline void
    _samAppendAlignment(
        TFragmentStore &fragStore,
        TReadSeq const &readSeq,
        TCigar &cigar,
        TPos &beginPos, TPos &endPos,
        TId &pairMatchId,
        FragStoreSAMContext<TFragmentStore> & contextSAM)
    {
        typedef typename TFragmentStore::TAlignedReadStore                      TAlignedReadStore;
        typedef typename Value<TAlignedReadStore>::Type                         TAlignedRead;
        typedef Gaps<TReadSeq, AnchorGaps<typename TAlignedRead::TGapAnchors> > TReadGaps;

        // insert alignment gaps
        clear(contextSAM.readGapAnchors);
        TReadGaps readGaps(readSeq, contextSAM.readGapAnchors);
        unsigned beginGaps = cigarToGapAnchorRead(cigar, readGaps);

        // adapt start or end (on reverse strand) position if alignment begins with gaps
        if (beginPos > endPos)
            endPos += beginGaps;
        else
            beginPos += beginGaps;

        // create a new entry in the aligned read store
        pairMatchId = appendAlignment(fragStore, contextSAM.readId, contextSAM.contigId, beginPos, endPos, contextSAM.readGapAnchors);
    }

    template <typename TCigar, typename TPos, typename TId, typename TFragmentStore>
    inline void
    _samAppendAlignmentWithoutSeq(
        TFragmentStore &fragStore,
        TCigar &cigar,
        TPos &beginPos, TPos &endPos,
        TId &pairMatchId,
        FragStoreSAMContext<TFragmentStore> & contextSAM)
    {
        typedef typename TFragmentStore::TAlignedReadStore                      TAlignedReadStore;
        typedef typename Value<TAlignedReadStore>::Type                         TAlignedRead;
        typedef Gaps<Nothing, AnchorGaps<typename TAlignedRead::TGapAnchors> >  TReadGaps;
        Nothing nothing;

        // insert alignment gaps
        clear(contextSAM.readGapAnchors);
        TReadGaps readGaps(nothing, contextSAM.readGapAnchors);
        unsigned beginGaps = cigarToGapAnchorRead(cigar, readGaps);

        // adapt start or end (on reverse strand) position if alignment begins with gaps
        if (beginPos > endPos)
            endPos += beginGaps;
        else
            beginPos += beginGaps;

        // create a new entry in the aligned read store
        pairMatchId = appendAlignment(fragStore, contextSAM.readId, contextSAM.contigId, beginPos, endPos, contextSAM.readGapAnchors);
    }

//////////////////////////////////////////////////////////////////////////////
// _readOneAlignment
//
// reads in one alignement section from a Sam file
    
    template <
        typename TFile,
        typename TSpec,
        typename TConfig,
        typename TContigAnchorGaps,
        typename TMatchMateInfos,
        typename TChar,
        typename TFragStore>
    inline void
    _readOneAlignment (
        TFile & file,
        FragmentStore<TSpec, TConfig> & fragStore,
        TContigAnchorGaps & contigAnchorGaps,
        TMatchMateInfos & matchMateInfos,
        TChar & c,
        FragStoreSAMContext<TFragStore> & contextSAM,
        FragStoreImportFlags const & importFlags)
    {
//IOREV _nodoc_
        // Basic types
        typedef FragmentStore<TSpec, TConfig>                                       TFragmentStore;
        typedef FragStoreSAMContext<TFragStore>                                     TSAMContext;
        typedef typename Id<TFragmentStore>::Type                                   TId;
        //typedef typename Size<TFragmentStore>::Type                                 TSize;

        // All fragment store element types
        //typedef typename Value<typename TFragmentStore::TContigStore>::Type         TContigElement;
        //typedef typename Value<typename TFragmentStore::TLibraryStore>::Type        TLibraryStoreElement;
        typedef typename Value<typename TFragmentStore::TMatePairStore>::Type       TMatePairElement;
        //typedef typename Value<typename TFragmentStore::TReadStore>::Type           TReadStoreElement;
        typedef typename Value<typename TFragmentStore::TAlignQualityStore>::Type   TAlignQualityElement;

        // Type for sequence in readstore
        typedef typename TFragmentStore::TReadSeq                                   TReadSeq;

        // Type for gap anchor
        typedef typename TFragmentStore::TContigPos                                 TContigPos;
        typedef Gaps<Nothing, AnchorGaps<typename TSAMContext::TContigAnchorGaps> > TContigGapsPW;

        // Type to temporarily store information about match mates
        typedef typename Value<TMatchMateInfos>::Type                               TMatchMateInfo;

        // read fields of alignments line
        _parseSkipWhitespace(file, c);

        // Read the query name.  The letters until the first
        // whitespace will be read into qname.  Then, we skip until we
        // hit the first tab character.
        String<char> qname;
        _parseReadSamIdentifier(file, qname, c);
        _parseSkipUntilChar(file, '\t', c);

        // read the flag
        int flag;
        flag = _parseReadNumber(file, c);
        _parseSkipWhitespace(file, c);
        bool reverse = (flag & (1 << 4)) == (1 << 4);

        // Read reference name.  Same behaviour as for query name:  Read up to
        // the first whitespace character and skip to next tab char.
        String<char> rname;
        _parseReadSamIdentifier(file, rname, c);
        _parseSkipUntilChar(file, '\t', c);

        // read begin position
        TContigPos beginPos;
        beginPos = _parseReadNumber(file, c);
        --beginPos; // Sam stores positions starting at 1 the fragment store starting at 0
        _parseSkipWhitespace(file, c);

        // read map quality
        TAlignQualityElement mapQ;
        mapQ.score = _parseReadNumber(file, c);
        _parseSkipWhitespace(file, c);

        // read CIGAR
        String<CigarElement<> > cigar;
        _parseReadCigar(file, cigar, c);
        _parseSkipWhitespace(file, c);
        
        // calculate the end position
        TContigPos endPos;
        _getLengthInRef(cigar, endPos);
        endPos = beginPos + endPos;

        // if the read is on the antisense strand switch begin and end position
        if (reverse)
        {
            TContigPos temp = beginPos;
            beginPos = endPos;
            endPos = temp;
        }

        // read mate reference name
        String<char> mrnm;
        _parseReadSamIdentifier(file, mrnm, c);
        _parseSkipWhitespace(file, c);

        // read mate position
        TContigPos mPos;
        mPos = _parseReadNumber(file, c);
        --mPos; // Sam stores positions starting at 1 the fragment store starting at 0
        _parseSkipWhitespace(file, c);

        // read iSize
        _parseReadNumber(file, c);
        _parseSkipWhitespace(file, c);

        // read in sequence
        TReadSeq readSeq;
        _parseReadDnaSeq(file, readSeq, c);
        _parseSkipWhitespace(file, c);

        // and associated qualities
        _parseReadSeqQual(file, readSeq, c);
        if (reverse)
            reverseComplement(readSeq);

        // read in Sam tags
        String<char> tags;
        _parseSkipSpace(file, c);
        _parseReadCharsUntilEndOfLine(file, tags, c);
        
        if (empty(qname) || empty(rname))
            return;
        
        // check if read sequence is already in the store.
        // if so get the ID, otherwise create new entries in the
        // read, read name and mate pair store
        
        contextSAM.readId = 0;
        if (importFlags.importRead)
        {
            if (!importFlags.importReadName)
                clear(qname);

            bool newRead = _storeAppendRead(fragStore, contextSAM.readId, qname, readSeq, flag, contextSAM);
            (void)newRead;
            
            SEQAN_ASSERT_NOT(newRead && empty(readSeq));
        }
        
        if (rname == "*")
            return; // stop here if read is not aligned
        
        // check if the contig is already in the store
        // get its ID or create a new one otherwise
        contextSAM.contigId = 0;
        _storeAppendContig(fragStore, contextSAM.contigId, rname);

        if (empty(cigar)) return;

        TId pairMatchId = 0;
        if (importFlags.importReadAlignment)
        {
            // generate gap anchor string for the read
            if (importFlags.importReadSeq)
                _samAppendAlignment(fragStore, fragStore.readSeqStore[contextSAM.readId], cigar, beginPos, endPos, pairMatchId, contextSAM);
            else
                _samAppendAlignmentWithoutSeq(fragStore, cigar, beginPos, endPos, pairMatchId, contextSAM);

            clear(contextSAM.contigGapAnchors);
            TContigGapsPW contigGaps(contextSAM.contigGapAnchors);
            cigarToGapAnchorContig(cigar, contigGaps);
            appendValue(contigAnchorGaps, contextSAM.contigGapAnchors);
        }
        
        if (importFlags.importReadAlignmentTags || importFlags.importReadAlignmentQuality)
        {
            // extract and delete some tags
            if (!empty(tags))
                for (unsigned pos = length(tags), right = length(tags); pos != 0; )
                {
                    --pos;
                    if (pos == 0 || tags[pos] == '\t')
                    {
                        unsigned left = pos;
                        if (tags[left] == '\t') ++left;                    

                        bool remove = false;
                        if (infix(tags, left, left + 2) == "MD")
                            remove = true;
                        
                        if (infix(tags, left, left + 2) == "NM")
                        {
                            std::string val;
                            int errors;
                            assign(val, infix(tags, left + 5, right)); // NM:i:x
                            std::istringstream stream(val);
                            stream >> errors;
                            mapQ.errors = errors;
                            remove = true;
                        }
                        
                        if (remove)
                            erase(tags, left, _min(right + 1, length(tags)));

                        right = pos;
                    }
                }
        }

        // create entries in Sam specific stores
        if (importFlags.importReadAlignmentQuality)
            appendValue(fragStore.alignQualityStore, mapQ, Generous());
        
        if (importFlags.importReadAlignmentTags)
            appendValue(fragStore.alignedReadTagStore, tags, Generous());

        // store additional data about match mate temporarily
        // used in the end of the read function to generate match mate IDs
        if (importFlags.importRead && importFlags.importReadAlignment && mrnm != "*")
        {
            TId mcontigId = contextSAM.contigId;
            if (mrnm != "=")
                _storeAppendContig(fragStore, mcontigId, mrnm);

            if (getMateNo(fragStore, contextSAM.readId) == 0)  // store mate info only for one mate
            {
                typename TMatePairElement::TId matePairId = TMatePairElement::INVALID_ID;
                if (contextSAM.readId < length(fragStore.readStore))
                    matePairId = fragStore.readStore[contextSAM.readId].matePairId;
                
                TMatchMateInfo matchMateInfo = {contextSAM.readId, mcontigId, pairMatchId, matePairId, (flag & 0x20) != 0, mPos};
                appendValue(matchMateInfos, matchMateInfo);
                back(fragStore.alignedReadStore).pairMatchId = pairMatchId;
            }
        }
    }


//////////////////////////////////////////////////////////////////////////////
// write functions for Sam
//////////////////////////////////////////////////////////////////////////////

// --------------------------------------------------------------------------
// Function _fillHeader()
// --------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void
_fillHeader(BamHeader & header,
            FragmentStore<TSpec, TConfig> & store)
{
    typedef FragmentStore<TSpec, TConfig>                       TFragmentStore;
    typedef typename TFragmentStore::TLibraryStore              TLibraryStore;
    typedef typename TFragmentStore::TContigStore               TContigStore;
    typedef typename TFragmentStore::TNameStore                 TNameStore;

    typedef typename Value<TContigStore>::Type                  TContig;
    typedef typename Iterator<TLibraryStore, Standard>::Type    TLibraryIter;
    typedef typename Iterator<TContigStore, Standard>::Type     TContigIter;
    typedef typename Iterator<TNameStore, Standard>::Type       TContigNameIter;
    typedef typename Id<TContig>::Type                          TId;

    typedef BamHeader::TSequenceInfo                            TSequenceInfo;
    typedef BamHeaderRecord::TTag                               TTag;

    // Fill first header line.
    BamHeaderRecord firstRecord;
    firstRecord.type = BAM_HEADER_FIRST;
    appendValue(firstRecord.tags, TTag("VN", "1.4"));
    appendValue(firstRecord.tags, TTag("SO", "unsorted"));
    appendValue(header.records, firstRecord);

    // Fill sequence info header line.
    TContigIter it          = begin(store.contigStore, Standard());
    TContigIter itEnd       = end(store.contigStore, Standard());
    TContigNameIter nit     = begin(store.contigNameStore, Standard());
    TContigNameIter nitEnd  = end(store.contigNameStore, Standard());

    for (; it != itEnd && nit != nitEnd; ++it, ++nit)
        appendValue(header.sequenceInfos, TSequenceInfo(*nit, length((*it).seq)));

    // Fill program header line.
    BamHeaderRecord pgRecord;
    pgRecord.type = BAM_HEADER_PROGRAM;
    appendValue(pgRecord.tags, TTag("ID", "SeqAn"));
    appendValue(header.records, pgRecord);

    // Fill library info header line.
    BamHeaderRecord rgRecord;
    rgRecord.type = BAM_HEADER_READ_GROUP;

    TLibraryIter lit    = begin(store.libraryStore, Standard());
    TLibraryIter litEnd = end(store.libraryStore, Standard());

    for (TId id = 0; lit != litEnd; ++lit, ++id)
    {
        appendValue(pgRecord.tags, TTag("ID", id + 1));
        appendValue(pgRecord.tags, TTag("LB", store.libraryNameStore[id]));
        appendValue(pgRecord.tags, TTag("PI", (int)store.libraryStore[id].mean));
        // Sample name needs to be included into fragment store.
        appendValue(pgRecord.tags, TTag("SM", "none"));
    }

}
    
//////////////////////////////////////////////////////////////////////////////
// _writeAlignments

    template<typename TFile, typename TSpec, typename TConfig, typename TAlignFunctor>
    inline void _writeAlignments(TFile & target,
                                 FragmentStore<TSpec, TConfig> & store,
                                 Sam,
                                 TAlignFunctor const & alignFunctor)
    {
//IOREV
        typedef FragmentStore<TSpec, TConfig>                           TFragmentStore;

        typedef typename TFragmentStore::TReadStore                     TReadStore;
        typedef typename TFragmentStore::TReadNameStore                 TReadNameStore;
        typedef typename TFragmentStore::TAlignedReadStore              TAlignedReadStore;
        typedef typename TFragmentStore::TContigStore                   TContigStore;
        typedef typename TFragmentStore::TReadSeq                       TReadSeq;

        typedef typename Value<TReadStore>::Type                        TRead;
        //typedef typename Size<TReadStore>::Type                         TSize;
        typedef typename Value<TContigStore>::Type                      TContig;
        typedef typename Value<TAlignedReadStore>::Type                 TAlignedRead;

        typedef typename Iterator<TAlignedReadStore, Standard>::Type    TAlignIter;
        typedef typename Iterator<TReadSeq, Standard>::Type             TReadSeqIter;
        typedef typename Id<TAlignedRead>::Type                         TId;

        typedef Gaps<Nothing, AnchorGaps<typename TContig::TGapAnchors> >   TContigGaps;

        String<int> mateIndex;  // store outer library size for each pair match (indexed by pairMatchId)
        calculateMateIndices(mateIndex, store);
        
        TAlignIter it = begin(store.alignedReadStore, Standard());
        TAlignIter itEnd = end(store.alignedReadStore, Standard());
        TAlignIter mit = it;
        CharString cigar, md;
        TReadSeq readSeq;
        
        typedef unsigned long TWord;
        const unsigned wordLen = BitsPerValue<TWord>::VALUE;
        String<TWord> readAligned;  // bitset to signal wether a read was aligned at least once
        resize(readAligned, (length(store.readStore) + wordLen - 1) / wordLen, (TWord)0);

        for(; it != itEnd; ++it)
        {
            TId alignedId = (*it).id;
            TId readId = (*it).readId;
            TId mateIdx = TRead::INVALID_ID;

            if ((*it).pairMatchId != TRead::INVALID_ID)
                mateIdx = mateIndex[2*(*it).pairMatchId + getMateNo(store, (*it).readId)];

            TContigGaps contigGaps(/*store.contigStore[(*it).contigId].seq, */store.contigStore[(*it).contigId].gaps);
            __int64 pos = positionGapToSeq(contigGaps, _min((*it).beginPos, (*it).endPos)) + 1;
            __int64 mpos = 0;
            int isize = 0;
            unsigned short flag = 0;

            if ((*it).beginPos > (*it).endPos)
                flag |= 0x0010;         

            // calculate flags, mpos, isize
            if (mateIdx < length(store.alignedReadStore))
            {
                mit = begin(store.alignedReadStore, Standard()) + mateIdx;
                if ((*it).contigId == (*mit).contigId)
                {
                    mpos = positionGapToSeq(contigGaps, _min((*mit).beginPos, (*mit).endPos)) + 1;
                    if ((*it).beginPos < (*mit).beginPos)
                        isize = positionGapToSeq(contigGaps, _max((*mit).beginPos, (*mit).endPos) - 1) + 2 - pos;
                    else
                        isize = mpos - positionGapToSeq(contigGaps, _max((*it).beginPos, (*it).endPos) - 1) - 2;
                }
                flag |= 0x0002;
                if ((*mit).beginPos > (*mit).endPos)
                    flag |= 0x0020;             
            }
            
            signed char mateNo = getMateNo(store, readId);
            if (mateNo == 0) flag |= 0x0040;    // this read is the first in the pair
            if (mateNo == 1) flag |= 0x0080;    // this read is the second in the pair
            
            // test for secondary alignment
            TWord mask = (TWord)1 << (readId % wordLen);
            bool secondary = (readAligned[readId / wordLen] & mask) != 0;
            readAligned[readId / wordLen] |= mask;
            if (secondary) flag |= 0x0100;      // we've already output an alignment for this read (this one is secondary)

            if (readId < length(store.readStore))
            {
                TRead &read = store.readStore[readId];
                if (read.matePairId != TRead::INVALID_ID)
                {
                    flag |= 0x0001;
                    if (mateIdx >= length(store.alignedReadStore))
                        flag |= 0x0008;         // mate is unmapped (actually we should check if the mate has no match at all)
                }
            }
            
            // <qname>
            if (readId < length(store.readNameStore))
            {
                typedef typename Iterator<typename Value<TReadNameStore>::Type const>::Type TReadNameIterator;
                if (empty(store.readNameStore[readId]))
                    continue;
                for (TReadNameIterator it = begin(store.readNameStore[readId]); it != end(store.readNameStore[readId]); ++it) {
                    if (*it == ' ' || *it == '\t' || *it == '\n' || *it == '\r')
                        break;
                    _streamPut(target, *it);
                }
            } else
                continue;
            _streamPut(target, '\t');
            
            // <flag>
            _streamPutInt(target, flag);
            _streamPut(target, '\t');
            
            // <rname>
            if ((*it).contigId < length(store.contigNameStore))
                _streamWrite(target, store.contigNameStore[(*it).contigId]);
            else
                _streamWrite(target, '.');  // No reference name given.  Standard says field must not be empty but gives no "NULL" value.
            _streamPut(target, '\t');
            
            // <pos>
            _streamPutInt(target, pos);
            _streamPut(target, '\t');
            
            // <mapq>
            if (alignedId < length(store.alignQualityStore) && store.alignQualityStore[alignedId].score > 0) // a value <= 0 equals -errors (forbidden in SAM)
                _streamPutInt(target, store.alignQualityStore[alignedId].score);
            else
                _streamPutInt(target, 255);
            _streamPut(target, '\t');
            
            // get read sequence
            if (readId < length(store.readSeqStore))
            {
                readSeq = store.readSeqStore[readId];
                if ((*it).beginPos > (*it).endPos) 
                    reverseComplement(readSeq);
            } else
                clear(readSeq);
            
            // <cigar>
            int errors = -1;
            if (alignedId < length(store.alignQualityStore))
                errors = store.alignQualityStore[alignedId].errors;
            // std::cout << store.readNameStore[readId] << "\n";
            alignAndGetCigarString(cigar, md, store.contigStore[(*it).contigId], readSeq, *it, errors, alignFunctor);
            _streamWrite(target, cigar);
			// _streamWrite(std::cout, cigar);
            // std::cout << "\n";
            _streamPut(target, '\t');
            
            // <mrnm>
            if ((mateIdx < length(store.alignedReadStore)))
            {
                if ((*it).contigId == (*mit).contigId)
                    _streamWrite(target, '=');
                else
                    if ((*mit).contigId < length(store.contigNameStore))
                        _streamWrite(target, store.contigNameStore[(*mit).contigId]);
            } else
                _streamWrite(target, '*');
                
            _streamPut(target, '\t');
            
            // <mpos>
            _streamPutInt(target, (int)mpos);
            _streamPut(target, '\t');
            
            // <isize>
            _streamPutInt(target, isize);
            _streamPut(target, '\t');

            // <seq>
            if (!secondary)
                _streamWrite(target, readSeq);
            else
                _streamPut(target, '*');
            _streamPut(target, '\t');
            
            // <qual>
            if (!secondary)
            {
                TReadSeqIter it = begin(readSeq, Standard());
                TReadSeqIter itEnd = end(readSeq, Standard());
                for (; it != itEnd; ++it)
                    _streamPut(target, (char)(getQualityValue(*it) + 33));
            }
            else
                _streamPut(target, '*');
            
            // <tags>
            
            if (errors != -1)
            {
                _streamWrite(target, "\tNM:i:");
                _streamPutInt(target, errors);
            }

            if (!empty(md))
            {
                _streamWrite(target, "\tMD:Z:");
                _streamWrite(target, md);
            }

            if (alignedId < length(store.alignedReadTagStore) && !empty(store.alignedReadTagStore[alignedId]))
            {
                _streamPut(target, '\t');
                _streamWrite(target, store.alignedReadTagStore[alignedId]);
            }
            
            _streamPut(target, '\n');
        }

#if 0  // DISABLE FOR NOW
        // Write out records for unaligned reads.
        TSize readCount = length(store.readSeqStore);
        for (unsigned readId = 0; readId < readCount; ++readId)
        {
            TWord mask = (TWord)1 << (readId % wordLen);
            if ((readAligned[readId / wordLen] & mask) != 0)
                continue;
            
            // <qname>
            if (readId < length(store.readNameStore)) {
                typedef typename Iterator<CharString, Standard>::Type TCharStringIterator;
                if (empty(store.readNameStore[readId]))
                    continue;
                for (TCharStringIterator it = begin(store.readNameStore[readId]); it != end(store.readNameStore[readId]); ++it) {
                    if (*it == ' ' || *it == '\t' || *it == '\n' || *it == '\r')
                        break;
                    _streamPut(target, *it);
                }
            } else
                continue;
            _streamPut(target, '\t');
            
            // <flag>
            unsigned short flag = 0x04;         // read is unaligned
            int mateNo = getMateNo(store, readId);
            if (mateNo == 0) flag |= 0x0041;    // this read is the first in the pair
            if (mateNo == 1) flag |= 0x0081;    // this read is the second in the pair

            _streamPutInt(target, flag);
            _streamPut(target, '\t');
            
            // <rname>
            _streamWrite(target, '*');
            _streamPut(target, '\t');
            
            // <pos>
            _streamWrite(target, '0');
            _streamPut(target, '\t');
            
            // <mapq>
            _streamWrite(target, '0');
            _streamPut(target, '\t');
            
            // <cigar>
            _streamWrite(target, '*');
            _streamPut(target, '\t');
            
            // <mrnm>
            _streamWrite(target, '*');
            _streamPut(target, '\t');
            
            // <mpos>
            _streamWrite(target, '0');
            _streamPut(target, '\t');
            
            // <isize>
            _streamWrite(target, '0');
            _streamPut(target, '\t');

            // <seq>
            readSeq = store.readSeqStore[readId];
            _streamWrite(target, readSeq);
            _streamPut(target, '\t');
            
            // <qual>
            TReadSeqIter it = begin(readSeq, Standard());
            TReadSeqIter itEnd = end(readSeq, Standard());
            for (; it != itEnd; ++it)
                _streamPut(target, (char)(getQualityValue(*it) + 33));

            // <tags>
            _streamWrite(target, "\tNH:i:0");
            
            _streamPut(target, '\n');
        }
#endif  // #if 0
    }
    
    template<typename TFile, typename TSpec, typename TConfig>
    inline void _writeAlignments(TFile & target,
                                 FragmentStore<TSpec, TConfig> & store,
                                 Sam)
    {
        _writeAlignments(target, store, Sam(), Nothing());
    }

//////////////////////////////////////////////////////////////////////////////
// write

///.Function.write.param.tag.type:Tag.File Format.tag.Sam
    
    template<typename TFile, typename TSpec, typename TConfig>
    inline void write(TFile & target,
                      FragmentStore<TSpec, TConfig> & store,
                      Sam)
    {
        typedef FragmentStore<TSpec, TConfig>               TFragmentStore;
        typedef typename TFragmentStore::TContigNameStore   TContigNameStore;
        typedef BamIOContext<TContigNameStore>              TBamIOContext;

        TBamIOContext context(store.contigNameStore, store.contigNameStoreCache);

        // 1. write header
        BamHeader header;

        // fill header with information from fragment store.
        _fillHeader(header, store);

        // write header to target.
        write2(target, header, context, Sam());

        // 2. write aligments
        _writeAlignments(target, store, Sam());
    }
    
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
