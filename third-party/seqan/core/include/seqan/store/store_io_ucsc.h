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

#ifndef SEQAN_HEADER_STORE_IO_UCSC_H
#define SEQAN_HEADER_STORE_IO_UCSC_H

/* IOREV
 *
 * _doc_
 *
 *
 * maybe move this to file/ because its a file format
 *
 */


namespace SEQAN_NAMESPACE_MAIN {

template <typename TSpec>
struct Ucsc_;

/**
.Tag.File Format.tag.Ucsc:
    Ucsc Genome Browser annotation file (a.k.a. knownGene format).
..include:seqan/store.h
*/

struct UcscKnownGene_;
typedef Tag<Ucsc_<UcscKnownGene_> > const Ucsc;

/**
.Tag.File Format.tag.UcscIsoforms:
    Ucsc Genome Browser isoforms file (a.k.a. knownIsoforms format).
..include:seqan/store.h
*/
struct UcscKnownIsoforms_;
typedef Tag<Ucsc_<UcscKnownIsoforms_> > const UcscIsoforms;

//////////////////////////////////////////////////////////////////////////////
// Read Ucsc
//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec = void>
struct IOContextUcsc_
{
    typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type      TAnnotation;
    typedef typename TAnnotation::TId                   TId;

    CharString      transName;
    CharString      contigName;
    __int64         cdsBegin;
    __int64         cdsEnd;
    String<__int64> exonBegin;
    String<__int64> exonEnd;
    CharString      proteinName;

    enum {KNOWN_GENE, KNOWN_ISOFORMS} format;
    TAnnotation annotation;
};

template <typename TFragmentStore, typename TSpec>
inline void clear(IOContextUcsc_<TFragmentStore, TSpec> & ctx)
{
    clear(ctx.transName);
    clear(ctx.contigName);
    clear(ctx.exonBegin);
    clear(ctx.exonEnd);
    clear(ctx.proteinName);
    clear(ctx.annotation.values);
}

//////////////////////////////////////////////////////////////////////////////
// _readOneAnnotation
//
// reads in one annotation line from a Gff file

// TODO(singer): return int instead of bool
template <typename TRecordReader, typename TFragmentStore, typename TSpec>
inline bool
_readOneAnnotation(
    TRecordReader & reader,
    IOContextUcsc_<TFragmentStore, TSpec> & ctx)
{
    typedef typename TFragmentStore::TContigPos         TContigPos;
    typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type      TAnnotation;

    clear(ctx);

    // read column 1: transcript name
    // The letters until the first whitespace will be read.
    // Then, we skip until we hit the first tab character.
    if (readUntilWhitespace(ctx.transName, reader))
        return false;

    if (!empty(ctx.transName) && ctx.transName[0] == '#')
    {
        if (skipLine(reader))
            return false;

        return false;
    }
    if (skipWhitespaces(reader))
        return false;

    // read column 2: contig name
    if (readUntilWhitespace(ctx.contigName, reader))
        return false;

    if (skipBlanks(reader))
        return false;

    // read column 3: orientation
    String<char> temp;
    if (readUntilWhitespace(temp, reader))
        return false;

    if (temp[0] != '+' && temp[0] != '-' && length(temp) == 1u)
    {
        ctx.format = ctx.KNOWN_ISOFORMS;
        insert(ctx.transName, 0, "GENE");
        if (skipLine(reader))
            return false;

        return true;
    }
    ctx.format = ctx.KNOWN_GENE;
    char orientation = temp[0];
    if (skipBlanks(reader))
        return false;


    // read column 4: transcript begin position
    clear(temp);
    if (readDigits(temp, reader))
        return false;

    if (length(temp) > 0u)
    {
        if (!lexicalCast2(ctx.annotation.beginPos, temp))
            return false;
    }
    else
    {
        ctx.annotation.beginPos = TAnnotation::INVALID_POS;
        if (skipUntilWhitespace(reader))
            return false;
    }
    if (skipBlanks(reader))
        return false;

    // read column 5: transcript end position
    clear(temp);
    if (readDigits(temp, reader))
        return false;

    if (length(temp) > 0u)
    {
        if (!lexicalCast2(ctx.annotation.endPos, temp))
            return false;
    }
    else
    {
        ctx.annotation.endPos = TAnnotation::INVALID_POS;
        if (skipUntilWhitespace(reader))
            return false;
    }
    if (skipBlanks(reader))
        return false;

    // read column 6: CDS begin position
    clear(temp);
    if (readDigits(temp, reader))
        return false;

    if (length(temp) > 0u)
    {
        if (!lexicalCast2(ctx.cdsBegin, temp))
            return false;
    }
    else
    {
        ctx.cdsBegin = TAnnotation::INVALID_POS;
        if (skipUntilWhitespace(reader))
            return false;
    }
    if (skipBlanks(reader))
        return false;

    // read column 7: CDS end position
    clear(temp);
    if (readDigits(temp, reader))
        return false;

    if (length(temp) > 0u)
    {
        if (!lexicalCast2(ctx.cdsEnd, temp))
            return false;
    }
    else
    {
        ctx.cdsEnd = TAnnotation::INVALID_POS;
        if (skipUntilWhitespace(reader))
            return false;
    }
    if (skipBlanks(reader))
        return false;

    // read column 8: exon count
    int exons = -1;
    clear(temp);
    if (readDigits(temp, reader))
        return false;

    if (length(temp) > 0u)
        if (!lexicalCast2(exons, temp))
            return false;

    if (skipBlanks(reader))
        return false;

    // read column 9: exon begin positions
    for (int i = 0; i < exons; ++i)
    {
        clear(temp);
        if (readDigits(temp, reader))
            return false;

        unsigned long long tempBegin;
        if (!lexicalCast2(tempBegin, temp))
            return false;

        appendValue(ctx.exonBegin, tempBegin, Generous());
        if (skipNCharsIgnoringWhitespace(reader, 1u))
            return false;

    }
    skipBlanks(reader);

    // read column 10: exon end positions
    for (int i = 0; i < exons; ++i)
    {
        clear(temp);
        if (readDigits(temp, reader))
            return false;

        unsigned long long tempEnd;
        if (!lexicalCast2(tempEnd, temp))
            return false;

        appendValue(ctx.exonEnd, tempEnd, Generous());
        if (skipNCharsIgnoringWhitespace(reader, 1u))
            return false;
    }
    if (skipUntilChar(reader, '\t'))
        return false;

    if (skipNChars(reader, 1u))
        return false;


    // read column 11: protein name
    if (readUntilWhitespace(ctx.proteinName, reader))
        return false;

    if (skipBlanks(reader))
        return false;

    // skip column 12
    if (skipLine(reader))
        return false;

    // adapt positions
    if (orientation == '-')
    {
        TContigPos tmp = ctx.annotation.beginPos;
        ctx.annotation.beginPos = ctx.annotation.endPos;
        ctx.annotation.endPos = tmp;
        tmp = ctx.cdsBegin;
        ctx.cdsBegin = ctx.cdsEnd;
        ctx.cdsEnd = tmp;
        for (int i = 0; i < exons; ++i)
        {
            tmp = ctx.exonBegin[i];
            ctx.exonBegin[i] = ctx.exonEnd[i];
            ctx.exonEnd[i] = tmp;
        }
    }

    return true;
}

template <typename TFragmentStore, typename TSpec>
inline void
_storeOneAnnotationKnownGene(
    TFragmentStore & fragStore,
    IOContextUcsc_<TFragmentStore, TSpec> & ctx)
{
    typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type      TAnnotation;
    typedef typename TAnnotation::TId                   TId;

    SEQAN_ASSERT_EQ(length(fragStore.annotationStore), length(fragStore.annotationNameStore));

    // add transcript and CDS
    TId transId = TAnnotation::INVALID_ID;
    _storeAppendAnnotationName(fragStore, transId, ctx.transName, (TId) TFragmentStore::ANNO_MRNA);
    TId cdsId = length(fragStore.annotationStore);
    appendName(fragStore.annotationNameStore, ctx.proteinName, fragStore.annotationNameStoreCache);

    resize(fragStore.annotationStore, cdsId + 1 + length(ctx.exonBegin), Generous());
    resize(fragStore.annotationNameStore, cdsId + 1 + length(ctx.exonBegin), Generous());

    // add contig name
    _storeAppendContig(fragStore, ctx.annotation.contigId, ctx.contigName);

    TAnnotation & transcript = fragStore.annotationStore[transId];
    TId geneId = transcript.parentId;
    if (geneId == TAnnotation::INVALID_ID)
        geneId = 0;
    transcript = ctx.annotation;
    transcript.parentId = geneId;
    transcript.typeId = TFragmentStore::ANNO_MRNA;

    TAnnotation & cds = fragStore.annotationStore[cdsId];
    cds = ctx.annotation;
    cds.parentId = transId;
    cds.typeId = TFragmentStore::ANNO_CDS;
    cds.beginPos = ctx.cdsBegin;
    cds.endPos = ctx.cdsEnd;
    _adjustParent(transcript, cds);

    // insert exons
    ctx.annotation.parentId = transId;
    ctx.annotation.typeId = TFragmentStore::ANNO_EXON;
    for (unsigned i = 0; i < length(ctx.exonBegin); ++i)
    {
        ctx.annotation.beginPos = ctx.exonBegin[i];
        ctx.annotation.endPos = ctx.exonEnd[i];
        fragStore.annotationStore[cdsId + 1 + i] = ctx.annotation;
        _adjustParent(transcript, ctx.annotation);
    }
    if (geneId != 0)
        _adjustParent(fragStore.annotationStore[geneId], transcript);
}

template <typename TFragmentStore, typename TSpec>
inline void
_storeOneAnnotationKnownIsoforms(
    TFragmentStore & fragStore,
    IOContextUcsc_<TFragmentStore, TSpec> & ctx)
{
    typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type      TAnnotation;
    typedef typename TAnnotation::TId                   TId;

    SEQAN_ASSERT_EQ(length(fragStore.annotationStore), length(fragStore.annotationNameStore));

    TId geneId = TAnnotation::INVALID_ID;
    TId transId = TAnnotation::INVALID_ID;

    // add transcript and CDS
    _storeAppendAnnotationName(fragStore, geneId, ctx.transName, (TId) TFragmentStore::ANNO_GENE);
    _storeAppendAnnotationName(fragStore, transId, ctx.contigName, (TId) TFragmentStore::ANNO_MRNA);

    // set parent link locus->root
    TAnnotation & locus = fragStore.annotationStore[geneId];
    locus.parentId = 0;
    locus.typeId = TFragmentStore::ANNO_GENE;

    // set parent link transcript->locus
    TAnnotation & transcript = fragStore.annotationStore[transId];
    transcript.parentId = geneId;
    transcript.typeId = TFragmentStore::ANNO_MRNA;

    _adjustParent(locus, transcript);
}

template <typename TFragmentStore, typename TSpec>
inline void
_storeOneAnnotation(
    TFragmentStore & fragStore,
    IOContextUcsc_<TFragmentStore, TSpec> & ctx)
{
    if (ctx.format == ctx.KNOWN_GENE)
        _storeOneAnnotationKnownGene(fragStore, ctx);
    else
        _storeOneAnnotationKnownIsoforms(fragStore, ctx);
}

// TODO(holtgrew): Change interface such that file is after store.
// TODO(singer): all the other read functions get a RecordReader not a file.
template <typename TFile, typename TSpec, typename TConfig, typename TFormatSpec>
inline void
read(
    TFile & file,
    FragmentStore<TSpec, TConfig> & fragStore,
    Tag<Ucsc_<TFormatSpec> > const)
{
//IOREV _nodoc_
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;

    if (streamEof(file))
        return;

    // get first character from the stream
    IOContextUcsc_<TFragmentStore> ctx;

    refresh(fragStore.contigNameStoreCache);
    refresh(fragStore.annotationNameStoreCache);
    refresh(fragStore.annotationTypeStoreCache);

    RecordReader<TFile, SinglePass<> > reader(file);

    while (!atEnd(reader))
    {
        if (_readOneAnnotation(reader, ctx))
            _storeOneAnnotation(fragStore, ctx);
    }
    _storeClearAnnoBackLinks(fragStore.annotationStore);
    _storeCreateAnnoBackLinks(fragStore.annotationStore);
    _storeRemoveTempAnnoNames(fragStore);
}

//////////////////////////////////////////////////////////////////////////////
// Write Ucsc
//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec, typename TAnnotation, typename TId>
inline bool
_retrieveOneAnnotation(
    TFragmentStore & fragStore,
    IOContextUcsc_<TFragmentStore, TSpec> & ctx,
    TAnnotation & annotation,
    TId id,
    Ucsc)
{
    if (annotation.typeId != TFragmentStore::ANNO_MRNA)
        return false;

    ctx.format = ctx.KNOWN_GENE;
    ctx.transName = getAnnoUniqueName(fragStore, id);
    if (annotation.contigId < length(fragStore.contigNameStore))
        ctx.contigName = fragStore.contigNameStore[annotation.contigId];
    else
        clear(ctx.contigName);

    ctx.annotation = annotation;
    clear(ctx.proteinName);
    clear(ctx.exonBegin);
    clear(ctx.exonEnd);

    TId lastChildId = annotation.lastChildId;
    TId i = lastChildId;
    do
    {
        i = fragStore.annotationStore[i].nextSiblingId;
        TAnnotation & anno = fragStore.annotationStore[i];
        if (anno.typeId == TFragmentStore::ANNO_CDS)
        {
            if (i < length(fragStore.annotationNameStore))
                ctx.proteinName = fragStore.annotationNameStore[i];
            ctx.cdsBegin = anno.beginPos;
            ctx.cdsEnd = anno.endPos;
        }
        if (anno.typeId == TFragmentStore::ANNO_EXON)
        {
            appendValue(ctx.exonBegin, anno.beginPos, Generous());
            appendValue(ctx.exonEnd, anno.endPos, Generous());
        }
    }
    while (i != lastChildId);
    return true;
}

template <typename TFragmentStore, typename TSpec, typename TAnnotation, typename TId>
inline bool
_retrieveOneAnnotation(
    TFragmentStore & fragStore,
    IOContextUcsc_<TFragmentStore, TSpec> & ctx,
    TAnnotation & annotation,
    TId id,
    UcscIsoforms)
{
    if (annotation.typeId != TFragmentStore::ANNO_MRNA)
        return false;

    if (annotation.parentId == TAnnotation::INVALID_ID || annotation.parentId == 0)
        return false;

    ctx.format = ctx.KNOWN_ISOFORMS;
    ctx.transName = getAnnoUniqueName(fragStore, annotation.parentId);
    ctx.contigName = getAnnoUniqueName(fragStore, id);
    return true;
}

template <typename TTargetStream, typename TFragmentStore, typename TSpec>
inline bool
_writeOneAnnotation(
    TTargetStream & file,
    IOContextUcsc_<TFragmentStore, TSpec> & ctx)
{
    typedef typename TFragmentStore::TContigPos         TContigPos;

    unsigned suf = 0;
    if (ctx.format == ctx.KNOWN_ISOFORMS && length(ctx.transName) >= 4 && prefix(ctx.transName, 4) == "GENE")
        suf = 4;

    // read column 1: transcript name
    // The letters until the first whitespace will be read.
    // Then, we skip until we hit the first tab character.
    if (length(suffix(ctx.transName, suf)) > 0u)
    {
        if (streamWriteBlock(file, &(suffix(ctx.transName, suf))[0], length(suffix(ctx.transName, suf))) < length(suffix(ctx.transName, suf)))
            return false;
    }

    if (streamPut(file, '\t'))
        return false;

    // read column 2: contig name
    if (length(ctx.contigName) > 0u)
    {
        if (streamWriteBlock(file, &ctx.contigName[0], length(ctx.contigName)) < length(ctx.contigName))
            return false;
    }

    if (ctx.format == ctx.KNOWN_ISOFORMS)
    {
        if (streamWriteChar(file, '\n'))
            return false;

        return true;
    }
    if (streamWriteChar(file, '\t'))
        return false;

    // read column 3: orientation
    TContigPos transBeginPos, transEndPos;
    TContigPos cdsBeginPos, cdsEndPos;
    if (ctx.annotation.beginPos < ctx.annotation.endPos)
    {
        if (streamWriteChar(file, '+'))
            return false;

        transBeginPos = ctx.annotation.beginPos;
        transEndPos = ctx.annotation.endPos;
        cdsBeginPos = ctx.cdsBegin;
        cdsEndPos = ctx.cdsEnd;
    }
    else
    {
        if (streamWriteChar(file, '-'))
            return false;

        transEndPos = ctx.annotation.beginPos;
        transBeginPos = ctx.annotation.endPos;
        cdsEndPos = ctx.cdsBegin;
        cdsBeginPos = ctx.cdsEnd;
    }
    if (streamWriteChar(file, '\t'))
        return false;

    // read column 4: transcript begin position
    if (streamPut(file, transBeginPos))
        return false;

    if (streamWriteChar(file, '\t'))
        return false;

    // read column 5: transcript end position
    if (streamPut(file, transEndPos))
        return false;

    if (streamWriteChar(file, '\t'))
        return false;

    // read column 6: CDS begin position
    if (streamPut(file, cdsBeginPos))
        return false;

    if (streamWriteChar(file, '\t'))
        return false;

    // read column 7: CDS end position
    if (streamPut(file, cdsEndPos))
        return false;

    if (streamWriteChar(file, '\t'))
        return false;

    // read column 8: exon count
    if (streamPut(file, length(ctx.exonBegin)))
        return false;

    if (streamWriteChar(file, '\t'))
        return false;

    // read column 9: exon begin positions
    for (unsigned i = 0; i < length(ctx.exonBegin); ++i)
    {
        if (streamPut(file, _min(ctx.exonBegin[i], ctx.exonEnd[i])))
            return false;

        if (streamWriteChar(file, ','))
            return false;
    }
    if (streamWriteChar(file, '\t'))
        return false;

    // read column 10: exon end positions
    for (unsigned i = 0; i < length(ctx.exonBegin); ++i)
    {
        if (streamPut(file, _max(ctx.exonBegin[i], ctx.exonEnd[i])))
            return false;

        if (streamWriteChar(file, ','))
            return false;
    }
    if (streamWriteChar(file, '\t'))
        return false;

    // read column 10: protein name
    if (length(ctx.proteinName) > 0u)
    {
        if (streamWriteBlock(file, &ctx.proteinName[0], length(ctx.proteinName)) < length(ctx.proteinName))
            return false;
    }

    if (streamWriteChar(file, '\t'))
        return false;

    // skip column 11
    if (length(ctx.transName) > 0u)
    {
        if (streamWriteBlock(file, &ctx.transName[0], length(ctx.transName)) < length(ctx.transName))
            return false;
    }

    if (streamWriteChar(file, '\n'))
        return false;

    return true;
}

template <typename TTargetStream, typename TSpec, typename TConfig, typename TFormatSpec>
inline void
write(
    TTargetStream & target,
    FragmentStore<TSpec, TConfig> & store,
    Tag<Ucsc_<TFormatSpec> > const format)
{
//IOREV _nodoc_
    typedef FragmentStore<TSpec, TConfig>                           TFragmentStore;
    typedef typename TFragmentStore::TAnnotationStore               TAnnotationStore;
    typedef typename Value<TAnnotationStore>::Type                  TAnnotation;
    typedef typename Iterator<TAnnotationStore, Standard>::Type     TAnnoIter;
    typedef typename Id<TAnnotation>::Type                          TId;

    IOContextUcsc_<TFragmentStore> ctx;

    TAnnoIter it = begin(store.annotationStore, Standard());
    TAnnoIter itEnd = end(store.annotationStore, Standard());

    for (TId id = 0; it != itEnd; ++it, ++id)
    {
        if (_retrieveOneAnnotation(store, ctx, *it, id, format))
            _writeOneAnnotation(target, ctx);
    }
}

} // namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
